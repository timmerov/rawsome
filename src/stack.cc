/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
register two source images.
**/

#include "canon.h"
#include "image.h"
#include "log.h"
#include "stack.h"

#include <cmath>
#include <sstream>

namespace {

/** everything **/
/*const int kFirstTag = 1227;
const int kLastTag = 1953;*/

/** first skip is 1251 -> 1252 **/
/*const int kFirstTag = 1250;
const int kLastTag = 1253;*/

/** up to the first skip **/
/*const int kFirstTag = 1227;
const int kLastTag = 1251;*/

/** first only **/
const int kFirstTag = 4117;
const int kLastTag = 4126;
//const int kLastTag = 4117;
//const int kFirstTag = 1227;
//const int kLastTag = 1227;

const char *kInputFile = "/home/timmer/Pictures/2022-03-29/stackable/IMG_";
const char *kOutputFile = "/home/timmer/Pictures/2022-03-29/stackable/stack.png";
const char *kBlackFile = "/home/timmer/Pictures/2022-03-29/stackable/black.rsm";
//const char *kInputFile = "/home/timmer/Pictures/2021-03-03/mars-pleiades/IMG_";
//const char *kOutputFile = "/home/timmer/Pictures/2021-03-03/mars-pleiades/stack.png";
//const char *kBlackFile = "/home/timmer/Pictures/2021-03-03/black/black.rsm";

const int kNoiseFloor = 500;
const int kMaxRegisterOffset = 300;

const int kStarFloor = 200;
const int kMaxStarRadius = 10;

class Luminage {
public:
    Luminage() = default;
    ~Luminage() = default;

    Image image_;
    Plane luminance_;
    std::vector<int> horz_;
    std::vector<int> vert_;
};

class StackImages {
public:
    StackImages() = default;
    ~StackImages() = default;

    Image black_;
    Luminage prev_;
    Luminage cur_;
    Image stack_;

    int noise_ = 0;
    int count_ = 0;
    int dx_ = 0;
    int dy_ = 0;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("stack images - work in progress.");

        /** load the black file. **/
        black_.load_rawsome(kBlackFile);

        /** load the first image. **/
        load_image(kFirstTag);
        cur_.image_.camera_.print();

        /** copy it to the stacked image. **/
        stack_ = cur_.image_;

        /** choose an arbitrary noise value. **/
        noise_ = kNoiseFloor;

        /** one image loaded so far. **/
        count_ = 1;

        /** do its sums **/
        sum_horz_vert_luminance();

        /** find stars **/
        find_stars();

        /** save current as previous. **/
        prev_ = std::move(cur_);

        /** stack the next images. **/
        for (int i = kFirstTag+1; i <= kLastTag; ++i) {
            /** load new current. **/
            load_image(i);
            if (cur_.image_.is_loaded_ == false) {
                return 1;
            }
            register_fast();
            stack_with_offset();

            /** save current as previous. **/
            prev_ = std::move(cur_);
        }

        /** scale the stacked image. **/
        if (count_ == 0) {
            LOG("no images.");
            return 1;
        }
        double brightness = 10.0;
        double factor = brightness / double(count_);
        stack_.planes_.multiply4(factor);

        /** save the stacked image. **/
        white_balance();
        combine_greens();
        convert_to_srgb();
        apply_user_gamma();
        apply_display_gamma();
        stack_.planes_.scale_to_8bits();
        stack_.save_png(kOutputFile);

        LOG("success.");
        return 0;
    }

    void load_image(
        int tag
    ) {
        std::stringstream ss;
        ss<<kInputFile<<tag<<".CR2";
        cur_.image_.load_raw(ss.str().c_str());
        if (cur_.image_.is_loaded_ == false) {
            return;
        }
        subtract_black();
        compute_luminance(cur_.image_.planes_, cur_.luminance_);
    }

    void subtract_black() {
        int sz = cur_.image_.planes_.r_.width_ * cur_.image_.planes_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            cur_.image_.planes_.r_.samples_[i] -= black_.planes_.r_.samples_[i];
            cur_.image_.planes_.g1_.samples_[i] -= black_.planes_.g1_.samples_[i];
            cur_.image_.planes_.g2_.samples_[i] -= black_.planes_.g2_.samples_[i];
            cur_.image_.planes_.b_.samples_[i] -= black_.planes_.b_.samples_[i];
        }
    }

    void register_fast() {
        LOG("fast registration...");
        sum_horz_vert_luminance();
        register_horz();
        register_vert();
    }

    void sum_horz_vert_luminance() {
        int wd = cur_.luminance_.width_;
        int ht = cur_.luminance_.height_;

        cur_.horz_.clear();
        cur_.vert_.clear();
        cur_.horz_.resize(wd, 0);
        cur_.vert_.resize(ht, 0);

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int val = cur_.luminance_.get(x, y);
                if (val >= noise_) {
                    cur_.horz_[x] += val;
                    cur_.vert_[y] += val;
                }
            }
        }

        /*std::stringstream ss;
        ss<<"horz: ["<<std::endl;
        for (int x = 0; x < wd; ++x) {
            ss<<" "<<cur_.horz_[x]<<std::endl;
        }
        ss<<" ]";
        LOG(ss.str());

        ss.str("");
        ss.clear();
        ss<<"vert: [";
        for (int y = 0; y < ht; ++y) {
            ss<<" "<<cur_.vert_[y];
        }
        ss<<" ]";
        LOG(ss.str());*/
    }

    class Star {
    public:
        Star() = default;
        ~Star() = default;

        double x_ = 0.0;
        double y_ = 0.0;
        int brightness_ = 0;
    };

    class FindStars {
    public:
        FindStars() = default;
        ~FindStars() = default;

        Plane luminance_;
        std::vector<std::int64_t> horz_;
        std::vector<std::int64_t> vert_;
        int xmax_ = 0;
        int ymax_ = 0;
        int radius_ = 0;
        std::vector<Star> stars_;

        void init(
            Plane& luminance
        ) {
            luminance_ = luminance;

            int wd = luminance_.width_;
            int ht = luminance_.height_;
            horz_.resize(wd);
            vert_.resize(ht);

            reinit_horz(0, wd);
            reinit_vert(0, ht);
        }

        void reinit_horz(
            int x0,
            int x1
        ) {
            int wd = luminance_.width_;
            int ht = luminance_.height_;
            int sz = wd * ht;

            for (int x = x0; x < x1; ++x) {
                horz_[x] = 0;
            }

            for (int y = 0; y < ht; ++y) {
                for (int x = x0; x < x1; ++x) {
                    int c16 = luminance_.get(x, y);
                    std::int64_t c = pin_to_16bits(c16);
                    /** unique-ify luminance. **/
                    c = c*sz + y*wd + x;
                    horz_[x] = std::max(horz_[x], c);
                }
            }
        }

        void reinit_vert(
            int y0,
            int y1
        ) {
            int wd = luminance_.width_;
            int ht = luminance_.height_;
            int sz = wd * ht;

            for (int y = y0; y < y1; ++y) {
                vert_[y] = 0;
            }

            for (int y = y0; y < y1; ++y) {
                for (int x = 0; x < wd; ++x) {
                    int c16 = luminance_.get(x, y);
                    std::int64_t c = pin_to_16bits(c16);
                    /** unique-ify luminance. **/
                    c = c*sz + y*wd + x;
                    vert_[y] = std::max(vert_[y], c);
                }
            }
        }

        int find_brightest() {
            int wd = luminance_.width_;
            int ht = luminance_.height_;
            int sz = wd * ht;

            xmax_ = 0;
            ymax_ = 0;

            std::int64_t cmax = 0;
            for (int x = 0; x < wd; ++x) {
                std::int64_t c = horz_[x];
                if (cmax < c) {
                    cmax = c;
                    xmax_ = x;
                }
            }
            cmax = 0;
            for (int y = 0; y < ht; ++y) {
                std::int64_t c = vert_[y];
                if (cmax < c) {
                    cmax = c;
                    ymax_ = y;
                }
            }
            LOG("max x,y="<<xmax_<<","<<ymax_);

            int brightness = cmax / sz;
            return brightness;
        }

        void brightest_centroid() {
            std::int64_t sumc = luminance_.get(xmax_, ymax_);
            std::int64_t sumx = sumc * xmax_;
            std::int64_t sumy = sumc * ymax_;
            int pixels = 1;
            radius_ = kMaxStarRadius;
            for (int r = 1; r <= kMaxStarRadius; ++r) {
                bool again = false;
                for (int y = ymax_ - r; y <= ymax_ + r; y += 2*r) {
                    for (int x = xmax_ - r; x <= xmax_ + r; ++x) {
                        int c = luminance_.get(x, y);
                        if (c >= kStarFloor) {
                            ++pixels;
                            sumc += c;
                            sumx += c * x;
                            sumy += c * y;
                            again = true;
                        }
                    }
                }
                for (int y = ymax_ - r + 1; y <= ymax_ + r - 1; ++y) {
                    for (int x = xmax_ - r; x <= xmax_ + r; x += 2*r) {
                        int c = luminance_.get(x, y);
                        if (c >= kStarFloor) {
                            ++pixels;
                            sumc += c;
                            sumx += c * x;
                            sumy += c * y;
                            again = true;
                        }
                    }
                }
                //LOG("r="<<r<<" sumc="<<sumc<<" sumx="<<sumx<<" sumy="<<sumy<<" pixels="<<pixels);
                if (again == false) {
                    radius_ = r;
                    break;
                }
            }
            double x = double(sumx) / double(sumc);
            double y = double(sumy) / double(sumc);
            LOG("c="<<sumc<<" x,y="<<x<<","<<y<<" pixels="<<pixels<<" radius="<<radius_);

            auto& star = stars_.emplace_back();
            star.x_ = x;
            star.y_ = y;
            star.brightness_ = sumc;
        }

        void erase_brightest() {
            int x0 = xmax_ - radius_;
            int x1 = xmax_ + radius_ + 1;
            int y0 = ymax_ - radius_;
            int y1 = ymax_ + radius_ + 1;

            for (int y = y0; y < y1; ++y) {
                for (int x = x0; x < x1; ++x) {
                    luminance_.set(x, y, 0);
                }
            }

            reinit_horz(x0, x1);
            reinit_vert(y0, y1);
        }
    };

    void find_stars() {
        FindStars fs;
        fs.init(cur_.luminance_);
        for (int i = 0; i < 30; ++i) {
            int brightness = fs.find_brightest();
            if (brightness < kStarFloor) {
                break;
            }
            fs.brightest_centroid();
            fs.erase_brightest();
        }
        int nstars = fs.stars_.size();
        LOG("nstars="<<nstars);
    }

    void register_horz() {
        int min_diff = 0;
        int max_diff = 0;
        std::int64_t sum_diff = 0;
        int min_dx = 0;

        int wd = cur_.horz_.size();
        for (int dx = - kMaxRegisterOffset; dx <= kMaxRegisterOffset; ++dx) {
            int sum = 0;
            for (int x1 = 0; x1 < wd; ++x1) {
                int x2 = x1 + dx;
                if (x2 < 0) {
                    continue;
                }
                if (x2 >= wd) {
                    break;
                }
                int c1 = cur_.horz_[x1];
                int c2 = prev_.horz_[x2];
                int dc = std::abs(c1 - c2);
                sum += dc;
            }

            if (min_diff == 0 || min_diff > sum) {
                min_diff = sum;
                min_dx = dx;
                //LOG("dx="<<dx<<" sum="<<sum);
            }
            max_diff = std::max(max_diff, sum);
            sum_diff += sum;
        }

        int n = kMaxRegisterOffset*2 + 1;
        double avg_diff = double(sum_diff) / double(n);
        LOG("dx="<<min_dx<<" min="<<min_diff<<" max="<<max_diff<<" avg="<<avg_diff);

        dx_ = min_dx;
    }

    void register_vert() {
        int min_diff = 0;
        int max_diff = 0;
        std::int64_t sum_diff = 0;
        int min_dy = 0;

        int ht = cur_.vert_.size();
        for (int dy = - kMaxRegisterOffset; dy <= kMaxRegisterOffset; ++dy) {
            int sum = 0;
            for (int y1 = 0; y1 < ht; ++y1) {
                int y2 = y1 + dy;
                if (y2 < 0) {
                    continue;
                }
                if (y2 >= ht) {
                    break;
                }
                int c1 = cur_.vert_[y1];
                int c2 = prev_.vert_[y2];
                int dc = std::abs(c1 - c2);
                sum += dc;
            }

            if (min_diff == 0 || min_diff > sum) {
                min_diff = sum;
                min_dy = dy;
                //LOG("dx="<<dx<<" sum="<<sum);
            }
            max_diff = std::max(max_diff, sum);
            sum_diff += sum;
        }

        int n = kMaxRegisterOffset*2 + 1;
        double avg_diff = double(sum_diff) / double(n);
        LOG("dy="<<min_dy<<" min="<<min_diff<<" max="<<max_diff<<" avg="<<avg_diff);

        dy_ = min_dy;
    }

    void stack_with_offset() {
        int wd = cur_.image_.planes_.r_.width_;
        int ht = cur_.image_.planes_.r_.height_;

        Planes dst;
        dst.init(wd, ht);

        for (int y1 = 0; y1 < ht; ++y1) {
            int y2 = y1 + dy_;
            if (y2 < 0) {
                continue;
            }
            if (y2 >= ht) {
                break;
            }
            for (int x1 = 0; x1 < wd; ++x1) {
                int x2 = x1 + dx_;
                if (x2 < 0) {
                    continue;
                }
                if (x2 >= wd) {
                    break;
                }
                int r = cur_.image_.planes_.r_.get(x1, y1);
                int g1 = cur_.image_.planes_.g1_.get(x1, y1);
                int g2 = cur_.image_.planes_.g2_.get(x1, y1);
                int b = cur_.image_.planes_.b_.get(x1, y1);
                r += stack_.planes_.r_.get(x2, y2);
                g1 += stack_.planes_.g1_.get(x2, y2);
                g2 += stack_.planes_.g2_.get(x2, y2);
                b += stack_.planes_.b_.get(x2, y2);
                dst.r_.set(x1, y1, r);
                dst.g1_.set(x1, y1, g1);
                dst.g2_.set(x1, y1, g2);
                dst.b_.set(x1, y1, b);
            }
        }

        stack_.planes_ = std::move(dst);
        ++count_;
    }

    void white_balance() {
        RggbDouble cam_mul;
        cam_mul.r_ = stack_.camera_.wb_r_;
        cam_mul.g1_ = 1.0;
        cam_mul.g2_ = 1.0;
        cam_mul.b_ = stack_.camera_.wb_b_;
        LOG("white balance: R="<<cam_mul.r_<<" B="<<cam_mul.b_);

        /** white balance. **/
        stack_.planes_.multiply(cam_mul);
    }

    void combine_greens() {
        auto& planes = stack_.planes_;
        int sz = planes.g1_.width_ * planes.g1_.height_;
        for (int i = 0; i < sz; ++i) {
            int g1 = planes.g1_.samples_[i];
            int g2 = planes.g2_.samples_[i];
            planes.g1_.samples_[i] = (g1 + g2 + 1)/2;
        }
    }

    void convert_to_srgb() {
        /**
        we are currently in canon bayer rgb space.
        convert to srgb space by matrix multiplication.
        the matrix is hard-coded.
        dcraw derives it roundaboutedly.
        **/
        double mat[3][3] = {
            {+1.901824, -0.972035, +0.070211},
            {-0.229410, +1.659384, -0.429974},
            {+0.042001, -0.519143, +1.477141}
        };

        double brightness = 1.3;
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                mat[i][k] *= brightness;
            }
        }

        auto& planes = stack_.planes_;
        int sz = planes.r_.width_ * planes.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = planes.r_.samples_[i];
            int in_g = planes.g1_.samples_[i];
            int in_b = planes.b_.samples_[i];

            /** transform by matrix multiplication. **/
            double out_r = mat[0][0]*in_r + mat[0][1]*in_g + mat[0][2]*in_b;
            double out_g = mat[1][0]*in_r + mat[1][1]*in_g + mat[1][2]*in_b;
            double out_b = mat[2][0]*in_r + mat[2][1]*in_g + mat[2][2]*in_b;

            /** overwrite old values. **/
            planes.r_.samples_[i] = out_r;
            planes.g1_.samples_[i] = out_g;
            planes.b_.samples_[i] = out_b;
        }
    }
    void apply_user_gamma() {
        double pwr = 0.6;
        stack_.planes_.apply_user_gamma(pwr);
    }

    void apply_display_gamma() {
        double pwr = 1.0 / stack_.camera_.gamma0_;
        double ts = stack_.camera_.gamma1_;
        int white = 0x10000;
        stack_.planes_.apply_display_gamma(pwr, ts, white);
    }
};
}

int stack_images(
    int argc,
    clo_argv_t argv
) {

    StackImages stack;
    int exit_code = stack.run(argc, argv);

    return exit_code;
}
