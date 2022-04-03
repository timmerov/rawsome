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

const int kFirstTag = 4154;
//const int kLastTag = 4156;
const int kLastTag = 4261;

const char *kInputFile = "/home/timmer/Pictures/2022-04-01/stackable/IMG_";
const char *kOutputFile = "/home/timmer/Pictures/2022-04-01/stackable/stack.png";

const int kNoiseFloor = 500;
const int kMaxRegisterOffset = 50;

class RegisterImage {
public:
    RegisterImage() = default;
    ~RegisterImage() = default;

    Image image_;
    std::vector<int> horz_;
    std::vector<int> vert_;
};

class StackImages {
public:
    StackImages() = default;
    ~StackImages() = default;

    RegisterImage cur_;
    RegisterImage stack_;
    Plane luminance_;

    int noise_ = 0;
    int count_ = 0;
    int dx_ = 0;
    int dy_ = 0;

    /**
    libraw gives us a maximum value=16383.
    however, the actual saturation point is less than that.
    it's around 13583.
    but we will need to determine the exact values.
    dump pixels into a histogram centered around 13583.
    if there are no saturated pixels then we should have a long tail.
    otherwise, there will be a spike at the end of the tail.
    ie we need to find the thagomizer.
    **/
    const int kSaturation = 13583 - 2046; /** we need to subtract black. **/

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("stack images - work in progress.");

        /** load the first image. **/
        load_image(kFirstTag);
        cur_.image_.camera_.print();

        /** choose an arbitrary noise value. **/
        noise_ = kNoiseFloor;

        /** one image loaded so far. **/
        count_ = 1;

        /** do its sums **/
        sum_horz_vert_luminance();

        /** move the current image to the stacked image. **/
        stack_ = std::move(cur_);

        /** find stars **/
        //find_stars();

        /** stack the next images. **/
        for (int i = kFirstTag+1; i <= kLastTag; ++i) {
            /** load new current. **/
            load_image(i);
            if (cur_.image_.is_loaded_ == false) {
                LOG("skipped.");
                continue;
            }
            register_fast();
            stack_with_offset();
            ++count_;
        }

        /** scale the stacked image. **/
        if (count_ == 0) {
            LOG("no images.");
            return 1;
        }

        LOG("stacked images:");
        analyze_image4(stack_.image_);

        /** scale the brightness sum. **/
        /** this is the desired full scale. **/
        double factor = 65535.0;
        /** this is the single image max value. **/
        factor /= double(kSaturation);
        /** average the stacked images. **/
        factor /= double(count_);
        stack_.image_.planes_.multiply4(factor);

        LOG("stacked scaled:");
        analyze_image4(stack_.image_);

        /** hack! remove more black. **/
        /*int sky_black = 5*256;
        RggbPixel pixel;
        pixel.r_ = sky_black;
        pixel.g1_ = sky_black;
        pixel.g2_ = sky_black;
        pixel.b_ = sky_black;
        stack_.image_.planes_.subtract(pixel);*/

        /** save the stacked image. **/
        combine_greens();
        LOG("stacked combine greens:");
        analyze_image3(stack_.image_);

        apply_user_gamma();
        LOG("stacked user gamma:");
        analyze_image3(stack_.image_);

        convert_to_srgb();
        LOG("stacked srgb:");
        analyze_image3(stack_.image_);

        apply_display_gamma();
        LOG("stacked display gamma:");
        analyze_image3(stack_.image_);

        stack_.image_.planes_.scale_to_8bits();
        stack_.image_.save_png(kOutputFile);

        LOG("success.");
        return 0;
    }

    void load_image(
        int tag
    ) {
        std::stringstream ss;
        ss<<kInputFile<<tag<<".CR2";
        LOG("loading image "<<ss.str()<<"...");
        cur_.image_.load_raw(ss.str().c_str());
        if (cur_.image_.is_loaded_ == false) {
            return;
        }
        analyze_image4(cur_.image_);
        interpolate();
        compute_luminance(cur_.image_.planes_, luminance_);
    }

    void analyze_image4(
        Image &image
    ) {
        analyze_plane("r ", image.planes_.r_);
        analyze_plane("g1", image.planes_.g1_);
        analyze_plane("g2", image.planes_.g2_);
        analyze_plane("b ", image.planes_.b_);
    }

    void analyze_image3(
        Image &image
    ) {
        analyze_plane("r", image.planes_.r_);
        analyze_plane("g", image.planes_.g1_);
        analyze_plane("b", image.planes_.b_);
    }

    void analyze_plane(
        const char *tag,
        Plane &plane
    ) {
        int min_s = 65536*256;
        int max_s = -65536*256;
        std::vector<int> histogram;
        histogram.resize(256, 0);

        int sz = plane.width_ * plane.height_;
        for (int i = 0; i < sz; ++i) {
            int s = plane.samples_[i];
            min_s = std::min(min_s, s);
            max_s = std::max(max_s, s);

            if (s < 0) {
                s = 0;
            }
            s /= 256;
            if (s >= 256) {
                s = 255;
            }
            ++histogram[s];
        }
        LOG("image:"<<tag<<": "<<min_s<<" "<<max_s);
        std::stringstream ss;
        for (int i = 0; i < 256; ++i) {
            ss << " " << histogram[i];
        }
        LOG("hist:"<<ss.str());

        int energy = 0;
        for (int i = 0; i < 256; ++i) {
            energy += i * histogram[i];
        }
        LOG("energy: "<<energy);
    }

    void interpolate() {
        LOG("interpolating...");
        cur_.image_.planes_.interpolate_1331();
        int wd = cur_.image_.planes_.r_.width_ - 1;
        int ht = cur_.image_.planes_.r_.height_ - 1;
        cur_.image_.planes_.r_.crop(1, 0, wd+1, ht);
        cur_.image_.planes_.g1_.crop(0, 0, wd, ht);
        cur_.image_.planes_.g2_.crop(1, 1, wd+1, ht+1);
        cur_.image_.planes_.b_.crop(0, 1, wd, ht+1);
    }

    void register_fast() {
        LOG("fast registration...");
        sum_horz_vert_luminance();
        register_horz();
        register_vert();
    }

    void sum_horz_vert_luminance() {
        int wd = luminance_.width_;
        int ht = luminance_.height_;

        cur_.horz_.clear();
        cur_.vert_.clear();
        cur_.horz_.resize(wd, 0);
        cur_.vert_.resize(ht, 0);

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int val = luminance_.get(x, y);
                /*if (val >= noise_)*/ {
                    cur_.horz_[x] = std::max(cur_.horz_[x], val);
                    cur_.vert_[y] = std::max(cur_.vert_[y], val);
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
                int c1 = stack_.horz_[x1];
                int c2 = cur_.horz_[x2];
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
                int c1 = stack_.vert_[y1];
                int c2 = cur_.vert_[y2];
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
        auto& stack_planes = stack_.image_.planes_;
        auto& cur_planes = cur_.image_.planes_;

        int wd = stack_planes.r_.width_;
        int ht = stack_planes.r_.height_;

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
                int r = stack_planes.r_.get(x1, y1);
                int g1 = stack_planes.g1_.get(x1, y1);
                int g2 = stack_planes.g2_.get(x1, y1);
                int b = stack_planes.b_.get(x1, y1);
                r += cur_planes.r_.get(x2, y2);
                g1 += cur_planes.g1_.get(x2, y2);
                g2 += cur_planes.g2_.get(x2, y2);
                b += cur_planes.b_.get(x2, y2);
                stack_planes.r_.set(x1, y1, r);
                stack_planes.g1_.set(x1, y1, g1);
                stack_planes.g2_.set(x1, y1, g2);
                stack_planes.b_.set(x1, y1, b);
            }
        }
    }

    void combine_greens() {
        auto& planes = stack_.image_.planes_;
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

        /*double brightness = 1.3;
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                mat[i][k] *= brightness;
            }
        }*/

        auto& planes = stack_.image_.planes_;
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
        stack_.image_.planes_.apply_user_gamma(pwr);
    }

    void apply_display_gamma() {
        double pwr = 1.0 / stack_.image_.camera_.gamma0_;
        double ts = stack_.image_.camera_.gamma1_;
        int white = 0x10000;
        stack_.image_.planes_.apply_display_gamma(pwr, ts, white);
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


#if 0
const int kStarFloor = 200;
const int kMaxStarRadius = 10;

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
#endif
