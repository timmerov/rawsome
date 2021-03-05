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

//const int kFirstTag = 1227;
//const int kLastTag = 1228;
//const int kLastTag = 1953;
const int kFirstTag = 1250;
const int kLastTag = 1253;
const char *kInputFile = "/home/timmer/Pictures/2021-03-03/mars-pleiades/IMG_";
const char *kOutputFile = "/home/timmer/Pictures/2021-02-27/pleiades800/stack.png";
const char *kBlackFile = "/home/timmer/Pictures/2021-03-03/black/black.rsm";

const int kNoiseFloor = 500;
const int kMaxRegisterOffset = 700;

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

    Luminage prev_;
    Luminage cur_;
    Luminage stack_;

    int noise_ = 0;
    int count_ = 0;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("stack images - work in progress.");

        /** load the black file. **/
        (void) kBlackFile;

        /** load the first image. **/
        load_image(kFirstTag);
        cur_.image_.camera_.print();

        /** copy it to the stacked image. **/
        stack_ = cur_;

        /** choose an arbitrary noise value. **/
        noise_ = kNoiseFloor;

        /** one image loaded so far. **/
        count_ = 1;

        /** do its sums **/
        sum_horz_vert_luminance();

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
            //register_with_stack();
            //stack_with_offset();

            /** save current as previous. **/
            prev_ = std::move(cur_);
        }

        white_balance();
        combine_greens();
        convert_to_srgb();
        apply_gamma();
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
        cur_.image_.load_raw(ss.str().c_str());
        if (cur_.image_.is_loaded_) {
            compute_luminance(cur_.image_.planes_, cur_.luminance_);
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
    }

    void white_balance() {
        RggbDouble cam_mul;
        cam_mul.r_ = stack_.image_.camera_.wb_r_;
        cam_mul.g1_ = 1.0;
        cam_mul.g2_ = 1.0;
        cam_mul.b_ = stack_.image_.camera_.wb_b_;
        LOG("white balance: R="<<cam_mul.r_<<" B="<<cam_mul.b_);

        /** white balance. **/
        stack_.image_.planes_.multiply_sat(cam_mul);
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

        double brightness = 1.3;
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                mat[i][k] *= brightness;
            }
        }

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

    void apply_gamma() {
        double pwr = 1.0 / stack_.image_.camera_.gamma0_;
        double ts = stack_.image_.camera_.gamma1_;
        int white = 0x10000;
        stack_.image_.planes_.apply_gamma(pwr, ts, white);
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
