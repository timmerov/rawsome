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

/*const int kFirstTag = 1179;
const int kLastTag = 1180;
const int kLastTag = 1198;
const char *kInputFile = "/home/timmer/Pictures/2021-02-27/pleiades800/IMG_";
const char *kOutputFile = "/home/timmer/Pictures/2021-02-27/pleiades800/stack.png";
const int kMaxRegisterOffset = 5;*/

const int kFirstTag = 1199;
const int kLastTag = 1209;
const char *kInputFile = "/home/timmer/Pictures/2021-02-27/black800/IMG_";
const char *kOutputFile = "/home/timmer/Pictures/2021-02-27/black800/stack.png";
const int kMaxRegisterOffset = 0;

class StackImages {
public:
    StackImages() = default;
    ~StackImages() = default;

    Image image_;
    Plane luminance_;
    int noise_ = 0;
    int dx_ = 0;
    int dy_ = 0;
    int count_ = 0;
    Plane stack_lum_;
    Image stack_image_;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("stack images - work in progress.");

        /** load the first image into luminance2. **/
        load_image(kFirstTag);
        image_.camera_.print();

        /** copy it to the stacked image. **/
        stack_lum_ = std::move(luminance_);
        stack_image_ = std::move(image_);

        /** perfectly aligned. **/
        dx_ = 0;
        dy_ = 0;

        /** one image loaded so far. **/
        count_ = 1;

        /** pick an arbitrary noise value. **/
        noise_ = 500;

        /** stack the next images. **/
        for (int i = kFirstTag+1; i <= kLastTag; ++i) {
            load_image(i);
            if (image_.is_loaded_ == false) {
                return 1;
            }
            register_with_stack();
            stack_with_offset();
        }

        white_balance();
        combine_greens();
        convert_to_srgb();
        apply_gamma();
        stack_image_.planes_.scale_to_8bits();
        stack_image_.save_png(kOutputFile);

        LOG("success.");
        return 0;
    }

    void load_image(
        int tag
    ) {
        std::stringstream ss;
        ss<<kInputFile<<tag<<".CR2";
        image_.load_raw(ss.str().c_str());
        if (image_.is_loaded_) {
            compute_luminance(image_.planes_, luminance_);
        }
    }

    void register_with_stack() {
        int mindx = - kMaxRegisterOffset;
        int maxdx = + kMaxRegisterOffset;
        int mindy = - kMaxRegisterOffset;
        int maxdy = + kMaxRegisterOffset;

        std::int64_t min_res = 0;
        dx_ = 0;
        dy_ = 0;

        for (int dy = mindy; dy <= maxdy; ++dy) {
            for (int dx = mindx; dx <= maxdx; ++dx) {
                std::int64_t res = residual(dx, dy);
                if (min_res == 0 || min_res > res) {
                    min_res = res;
                    dx_ = dx;
                    dy_ = dy;
                    LOG("min_res="<<min_res<<" dx,dy="<<dx<<","<<dy);
                }
            }
        }
    }

    std::int64_t residual(
        int dx,
        int dy
    ) {
        std::int64_t res = 0;

        int wd = luminance_.width_;
        int ht = luminance_.height_;
        for (int y1 = 0; y1 < ht; ++y1) {
            int y2 = y1 + dy;
            if (y2 < 0 || y2 >= ht) {
                continue;
            }
            for (int x1 = 0; x1 < wd; ++x1) {
                int x2 = x1 + dx;
                if (x2 < 0 || x2 >= wd) {
                    continue;
                }
                int cl = luminance_.get(x1, y1);
                int cs = stack_lum_.get(x2, y2);

                cs = (cs + count_/2) / count_;
                if (cs <= noise_) {
                    cs = 0;
                }
                if (cl <= noise_) {
                    cl = 0;
                }

                int df = std::abs(cs - cl);
                res += df;
            }
        }

        return res;
    }

    void stack_with_offset() {
        Plane dst_lum;
        Planes dst_rggb;

        int wd = stack_lum_.width_;
        int ht = stack_lum_.height_;
        dst_lum.init(wd, ht);
        dst_rggb.init(wd, ht);

        for (int y1 = 0; y1 < ht; ++y1) {
            int y2 = y1 + dy_;
            if (y2 < 0 || y2 >= ht) {
                continue;
            }
            for (int x1 = 0; x1 < wd; ++x1) {
                int x2 = x1 + dx_;
                if (x2 < 0 || x2 >= wd) {
                    continue;
                }
                int lum = luminance_.get(x1, y1);
                lum += stack_lum_.get(x2, y2);
                dst_lum.set(x1, y1, lum);

                int r = image_.planes_.r_.get(x1, y1);
                int g1 = image_.planes_.g1_.get(x1, y1);
                int g2 = image_.planes_.g2_.get(x1, y1);
                int b = image_.planes_.b_.get(x1, y1);
                r += stack_image_.planes_.r_.get(x2, y2);
                g1 += stack_image_.planes_.g1_.get(x2, y2);
                g2 += stack_image_.planes_.g2_.get(x2, y2);
                b += stack_image_.planes_.b_.get(x2, y2);
                dst_rggb.r_.set(x1, y1, r);
                dst_rggb.g1_.set(x1, y1, g1);
                dst_rggb.g2_.set(x1, y1, g2);
                dst_rggb.b_.set(x1, y1, b);
            }
        }

        stack_lum_ = std::move(dst_lum);
        stack_image_.planes_ = std::move(dst_rggb);
        ++count_;
    }

    void white_balance() {
        /**
        normally, we would get the rggb camera multipliers from the raw_image.
        note order permutation.
        **/
        RggbDouble cam_mul;
        cam_mul.r_ = stack_image_.camera_.wb_r_;
        cam_mul.g1_ = 1.0;
        cam_mul.g2_ = 1.0;
        cam_mul.b_ = stack_image_.camera_.wb_b_;
        LOG("white balance: R="<<cam_mul.r_<<" B="<<cam_mul.b_);

        /** white balance. **/
        stack_image_.planes_.multiply_sat(cam_mul);
    }

    void combine_greens() {
        int sz = stack_image_.planes_.g1_.width_ * stack_image_.planes_.g1_.height_;
        for (int i = 0; i < sz; ++i) {
            int g1 = stack_image_.planes_.g1_.samples_[i];
            int g2 = stack_image_.planes_.g2_.samples_[i];
            stack_image_.planes_.g1_.samples_[i] = (g1 + g2 + 1)/2;
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

        int sz = stack_image_.planes_.r_.width_ * stack_image_.planes_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = stack_image_.planes_.r_.samples_[i];
            int in_g = stack_image_.planes_.g1_.samples_[i];
            int in_b = stack_image_.planes_.b_.samples_[i];

            /** transform by matrix multiplication. **/
            double out_r = mat[0][0]*in_r + mat[0][1]*in_g + mat[0][2]*in_b;
            double out_g = mat[1][0]*in_r + mat[1][1]*in_g + mat[1][2]*in_b;
            double out_b = mat[2][0]*in_r + mat[2][1]*in_g + mat[2][2]*in_b;

            /** overwrite old values. **/
            stack_image_.planes_.r_.samples_[i] = out_r;
            stack_image_.planes_.g1_.samples_[i] = out_g;
            stack_image_.planes_.b_.samples_[i] = out_b;
        }
    }

    void apply_gamma() {
        double pwr = 1.0 / stack_image_.camera_.gamma0_;
        double ts = stack_image_.camera_.gamma1_;
        int white = 0x10000;
        stack_image_.planes_.apply_gamma(pwr, ts, white);
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
