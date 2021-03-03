/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
register two source images.
**/

#include "canon.h"
#include "image.h"
#include "log.h"
#include "register.h"

#include <cmath>
#include <sstream>

namespace {

const char *kInputFile1 = "/home/timmer/Pictures/2021-02-27/pleiades800/IMG_1179.CR2";
const char *kInputFile2 = "/home/timmer/Pictures/2021-02-27/pleiades800/IMG_1180.CR2";

const char *kOutputFileWip1 = "/home/timmer/Pictures/2021-02-27/pleiades800/wip1.png";
const char *kOutputFileWip2 = "/home/timmer/Pictures/2021-02-27/pleiades800/wip2.png";
const char *kOutputFileDiff = "/home/timmer/Pictures/2021-02-27/pleiades800/diff.png";
const char *kOutputFileStack = "/home/timmer/Pictures/2021-02-27/pleiades800/stack.png";

const int kMaxRegisterOffset = 5;

class RegisterTwoImages {
public:
    RegisterTwoImages() = default;
    ~RegisterTwoImages() = default;

    Image image1_;
    Image image2_;
    Plane luminance1_;
    Plane luminance2_;
    Plane diff_;
    int noise_ = 0;
    int dx_ = 0;
    int dy_ = 0;
    Plane stack_;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("register two images - work in progress.");

        image1_.load_raw(kInputFile1);
        if (image1_.is_loaded_ == false) {
            return 1;
        }
        image2_.load_raw(kInputFile2);
        if (image2_.is_loaded_ == false) {
            return 1;
        }

        compute_luminance(image1_.planes_, luminance1_);
        compute_luminance(image2_.planes_, luminance2_);
        difference();
        residual();
        stack_with_offset();

        find_brightest(luminance1_);
        find_brightest(luminance2_);

        /*histogram(luminance1_);
        histogram(luminance2_);*/

        threshold(luminance1_);
        threshold(luminance2_);

        //luminance1_.multiply(2.3);
        //luminance2_.multiply(2.3);

        save_luminance(image1_, luminance1_, kOutputFileWip1);
        save_luminance(image2_, luminance2_, kOutputFileWip2);
        save_plane(diff_, kOutputFileDiff);
        save_plane(stack_, kOutputFileStack);

        LOG("success.");
        return 0;
    }

    void difference() {
        int wd = luminance1_.width_;
        int ht = luminance1_.height_;
        diff_.init(wd, ht);

        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            int c1 = luminance1_.samples_[i];
            int c2 = luminance2_.samples_[i];
            int df = c1 - c2 + 32767;
            diff_.samples_[i] = df;
        }
    }

    void residual() {
        noise_ = std::max(image1_.noise_, image2_.noise_);
        noise_ *= 3;

        std::int64_t min_res = 0;
        min_res *= luminance1_.width_ * luminance1_.height_;
        min_res *= 100;

        for (int dy = - kMaxRegisterOffset; dy <= kMaxRegisterOffset; ++dy) {
            for (int dx = - kMaxRegisterOffset; dx <= kMaxRegisterOffset; ++dx) {
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

        int wd = luminance1_.width_ - kMaxRegisterOffset - 1;
        int ht = luminance1_.height_- kMaxRegisterOffset - 1;
        for (int y = kMaxRegisterOffset; y < ht; ++y) {
            for (int x = kMaxRegisterOffset; x < wd; ++x) {
                int c1 = luminance1_.get(x, y);
                int c2 = luminance2_.get(x + dx, y + dy);
                if (c1 <= noise_) {
                    c1 = 0;
                }
                if (c2 <= noise_) {
                    c2 = 0;
                }
                int df = std::abs(c1 - c2);
                df *= df;
                res += df;
            }
        }

        return std::sqrt(res);
    }

    void stack_with_offset() {
        int wd = luminance1_.width_;
        int ht = luminance1_.height_;
        stack_.init(wd, ht);

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
                int c1 = luminance1_.get(x1, y1);
                int c2 = luminance2_.get(x2, y2);
                int c = c1 + c2;
                stack_.set(x1, y1, c);
            }
        }
    }

    void find_brightest(
        Plane &luminance
    ) {
        int maxc = 0;
        int sz = luminance.width_ * luminance.height_;
        for (int i = 0; i < sz; ++i) {
            int c = luminance.samples_[i];
            maxc = std::max(maxc, c);
        }
        LOG("brightest: "<<maxc);
    }

    void histogram(
        Plane &luminance
    ) {
        std::vector<int> histogram;
        histogram.resize(6000, 0);

        int sz = luminance.width_ * luminance.height_;
        for (int i = 0; i < sz; ++i) {
            int c = luminance.samples_[i];
            if (c < 0) {
                c = 0;
            } else if (c >= 6000) {
                c = 6000-1;
            }
            ++histogram[c];
        }

        std::stringstream ss;
        ss<<"histogram: [";
        for (int i = 200; i < 6000; ++i) {
            ss<<" "<<histogram[i];
        }
        ss<<" ]";
        LOG(ss.str());
    }

    void threshold(
        Plane &luminance
    ) {
        int count = 0;
        int sz = luminance.width_ * luminance.height_;
        for (int i = 0; i < sz; ++i) {
            int c = luminance.samples_[i];
            if (c < 180) {
                c = 0;
            } else {
                c = 65535;
                ++count;
            }
            luminance.samples_[i] = c;
        }
        LOG("threshold: "<<count);
    }

    void save_luminance(
        Image &image,
        Plane &luminance,
        const char *filename
    ) {
        /*double pwr = 1.0 / image.camera_.gamma0_;
        double ts = image.camera_.gamma1_;
        int white = 0x10000;
        luminance.apply_gamma(pwr, ts, white);*/

        luminance.scale_to_8bits();

        /** overwrite rgb with luminance. **/
        image.planes_.r_ = luminance;
        image.planes_.g1_ = luminance;
        image.planes_.b_ = luminance;

        image.save_png(filename);
    }

    void save_plane(
        Plane &plane,
        const char *filename
    ) {
        plane.scale_to_8bits();

        /** overwrite rgb with difference. **/
        Image image;
        image.planes_.init(plane.width_, plane.height_);
        image.planes_.r_ = plane;
        image.planes_.g1_ = plane;
        image.planes_.b_ = plane;

        image.save_png(filename);
    }
};
}

int register_two_images(
    int argc,
    clo_argv_t argv
) {

    RegisterTwoImages reg;
    int exit_code = reg.run(argc, argv);

    return exit_code;
}
