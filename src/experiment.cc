/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
run the experiment of the day.
**/

#include "experiment.h"
#include "log.h"
#include "rs_png.h"

#include <vector>

namespace {

const char *kSourceFile = "moon-venus.png";
const char *kBlurFile = "simple-blur.png";
//const char *kBlurFile = "dual-blur.png";
//const char *kBlurFile = "complex-blur.png";

class Buffer {
public:
    int width_ = 0;
    int height_ = 0;
    std::vector<double> samples_;

    void init(
        int wd,
        int ht
    ) {
        width_ = wd;
        height_ = ht;
        int sz = wd * ht;
        samples_.resize(sz);
    }

    double get(
        int x,
        int y
    ) {
        int idx = y * width_ + x;
        double s = samples_[idx];
        return s;
    }

    void set(
        int x,
        int y,
        double s
    ) {
        int idx = y * width_ + x;
        samples_[idx] = s;
    }
};

class Experiment {
public:
    Experiment() = default;
    ~Experiment() = default;

    Buffer source;
    Buffer kernel;
    Buffer observed;
    Buffer estimate;
    Buffer denom;
    Buffer error;
    Buffer temp;
    Buffer lenrek;
    Buffer scale;

    void run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;
        LOG("experiment of the day.");

        /** load the source image. **/
        load_png(kSourceFile, source);
        save_png("x-real-source.png", source);

        /** load the actual blur kernel. **/
        load_png(kBlurFile, kernel);
        save_png("x-real-kernel.png", kernel);

        /** blur the source with the kernel. **/
        int swd = source.width_;
        int sht = source.height_;
        int ssz = swd * sht;
        observed.init(swd, sht);
        convolve(source, kernel, observed);
        save_png("x-observed.png", observed);

        /** "estimate" the blur kernel. **/
        int kwd = kernel.width_;
        int kht = kernel.height_;
        int ksz = kwd * kht;
        for (int i = 0; i < ksz; ++i) {
            double s = kernel.samples_[i];
            if (s < 0.2) {
                s = 0.0;
            } else {
                s = 0.88;
            }
            kernel.samples_[i] = s;
        }
        save_png("x-est-kernel.png", kernel);

        /**
        the algorithm doesn't like black.
        ensure the source has no black.
        **/
        for (int i = 0; i < ssz; ++i) {
            double s = observed.samples_[i];
            s *= 0.9;
            s += 0.1;
            observed.samples_[i] = s;
        }
        save_png("x-deblacked.png", observed);

        /** start with the source image. **/
        estimate = observed;

        /** initialize temporary space. **/
        denom.init(swd, sht);
        error.init(swd, sht);
        temp.init(swd, sht);
        lenrek.init(kwd, kht);
        scale.init(swd, sht);

        /** convolve the estimated image with the estimated blur kernel. **/
        convolve(estimate, kernel, denom);
        save_png("x-denom.png", denom);

        /** divide the observed image by the denom. **/
        divide(observed, denom, error);
        for (int i = 0; i < ssz; ++i) {
            double s = error.samples_[i];
            s /= 2.0;
            temp.samples_[i] = s;
        }
        save_png("x-error.png", temp);

        /** flip the kernel. **/
        rotate_180(kernel, lenrek);
        save_png("x-flipped-est-kernel.png", lenrek);

        /** convolve the error with the flipped kernel. **/
        convolve(error, lenrek, scale);
        for (int i = 0; i < ssz; ++i) {
            double s = scale.samples_[i];
            s /= 2.0;
            temp.samples_[i] = s;
        }
        save_png("x-scale.png", temp);

        /** multiply the estimate by the scale factor in place. **/
        multiply(estimate, scale);
        save_png("x-estimate-1.png", estimate);
    }

    void load_png(
        const char *filename,
        Buffer &buffer
    ) {
        LOG("loading png: "<<filename);

        /** load the png. **/
        Png png;
        png.read(filename);

        /** copy it to the buffer. **/
        int wd = png.wd_;
        int ht = png.ht_;
        int stride = png.stride_;
        buffer.init(wd, ht);
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int png_idx = y * stride + 3 * x;
                /** assume grayscale **/
                int png_s = png.data_[png_idx];
                int buf_idx = y * wd + x;
                double buf_s = double(png_s) / 255.0;
                buffer.samples_[buf_idx] = buf_s;
            }
        }
    }

    void save_png(
        const char *filename,
        Buffer &buffer
    ) {
        LOG("saving png: "<<filename);

        /** copy the buffer to the png. **/
        Png png;
        int wd = buffer.width_;
        int ht = buffer.height_;
        png.init(wd, ht);
        int stride = png.stride_;
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int buf_idx = y * wd + x;
                double buf_s = buffer.samples_[buf_idx];
                int png_idx = y * stride + 3* x;
                int png_s = buf_s * 255.0;
                /** rgb **/
                png.data_[png_idx+0] = png_s;
                png.data_[png_idx+1] = png_s;
                png.data_[png_idx+2] = png_s;
            }
        }

        /** save the png. p**/
        png.write(filename);
    }

    void convolve(
        Buffer &source,
        Buffer &kernel,
        Buffer &dest
    ) {
        /** assume source and result are the same size. **/
        int swd = source.width_;
        int sht = source.height_;
        /** assume kernel is not. **/
        int kwd = kernel.width_;
        int kht = kernel.height_;
        for (int dy = 0; dy < sht; ++dy) {
            for (int dx = 0; dx < swd; ++dx) {
                double sums = 0.0;
                double sumk = 0.0;
                for (int ky = 0; ky < kht; ++ky) {
                    int sy = dy + ky - kht/2;
                    if (sy < 0) {
                        continue;
                    }
                    if (sy >= sht) {
                        break;
                    }
                    for (int kx = 0; kx < kwd; ++kx) {
                        int sx = dx + kx - kwd/2;
                        if (sx < 0) {
                            continue;
                        }
                        if (sx >= swd) {
                            break;
                        }
                        double k = kernel.get(kx, ky);
                        double ss = source.get(sx, sy);
                        sums += k * ss;
                        sumk += k;
                    }
                }
                double ds = 0.0;
                if (sumk > 0.0) {
                    ds = sums / sumk;
                }
                dest.set(dx, dy, ds);
            }
        }
    }

    void divide(
        Buffer &numer,
        Buffer &denom,
        Buffer &result
    ) {
        int sz = numer.width_ * numer.height_;
        for (int i = 0; i < sz; ++i) {
            double n = numer.samples_[i];
            double d = denom.samples_[i];
            double r = 1.0;
            if (d != 0.0) {
                r = n / d;
            }
            r = std::max(r, 0.0);
            r = std::min(r, 2.0);
            result.samples_[i] = r;
        }
    }

    void rotate_180(
        Buffer &src,
        Buffer &dst
    ) {
        int wd = src.width_;
        int ht = src.height_;
        for (int dy = 0; dy < ht; ++dy) {
            int sy = ht - 1 - dy;
            for (int dx = 0; dx < wd; ++dx) {
                int sx = wd - 1 - dx;
                double s = src.get(sx, sy);
                dst.set(dx, dy, s);
            }
        }
    }

    void multiply(
        Buffer &buffer,
        Buffer &scale
    ) {
        int sz = buffer.width_ * buffer.height_;
        for (int i = 0; i < sz; ++i) {
            double a = buffer.samples_[i];
            double s = scale.samples_[i];
            double b = a * s;
            b = std::max(b, 0.1);
            b = std::min(b, 1.0);
            buffer.samples_[i] = b;
        }
    }
};

}

int experiment_of_the_day(
    int argc,
    clo_argv_t argv
) {
    Experiment x;
    x.run(argc, argv);

    return 0;
}
