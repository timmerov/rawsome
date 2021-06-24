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
    Buffer blurred;

    void run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;
        LOG("experiment of the day.");

        load_png(kSourceFile, source);
        load_png(kBlurFile, kernel);

        int wd = source.width_;
        int ht = source.height_;
        blurred.init(wd, ht);
        convolve(source, kernel, blurred);

        save_png("experiment-blurred.png", blurred);
        save_png("experiment-kernel.png", kernel);
        save_png("experiment-source.png", source);
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
                    int sy = dy - ky + kht/2;
                    if (sy < 0) {
                        continue;
                    }
                    if (sy >= sht) {
                        break;
                    }
                    for (int kx = 0; kx < kwd; ++kx) {
                        int sx = dx - kx + kwd/2;
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
