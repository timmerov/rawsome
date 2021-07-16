/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
run the experiment of the day.

notes for future experiments:
estimate sharpenss using image gradient.
some info from here:
https://feichtenhofer.github.io/pubs/BSc_thesis_feichtenhofer_sharpness_10-2011.pdf
gradient convolution:
    +1 0 -1
    +2 0 -2
    +1 0 -1
for horizontal Gx.
and transposed for vertical Gy.
G^2 = Gx^2 + Gy^2
threshold = alpha sqrt(average(G))
where alpha = 2.
the logic in the paper here is a bit sketchy.
irregardless...
there's a threshold. ;->

so i took a picture of the moon and venus reflected in a pond with silhouettes of
people on a bridge.
there was lots of camera shake.
the problem is how to undo the shake and recover a clean image.

B = I * C + N

B: observed image.
I: clean image.
C: blur kernel.
N: noise.
*: convolution operation.

the problem is under-determined.
cause many I's and C's will give the correct B.

the user must estimate the size of the blur kernel.

the literature is full of people trying to recover the clean image and the kernel
all in one go.
it's slow.
understandably.
and results not guaranteed.

i propose an iterative solution.

shrink the images and the kernel so the kernel is 1x1.
at this resolution the clean image and the observed image are the same.
this is the initial estimate of the clean image.

double the width (but not the height) of the estimate and the blur kernel.
shrink the observed image to the new dimensions.
*waving hands* figure out the new estimate and the new kernel.

transpose images and kernel.
repeat until a full size image is acheived.

great we still have the same basic problem.
but we have a decent place to start at each step.
ie we have a reasonable estimate image and a reasonable kernel.

for updating the estimated image the question becomes:
how do we create two pixels from one pixel?
natural images have a particular distribution of gradients.
specifically, a heavy tail.
ie most gradients are really small but some are really big.
many methods appear to be kinda cheap to me.
lucy-richardson can be pretty fast assuming the kernel is known.

for updating the estimated kernel... ?
the literature is really weak here.
and worse, the options seem to be conputationally expensive.
there are some restrictions on the kernel that can help.
the blur kernel is likely to be sparse.
elements must be non-negative.
and the non-zeros must be contiguous.

**/

#include "experiment.h"
#include "log.h"
#include "rs_png.h"

#include <vector>

namespace {

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

    Buffer estimate0;
    Buffer estimate1;
    Buffer estimate2;
    Buffer kernel0;
    Buffer kernel1;
    Buffer observed1;
    Buffer estimate1_real;
    Buffer kernel1_real;
    Buffer blurred1;
    Buffer lenrek1;
    Buffer error1;
    Buffer scaling1;
    Buffer temp1;

    void run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;
        LOG("experiment of the day.");

        /**
        these files were downsampled from fullsize by gimp.
        **/
        const char *kEstimate128x128 = "sunset/sunset-gray-128x128.png";
        const char *kEstimate256x128 = "sunset/sunset-gray-256x128.png";
        const char *kKernel4x4       = "sunset/sunset-kernel-4x4.png";
        const char *kKernel8x4       = "sunset/sunset-kernel-8x4.png";
        const char *kShook256x128    = "sunset/sunset-shook-256x128.png";

        /**
        load the half-sized estimate image,
        half-sized kernel,
        and full-sized observed image.
        **/
        load_png(kEstimate128x128, estimate0);
        load_png(kKernel4x4, kernel0);
        load_png(kShook256x128, observed1);

        /**
        double the width of the estimate and kernel by replicating pixels.
        save them.
        **/
        replicate_x(estimate0, estimate1);
        replicate_x(kernel0, kernel1);
        save_png("sunset/x-estimate1.png", estimate1);
        save_png("sunset/x-kernel1.png", kernel1);

        /**
        load the full-sized clean image and full-sized kernel
        for comparison with the expanded estimate and kernel.
        **/
        load_png(kEstimate256x128, estimate1_real);
        load_png(kKernel8x4, kernel1_real);
        save_diff("sunset/x-estimate1-diff.png", estimate1, estimate1_real);
        save_diff("sunset/x-kernel-diff.png", kernel1, kernel1_real);

        /**
        blur the estimate.
        compare to the observed.
        save them.
        */
        convolve(estimate1, kernel1, blurred1);
        save_png("sunset/x-blurred1.png", blurred1);
        save_diff("sunset/x-blurred1-diff.png", blurred1, observed1);

        /**
        lucy-richardson algorithm.
        **/
        divide(observed1, blurred1, error1);
        scale(error1, temp1, 0.5);
        save_png("sunset/x-error1.png", temp1);
        rotate_180(kernel1, lenrek1);
        save_png("sunset/x-lenrek1.png", lenrek1);
        convolve(error1, lenrek1, scaling1);
        scale(scaling1, temp1, 0.5);
        save_png("sunset/x-scaling1.png", temp1);
        multiply(estimate1, scaling1, estimate2);
        save_png("sunset/x-estimate2.png", estimate2);
        save_diff("sunset/x-estimate2-diff.png", estimate2, observed1);
    }

    void replicate_x(
        Buffer &src,
        Buffer &dst
    ) {
        int wd = src.width_;
        int ht = src.height_;
        dst.init(2 * wd, ht);

        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            double s = src.samples_[i];
            int k = 2 * i;
            dst.samples_[k+0] = s;
            dst.samples_[k+1] = s;
        }
    }

    void save_diff(
        const char *filename,
        Buffer &abuf,
        Buffer &bbuf
    ) {
        int wd = abuf.width_;
        int ht = bbuf.height_;
        Buffer diff;
        diff.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            double a = abuf.samples_[i];
            double b = bbuf.samples_[i];
            double d = a - b;
            d = d / 2.0 + 0.5;
            d = std::max(d, 0.0);
            d = std::min(d, 1.0);
            diff.samples_[i] = d;
        }
        save_png(filename, diff);
    }

    void generate_fullsize_shook_image() {
        const char *kSourceFile = "sunset/sunset-gray-1024.png";
        const char *kKernelFile = "sunset/sunset-kernel-32.png";
        const char *kShookFile = "sunset/sunset-shook-1024.png";

        Buffer source;
        Buffer kernel;
        Buffer observed;

        /** load the fullsize source and kernel. **/
        load_png(kSourceFile, source);
        load_png(kKernelFile, kernel);

        /** blur the source with the kernel. **/
        int swd = source.width_;
        int sht = source.height_;
        observed.init(swd, sht);
        convolve(source, kernel, observed);

        /** save it. **/
        save_png(kShookFile, observed);
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
        dest.init(swd, sht);
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
        int wd = numer.width_;
        int ht = numer.height_;
        result.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            double n = numer.samples_[i];
            double d = denom.samples_[i];
            double r = 1.0;
            if (std::abs(d) > 1e-9) {
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
        dst.init(wd, ht);
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
        Buffer &abuf,
        Buffer &bbuf,
        Buffer &result
    ) {
        int wd = abuf.width_;
        int ht = abuf.height_;
        result.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            double a = abuf.samples_[i];
            double b = bbuf.samples_[i];
            double r = a * b;
            r = std::max(r, 0.1);
            r = std::min(r, 1.0);
            result.samples_[i] = r;
        }
    }

    void scale(
        Buffer &src,
        Buffer &dst,
        double factor
    ) {
        int wd = src.width_;
        int ht = src.height_;
        dst.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            double s = src.samples_[i];
            s *= factor;
            s = std::max(s, 0.0);
            s = std::min(s, 1.0);
            dst.samples_[i] = s;
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
