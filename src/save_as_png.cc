/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
process pictures of the moon.

some things we've learned about libraw and/or this camera.
open_file(filename)
unpack()
raw2image();
imgdata.sizes.width
imgdata.sizes.height
idx = x + y * width
imgdata.image[idx][c]
where c is the color component index is:
0: red
1: green 1
2: blue
3: green 2
note: ONE of those 4 pixels will have a value.
the other THREE will be zero.
pretty sure the bayer pattern is:
     x     x+1
y  : R/0   G1/1
y+1: G2/3  B/2
where the /# is the above color component index as per above.
one assumes one is expected to interpolate rgb values for everywhere else.
the rgb component values range from about 2k to 14k.
it's not clear how to get these numbers from the imgdata.
each r g1 b g2 component has its own slightly different black value.
r and b seem to correlate very tightly.
but g1 g2 seem to be out to lunch.
they seem to need need to be scaled (after subracting black) by about 1.9.

this camera seems to have some unused pixels:
the top 36 rows.
and the left 76 columns.
these numbers are hard coded in dcraw.
weird.

the height of this camera is 4051.
which is literally and figuratively odd.
should probably ignore the last row.

we should probably crop to:
left   = 76 unused + 3 interpolation + 1 round = 80
top    = 36 unused + 3 interpolation + 1 round = 40
right  = 6096 width - 3 interpolation - 1 round = 6092
bottom = 4051 height - 3 interpolation - 2 round = 4046
width  = 6092 - 80 = 6012
height = 4046 - 40 = 4006

dcraw gets the black value from the unused pixels.
that's a LOT of black pixels.
be careful of overflow.

we need to subtract black first.

second we need to do white balance and pixel scaling.
we're given cam_mul.
which gives the relative magnitudes of rgbg.
relative is the key word.
find the smallest: cam_mul_min.
range = 16384-1 - black
the scaling values are: cam_mul[i] / cam_mul_min * 65536-1 / range
multiply every pixel by the scaling values.
and now the pixel values span 0 .. 65535.
yay.

third we need to gamma correct.
we're given imgdata.params.gamm[4].
steal the gamma_curve function from dcraw.c.
dcraw does gamma correction when converting to rgb.
i don't know the values of the gamma curve.

note: if we don't care about component...
like say for to compute black.
then we can use imgdata.rawdata.raw_image.
it's the single component version of imgdata.image.

re interpolation...
currently we interpolate rgbg values for every pixel.
using a 1331 filter.
there are probably other methods.

note: there seems to be a black band top and left of some images.
perhaps there's padding in the image data.

note: there seem to be more r&b pixels than g pixels.
strange. why?

note: the base number seems to be consistent across shots.
the factor number is not.
some shots the g1 and g2 histograms are somewhat different.
it almost looks like the greens need to be offset.
maybe by the difference in the averages excluding the near blacks?

note: where are the pixels with base values?
are they near the borders?
like maybe not illuminated?
are they maybe literally camera noise?
**/

#include "canon.h"
#include "dump.h"
#include "image.h"
#include "log.h"
#include "options.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <sstream>
#include <thread>

namespace {

const int kFullySaturated = (1<<30)-1;

void convolve_plane_thread(
    Plane *src,
    Plane *kernel,
    Plane *dst,
    int ystart,
    int dy
) {
    int swd = src->width_;
    int sht = src->height_;
    int kwd = kernel->width_;
    int kht = kernel->height_;

    for (int y1 = ystart; y1 < sht; y1 += dy) {
        if (y1 % 100 == 0) {
            LOG("y1="<<y1);
        }
        for (int x1 = 0; x1 < swd; ++x1) {
            double s = 0;
            double sumk = 0;
            for (int yk = 0; yk < kht; ++yk) {
                int y0 = y1 + yk - kht / 2;
                if (y0 < 0) {
                    continue;
                }
                if (y0 >= sht) {
                    break;
                }
                for (int xk = 0; xk < kwd; ++xk) {
                    int x0 = x1 + xk - kwd / 2;
                    if (x0 < 0) {
                        continue;
                    }
                    if (x0 >= swd) {
                        break;
                    }
                    double k = kernel->get(xk, yk);
                    if (k != 0) {
                        s += k * src->get(x0, y0);
                        sumk += k;
                    }
                }
            }
            if (sumk > 0) {
                s /= sumk;
            } else {
                s = 32768;
            }
            dst->set(x1, y1, s);
        }
    }
}

#if 0
void convolve_planes_thread(
    Planes *src,
    Plane *kernel,
    Planes *dst,
    int ystart,
    int dy
) {
    int swd = src->r_.width_;
    int sht = src->r_.height_;
    int kwd = kernel->width_;
    int kht = kernel->height_;

    for (int y1 = ystart; y1 < sht; y1 += dy) {
        if (y1 % 100 == 0) {
            LOG("y1="<<y1);
        }
        for (int x1 = 0; x1 < swd; ++x1) {
            double r = 0;
            double g1 = 0;
            double g2 = 0;
            double b = 0;
            double sumk = 0;
            for (int yk = 0; yk < kht; ++yk) {
                int y0 = y1 + yk - kht / 2;
                if (y0 < 0) {
                    continue;
                }
                if (y0 >= sht) {
                    break;
                }
                for (int xk = 0; xk < kwd; ++xk) {
                    int x0 = x1 + xk - kwd / 2;
                    if (x0 < 0) {
                        continue;
                    }
                    if (x0 >= swd) {
                        break;
                    }
                    double k = kernel->get(xk, yk);
                    if (k != 0) {
                        r += k * src->r_.get(x0, y0);
                        g1 += k * src->g1_.get(x0, y0);
                        g2 += k * src->g2_.get(x0, y0);
                        b += k * src->b_.get(x0, y0);
                        sumk += k;
                    }
                }
            }
            if (sumk > 0) {
                r /= sumk;
                g1 /= sumk;
                g2 /= sumk;
                b /= sumk;
            }
            dst->r_.set(x1, y1, r);
            dst->g1_.set(x1, y1, g1);
            dst->g2_.set(x1, y1, g2);
            dst->b_.set(x1, y1, b);
        }
    }
}
#endif

void convolve_mt(
    Plane &src,
    Plane &kernel,
    Plane &dst
) {
    std::thread th0(convolve_plane_thread, &src, &kernel, &dst, 0, 8);
    std::thread th1(convolve_plane_thread, &src, &kernel, &dst, 1, 8);
    std::thread th2(convolve_plane_thread, &src, &kernel, &dst, 2, 8);
    std::thread th3(convolve_plane_thread, &src, &kernel, &dst, 3, 8);
    std::thread th4(convolve_plane_thread, &src, &kernel, &dst, 4, 8);
    std::thread th5(convolve_plane_thread, &src, &kernel, &dst, 5, 8);
    std::thread th6(convolve_plane_thread, &src, &kernel, &dst, 6, 8);
    std::thread th7(convolve_plane_thread, &src, &kernel, &dst, 7, 8);
    th0.join();
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
}

#if 0
void convolve_mt(
    Planes &src,
    Plane &kernel,
    Planes &dst
) {
    std::thread th0(convolve_planes_thread, &src, &kernel, &dst, 0, 8);
    std::thread th1(convolve_planes_thread, &src, &kernel, &dst, 1, 8);
    std::thread th2(convolve_planes_thread, &src, &kernel, &dst, 2, 8);
    std::thread th3(convolve_planes_thread, &src, &kernel, &dst, 3, 8);
    std::thread th4(convolve_planes_thread, &src, &kernel, &dst, 4, 8);
    std::thread th5(convolve_planes_thread, &src, &kernel, &dst, 5, 8);
    std::thread th6(convolve_planes_thread, &src, &kernel, &dst, 6, 8);
    std::thread th7(convolve_planes_thread, &src, &kernel, &dst, 7, 8);
    th0.join();
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
}
#endif

class SaveAsPng {
public:
    SaveAsPng() = default;
    ~SaveAsPng() = default;

    Options opt_;
    Image image_;
    int saturation_;
    std::vector<int> histogram_;
    Plane luminance_;
    double lum_rx_ = 0.0;
    double lum_gx_ = 0.0;
    double lum_bx_ = 0.0;
    double auto_brightness_ = 0.0;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        bool result = opt_.parse(argc, argv);
        if (result == false) {
            Options::print_usage();
            return 1;
        }

        image_.load_raw(opt_.in_filename_.c_str());
        if (image_.is_loaded_ == false) {
            return 1;
        }

        image_.camera_.print();
        show_special_pixel();
        fix_bad_pixels();
        determine_saturation();
        desaturate_pixels();
        show_special_pixel();
        scale_image();
        show_special_pixel();
        deblur_image();
        show_special_pixel();
        adjust_dynamic_range();
        show_special_pixel();
        interpolate();
        combine_greens();
        show_special_pixel();
        determine_auto_brightness();
        convert_to_srgb();
        show_special_pixel();
        enhance_colors();
        show_special_pixel();
        apply_gamma();
        show_special_pixel();
        scale_to_8bits();
        image_.save_png(opt_.out_filename_.c_str());
        return 0;
    }

    void show_special_pixel() {
        /*int x = 2346;
        int y = 1186;
        int r = planes_.r_.get(x, y);
        int g1 = planes_.g1_.get(x, y);
        int g2 = planes_.g2_.get(x, y);
        int b = planes_.b_.get(x, y);
        LOG("special pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);*/
    }

    void fix_bad_pixels() {
        /** hard coded list of bad pixels. **/
        fix_bad_pixel(4011, 2495, 0);
    }

    void fix_bad_pixel(
        int x,
        int y,
        int which
    ) {
        const char *what = "";
        Plane *plane = nullptr;
        switch (which) {
            case 0:
                plane = &image_.planes_.r_;
                what = "red";
                break;
            case 1:
                plane = &image_.planes_.g1_;
                what = "green1";
                break;
            case 2:
                plane = &image_.planes_.g2_;
                what = "green2";
                break;
            case 3:
                plane = &image_.planes_.b_;
                what = "blue";
                break;
        }

        LOG("fixing bad "<<what<<" pixel at "<<x<<","<<y);

        /**
        we found the bad pixel by inspection.
        x,y are in final coordinates.
        ie after interpolation. hence the /2.
        also after cropping border pixels. hence the +2.
        **/
        x = x / 2 + 2;
        y = y / 2 + 2;

        /** overwrite the bad pixel with the average of its neighbors. **/
        int c0 = plane->get(x, y-1);
        int c1 = plane->get(x-1, y);
        int c2 = plane->get(x+1, y);
        int c3 = plane->get(x, y+1);
        int c = (c0 + c1 + c2 + c3 + 2) / 4;
        plane->set(x, y, c);
    }

    void determine_saturation() {
        LOG("determining saturation...");

        if (opt_.saturation_ > 0) {
            LOG("using user saturation...");
            saturation_ = opt_.saturation_;
        } else {
            /**
            the camera sensor saturates at less than the maximum value.
            **/
            int satr = determine_saturation(image_.planes_.r_);
            int satg1 = determine_saturation(image_.planes_.g1_);
            int satg2 = determine_saturation(image_.planes_.g2_);
            int satb = determine_saturation(image_.planes_.b_);
            saturation_ = std::max(satr, satg1);
            saturation_ = std::max(saturation_, satg2);
            saturation_ = std::max(saturation_, satb);
        }

        LOG("saturation is: "<<saturation_);
    }

    int determine_saturation(
        Plane &plane
    ) {
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
        const int kSatMax = 16383;
        const int kSatExpected = 13583 - 2046; /** we need to subtract black. **/
        const int kSatDiff = kSatMax - kSatExpected;
        const int kSatSize = 2*kSatDiff;
        const int kSatThreshold = kSatMax - kSatSize;

        histogram_.clear();
        histogram_.resize(kSatSize, 0);
        int sz = plane.width_ * plane.height_;
        for (int i = 0; i < sz; ++i) {
            int c = plane.samples_[i];
            c -= kSatThreshold;
            if (c >= 0 && c < kSatSize) {
                ++histogram_[c ];
            }
        }

        /** skip trailing zeros. **/
        int saturation = kSatSize-1;
        for (; saturation >= 0; --saturation) {
            if (histogram_[saturation] > 0) {
                break;
            }
        }

        /** for very black pictures, use the expected value. **/
        if (saturation < 12) {
            return kSatExpected;
        }

        /**
        look at the previous 12 values.

        if they kinda tail away then no pixels are saturated.
        use max of saturation+1 and expectation.

        if they're mostly small values but there's a big jump.
        then set the saturation before the big jump.
        **/
        int idx = saturation - 12;
        int c0 = histogram_[idx];
        c0 = std::min(c0, 100);
        for (++idx; idx <= saturation; ++idx) {
            int c1 = histogram_[idx];
            if (c1 > 3*c0) {
                /**
                we found a big jump.
                idx is the saturation value.
                might as well back it off a bit.
                **/
                saturation = idx - 2;
                saturation += kSatThreshold;
                return saturation;
            }
        }

        /** no big jump. **/
        saturation += 2;
        saturation += kSatThreshold;
        saturation = std::max(saturation, kSatExpected);
        return saturation;
    }

    void desaturate_pixels() {
        LOG("setting saturated pixels to white...");
        /** if any component is saturated then saturate all components. **/
        int count = 0;
        int wd = image_.planes_.r_.width_;
        int ht = image_.planes_.r_.height_;
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            int r = image_.planes_.r_.samples_[i];
            int g1 = image_.planes_.g1_.samples_[i];
            int g2 = image_.planes_.g2_.samples_[i];
            int b = image_.planes_.b_.samples_[i];
            if (r >= kSaturated
            ||  g1 >= kSaturated
            ||  g2 >= kSaturated
            ||  b >= kSaturated) {
                image_.planes_.r_.samples_[i] = kSaturated;
                image_.planes_.g1_.samples_[i] = kSaturated;
                image_.planes_.g2_.samples_[i] = kSaturated;
                image_.planes_.b_.samples_[i] = kSaturated;
                ++count;
            }
        }
        double pct = 100.0 * double(count) / double(sz);
        LOG("saturated pixels: "<<count<<" "<<pct<<"%");
    }

    void scale_image() {
        LOG("scaling image...");
        /**
        combine three things:
        1. subtract black
        2. white balance
        3. scale to full 32 bits.
        **/

        /**
        normally, we would get the rggb camera multipliers from the raw_image.
        note order permutation.
        **/
        RggbDouble cam_mul;
        if (opt_.wb_r_ > 0.0 && opt_.wb_b_ >= 0.0) {
            LOG("using user white balance...");
            cam_mul.r_ = opt_.wb_r_;
            cam_mul.g1_ = 1.0;
            cam_mul.g2_ = 1.0;
            cam_mul.b_ = opt_.wb_b_;
        } else {
            cam_mul.r_ = image_.camera_.wb_r_;
            cam_mul.g1_ = 1.0;
            cam_mul.g2_ = 1.0;
            cam_mul.b_ = image_.camera_.wb_b_;
        }

        /** find the smallest multiplier. **/
        double f = std::min(cam_mul.r_, cam_mul.g1_);
        f = std::min(f, cam_mul.g2_);
        f = std::min(f, cam_mul.b_);

        /** scale so smallest is 1.0 **/
        cam_mul.r_ /= f;
        cam_mul.g1_ /= f;
        cam_mul.g2_ /= f;
        cam_mul.b_ /= f;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);
        LOG("white balance R="<<cam_mul.r_<<" B="<<cam_mul.b_);

        /**
        adjust so fully saturated scales to 1.0 when cam_mul is 1.0.
        (16383 - black) * cam_mul_new = 1.0 = cam_mul_org
        cam_mul_new = cam_mul_org / (16383 - black)
        **/
        cam_mul.r_ /= saturation_;
        cam_mul.g1_ /= saturation_;
        cam_mul.g2_ /= saturation_;
        cam_mul.b_ /= saturation_;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** adjust to span full 32 bit range. **/
        cam_mul.r_ *= 65535.0;
        cam_mul.g1_ *= 65535.0;
        cam_mul.g2_ *= 65535.0;
        cam_mul.b_ *= 65535.0;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** white balance. **/
        image_.planes_.multiply_sat(cam_mul);
    }

    void deblur_image() {
        LOG("deblurring image...");
        int l = opt_.deblur_left_;
        int t = opt_.deblur_top_;
        int r = opt_.deblur_right_;
        int b = opt_.deblur_bottom_;
        if (l <= 0 || t <= 0 || r <= 0 || b <= 0 || l >= r || t >= b) {
            LOG("deblurring is not enabled.");
            return;
        }
        LOG("deblurring patch: "<<l<<","<<t<<","<<r<<","<<b);

        /**
        we are going to be working in grayscale.
        so let's just start there.
        **/
        Plane cropped;
        compute_luminance(image_.planes_, cropped);

        /**
        we found the patch coordinages by inspection.
        t,l,b,r are in final coordinates.
        ie after interpolation. hence the /2.
        also after cropping border pixels. hence the +2.
        **/
        t = t / 2 + 2;
        l = l / 2 + 2;
        b = b / 2 + 2;
        r = r / 2 + 2;
        int kwd = r - l;
        int kht = b - t;
        /**
        we want the patch to have a center.
        so we ensure the width and height are odd.
        **/
        kwd |= 1;
        kht |= 1;
        r = l + kwd;
        b = t + kht;

        /**
        estimate the point spread function (the blur kernel)
        from the luninance of the patch.
        threshold it so most of the values are zero.
        **/
        Plane kernel = cropped;
        kernel.crop(l, t, r, b);

        /** this is approximately the number of kernel pixels we want lit up. **/
        int litup = 1 * (kwd + kht);

        /** find a threshold that gives us that many pixels. **/
        Plane sorted = kernel;
        std::sort(sorted.samples_.begin(), sorted.samples_.end());
        int ksz = kwd * kht;
        int idx = ksz - litup - 1;
        int threshold = sorted.samples_[idx];
        int maxs = sorted.samples_[ksz-1];
        LOG("deblur threshold: "<<threshold<<" max: "<<maxs);

        /** obliterate pixels below the threshold. **/
        for (int i = 0; i < ksz; ++i) {
            int s = kernel.samples_[i];
            s -= threshold;
            s = std::max(s, 0);
            /** arbitrary hard coded value to make the kernel png look pretty. **/
            s *= 2;
            kernel.samples_[i] = s;
        }

        /**
        hack: overwrite the image to ensure we have the right patch.
        for venus and crescent moon with people on bridge: 4520,1470,4642,1554
        **/
        /*image_.planes_.r_ = kernel;
        image_.planes_.g1_ = kernel;
        image_.planes_.g2_ = kernel;
        image_.planes_.b_ = kernel;
        return;*/

        /**
        apparently this method generates a kernel that's flipped x and y.
        ie rotated 180 degrees.
        so we need to correct that.
        **/
        rotate_180(kernel);

        /**
        hack: the image is really freaking big compared to the patch.
        crop it.
        **/
        /*LOG("deblur hack! cropping image");
        const int margin = 2 * std::max(kwd, kht);
        cropped.crop(l - margin, t - margin, r + margin, b + margin);*/

        /**
        okay now we're going to do the Lucy-Richardson deconvolution algorithm.
        this guy explains: https://stargazerslounge.com/topic/228147-lucy-richardson-deconvolution-so-what-is-it/
        also on wikipedia:

        estimate(i+1) = estimate(i) * ( ( observed / estimate(i) ** kernel ) ** flipped(kernel) )

        1. convolve estimate(i) with kernel.
        2. elementwise divide observed imaged by 1.
        3. convolve 2 with flipped kernel.
        4. elementwise multiply estimate(i) by 3.

        okay so there seems to be a kinda stupid problem.
        what do you do where the observed pixels are black?
        quite likely the estimate(i) pixels are also black.
        now you're dividing 0/0.
        **/

        /** allocate temporary storage. **/
        Plane estimate;
        Plane blurred;
        Plane error;
        Plane correction;
        int swd = cropped.width_;
        int sht = cropped.height_;
        estimate.init(swd, sht);
        blurred.init(swd, sht);
        error.init(swd, sht);
        correction.init(swd, sht);

        /** start with 50% gray. **/
        int ssz = swd * sht;
        for (int i = 0; i < ssz; ++i) {
            estimate.samples_[i] = 32768;
        }

        /**
        offset and scale the source image to be gray to white instead of black to white.
        this avoids the 0/0 problem above.
        **/
        LOG("deblur hack! changing image range to gray to white.");
        black_to_gray(cropped);

        int niterations = 1;
        for (int n = 1; ; ++n) {
            LOG("iteration: "<<n<<" of "<<niterations);

            /** convolve the first estimate with blur kernel **/
            if (n == 1) {
                blurred = estimate;
            } else {
                convolve_mt(estimate, kernel, blurred);
            }

            /** divide the observed image by the blurred estimate. **/
            divide(cropped, blurred, error);

            /** convolve the scaling factor with the flipped blur kernel **/
            rotate_180(kernel);
            convolve_mt(error, kernel, correction);
            rotate_180(kernel);

            /** multiply the estimate by the correction to get the new estimate. **/
            multiply(estimate, correction, 1.0);

            /** stop **/
            if (n == niterations) {
                break;
            }

            #if 0
            /**
            now we need to refine the kernel.
            it's exactly the same procedure.
            except the roles of the kernel and the estimated image are reversed.
            and we want to work in luminance space instead of rggb space.
            make use of the fact that convolution is commutative.
            **/
            /** these are the size of the source. **/
            convolve_mt(estimate, kernel, blurred);

            /** hack! save the blurred. **/
            /*image_.planes_.r_ = blurred;
            image_.planes_.g1_ = blurred;
            image_.planes_.g2_ = blurred;
            image_.planes_.b_ = blurred;
            return;*/

            divide(cropped, blurred, error);

            /** hack! save the error. **/
            for (int i = 0; i < ssz; ++i) {
                int s = error.samples_[i];
                s /= 2;
                error.samples_[i] = s;
            }
            image_.planes_.r_ = error;
            image_.planes_.g1_ = error;
            image_.planes_.g2_ = error;
            image_.planes_.b_ = error;
            return;

            rotate_180(estimate);
            convolve_mt(error, estimate, correction);

            /** hack! save the correction. **/
            image_.planes_.r_ = correction;
            image_.planes_.g1_ = correction;
            image_.planes_.g2_ = correction;
            image_.planes_.b_ = correction;
            return;

            rotate_180(estimate);
            /**
            size mismatch here.
            mathematically we can expand the kernel by padding with zeros
            on the top, left, bottom, and right.
            which means we could crop the middle of the correction.
            but need to verify.
            **/
            correction.crop(margin, margin, swd - margin, sht - margin);
            multiply(kernel, correction, 1.0);
            #endif
        }

        /** hack! save the updated kernel. **/
        /*image_.planes_.r_ = kernel;
        image_.planes_.g1_ = kernel;
        image_.planes_.g2_ = kernel;
        image_.planes_.b_ = kernel;
        return;*/

        LOG("deblur hack! changing image range to black to white.");
        gray_to_black(estimate);

        /** save it. **/
        image_.planes_.r_ = estimate;
        image_.planes_.g1_ = estimate;
        image_.planes_.g2_ = estimate;
        image_.planes_.b_ = estimate;
    }

    void rotate_180(
        Plane &src
    ) {
        flip_x(src);
        src.transpose();
        flip_x(src);
        src.transpose();
    }

    void flip_x(
        Plane &src
    ) {
        int wd = src.width_;
        int ht = src.height_;
        Plane dst;
        for (int y = 0; y < ht; ++y) {
            int xl = 0;
            int xr = wd - 1;
            while (xl < xr) {
                int sl = src.get(xl, y);
                int sr = src.get(xr, y);
                src.set(xl, y, sr);
                src.set(xr, y, sl);
                ++xl;
                --xr;
            }
        }
    }

    void black_to_gray(
        Plane &plane
    ) {
        int ssz = plane.width_ * plane.height_;
        for (int i = 0; i < ssz; ++i) {
            int s = plane.samples_[i];
            s = s / 2 + 32768;
            plane.samples_[i] = s;
        }
    }

    void gray_to_black(
        Planes &planes
    ) {
        gray_to_black(planes.r_);
        gray_to_black(planes.g1_);
        gray_to_black(planes.g2_);
        gray_to_black(planes.b_);
    }

    void gray_to_black(
        Plane &plane
    ) {
        int ssz = plane.width_ * plane.height_;
        for (int i = 0; i < ssz; ++i) {
            int s = plane.samples_[i];
            s = (s - 32768) * 2;
            plane.samples_[i] = s;
        }
    }

    void divide(
        Planes &numer,
        Planes &denom,
        Planes &result
    ) {
        divide(numer.r_, denom.r_, result.r_);
        divide(numer.g1_, denom.g1_, result.g1_);
        divide(numer.g2_, denom.g2_, result.g2_);
        divide(numer.b_, denom.b_, result.b_);
    }

    void divide(
        Plane &numer,
        Plane &denom,
        Plane &result
    ) {
        int ssz = numer.width_ * numer.height_;
        for (int i = 0; i < ssz; ++i) {
            int64_t n = numer.samples_[i];
            n *= 65536;
            int64_t d = denom.samples_[i];
            int64_t r = 0;
            if (d != 0) {
                r = n / d;
            }
            result.samples_[i] = r;
        }
    }

    void multiply(
        Planes &planes,
        Planes &correction,
        double reduction
    ) {
        multiply(planes.r_, correction.r_, reduction);
        multiply(planes.g1_, correction.g1_, reduction);
        multiply(planes.g2_, correction.g2_, reduction);
        multiply(planes.b_, correction.b_, reduction);
    }

    void multiply(
        Plane &plane,
        Plane &correction,
        double reduction
    ) {
        int ssz = plane.width_ * plane.height_;
        for (int i = 0; i < ssz; ++i) {
            double s = correction.samples_[i];
            s /= 65536.0;

            /**
            the scaling factor should be about 1.0.
            make it closer to 1.0 by the given reduction factor.
            **/
            s = (s - 1.0) * reduction + 1.0;

            s *= (double) plane.samples_[i];
            plane.samples_[i] = s;
        }
    }

    void decimate(
        Plane &plane
    ) {
        decimate_x(plane);
        plane.transpose();
        decimate_x(plane);
        plane.transpose();
    }

    void decimate_x(
        Plane &src
    ) {
        int wds = src.width_;
        int ht = src.height_;
        int wdd = wds / 2;
        Plane dst;
        dst.init(wdd, ht);
        for (int y = 0; y < ht; ++y) {
            for (int xd = 0; xd < wdd; ++xd) {
                int xs = xd * 2;
                int s0 = src.get(xs, y);
                int s1 = src.get(xs + 1, y);
                int d = (s0 + s1 + 1) / 2;
                dst.set(xd, y, d);
            }
        }
        src = std::move(dst);
    }

    void adjust_dynamic_range() {
        LOG("adjusting dynamic range...");
        if (opt_.drama_ <= 0.0) {
            LOG("dynamic range is not enabled.");
            return;
        }

        /**
        this technique has promise.
        but the algorithm needs attention.
        the shift and drama parameters work.
        but small values can make images look ridiculous.
        the downsampling technique has severe issues with macro blocking.
        the shifting technique doesn't play nicely with gamma correction.
        expanding mid range luminances works nicely for the moon.
        the shifting technique tends to just amplify noise in the dark areas.

        need to handle areas that overflow one component.
        that cap will shift the color.

        am looking at this paper:
        Image Display Algorithms for High and Low Dynamic Range Display Devices
        by Erik Reinhard Timo Kunkel Yoann Marion Kadi Bouatouch Jonathan Brouillat
        **/

        /** compute luminance for all pixels. **/
        maybe_compute_luminance();

        int window = 32;
        if (opt_.window_ > 0) {
            LOG("using user dynamic range window...");
            window = opt_.window_;
        }
        LOG("dynamic range window: "<<window);

        /** apply gaussian blur to get average lumance. **/
        Plane average(luminance_);
        average.gaussian(window);

        /** expand the noise floor a bit. **/
        int noise = image_.noise_;
        if (opt_.noise_ > 0) {
            LOG("using user noise floor...");
            noise = opt_.noise_;
        }
        LOG("noise floor: "<<noise);
        noise = noise * 150 / 100;

        /** dramatically move a pixel from its average luminance. **/
        double drama = opt_.drama_;
        LOG("dynamic range expansion factor: "<<drama);

        int wd = luminance_.width_;
        int ht = luminance_.height_;
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int lum = luminance_.get(x, y);

                /** ignore pixels below the noise floor. **/
                if (lum <= noise) {
                    continue;
                }
                /** ignore pixels below the noise floor and above the saturation level. **/
                if (lum >= kSaturated) {
                    continue;
                }

                /** rescale the average luminance. **/
                double avg_lum = average.get(x, y);

                /** move the colors of this pixel away from the average luminance. **/
                double target_lum = avg_lum + drama*(lum - avg_lum);

                /** get the components. **/
                double r = image_.planes_.r_.get(x, y);
                double g1 = image_.planes_.g1_.get(x, y);
                double g2 = image_.planes_.g2_.get(x, y);
                double b = image_.planes_.b_.get(x, y);

                /** scale the components. **/
                double factor = target_lum / lum;
                r *= factor;
                g1 *= factor;
                g2 *= factor;
                b *= factor;

                /** store the dynamically enhanced pixel **/
                image_.planes_.r_.set(x, y, r);
                image_.planes_.g1_.set(x, y, g1);
                image_.planes_.g2_.set(x, y, g2);
                image_.planes_.b_.set(x, y, b);
            }
        }
    }

    void maybe_compute_luminance() {
        /**
        optionally used multiple places.
        no sense computing luminance more than once.
        **/
        int wd = luminance_.width_;
        int ht = luminance_.height_;
        if (wd > 0 && ht > 0) {
            return;
        }

        compute_luminance(image_.planes_, luminance_);
    }

    void interpolate() {
        LOG("interpolating pixels...");

        /** halfsize option disables demosaicing. **/
        if (opt_.halfsize_) {
            LOG("interpolation disabled by --halfsize option.");
            return;
        }

        /**
        how do we keep saturated pixels from bleeding into other pixels
        and becoming not-saturated hot pink?

        n n n s s s s n n n
        1 3 3 1     1 3 3 1
           ?           ?

        n n s s s s s s n n
        1 3 3 1     1 3 3 1
           ?           ?

        n s s s s s s s s n
        1 3 3 1     1 3 3 1
           ?           ?
        **/

        /**
        applying the 1331 filter is really fast for horizontal pixels.
        so we transpose while small.
        add pixels horizontally.
        tanspose again while medium sized.
        add pixels horizontally again.
        **/
        image_.planes_.interpolate_1331_sat();

        /**
        at this point we have alignment issues.
        every plane has 1 too many rows and 1 too many columns.
        but every plane needs to chop a different set.
        right now every plane looks like this:
        r i r i r i
        i i i i i i
        r i r i r i
        i i i i i i
        where r is a "real" pixel and i is an interpolated pixels.
        we need to make them like up with the bayer pattern.
        r g r g r g
        g b g b g b
        r g r g r g
        g b g b g b
        **/
        int wd = image_.planes_.r_.width_ - 1;
        int ht = image_.planes_.r_.height_ - 1;
        image_.planes_.r_.crop(1, 1, wd+1, ht+1);
        image_.planes_.g1_.crop(0, 1, wd, ht+1);
        image_.planes_.g2_.crop(1, 0, wd+1, ht);
        image_.planes_.b_.crop(0, 0, wd, ht);

        LOG("interpolated width ="<<image_.planes_.r_.width_);
        LOG("interpolated height="<<image_.planes_.r_.height_);

        /** redo-desaturation after crop-shifting. **/
        desaturate_pixels();
    }

    void combine_greens() {
        LOG("combining greens...");
        int sz = image_.planes_.g1_.width_ * image_.planes_.g1_.height_;
        for (int i = 0; i < sz; ++i) {
            int g1 = image_.planes_.g1_.samples_[i];
            int g2 = image_.planes_.g2_.samples_[i];
            image_.planes_.g1_.samples_[i] = (g1 + g2 + 1)/2;
        }
    }

    void determine_auto_brightness() {
        LOG("determining auto brightness...");

        std::vector<int> histogram;
        histogram.resize(65536, 0);

        maybe_compute_luminance();

        int wd = luminance_.width_;
        int ht = luminance_.height_;
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int lum = luminance_.get(x, y);
                int idx = pin_to_16bits(lum);
                ++histogram[idx];
            }
        }

        for (int i = 65535; i >= 0; --i) {
            if (histogram[i] > 0) {
                double brightest = double(i) / 65535.0;
                double brightness = 1.0 / brightest;
                LOG("current brightness: "<<brightest);
                LOG("proper brightness: "<<brightness);
                break;
            }
        }

        auto_brightness_ = 0.0;
        if (opt_.auto_brightness_ < 0.0) {
            return;
        }

        int sz = wd * ht;
        int target = sz * opt_.auto_brightness_;
        int count = 0;
        for (int i = 65535; i >= 0; --i) {
            count += histogram[i];
            if (count >= target) {
                auto_brightness_ = 65535.0 / double(i);
                break;
            }
        }

    }

    void convert_to_srgb() {
        LOG("converting to sRGB...");
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

        /** account for user specified auto or linear brightness. **/
        double brightness = 1.0;
        if (opt_.auto_brightness_ >= 0.0) {
            LOG("using user auto brightness: "<<opt_.auto_brightness_<<" "<<auto_brightness_);
            brightness = auto_brightness_;
        } else if (opt_.linear_brightness_ > 0.0) {
            LOG("using user linear brightness: "<<opt_.linear_brightness_);
            brightness = opt_.linear_brightness_;
        }
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                mat[i][k] *= brightness;
            }
        }

        int sz = image_.planes_.r_.width_ * image_.planes_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = image_.planes_.r_.samples_[i];
            int in_g = image_.planes_.g1_.samples_[i];
            int in_b = image_.planes_.b_.samples_[i];

            double out_r;
            double out_g;
            double out_b;

            /** ensure saturated pixels stay saturated. **/
            if (in_r >= kSaturated
            &&  in_g >= kSaturated
            &&  in_b >= kSaturated) {
                out_r = 65535.0;
                out_g = 65535.0;
                out_b = 65535.0;
            } else {
                /** transform by matrix multiplication. **/
                out_r = mat[0][0]*in_r + mat[0][1]*in_g + mat[0][2]*in_b;
                out_g = mat[1][0]*in_r + mat[1][1]*in_g + mat[1][2]*in_b;
                out_b = mat[2][0]*in_r + mat[2][1]*in_g + mat[2][2]*in_b;

                /** ensure we don't change color on overflow. **/
                int maxc = std::max(std::max(out_r, out_g), out_b);
                if (maxc > 65535.0) {
                    double factor = 65535.0 / maxc;
                    out_r *= factor;
                    out_g *= factor;
                    out_b *= factor;
                }
            }

            /** overwrite old values. **/
            image_.planes_.r_.samples_[i] = out_r;
            image_.planes_.g1_.samples_[i] = out_g;
            image_.planes_.b_.samples_[i] = out_b;
        }
    }

    void enhance_colors() {
        LOG("enhancing colors...");
        double factor = opt_.color_enhancement_;
        if (factor <= 0.0 || factor == 1.0) {
            LOG("color enhancement is disabled.");
            return;
        }
        LOG("color enhancement factor: "<<factor);

        /**
        rgb to yuv
        Y =  0.257 * R + 0.504 * G + 0.098 * B +  16;
        U = -0.148 * R - 0.291 * G + 0.439 * B + 128;
        V =  0.439 * R - 0.368 * G - 0.071 * B + 128;
        yuv to rgb
        Y -= 16;
        U -= 128;
        V -= 128;
        R = 1.164 * Y             + 1.596 * V;
        G = 1.164 * Y - 0.392 * U - 0.813 * V;
        B = 1.164 * Y + 2.017 * U;
        re-arranging:
        Y =  0.257*1.164 * R + 0.504*1.164 * G + 0.098*1.164 * B;
        U = -0.148 * R - 0.291 * G + 0.439 * B;
        V =  0.439 * R - 0.368 * G - 0.071 * B;
        R = Y             + 1.596 * V;
        G = Y - 0.392 * U - 0.813 * V;
        B = Y + 2.017 * U;
        **/

        int sz = image_.planes_.r_.width_ * image_.planes_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = image_.planes_.r_.samples_[i];
            int in_g = image_.planes_.g1_.samples_[i];
            int in_b = image_.planes_.b_.samples_[i];

            double out_r;
            double out_g;
            double out_b;

            /** ensure saturated pixels stay saturated. **/
            if (in_r >= kSaturated
            &&  in_g >= kSaturated
            &&  in_b >= kSaturated) {
                out_r = 65535.0;
                out_g = 65535.0;
                out_b = 65535.0;
            } else {
                /** transform to yuv. **/
                double y = +0.257*in_r + 0.504*in_g + 0.098*in_b;
                double u = -0.148*in_r - 0.291*in_g + 0.439*in_b;
                double v = +0.439*in_r - 0.368*in_g - 0.071*in_b;

                /** enhance colors. **/
                y *= 1.164;
                u *= factor;
                v *= factor;

                /** transform back to rgb. **/
                out_r = y + 1.596*v;
                out_g = y - 0.392*u - 0.813*v;
                out_b = y + 2.017*u;

                /** ensure we don't change color on overflow. **/
                int maxc = std::max(std::max(out_r, out_g), out_b);
                if (maxc > 65535.0) {
                    double factor = 65535.0 / maxc;
                    out_r *= factor;
                    out_g *= factor;
                    out_b *= factor;
                }
            }

            /** overwrite old values. **/
            image_.planes_.r_.samples_[i] = out_r;
            image_.planes_.g1_.samples_[i] = out_g;
            image_.planes_.b_.samples_[i] = out_b;
        }
    }

    void apply_gamma() {
        LOG("applying gamma correction...");

        double gamma0;
        double gamma1;
        if (opt_.gamma0_ > 0.0 && opt_.gamma1_ > 0.0) {
            LOG("using user gamma correction...");
            gamma0 = opt_.gamma0_;
            gamma1 = opt_.gamma1_;
        } else {
            /** using camera's gamma. **/
            gamma0 = image_.camera_.gamma0_;
            gamma1 = image_.camera_.gamma1_;
        }
        double pwr = 1.0 / gamma0;
        double ts = gamma1;

        if (pwr == 1.0 && ts == 1.0) {
            LOG("gamma correction is disabled.");
            return;
        }

        double ipwr = 1.0 / pwr;
        LOG("gamma correction: "<<ipwr<<" "<<ts);

        int white = 0x10000;
        image_.planes_.apply_gamma(pwr, ts, white);
    }

    void scale_to_8bits() {
        LOG("scaling to 8 bits per sample..");
        image_.planes_.scale_to_8bits();
    }
};

}

int save_as_png(
    int argc,
    clo_argv_t argv
) {
    SaveAsPng sap;
    int exit_code = sap.run(argc, argv);

    return exit_code;
}
