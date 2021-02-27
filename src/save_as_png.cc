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

#include "dump.h"
#include "image.h"
#include "log.h"
#include "options.h"
#include "rs_png.h"

#include <libraw/libraw.h>

#include <cmath>
#include <ctime>
#include <iomanip>
#include <vector>
#include <sstream>

namespace {

/**
moon:   "/home/timmer/Pictures/2020-05-12/moon/IMG_0393.CR2"
santa:  "/home/timmer/Pictures/2021-01-29/santa.CR2"
comet1: "/home/timmer/Pictures/2020-07-11/IMG_0480.CR2"
comet2: "/home/timmer/Pictures/2020-07-11/IMG_0481.CR2"
***/

const int kFullySaturated = (1<<30)-1;

class SaveAsPng {
public:
    SaveAsPng() = default;
    ~SaveAsPng() = default;

    Options opt_;
    Image image_;
    std::vector<int> gamma_curve_;
    RggbPixel saturation_;
    std::vector<int> histogram_;
    RggbPixel saturated_;
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
        determine_saturation();
        desaturate_pixels();
        show_special_pixel();
        scale_image();
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
        write_png();
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

    void determine_saturation() {
        LOG("determining saturation...");

        if (opt_.saturation_ > 0) {
            LOG("using user saturation...");
            saturation_.r_ = opt_.saturation_;
            saturation_.g1_ = opt_.saturation_;
            saturation_.g2_ = opt_.saturation_;
            saturation_.b_ = opt_.saturation_;
        } else {
            /**
            the camera sensor saturates at less than the maximum value.
            **/
            saturation_.r_ = determine_saturation(image_.planes_.r_);
            saturation_.g1_ = determine_saturation(image_.planes_.g1_);
            saturation_.g2_ = determine_saturation(image_.planes_.g2_);
            saturation_.b_ = determine_saturation(image_.planes_.b_);
        }

        LOG("saturation is: "<<saturation_.r_<<" "<<saturation_.g1_<<" "<<saturation_.g2_<<" "<<saturation_.b_);

        /** the fully saturated values will be modified. **/
        saturated_ = saturation_;
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
            if (r >= saturated_.r_
            ||  g1 >= saturated_.g1_
            ||  g2 >= saturated_.g2_
            ||  b >= saturated_.b_) {
                image_.planes_.r_.samples_[i] = saturated_.r_;
                image_.planes_.g1_.samples_[i] = saturated_.g1_;
                image_.planes_.g2_.samples_[i] = saturated_.g2_;
                image_.planes_.b_.samples_[i] = saturated_.b_;
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
        cam_mul.r_ /= saturation_.r_;
        cam_mul.g1_ /= saturation_.g1_;
        cam_mul.g2_ /= saturation_.g2_;
        cam_mul.b_ /= saturation_.b_;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** adjust to span full 32 bit range. **/
        cam_mul.r_ *= 65535.0;
        cam_mul.g1_ *= 65535.0;
        cam_mul.g2_ *= 65535.0;
        cam_mul.b_ *= 65535.0;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** white balance. **/
        image_.planes_.multiply(cam_mul);

        /** update the saturated values. **/
        saturated_.r_ *= cam_mul.r_;
        saturated_.g1_ *= cam_mul.g1_;
        saturated_.g2_ *= cam_mul.g2_;
        saturated_.b_ *= cam_mul.b_;
        //LOG("saturated="<<saturated_.r_<<" "<<saturated_.g1_<<" "<<saturated_.g2_<<" "<<saturated_.b_);
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
        compute_luminance();

        /** compute the saturated luminance. **/
        int sat_lum = saturated_.r_*lum_rx_ + saturated_.g1_*lum_gx_ + saturated_.g2_*lum_gx_ + saturated_.b_*lum_bx_;

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
                if (lum >= sat_lum) {
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

    void compute_luminance() {
        /**
        optionally used multiple places.
        no sense computing luminance more than once.
        **/
        int wd = luminance_.width_;
        int ht = luminance_.height_;
        if (wd > 0 && ht > 0) {
            return;
        }

        /** convert canon rggb to srgb by multiplying by this matrix. **/
        static double mat[3][3] = {
            {+1.901824, -0.972035, +0.070211},
            {-0.229410, +1.659384, -0.429974},
            {+0.042001, -0.519143, +1.477141}
        };

        /** convert srgb to luminance by multiplying by this vector. **/
        static double lum_vec[3] = { 0.2125, 0.7154, 0.0721};

        /** combine them **/
        lum_rx_ = mat[0][0]*lum_vec[0] + mat[0][1]*lum_vec[1] + mat[0][2]*lum_vec[2];
        lum_gx_ = mat[1][0]*lum_vec[0] + mat[1][1]*lum_vec[1] + mat[1][2]*lum_vec[2];
        lum_bx_ = mat[2][0]*lum_vec[0] + mat[2][1]*lum_vec[1] + mat[2][2]*lum_vec[2];

        /** we still have two greens. **/
        lum_gx_ /= 2.0;

        wd = image_.planes_.r_.width_;
        ht = image_.planes_.r_.height_;
        luminance_.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            int r = image_.planes_.r_.samples_[i];
            int g1 = image_.planes_.g1_.samples_[i];
            int g2 = image_.planes_.g2_.samples_[i];
            int b = image_.planes_.b_.samples_[i];
            int lum = r*lum_rx_ + g1*lum_gx_ + g2*lum_gx_ + b*lum_bx_;
            luminance_.samples_[i] = lum;
        }
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
        image_.planes_.transpose();
        interpolate_horz_1331();
        image_.planes_.transpose();
        interpolate_horz_1331();

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

    void interpolate_horz_1331() {
        image_.planes_.r_.interpolate_horz_1331_mt(saturated_.r_);
        image_.planes_.g1_.interpolate_horz_1331_mt(saturated_.g1_);
        image_.planes_.g2_.interpolate_horz_1331_mt(saturated_.g2_);
        image_.planes_.b_.interpolate_horz_1331_mt(saturated_.b_);
    }

    void combine_greens() {
        LOG("combining greens...");
        int sz = image_.planes_.g1_.width_ * image_.planes_.g1_.height_;
        for (int i = 0; i < sz; ++i) {
            int g1 = image_.planes_.g1_.samples_[i];
            int g2 = image_.planes_.g2_.samples_[i];
            image_.planes_.g1_.samples_[i] = (g1 + g2 + 1)/2;
        }

        /** update the saturated values. **/
        saturated_.g1_ = (saturated_.g1_ + saturated_.g2_ + 1)/2;
        //LOG("saturated="<<saturated_.r_<<" "<<saturated_.g1_<<" "<<saturated_.g2_<<" "<<saturated_.b_);
    }

    void determine_auto_brightness() {
        LOG("determining auto brightness...");

        std::vector<int> histogram;
        histogram.resize(65536, 0);

        compute_luminance();

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
            if (in_r >= saturated_.r_
            &&  in_g >= saturated_.g1_
            &&  in_b >= saturated_.b_) {
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
            if (in_r >= saturated_.r_
            &&  in_g >= saturated_.g1_
            &&  in_b >= saturated_.b_) {
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
        gamma_curve(pwr, ts, white);

        apply_gamma(image_.planes_.r_);
        apply_gamma(image_.planes_.g1_);
        apply_gamma(image_.planes_.b_);
    }

    /** derived from dcraw. **/
    void gamma_curve(
        double pwr,
        double ts,
        int imax
    ){
        gamma_curve_.resize(0x10000);

        #define SQR(x) ((x)*(x))
        double g2 = 0.0;
        double g3 = 0.0;
        double g4 = 0.0;

        double bnd[2] = {0.0, 0.0};
        double r;

        pwr = pwr;
        ts = ts;
        g2 = g3 = g4 = 0;
        bnd[ts >= 1] = 1;
        if ((ts-1)*(pwr-1) <= 0) {
            for (int i = 0; i < 48; ++i) {
                g2 = (bnd[0] + bnd[1])/2;
                bnd[(std::pow(g2/ts,-pwr) - 1)/pwr - 1/g2 > -1] = g2;
            }
            g3 = g2 / ts;
            g4 = g2 * (1/pwr - 1);
        }
        for (int i = 0; i < 0x10000; ++i) {
            gamma_curve_[i] = 0xffff;
            r = (double) i / imax;
            if (r < 1) {
                gamma_curve_[i] = 0x10000 * (r < g3 ? r*ts : std::pow(r,pwr)*(1+g4)-g4);
            }
        }

        /*std::stringstream ss;
        for (int i = 0; i < 0x10000; i += 0x100) {
            ss<<gamma_curve_[i]<<" ";
        }
        LOG("gamma_curve=[ "<<ss.str()<<"]");*/
    }

    void apply_gamma(
        Plane& plane
    ) {
        int sz = plane.width_ * plane.height_;
        for (int i = 0; i < sz; ++i) {
            int c = plane.samples_[i];
            c = pin_to_16bits(c);
            c = gamma_curve_[c];
            plane.samples_[i] = c;
        }
    }

    void scale_to_8bits() {
        LOG("scaling to 8 bits per sample..");
        double factor = 255.0/65535.0;
        image_.planes_.r_.multiply(factor);
        image_.planes_.g1_.multiply(factor);
        image_.planes_.b_.multiply(factor);
    }

    void write_png() {
        LOG("writing to: \""<<opt_.out_filename_<<"\"...");
        int wd = image_.planes_.r_.width_;
        int ht = image_.planes_.r_.height_;
        Png png;
        png.init(wd, ht);
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int r = image_.planes_.r_.get(x, y);
                int g = image_.planes_.g1_.get(x, y);
                int b = image_.planes_.b_.get(x, y);
                int idx = 3*x + y*png.stride_;
                png.data_[idx] = pin_to_8bits(r);
                png.data_[idx+1] = pin_to_8bits(g);
                png.data_[idx+2] = pin_to_8bits(b);
            }
        }
        png.write(opt_.out_filename_.c_str());
    }

    int pin_to_8bits(
        int x
    ) {
        if (x < 0) {
            return 0;
        }
        if (x > 255) {
            return 255;
        }
        return x;
    }

    int pin_to_16bits(
        int x
    ) {
        if (x < 0) {
            return 0;
        }
        if (x > 65535) {
            return 65535;
        }
        return x;
    }
};

}

int save_as_png(
    int argc,
    clo_argv_t argv
) {

    SaveAsPng rawsome;
    int exit_code = rawsome.run(argc, argv);

    return exit_code;
}
