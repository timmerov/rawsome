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

const int kWidth = 640/2;
const int kHeight = 360/2;
const int kSize = kWidth * kHeight;
const int kThreshold = 3000;
const double kJupiterX = 103.5;
const double kJupiterY = 76.5;

class RggbPixel {
public:
    int r_ = 0;
    int g1_ = 0;
    int g2_ = 0;
    int b_ = 0;
};

class Balance {
public:
    double r_ = 1.0;
    double g1_ = 1.0;
    double g2_ = 1.0;
    double b_ = 1.0;
};

class Rawsome {
public:
    Rawsome() = default;
    ~Rawsome() {
        if (is_loaded_) {
            raw_image_.recycle();
        }
    }

    Options opt_;
    LibRaw raw_image_;
    bool is_loaded_ = false;
    Image image_;
    RggbPixel black_;
    int noise_;
    std::vector<int> gamma_curve_;
    RggbPixel saturation_;
    std::vector<int> histogram_;
    RggbPixel saturated_;
    Plane luminance_;
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

        load_raw_image();
        if (is_loaded_ == false) {
            return 1;
        }

        show_camera_parameters();
        copy_raw_to_image();
        determine_black();
        crop_black();
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
        int r = image_.r_.get(x, y);
        int g1 = image_.g1_.get(x, y);
        int g2 = image_.g2_.get(x, y);
        int b = image_.b_.get(x, y);
        LOG("special pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);*/
    }

    void load_raw_image() {
        LOG("loading raw image from: \""<<opt_.in_filename_<<"\"...");
        is_loaded_ = false;
        raw_image_.open_file(opt_.in_filename_.c_str());
        int raw_wd = raw_image_.imgdata.sizes.width;
        int raw_ht = raw_image_.imgdata.sizes.height;
        LOG("raw_wd="<<raw_wd);
        LOG("raw_ht="<<raw_ht);
        if (raw_wd <= 0 || raw_ht <= 0) {
            LOG("File not opened: "<<opt_.in_filename_);
            return;
        }
        raw_image_.unpack();
        raw_image_.raw2image();
        //dump(raw_image_);
        is_loaded_ = true;
    }

    void show_camera_parameters() {
        double gamma0 = raw_image_.imgdata.params.gamm[0];
        double gamma1 = raw_image_.imgdata.params.gamm[1];
        gamma0 = 1.0 / gamma0;

        double white_balance_r = raw_image_.imgdata.color.cam_mul[0];
        double white_balance_g = raw_image_.imgdata.color.cam_mul[1];
        double white_balance_b = raw_image_.imgdata.color.cam_mul[2];
        white_balance_r /= white_balance_g;
        white_balance_b /= white_balance_g;

        auto tm = std::localtime(&raw_image_.imgdata.other.timestamp);

        LOG("make         : "<<raw_image_.imgdata.idata.make);
        LOG("model        : "<<raw_image_.imgdata.idata.model);
        LOG("lens         : "<<raw_image_.imgdata.lens.Lens);
        LOG("gamma        : "<<gamma0<<" "<<gamma1);
        LOG("white balance: R="<<white_balance_r<<" B="<<white_balance_b);
        LOG("iso          : "<<raw_image_.imgdata.other.iso_speed);
        LOG("shutter      : "<<raw_image_.imgdata.other.shutter<<"s");
        LOG("aperture     : f/"<<raw_image_.imgdata.other.aperture);
        LOG("focal_len    : "<<raw_image_.imgdata.other.focal_len);
        LOG("timestamp    : "<<std::put_time(tm, "%c %Z"));
        LOG("temp         : "<<raw_image_.imgdata.other.CameraTemperature);
    }

    void copy_raw_to_image() {
        LOG("copying raw to image...");
        /**
        raw_image_.imgdata.rawdata.raw_image is the raw samples in bayer format.
        evem rows: G B G B G B ...
        odd  rows: R G R G R G ...
        we want this pattern:
        every row: R G G B R G G B ...

        i don't know why the bayer pattern is GB/RG.
        i was expecting it to be RG/GB.
        but experimentation shows otherwise.
        maybe it has to do with the odd number of rows?
        and the pattern is bottom justified?
        dunno.
        weird.
        **/
        int raw_wd = raw_image_.imgdata.sizes.width;
        int raw_ht = raw_image_.imgdata.sizes.height;
        int wd = raw_wd / 2;
        int ht = raw_ht / 2;
        image_.init(wd, ht);

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int raw_idx = 2*x + 2*y*raw_wd;

                /** extract rggb from the GB/RG bayer pattern **/
                int g2 = raw_image_.imgdata.rawdata.raw_image[raw_idx];
                int b  = raw_image_.imgdata.rawdata.raw_image[raw_idx+1];
                int r  = raw_image_.imgdata.rawdata.raw_image[raw_idx+raw_wd];
                int g1 = raw_image_.imgdata.rawdata.raw_image[raw_idx+raw_wd+1];

                image_.r_.set(x, y, r);
                image_.g1_.set(x, y, g1);
                image_.g2_.set(x, y, g2);
                image_.b_.set(x, y, b);
            }
        }
    }

    void determine_black() {
        LOG("determining black...");
        /**
        the top rows and left columns of pixels are black.
        also set the black noise level.
        **/
        noise_ = 0;
        black_.r_ = determine_black(image_.r_);
        black_.g1_ = determine_black(image_.g1_);
        black_.g2_ = determine_black(image_.g2_);
        black_.b_ = determine_black(image_.b_);
        LOG("black is: "<<black_.r_<<" "<<black_.g1_<<" "<<black_.g2_<<" "<<black_.b_);

        /** adjust noise for black levels. **/
        int min_black = std::min(black_.r_, black_.g1_);
        min_black = std::min(min_black, black_.g2_);
        min_black = std::min(min_black, black_.b_);
        noise_ -= min_black;
        LOG("noise is: "<<noise_);
    }

    int determine_black(
        Plane &plane
    ) {
        std::int64_t sum = 0;
        /**
        the left 37 columns are black.
        the top 16 rows are black.
        row 17 is garbage.
        the rest of the pixels are the actual image.
        **/
        int top_count = 16 * plane.width_;
        for (int y = 0; y < 16; ++y) {
            for (int x = 0; x < plane.width_; ++x) {
                int c = plane.get(x, y);
                sum += c;
                noise_ = std::max(noise_, c);
            }
        }

        int left_count = 37 * (plane.height_ - 16);
        for (int y = 16; y < plane.height_; ++y) {
            for (int x = 0; x < 37; ++x) {
                sum += plane.get(x, y);
            }
        }

        int count = top_count + left_count;
        int black = sum / count;
        return black;
    }

    void crop_black() {
        LOG("cropping black pixels...");
        /**
        the left 37 columns are black.
        the top 16 rows are black.
        row 17 is garbage.
        the rest of the pixels are the actual image.
        **/
        /** left, top, right, bottom **/
        image_.crop(38, 18, image_.r_.width_, image_.r_.height_);
        LOG("cropped width ="<<image_.r_.width_);
        LOG("cropped height="<<image_.r_.height_);
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
            saturation_.r_ = determine_saturation(image_.r_);
            saturation_.g1_ = determine_saturation(image_.g1_);
            saturation_.g2_ = determine_saturation(image_.g2_);
            saturation_.b_ = determine_saturation(image_.b_);
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
        const int kSatExpected = 13583;
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
        int wd = image_.r_.width_;
        int ht = image_.r_.height_;
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            int r = image_.r_.samples_[i];
            int g1 = image_.g1_.samples_[i];
            int g2 = image_.g2_.samples_[i];
            int b = image_.b_.samples_[i];
            if (r >= saturated_.r_
            ||  g1 >= saturated_.g1_
            ||  g2 >= saturated_.g2_
            ||  b >= saturated_.b_) {
                image_.r_.samples_[i] = saturated_.r_;
                image_.g1_.samples_[i] = saturated_.g1_;
                image_.g2_.samples_[i] = saturated_.g2_;
                image_.b_.samples_[i] = saturated_.b_;
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
        Balance cam_mul;
        if (opt_.wb_r_ > 0.0 && opt_.wb_b_ >= 0.0) {
            LOG("using user white balance...");
            cam_mul.r_ = opt_.wb_r_;
            cam_mul.g1_ = 1.0;
            cam_mul.g2_ = 1.0;
            cam_mul.b_ = opt_.wb_b_;
        } else {
            auto& raw_cam_mul = raw_image_.imgdata.rawdata.color.cam_mul;
            cam_mul.r_ = raw_cam_mul[0];
            cam_mul.g1_ = raw_cam_mul[1];
            cam_mul.g2_ = raw_cam_mul[3];
            cam_mul.b_ = raw_cam_mul[2];
            //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);
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
        cam_mul.r_ /= saturation_.r_ - black_.r_;
        cam_mul.g1_ /= saturation_.g1_ - black_.g1_;
        cam_mul.g2_ /= saturation_.g2_ - black_.g2_;
        cam_mul.b_ /= saturation_.b_ - black_.b_;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** adjust to span full 32 bit range. **/
        cam_mul.r_ *= 65535.0;
        cam_mul.g1_ *= 65535.0;
        cam_mul.g2_ *= 65535.0;
        cam_mul.b_ *= 65535.0;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** subtract black and white balance. **/
        image_.r_.scale(black_.r_, cam_mul.r_);
        image_.g1_.scale(black_.g1_, cam_mul.g1_);
        image_.g2_.scale(black_.g2_, cam_mul.g2_);
        image_.b_.scale(black_.b_, cam_mul.b_);

        /** update the saturated values. **/
        saturated_.r_ = (saturated_.r_ - black_.r_) * cam_mul.r_;
        saturated_.g1_ = (saturated_.g1_ - black_.g1_) * cam_mul.g1_;
        saturated_.g2_ = (saturated_.g2_ - black_.g2_) * cam_mul.g2_;
        saturated_.b_ = (saturated_.b_ - black_.b_) * cam_mul.b_;
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

        /** convert canon rggb to srgb by multiplying by this matrix. **/
        static double mat[3][3] = {
            {+1.901824, -0.972035, +0.070211},
            {-0.229410, +1.659384, -0.429974},
            {+0.042001, -0.519143, +1.477141}
        };

        /** convert srgb to luminance by multiplying by this vector. **/
        static double lum_vec[3] = { 0.2125, 0.7154, 0.0721};

        /** combine them **/
        double rx = mat[0][0]*lum_vec[0] + mat[0][1]*lum_vec[1] + mat[0][2]*lum_vec[2];
        double gx = mat[1][0]*lum_vec[0] + mat[1][1]*lum_vec[1] + mat[1][2]*lum_vec[2];
        double bx = mat[2][0]*lum_vec[0] + mat[2][1]*lum_vec[1] + mat[2][2]*lum_vec[2];

        /** we still have two greens. **/
        gx /= 2.0;

        /** compute luminance for all pixels. **/
        int wd = image_.r_.width_;
        int ht = image_.r_.height_;
        luminance_.init(wd, ht);
        int sz = wd * ht;
        for (int i = 0; i < sz; ++i) {
            int r = image_.r_.samples_[i];
            int g1 = image_.g1_.samples_[i];
            int g2 = image_.g2_.samples_[i];
            int b = image_.b_.samples_[i];
            int lum = r*rx + g1*gx + g2*gx + b*bx;
            luminance_.samples_[i] = lum;
        }

        /** compute the saturated luminance. **/
        int sat_lum = saturated_.r_*rx + saturated_.g1_*gx + saturated_.g2_*gx + saturated_.b_*bx;

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
        int noise = noise_;
        if (opt_.noise_ > 0) {
            LOG("using user noise floor...");
            noise = opt_.noise_;
        }
        LOG("noise floor: "<<noise);
        noise = noise * 150 / 100;

        /** dramatically move a pixel from its average luminance. **/
        double drama = opt_.drama_;
        LOG("dynamic range expansion factor: "<<drama);

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
                double r = image_.r_.get(x, y);
                double g1 = image_.g1_.get(x, y);
                double g2 = image_.g2_.get(x, y);
                double b = image_.b_.get(x, y);

                /** scale the components. **/
                double factor = target_lum / lum;
                r *= factor;
                g1 *= factor;
                g2 *= factor;
                b *= factor;

                /** store the dynamically enhanced pixel **/
                image_.r_.set(x, y, r);
                image_.g1_.set(x, y, g1);
                image_.g2_.set(x, y, g2);
                image_.b_.set(x, y, b);
            }
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
        image_.transpose();
        interpolate_horz_1331();
        image_.transpose();
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
        int wd = image_.r_.width_ - 1;
        int ht = image_.r_.height_ - 1;
        image_.r_.crop(1, 1, wd+1, ht+1);
        image_.g1_.crop(0, 1, wd, ht+1);
        image_.g2_.crop(1, 0, wd+1, ht);
        image_.b_.crop(0, 0, wd, ht);

        LOG("interpolated width ="<<image_.r_.width_);
        LOG("interpolated height="<<image_.r_.height_);

        /** redo-desaturation after crop-shifting. **/
        desaturate_pixels();
    }

    void interpolate_horz_1331() {
        image_.r_.interpolate_horz_1331_mt(saturated_.r_);
        image_.g1_.interpolate_horz_1331_mt(saturated_.g1_);
        image_.g2_.interpolate_horz_1331_mt(saturated_.g2_);
        image_.b_.interpolate_horz_1331_mt(saturated_.b_);
    }

    void combine_greens() {
        LOG("combining greens...");
        int sz = image_.g1_.width_ * image_.g1_.height_;
        for (int i = 0; i < sz; ++i) {
            int g1 = image_.g1_.samples_[i];
            int g2 = image_.g2_.samples_[i];
            image_.g1_.samples_[i] = (g1 + g2 + 1)/2;
        }

        /** update the saturated values. **/
        saturated_.g1_ = (saturated_.g1_ + saturated_.g2_ + 1)/2;
        //LOG("saturated="<<saturated_.r_<<" "<<saturated_.g1_<<" "<<saturated_.g2_<<" "<<saturated_.b_);
    }

    void determine_auto_brightness() {
        LOG("determining auto brightness...");

        std::vector<int> histogram;
        histogram.resize(65536, 0);

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

        int sz = image_.r_.width_ * image_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = image_.r_.samples_[i];
            int in_g = image_.g1_.samples_[i];
            int in_b = image_.b_.samples_[i];

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
            image_.r_.samples_[i] = out_r;
            image_.g1_.samples_[i] = out_g;
            image_.b_.samples_[i] = out_b;
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

        int sz = image_.r_.width_ * image_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = image_.r_.samples_[i];
            int in_g = image_.g1_.samples_[i];
            int in_b = image_.b_.samples_[i];

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
            image_.r_.samples_[i] = out_r;
            image_.g1_.samples_[i] = out_g;
            image_.b_.samples_[i] = out_b;
        }
    }

    void apply_gamma() {
        LOG("applying gamma correction...");

        double pwr;
        double ts;
        if (opt_.gamma0_ > 0.0 && opt_.gamma1_ > 0.0) {
            LOG("using user gamma correction...");
            pwr =  1.0 / opt_.gamma0_;
            ts = opt_.gamma1_;
        } else {
            /** using camera's gamma. **/
            pwr = raw_image_.imgdata.params.gamm[0];
            ts = raw_image_.imgdata.params.gamm[1];
        }

        if (pwr == 1.0 && ts == 1.0) {
            LOG("gamma correction is disabled.");
            return;
        }

        double ipwr = 1.0 / pwr;
        LOG("gamma correction: "<<ipwr<<" "<<ts);

        int white = 0x10000;
        gamma_curve(pwr, ts, white);

        apply_gamma(image_.r_);
        apply_gamma(image_.g1_);
        apply_gamma(image_.b_);
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
        image_.r_.scale(0, factor);
        image_.g1_.scale(0, factor);
        image_.b_.scale(0, factor);
    }

    void write_png() {
        LOG("writing to: \""<<opt_.out_filename_<<"\"...");
        int wd = image_.r_.width_;
        int ht = image_.r_.height_;
        Png png;
        png.init(wd, ht);
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int r = image_.r_.get(x, y);
                int g = image_.g1_.get(x, y);
                int b = image_.b_.get(x, y);
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

    Image stack_;
    double jupiterx_ = 0.0;
    double jupitery_ = 0.0;

    int special_stack() {

        stack_.init(kWidth, kHeight);

        //const int kFirst = 873;
        //const int kLast = 882;
        const int kFirst = 884;
        const int kLast = 933;
        //const int kLast = 886;
        //const int kLast = 884;
        int count = 0;
        for (int i = kFirst; i <= kLast; ++i) {
            if (i == 890) {
                continue;
            }
            std::stringstream ss;
            ss<<"/home/timmer/Pictures/2020-12-22/stackable/IMG_0"<<i<<".CR2";
            opt_.in_filename_ = ss.str();
            ss.clear();
            ss.str("");
            ss<<"/home/timmer/Pictures/2020-12-22/stackable/IMG_0"<<i<<".png";
            opt_.out_filename_ = ss.str();

            load_raw_image();
            if (is_loaded_ == false) {
                continue;
            }
            ++count;

            copy_raw_to_image();
            determine_black();
            crop_black();

            //histogram();
            hack_crop();
            find_jupiter();

            determine_saturation();
            //desaturate_pixels();
            scale_image();

            shift_image();
            stack_image();

            //adjust_dynamic_range();
            interpolate();
            combine_greens();
            convert_to_srgb();
            //enhance_colors();
            apply_gamma();
            scale_to_8bits();
            write_png();
        }

        opt_.out_filename_  = "/home/timmer/Pictures/2020-12-22/stackable/stack.png";
        image_.init(kWidth, kHeight);
        for (int k = 0; k < kSize; ++k) {
            image_.r_.samples_[k] = (stack_.r_.samples_[k] + count/2) / count;
            image_.g1_.samples_[k] = (stack_.g1_.samples_[k] + count/2) / count;
            image_.g2_.samples_[k] = (stack_.g2_.samples_[k] + count/2) / count;
            image_.b_.samples_[k] = (stack_.b_.samples_[k] + count/2) / count;
        }

        image_.crop(2, 2, kWidth-2, kHeight-2);

        //opt_.drama_ = 1.1;
        //opt_.window_ = 4;
        //adjust_dynamic_range();

        //interpolate();
        combine_greens();

        //opt_.linear_brightness_ = 40.0;
        convert_to_srgb();

        //enhance_colors();

        opt_.gamma0_ = 10;
        opt_.gamma1_ = 20.0;
        apply_gamma();

        scale_to_8bits();
        write_png();

        return 0;
    }

    void hack_crop() {
        int wd = image_.r_.width_;
        int ht = image_.r_.height_;

        int maxx = 0;
        int maxy = 0;
        int minx = wd;
        int miny = ht;

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int r = image_.r_.get(x, y);
                int g1 = image_.g1_.get(x, y);
                int g2 = image_.g2_.get(x, y);
                int b = image_.b_.get(x, y);
                if (r >= kThreshold
                ||  g1 >= kThreshold
                ||  g2 >= kThreshold
                ||  b >= kThreshold) {
                    maxx = std::max(maxx, x);
                    maxy = std::max(maxy, y);
                    minx = std::min(minx, x);
                    miny = std::min(miny, y);
                }
            }
        }

        int centerx = (minx + maxx + 1) / 2;
        int centery = (miny + maxy + 1) / 2;
        LOG("center="<<centerx<<","<<centery);

        int left = centerx - kWidth/2;
        int top = centery - kHeight/2;
        image_.crop(left, top, left+kWidth, top+kHeight);
    }

    void histogram() {
        std::vector<int> hist;
        hist.resize(16384, 0);

        int sz = image_.r_.width_ * image_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            int r = image_.r_.samples_[i];
            int g1 = image_.g1_.samples_[i];
            int g2 = image_.g2_.samples_[i];
            int b = image_.b_.samples_[i];
            r = pin_to_14bits(r);
            g1 = pin_to_14bits(g1);
            g2 = pin_to_14bits(g2);
            b = pin_to_14bits(b);
            ++hist[r];
            ++hist[g1];
            ++hist[g2];
            ++hist[b];
        }

        LOG("histogram:");
        for (int i = 0; i < 16384;i += 32) {
            std::stringstream ss;
            ss<<i<<":";
            for (int k = 0; k < 32; ++k) {
                ss<<" "<<hist[i+k];
            }
            LOG(ss.str());
        }
    }

    int pin_to_14bits(
        int x
    ) {
        if (x < 0) {
            return 0;
        }
        if (x > 16384) {
            return 16384;
        }
        return x;
    }

    void find_jupiter() {
        int wd = image_.r_.width_ / 2;
        int ht = image_.r_.height_;
        double sumx = 0;
        double sumy = 0;
        int sumw = 0;
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int r = image_.r_.get(x, y);
                int g1 = image_.g1_.get(x, y);
                int g2 = image_.g2_.get(x, y);
                int b = image_.b_.get(x, y);
                r -= kThreshold;
                g1 -= kThreshold;
                g2 -= kThreshold;
                b -= kThreshold;
                if (r > 0) {
                    sumx += (x + 0.25) * r;
                    sumy += (y + 0.25) * r;
                    sumw += r;
                }
                if (g1 > 0) {
                    sumx += (x + 0.75) * g1;
                    sumy += (y + 0.25) * g1;
                    sumw += g1;
                }
                if (g2 > 0) {
                    sumx += (x + 0.25) * g2;
                    sumy += (y + 0.75) * g2;
                    sumw += g2;
                }
                if (b > 0) {
                    sumx += (x + 0.75) * b;
                    sumy += (y + 0.75) * b;
                    sumw += b;
                }
            }
        }
        jupiterx_ = sumx / double(sumw);
        jupitery_ = sumy / double(sumw);
        LOG("jupiter: "<<jupiterx_<<","<<jupitery_);
    }

    void shift_image() {
        double dx = jupiterx_ - kJupiterX;
        double dy = jupitery_ - kJupiterY;
        LOG("dx="<<dx<<" dy="<<dy);

        /**
        interpolate the pixel at destination kJupiterX=103.5, kJupiterY=76.5.
        the first problem is the coordinates of the source pixels are offset.
        different amounts depending on which component.
        the coordinates for the r samples are x+0.25, y+0.25 where x,y are integers.
        the coordinates for the b samples are x+0.75, y+0.75.

        dx is the amount we need to add to kJupiterX to get the center of the source pixel.
        **/

        shift_image(image_.r_,  dx-0.25);
        shift_image(image_.g1_, dx+0.25);
        shift_image(image_.g2_, dx+0.25);
        shift_image(image_.b_,  dx-0.25);
        image_.transpose();
        shift_image(image_.r_,  dy-0.25);
        shift_image(image_.g1_, dy-0.25);
        shift_image(image_.g2_, dy+0.25);
        shift_image(image_.b_,  dy+0.25);
        image_.transpose();
    }

    void shift_image(
        Plane &src,
        double dx
    ) {
        int wd = src.width_;
        int ht = src.height_;

        Plane dst;
        dst.init(wd, ht);

        double pi = std::atan(1.0)*4.0;

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int srcx = x + 0.5 + dx;
                if (srcx < 2 || srcx >= wd-2) {
                    continue;
                }
                int s0 = src.get(srcx - 1, y);
                int s1 = src.get(srcx + 0, y);
                int s2 = src.get(srcx + 1, y);
                int s3 = src.get(srcx + 2, y);
                double sa = (- s0 + 3*s1 + 3*s2 - s3 + 2.0) / 4.0;
                double sb = s1;
                double fraction = srcx - int(srcx);
                if (fraction > 0.5) {
                    sb = s2;
                    fraction = 1.0 - fraction;
                }
                double w = std::sin(pi*fraction)*0.5 + 0.5;
                double sc = w*sa + (1.0 - w)*sb;
                dst.set(x, y, sc);

                /*if (y == 103) {
                    if (x >= 125-4 && x <= 125+5) {
                        LOG("x="<<x<<" srcx="<<srcx<<" s0="<<s0<<" s1="<<s1<<" s2="<<s2<<" sa="<<sa<<" f="<<fraction<<" sb="<<sb<<" w="<<w<<" sc="<<sc);
                    }
                }*/
            }
        }

        src = std::move(dst);
    }

    void stack_image() {
        stack_image(stack_.r_, image_.r_);
        stack_image(stack_.g1_, image_.g1_);
        stack_image(stack_.g2_, image_.g2_);
        stack_image(stack_.b_, image_.b_);
    }

    void stack_image(
        Plane &dst,
        Plane &src
    ) {
        int sz = dst.width_ * dst.height_;
        for (int i = 0; i < sz; ++i) {
            dst.samples_[i] += src.samples_[i];
        }
    }
};

}

int main(
    int argc,
    clo_argv_t argv
) noexcept {
    rs_log::init("rawsome.log");

    Rawsome rawsome;
    //int exit_code = rawsome.run(argc, argv);
    int exit_code = rawsome.special_stack();
    (void) argc;
    (void) argv;

    return exit_code;
}
