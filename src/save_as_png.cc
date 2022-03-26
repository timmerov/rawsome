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
y  : G2/3  B/2
y+1: R/0   G1/1
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
#include <random>
#include <sstream>
#include <thread>
#include <vector>


namespace {

const int kFullySaturated = (1<<30)-1;

class SaveAsPng {
public:
    SaveAsPng() = default;
    ~SaveAsPng() = default;

    Options opt_;
    Image image_;
    std::vector<int> histogram_;
    Plane luminance_;
    int black_ = 0;
    int white_ = 65535;

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
        desaturate_pixels();
        show_special_pixel();
        scale_image();
        show_special_pixel();
        interpolate();
        combine_greens();
        show_special_pixel();
        auto_white_balance();
        show_special_pixel();
        determine_black_and_white();
        convert_to_srgb();
        show_special_pixel();
        enhance_colors();
        show_special_pixel();
        apply_user_gamma();
        show_special_pixel();
        apply_display_gamma();
        show_special_pixel();
        scale_to_8bits();
        image_.save_png(opt_.out_filename_.c_str());
        return 0;
    }

    void show_special_pixel() {
        /** show the special pixel. **/
        /*
        int x = 1924;
        int y = 906;
        int r = image_.planes_.r_.get(x, y);
        int g1 = image_.planes_.g1_.get(x, y);
        int g2 = image_.planes_.g2_.get(x, y);
        int b = image_.planes_.b_.get(x, y);
        LOG("special pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);
        */
        /** show its neighbors. **/
        /*
        --y;
        r = image_.planes_.r_.get(x, y);
        g1 = image_.planes_.g1_.get(x, y);
        g2 = image_.planes_.g2_.get(x, y);
        b = image_.planes_.b_.get(x, y);
        LOG("neighbor pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);
        ++y;
        --x;
        r = image_.planes_.r_.get(x, y);
        g1 = image_.planes_.g1_.get(x, y);
        g2 = image_.planes_.g2_.get(x, y);
        b = image_.planes_.b_.get(x, y);
        LOG("neighbor pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);
        x += 2;
        r = image_.planes_.r_.get(x, y);
        g1 = image_.planes_.g1_.get(x, y);
        g2 = image_.planes_.g2_.get(x, y);
        b = image_.planes_.b_.get(x, y);
        LOG("neighbor pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);
        --x;
        ++y;
        r = image_.planes_.r_.get(x, y);
        g1 = image_.planes_.g1_.get(x, y);
        g2 = image_.planes_.g2_.get(x, y);
        b = image_.planes_.b_.get(x, y);
        LOG("neighbor pixel x,y="<<x<<","<<y<<" rggb="<<r<<" "<<g1<<" "<<g2<<" "<<b);
        */
    }

    void desaturate_pixels() {
        LOG("desaturating pixels...");

        /** by component **/
        desaturate_pixels(image_.planes_.r_);
        desaturate_pixels(image_.planes_.g1_);
        desaturate_pixels(image_.planes_.g2_);
        desaturate_pixels(image_.planes_.b_);
    }

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

    void desaturate_pixels(
        Plane &plane
    ) {
        const int kSatMax = 16383;
        const int kSatHalf = kSatMax - kSaturation;
        const int kSatSize = 2*kSatHalf;
        const int kSatThreshold = kSatMax - kSatSize;

        histogram_.clear();
        histogram_.resize(kSatSize, 0);
        int sz = plane.width_ * plane.height_;
        for (int i = 0; i < sz; ++i) {
            int c = plane.samples_[i];
            c -= kSatThreshold;
            if (c >= kSatSize) {
                c = kSatSize-1;
            }
            if (c >= 0) {
                ++histogram_[c];
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
            LOG("saturation: "<<kSaturation<<" (default)");
            return;
        }

        /**
        look at the previous 12 values.

        if they kinda tail away then no pixels are saturated.
        use max of saturation+1 and expectation.

        if they're mostly small values but there's a big jump.
        then set the saturation before the big jump.
        **/
        bool big_jump = false;
        int idx = saturation - 12;
        int c0 = histogram_[idx];
        c0 = std::max(c0, 100);
        int threshold = 3*c0;
        for (++idx; idx <= saturation; ++idx) {
            int c1 = histogram_[idx];
            if (c1 > threshold) {
                /**
                we found a big jump.
                idx is the saturation value.
                might as well back it off a bit.
                **/
                saturation = idx - 2;
                saturation += kSatThreshold;
                LOG("saturation: "<<saturation<<" (saturated)");
                big_jump = true;
                break;
            }
        }

        /** no big jump. **/
        if (big_jump == false) {
            saturation += 2;
            saturation += kSatThreshold;
            saturation = std::max(saturation, kSaturation);
            LOG("saturation: "<<saturation<<" (small tail)");
            return;
        }

        /**
        replace saturated pixels with a guess.
        estimate the super-saturated range.
        **/
        int limit = saturation - kSatThreshold;
        int nbright = 0;
        for (int i = 0; i < limit; ++i) {
            nbright += histogram_[i];
        }
        int nsaturated = 0;
        for (int i = limit; i < kSatSize; ++i) {
            nsaturated += histogram_[i];
        }
        /**
        make the super saturated range somewhat less dense
        than the density of the bright pixels.
        **/
        int range = 3.0 * double(nsaturated) * double(limit) / double(nbright);
        LOG("supersaturated range="<<range);

        /**
        initialize the stupid c++ stupid random number generator.
        use a constant seed so we get consistent results every run.
        **/
        std::seed_seq ss{int(0x12345678), int(0x9ABCDEF0)};
        std::mt19937_64 rng;
        rng.seed(ss);
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        /**
        find satured pixels.
        replace them with a random value in the super-saturated range.
        **/
        for (int i = 0; i < sz; ++i) {
            int p = plane.samples_[i];
            if (p < saturation) {
                continue;
            }
            int rnd = range * unif(rng);
            plane.samples_[i] = saturation + rnd;
        }
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
        cam_mul.r_ = image_.camera_.wb_r_;
        cam_mul.g1_ = 1.0;
        cam_mul.g2_ = 1.0;
        cam_mul.b_ = image_.camera_.wb_b_;

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
        cam_mul.r_ /= kSaturation;
        cam_mul.g1_ /= kSaturation;
        cam_mul.g2_ /= kSaturation;
        cam_mul.b_ /= kSaturation;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** adjust to span full 16 bit range. **/
        cam_mul.r_ *= 65535.0;
        cam_mul.g1_ *= 65535.0;
        cam_mul.g2_ *= 65535.0;
        cam_mul.b_ *= 65535.0;
        //LOG("cam_mul="<<cam_mul.r_<<" "<<cam_mul.g1_<<" "<<cam_mul.g2_<<" "<<cam_mul.b_);

        /** white balance. **/
        image_.planes_.multiply(cam_mul);
    }

    void interpolate() {
        LOG("interpolating pixels...");

        /** halfsize option disables demosaicing. **/
        if (opt_.half_size_) {
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
        image_.planes_.interpolate_1331();

        /**
        at this point we have alignment issues.
        every plane has 1 too many rows and 1 too many columns.
        but every plane needs to chop a different set.
        right now every plane looks like this:
        r i r i r i r
        i i i i i i i
        r i r i r i r
        i i i i i i i
        r i r i r i r
        where r is a "real" pixel and i is an interpolated pixels.
        we need to make them like up with the bayer pattern.
        g2 b g2 b g2 b g2 b
        r g1 r g1 r g1 r g1
        g2 b g2 b g2 b g2 b
        r g1 r g1 r g1 r g1

        yes the actual bayer pattern is g2b/rg1.
        see notes in copy_raw_to_planes() in image.cc.

        real r pixels expect to be in column+0 and row+1.
        x x x x x x x
        r i r i r i r
        i i i i i i i
        r i r i r i r
        i i i i i i i
        r i r i r i r

        real g1 pixels expect to be in column+1 and row+1.
        x x x x x x x x
        x g i g i g i g
        x i i i i i i i
        x g i g i g i g
        x i i i i i i i
        x g i g i g i g

        real g2 pixels expect to be in column+0 and row+0.
        G i G i G i G
        i i i i i i i
        G i G i G i G
        i i i i i i i
        G i G i G i G

        real b pixels expect to be in column+1 and row+0.
        x b i b i b i b
        x i i i i i i i
        x b i b i b i b
        x i i i i i i i
        x b i b i b i b

        we can't use row 0 or column 0 as shown above.
        cause we don't have data for some of the components.
        therefore we must crop row 0 and column 0.
        where each component may have already cropped one or both.

        r must crop only column 0.
        g1 must crop neither.
        g2 must crop both row 0 and column 0.
        b must crop only row 0.
        **/
        int wd = image_.planes_.r_.width_ - 1;
        int ht = image_.planes_.r_.height_ - 1;
        image_.planes_.r_.crop(1, 0, wd+1, ht);
        image_.planes_.g1_.crop(0, 0, wd, ht);
        image_.planes_.g2_.crop(1, 1, wd+1, ht+1);
        image_.planes_.b_.crop(0, 1, wd, ht+1);

        LOG("interpolated width : "<<image_.planes_.r_.width_);
        LOG("interpolated height: "<<image_.planes_.r_.height_);
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

    void auto_white_balance() {
        LOG("auto balancing white...");
        if (opt_.auto_white_balance_ == false) {
            LOG("auto balancing white is disabled.");
            return;
        }

        int sz = image_.planes_.r_.width_ * image_.planes_.r_.height_;
        std::int64_t energy_r = 0;
        std::int64_t energy_g = 0;
        std::int64_t energy_b = 0;
        for (int i = 0; i < sz; ++i) {
            energy_r += image_.planes_.r_.samples_[i];
            energy_g += image_.planes_.g1_.samples_[i];
            energy_b += image_.planes_.b_.samples_[i];
        }

        double energy = energy_r + energy_g + energy_b;
        energy /= 3.0;
        double factor_r = energy / double(energy_r);
        double factor_g = energy / double(energy_g);
        double factor_b = energy / double(energy_b);
        LOG("factors: "<<factor_r<<" "<<factor_g<<" "<<factor_b);

        for (int i = 0; i < sz; ++i) {
            image_.planes_.r_.samples_[i]  *= factor_r;
            image_.planes_.g1_.samples_[i] *= factor_g;
            image_.planes_.b_.samples_[i]  *= factor_b;
        }
    }

    void determine_black_and_white() {
        LOG("determining black and white levels...");

        std::vector<int> histogram;
        const int kHistSize = 2 * 65536;
        histogram.resize(kHistSize, 0);

        compute_luminance(image_.planes_, luminance_);

        int wd = luminance_.width_;
        int ht = luminance_.height_;
        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int lum = luminance_.get(x, y);
                if (lum < 0) {
                    lum = 0;
                } else if (lum >= kHistSize) {
                    lum = kHistSize-1;
                }
                ++histogram[lum];
            }
        }

        /**
        find the brightest and darkest pixels.
        guides the user's non-auto choices.
        **/
        int cur_black = 0;
        int cur_white = 0;
        for (int i = kHistSize-1; i >= 0; --i) {
            if (histogram[i] > 0) {
                cur_white = i;
                break;
            }
        }
        for (int i = 0; i < kHistSize; ++i) {
            if (histogram[i] > 0) {
                cur_black = i;
                break;
            }
        }
        double obs_black = double(cur_black) / 65535.0;
        double obs_white = double(cur_white) / 65535.0;
        LOG("observed black: "<<obs_black);
        LOG("observed white: "<<obs_white);

        double auto_black = opt_.auto_black_;
        double auto_white = opt_.auto_white_;
        LOG("auto black: "<<auto_black);
        LOG("auto white: "<<auto_white);

        black_ = 0;
        white_ = 65535;
        const char *which_black = "default";
        const char *which_white = "default";

        /** user black has priority over user auto-black. **/
        if (opt_.black_ >= 0.0) {
            which_black = "user";
            black_ = opt_.black_ * 65535.0;
        } else if (auto_black >= 0.0) {
            which_black = "auto";
            int sz = wd * ht;
            int target = sz * auto_black;
            int count = 0;
            for (int i = 0; i < kHistSize; ++i) {
                int c = histogram[i];
                if (c > 0) {
                    count += c;
                    if (count >= target) {
                        black_ = i;
                        break;
                    }
                }
            }
        }
        /** user white has priority over user auto-white. **/
        if (opt_.white_ >= 0.0) {
            which_white = "user";
            white_ = opt_.white_ * 65535.0;
        } else if (auto_white >= 0.0) {
            which_white = "auto";
            int sz = wd * ht;
            int target = sz * auto_white;
            int count = 0;
            for (int i = kHistSize-1; i >= 0; --i) {
                int c = histogram[i];
                if (c > 0) {
                    count += c;
                    if (count >= target) {
                        white_ = i;
                        break;
                    }
                }
            }
        }

        double black = double(black_) / 65535.0;
        double white = double(white_) / 65535.0;
        LOG("black: "<<black<<" ("<<which_black<<")");
        LOG("white: "<<white<<" ("<<which_white<<")");
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

        /** account for user specified black and white levels. **/
        int range = white_ - black_;
        double scale = 65535.0 / double(range);
        for (int k = 0; k < 3; ++k) {
            for (int i = 0; i < 3; ++i) {
                mat[i][k] *= scale;
            }
        }

        int sz = image_.planes_.r_.width_ * image_.planes_.r_.height_;
        for (int i = 0; i < sz; ++i) {
            /** start with the original rgb. **/
            int in_r = image_.planes_.r_.samples_[i];
            int in_g = image_.planes_.g1_.samples_[i];
            int in_b = image_.planes_.b_.samples_[i];

            /** subtrack user black. **/
            in_r -= black_;
            in_g -= black_;
            in_b -= black_;

            double out_r;
            double out_g;
            double out_b;

            /** transform by matrix multiplication. **/
            out_r = mat[0][0]*in_r + mat[0][1]*in_g + mat[0][2]*in_b;
            out_g = mat[1][0]*in_r + mat[1][1]*in_g + mat[1][2]*in_b;
            out_b = mat[2][0]*in_r + mat[2][1]*in_g + mat[2][2]*in_b;

            /** ensure we don't change color on overflow. **/
            /*int maxc = std::max(std::max(out_r, out_g), out_b);
            if (maxc > 65535.0) {
                double factor = 65535.0 / maxc;
                out_r *= factor;
                out_g *= factor;
                out_b *= factor;
            }*/

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

            /** overwrite old values. **/
            image_.planes_.r_.samples_[i] = out_r;
            image_.planes_.g1_.samples_[i] = out_g;
            image_.planes_.b_.samples_[i] = out_b;
        }
    }

    void apply_user_gamma() {
        LOG("applying user gamma...");

        double pwr = opt_.gamma_;
        if (pwr == 1.0) {
            LOG("user gamma is disabled.");
            return;
        }

        image_.planes_.apply_user_gamma(pwr);
    }

    void apply_display_gamma() {
        LOG("applying gamma correction for displays...");

        /** use camera's gamma. **/
        double gamma0 = image_.camera_.gamma0_;
        double gamma1 = image_.camera_.gamma1_;
        double pwr = 1.0 / gamma0;
        double ts = gamma1;

        if (pwr == 1.0 && ts == 1.0) {
            LOG("gamma correction is disabled.");
            return;
        }

        double ipwr = 1.0 / pwr;
        LOG("gamma correction: "<<ipwr<<" "<<ts);

        int white = 0x10000;
        image_.planes_.apply_display_gamma(pwr, ts, white);
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
