/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
functions specific to canon.
**/

#include "canon.h"
#include "log.h"

#include <cmath>

namespace {

int determine_black(
    Plane &plane,
    int &noise
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
            noise = std::max(noise, c);
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

void fix_bad_pixel(
    Plane &plane,
    int x,
    int y
) {
    LOG("fixing bad pixel at "<<x<<","<<y);

    /**
    we found the bad pixel by inspection.
    x,y are in final coordinates.
    ie after interpolation. hence the /2.
    also after cropping border pixels. hence the +2.
    **/
    x = x / 2 + 2;
    y = y / 2 + 2;

    /** overwrite the bad pixel with the average of its neighbors. **/
    int c0 = plane.get(x, y-1);
    int c1 = plane.get(x-1, y);
    int c2 = plane.get(x+1, y);
    int c3 = plane.get(x, y+1);
    int c = (c0 + c1 + c2 + c3 + 2) / 4;
    plane.set(x, y, c);
}
} // anonymous namespace

void determine_black(
    Planes &planes,
    RggbPixel &black,
    int &noise
) {
    LOG("determining black...");
    /**
    the top rows and left columns of pixels are black.
    also set the black noise level.
    **/
    noise = 0;
    black.r_ = determine_black(planes.r_, noise);
    black.g1_ = determine_black(planes.g1_, noise);
    black.g2_ = determine_black(planes.g2_, noise);
    black.b_ = determine_black(planes.b_, noise);
    LOG("black is: "<<black.r_<<" "<<black.g1_<<" "<<black.g2_<<" "<<black.b_);

    /** adjust noise for black levels. **/
    int min_black = std::min(black.r_, black.g1_);
    min_black = std::min(min_black, black.g2_);
    min_black = std::min(min_black, black.b_);
    noise -= min_black;
    LOG("noise is: "<<noise);
}

void crop_black(
    Planes &planes
) {
    LOG("cropping black pixels...");
    /**
    the left 37 columns are black.
    the top 16 rows are black.
    row 17 is garbage.
    the rest of the pixels are the actual image.
    **/
    /** left, top, right, bottom **/
    planes.crop(38, 18, planes.r_.width_, planes.r_.height_);
    LOG("cropped width ="<<planes.r_.width_);
    LOG("cropped height="<<planes.r_.height_);
}

void fix_bad_pixels(
    Planes &planes
) {
    /** hard coded list of bad pixels. **/
    fix_bad_pixel(planes.r_, 4011, 2495);
}

void compute_luminance(
    Planes &planes,
    Plane &luminance
) {
    /** convert canon rggb to srgb by multiplying by this matrix. **/
    static double mat[3][3] = {
        {+1.901824, -0.972035, +0.070211},
        {-0.229410, +1.659384, -0.429974},
        {+0.042001, -0.519143, +1.477141}
    };

    /** convert srgb to luminance by multiplying by this vector. **/
    static double lum_vec[3] = { 0.2125, 0.7154, 0.0721};

    /** combine them **/
    double lum_rx = mat[0][0]*lum_vec[0] + mat[0][1]*lum_vec[1] + mat[0][2]*lum_vec[2];
    double lum_gx = mat[1][0]*lum_vec[0] + mat[1][1]*lum_vec[1] + mat[1][2]*lum_vec[2];
    double lum_bx = mat[2][0]*lum_vec[0] + mat[2][1]*lum_vec[1] + mat[2][2]*lum_vec[2];

    /** we still have two greens. **/
    lum_gx /= 2.0;

    int wd = planes.r_.width_;
    int ht = planes.r_.height_;
    luminance.init(wd, ht);
    int sz = wd * ht;
    for (int i = 0; i < sz; ++i) {
        int r = planes.r_.samples_[i];
        int g1 = planes.g1_.samples_[i];
        int g2 = planes.g2_.samples_[i];
        int b = planes.b_.samples_[i];
        int lum = std::round(r*lum_rx + (g1 + g2)*lum_gx + b*lum_bx);
        luminance.samples_[i] = lum;
    }
}
