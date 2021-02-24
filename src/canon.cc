/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
functions specific to canon.
**/

#include "canon.h"
#include "log.h"

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
} // anonymous namespace

void determine_black(
    Image &image,
    RggbPixel &black,
    int &noise
) {
    LOG("determining black...");
    /**
    the top rows and left columns of pixels are black.
    also set the black noise level.
    **/
    noise = 0;
    black.r_ = determine_black(image.planes_.r_, noise);
    black.g1_ = determine_black(image.planes_.g1_, noise);
    black.g2_ = determine_black(image.planes_.g2_, noise);
    black.b_ = determine_black(image.planes_.b_, noise);
    LOG("black is: "<<black.r_<<" "<<black.g1_<<" "<<black.g2_<<" "<<black.b_);

    /** adjust noise for black levels. **/
    int min_black = std::min(black.r_, black.g1_);
    min_black = std::min(min_black, black.g2_);
    min_black = std::min(min_black, black.b_);
    noise -= min_black;
    LOG("noise is: "<<noise);
}
