/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
data structures for image processing.
**/

#include "common.h"

#include <vector>

class Plane {
public:
    Plane();
    ~Plane();

    int width_ = 0;
    int height_ = 0;
    std::vector<int> samples_;

    /** allocate space for the samples. **/
    void init(int wd, int ht);

    /** set and get the sample value at the location. **/
    void set(int x, int y, int value);
    int get(int x, int y);

    /** reduce every sample by delta. **/
    void subtract(int delta);

    /** multiply every sample by factor. **/
    void multiply(double factor);

    /** crop to rectangle **/
    void crop(int left, int top, int right, int bottom);

    /** transpose the image. **/
    void transpose();
    void transpose_mt();

    /** interpolate intermediate pixels. **/
    void interpolate_1331();
    void interpolate_1331_mt();

    /** interpolate intermediate horizontal pixels. **/
    void interpolate_horz_1331();
    void interpolate_horz_1331_mt();

    /** downsample by 2x. **/
    void downsample(Plane &src);

    /** apply 121 gaussian N times. **/
    void gaussian(int n);
    void gaussian_horz(int n);
    void gaussian_horz_mt(int n);

    /** apply gamma curve to compensate for non-linearity of displays/eyeballs. **/
    void apply_display_gamma(double pwr, double ts, int white);
    void apply_display_gamma(std::vector<int> &curve);

    /** scale to 8 bit values. **/
    void scale_to_8bits();
};

class Planes {
public:
    Planes();
    ~Planes();

    Plane r_;
    Plane g1_;
    Plane g2_;
    Plane b_;

    /** initialize all of the planes. **/
    void init(int wd, int ht);

    /** set the pixel **/
    void set(int x, int y, RggbPixel p);
    void set(int x, int y, int value);

    /** reduce every pixel by delta. **/
    void subtract(RggbPixel &delta);

    /** multiply every pixel by factor. **/
    void multiply(RggbDouble &factor);
    void multiply3(double factor);
    void multiply4(double factor);

    /** crop to rectangle **/
    void crop(int left, int top, int right, int bottom);

    /** transpose the image. **/
    void transpose();

    /** interpolate intermediate pixels. **/
    void interpolate_1331();

    /** apply user gamma curve. **/
    void apply_user_gamma(double pwr);

    /** apply gamma. **/
    void apply_display_gamma(double pwr, double ts, int white);

    /** scale to 8 bit values. **/
    void scale_to_8bits();
};
