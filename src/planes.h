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

    /** interpolate intermediate horizontal pixels. **/
    void interpolate_horz_1331();
    void interpolate_horz_1331(int sat);
    void interpolate_horz_1331_mt();
    void interpolate_horz_1331_mt(int sat);

    /** downsample by 2x. **/
    void downsample(Plane &src);

    /** apply 121 gaussian N times. **/
    void gaussian(int n);
    void gaussian_horz(int n);
    void gaussian_horz_mt(int n);
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

    /** reduce every pixel by delta. **/
    void subtract(RggbPixel &delta);

    /** multiply every pixel by factor. **/
    void multiply(RggbDouble &factor);

    /** crop to rectangle **/
    void crop(int left, int top, int right, int bottom);

    /** transpose the image. **/
    void transpose();

    /** interpolate intermediate horizontal pixels. **/
    void interpolate_horz_1331();
};