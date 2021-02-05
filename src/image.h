/*
Copyright (C) 2012-2020 tim cotter. All rights reserved.
*/

/**
data structures for image processing.
**/

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

    /** scale low to 0 and multiply by the factor. **/
    void scale(int low, double factor);

    /** crop to rectangle **/
    void crop(int left, int top, int right, int bottom);

    /** transpose the image. **/
    void transpose();

    /** interpolate intermediate horizontal pixels. **/
    void interpolate_horz_1331();
    void interpolate_horz_1331(int sat);

    /** downsample by 2x. **/
    void downsample(Plane &src);

    /** apply 121 gaussian N times. **/
    void gaussian(int n);
    void gaussian_horz(int n);
};

class Image {
public:
    Image();
    ~Image();

    Plane r_;
    Plane g1_;
    Plane g2_;
    Plane b_;

    /** initialize all of the planes. **/
    void init(int wd, int ht);

    /** crop to rectangle **/
    void crop(int left, int top, int right, int bottom);

    /** transpose the image. **/
    void transpose();

    /** interpolate intermediate horizontal pixels. **/
    void interpolate_horz_1331();
};
