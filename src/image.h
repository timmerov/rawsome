/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
data structures for image processing.
**/

#include "planes.h"

#include <ctime>
#include <string>

class CameraParams {
public:
    CameraParams() = default;
    ~CameraParams() = default;

    std::string make_;
    std::string model_;
    std::string lens_;
    double gamma0_ = 0.0;
    double gamma1_ = 0.0;
    double wb_r_ = 0.0;
    double wb_b_ = 0.0;
    int iso_ = 0;
    double shutter_ = 0.0;
    double aperture_ = 0.0;
    double focal_length_ = 0.0;
    std::time_t timestamp_ = 0;
    double temperature_ = 0.0;

    void print();
};

class UserParams {
public:
    std::string in_filename_;
    std::string out_filename_;
    bool halfsize_ = false;
    int saturation_ = 0;
    double wb_r_ = 0.0;
    double wb_b_ = 0.0;
    double gamma0_ = 0.0;
    double gamma1_ = 0.0;
    int noise_ = 0;
    double drama_ = 0.0;
    int window_ = 0;
    double auto_brightness_ = -1.0;
    double linear_brightness_ = 0.0;
    double color_enhancement_ = 0.0;
};

class Image {
public:
    Image() = default;
    ~Image() = default;

    bool is_loaded_ = false;
    CameraParams camera_;
    Planes planes_;

    void load_raw(const char *filename);
};
