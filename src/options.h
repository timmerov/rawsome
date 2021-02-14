/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
rawsome usage and command line options.
**/

#include "cmd_line.h"

#include <string>

class Options {
public:
    Options() = default;
    ~Options() = default;

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

    bool parse(int argc, clo_argv_t argv);

    static void print_usage();
};
