/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

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
    double black_ = -1.0;
    double white_ = -1.0;
    double auto_black_ = -1.0;
    double auto_white_ = -1.0;
    double color_enhancement_ = 0.0;
    double gamma_ = 1.0;

    bool parse(int argc, clo_argv_t argv);

    static void print_usage();
};
