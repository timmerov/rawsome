/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
rawsome usage and command line options.
**/

/**
the moon is about 724 pixels wide.
it's area is about 400k pixels.
the image is 6000*4000
so the moon is about 1.7% of the pixels.
count pixels from brightest to get to 1%.
above that threshold we're definitely in the moon.
90% of the pixels are black sky pixels.
count pixels from the darkest to get to 90%.
below that threshold we're definitely in the black.
average the two thresholds for moon/not-moon pixels.

or...
use joint index.
where X% of the pixels consume 1-X% of the range.
**/

#include "cmd_line.h"
#include "options.h"
#include "log.h"

#include <string>
#include <sstream>

namespace {
/** we may want this again the future. **/
#if 0
bool comma_separated_doubles_2(
    const char *s,
    double &x,
    double &y
) {
    /**
    wtf can c++ not just do this:
    stringstream ss(s);
    ss>>x>>",">>y;
    return ss.good();
    ?
    sigh.
    **/
    std::stringstream ss(s);
    ss >> x;
    if (ss.fail()) {
        return false;
    }
    char ch = 0;
    ss >> ch;
    if (ss.fail()) {
        return false;
    }
    if (ch != ',') {
        return false;
    }
    ss >> y;
    if (ss.fail()) {
        return false;
    }
    return true;
}
#endif

bool comma_separated_ints_4(
    const char *s,
    int &a,
    int &b,
    int &c,
    int &d
) {
    /**
    wtf can c++ not just do this:
    stringstream ss(s);
    ss>>a>>",">>b>>",">>c>>",">>d;
    return ss.good();
    ?
    sigh.
    **/
    std::stringstream ss(s);
    ss >> a;
    if (ss.fail()) {
        return false;
    }
    char ch = 0;
    ss >> ch;
    if (ss.fail()) {
        return false;
    }
    if (ch != ',') {
        return false;
    }
    ss >> b;
    if (ss.fail()) {
        return false;
    }
    ss >> ch;
    if (ss.fail()) {
        return false;
    }
    if (ch != ',') {
        return false;
    }
    ss >> c;
    if (ss.fail()) {
        return false;
    }
    ss >> ch;
    if (ss.fail()) {
        return false;
    }
    if (ch != ',') {
        return false;
    }
    ss >> d;
    if (ss.fail()) {
        return false;
    }
    return true;
}

void set_out_filename(
    Options &opt
) {
    /**
    create an output filename from the input filename.
    unlesse the user supplied one.
    **/
    if (opt.out_filename_.empty() == false) {
        return;
    }

    /** split the path from the filename plus extension. **/
    std::string path;
    std::string fn_ext;
    auto len = opt.in_filename_.size();
    auto found = opt.in_filename_.find_last_of("/");
    if (found >= len) {
        path = "";
        fn_ext = opt.in_filename_;
    } else {
        path = opt.in_filename_.substr(0, found);
        fn_ext = opt.in_filename_.substr(found + 1);
    }

    /** split the filename from the extension. **/
    std::string fn;
    std::string ext;
    len = fn_ext.size();
    found = fn_ext.find_last_of(".");
    if (found >= len) {
        fn = fn_ext;
        ext = "";
    } else {
        fn = fn_ext.substr(0, found);
        ext = fn_ext.substr(found + 1);
    }

    /** replace the extension in the output filename. **/
    opt.out_filename_ = path;
    if (opt.out_filename_.empty() == false) {
        opt.out_filename_ += "/";
    }
    opt.out_filename_ += fn;
    opt.out_filename_ += ".png";
}
} // anonymous namespace

bool Options::parse(
    int argc,
    clo_argv_t argv
) {
    LOG("parsing command line options...");

    /** set default values for all options. **/
    in_filename_.clear();
    out_filename_.clear();
    halfsize_ = false;
    noise_ = 0.0;
    drama_ = 0.0;
    window_ = 0;
    auto_brightness_ = -1.0;
    linear_brightness_ = 0.0;
    color_enhancement_ = 0.0;
    deblur_left_ = 0;
    deblur_top_ = 0;
    deblur_right_ = 0;
    deblur_bottom_ = 0;

    const char *options_short = "?i:o:n:d:W:ha:b:c:D:";
    CmdLineOptions::LongFormat options_long[] = {
        {'?', "help"},
        {'i', "input"},
        {'o', "output"},
        {'n', "noise"},
        {'d', "dynamic"},
        {'W', "window"},
        {'h', "halfsize"},
        {'a', "auto-brightness"},
        {'b', "linear-brightness"},
        {'c', "color-enhancement"},
        {'D', "deblur"},
        {0, nullptr}
    };
    CmdLineOptions clo(argc, argv, options_short, options_long);
    for(;;) {
        bool success = clo.get();
        if (success == false) {
            if (clo.error_) {
                return false;
            }
            break;
        }

        switch (clo.option_) {
        case '?':
            /** return a parsing error so usage is shown. **/
            return false;
        case 'i':
            in_filename_ = clo.value_;
            break;
        case 'o':
            out_filename_ = clo.value_;
            break;
        case 'n':
            noise_ = std::atof(clo.value_);
            break;
        case 'd':
            drama_ = std::atof(clo.value_);
            break;
        case 'W':
            window_ = std::atoi(clo.value_);
            break;
        case 'h':
            halfsize_ = true;
            break;
        case 'a':
            auto_brightness_ = std::atof(clo.value_);
            break;
        case 'b':
            linear_brightness_ = std::atof(clo.value_);
            break;
        case 'c':
            color_enhancement_ = std::atof(clo.value_);
            break;
        case 'D': {
            bool good = comma_separated_ints_4(clo.value_, deblur_left_, deblur_top_, deblur_right_, deblur_bottom_);
            if (good == false) {
                return false;
            }
            break;
        }
        }
    }

    if (in_filename_.empty()) {
        return false;
    }
    set_out_filename(*this);

    return true;
}


void Options::print_usage() {
    LOG("usage:");
    LOG("rawsome [options]...");
    LOG("  -? --help        : prints usage.");
    LOG("  -i --input file  : input filename - required.");
    LOG("  -o --output file : output filename.");
    LOG("     derived from input filename if not specified.");
    LOG("     the extension is replaced with \".png\".");
    LOG("");
    LOG("these operations are applied in this order:");
    LOG("");
    LOG("  -n --noise noise    : override noise floor when expanding dynamic range.");
    LOG("  -d --dynamic dynamic: sets the dynamic range expansion factor. default 1.5");
    LOG("  -W --window window  : set the dynamic range gaussian window. default 32");
    LOG("");
    LOG("  -h --halfsize       : disables demosaicing.");
    LOG("     the bayer block of RGGB is treated as a single pixel.");
    LOG("");
    LOG("  -a --auto-brightness fraction:");
    LOG("     this fraction of the brightest pixels will be forced to full white.");
    LOG("     set to -1 to disable. default 0.");
    LOG("  -b --linear-brightness brightness:");
    LOG("     uniformly scale pixel values. default 1.0.");
    LOG("");
    LOG("  -c --color-enhancement factor:");
    LOG("     convert to yuv. scale uv by factor. convert back to rgb.");
    LOG("     default 0.0.");
    LOG("");
    LOG("  -D --deblur t,l,b,r : deblur image. default 0,0,0,0");
    LOG("     derive the blur kernel from the patch.");
}
