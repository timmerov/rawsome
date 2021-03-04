/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
load and stack black images.
save in rawsome format.
**/

#include "black.h"
#include "log.h"

#include <sstream>

namespace {

const char kDefaultOutputFilename[] = "black.rsm";

class StackBlacks {
public:
    StackBlacks() = default;
    ~StackBlacks() = default;

    std::string in_filename_;
    std::string out_filename_;
    int first_ = 0;
    int last_ = 0;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        LOG("stack blacks - not yet implemented.");
        int result = set_options(argc, argv);
        if (result == false) {
            print_usage();
            return 1;
        }

        LOG("in : \""<<in_filename_<<"\"");
        LOG("first: "<<first_);
        LOG("last : "<<last_);
        LOG("out: \""<<out_filename_<<"\"");

        for (int i = first_; i <= last_; ++i) {
            std::string fn = get_filename(i);
            LOG("fn: "<<fn);
        }

        return 0;
    }

    int set_options(
        int argc,
        clo_argv_t argv
    ) {
        /** define the command line options. **/
        const char *options_short = "?hi:o:";
        CmdLineOptions::LongFormat options_long[] = {
            {'?', "help"},
            {'h', "help"},
            {'i', "--input"},
            {'o', "--output"},
            {0, nullptr}
        };
        CmdLineOptions clo(argc, argv, options_short, options_long);

        /** the first option is required. **/
        while (clo.get()) {
            switch (clo.option_) {
            default:
            case '?':
            case 'h':
                return false;
            case 'i':
                in_filename_ = clo.value_;
                break;
            case 'o':
                out_filename_ = clo.value_;
                break;
            }
        }
        if (clo.error_) {
            return false;
        }
        if (in_filename_.empty()) {
            return false;
        }
        bool result = set_first_last();
        if (result == false) {
            return false;
        }
        set_out_filename();

        return true;
    }

    void print_usage() {
        LOG("usage: rawsome --stack-blacks -i input[,first,last] [-o ouput]");
        LOG("  -h, -?, --help   : show usage.");
        LOG("  -i, --input file[,first,last] :");
        LOG("    example: ~/blacks/IMG_#.CR2,1101,1121");
        LOG("    the # is optional.");
        LOG("    if included then ,first,last is required.");
        LOG("    for i = first to last");
        LOG("    replace # with i in the input filename.");
        LOG("    load and stack the black image.");
        LOG("  -o, --ouput file : save output file in rawsome format.");
        LOG("    default: black.rsm");
        LOG("load and stack a number of images captures with the lens cap on.");
        LOG("this saved black image can then be subtracted from source images.");
    }

    bool set_first_last() {
        /** split the path from the filename plus extension. **/
        std::string path;
        std::string fn_tags;
        auto len = in_filename_.size();
        auto found = in_filename_.find_last_of("/");
        if (found >= len) {
            path = "";
            fn_tags = in_filename_;
        } else {
            path = in_filename_.substr(0, found);
            fn_tags = in_filename_.substr(found + 1);
        }

        /** find the optional # in the filename. **/
        len = fn_tags.size();
        found = fn_tags.find("#");
        if (found >= len) {
            return true;
        }

        /** strip and extract the first,last tags. **/
        found = fn_tags.find(",");
        if (found >= len) {
            return false;
        }
        std::string fn = fn_tags.substr(0, found);
        std::string tags = fn_tags.substr(found + 1);
        in_filename_ = path + fn;

        /** separate first and last **/
        len = tags.size();
        found = tags.find(",");
        if (found >= len) {
            return false;
        }
        std::string first = tags.substr(0, found);
        std::string last = tags.substr(found + 1);
        first_ = std::atoi(first.c_str());
        last_ = std::atoi(last.c_str());

        return true;
    }

    void set_out_filename() {
        /**
        create an output filename from the input filename.
        unlesse the user supplied one.
        **/
        if (out_filename_.empty() == false) {
            return;
        }

        /** split the path from the filename plus extension. **/
        std::string path;
        auto len = in_filename_.size();
        auto found = in_filename_.find_last_of("/");
        if (found < len) {
            path = in_filename_.substr(0, found+1);
        }

        /** add the default name to the path. **/
        out_filename_ = path;
        out_filename_ += kDefaultOutputFilename;
    }

    std::string get_filename(
        int i
    ) {
        auto len = in_filename_.size();
        auto found = in_filename_.find("#");
        if (found >= len) {
            return in_filename_;
        }

        std::string prefix = in_filename_.substr(0, found);
        std::string suffix = in_filename_.substr(found + 1);
        std::stringstream ss;
        ss<<prefix<<i<<suffix;
        return ss.str();
    }
};

}

int stack_blacks(
    int argc,
    clo_argv_t argv
) {
    StackBlacks sb;
    int exit_code = sb.run(argc, argv);

    return exit_code;
}
