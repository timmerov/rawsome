/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
process images.
**/

#include "cmd_line.h"
#include "log.h"
#include "save_as_png.h"

namespace {
class Rawsome {
public:
    Rawsome() = default;
    ~Rawsome() = default;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        /** define the command line options. **/
        const char *options_short = "?hp";
        CmdLineOptions::LongFormat options_long[] = {
            {'?', "help"},
            {'h', "help"},
            {'p', "--save-as-png"},
            {0, nullptr}
        };
        CmdLineOptions clo(argc, argv, options_short, options_long);

        /** the first option is required. **/
        bool success = clo.get();
        if (success == false) {
            show_usage();
            return 1;
        }

        /** dispatch to the command. **/
        int exit_code = 0;
        switch (clo.option_) {
        default:
        case '?':
        case 'h':
            show_usage();
            break;
        case 'p':
            exit_code = save_as_png(argc-1, argv+1);
            break;
        }

        return exit_code;
    }

    void show_usage() {
        LOG("usage: rawsome command [...]");
        LOG("  -h, -?, --help    : show usage.");
        LOG("  -p, --save-as-png : save image file as png.");
        LOG("further arguments are command dependent.");
    }
};
} // anonymouse namespace

int main(
    int argc,
    clo_argv_t argv
) noexcept {
    rs_log::init("rawsome.log");

    Rawsome rawsome;
    int exit_code = rawsome.run(argc, argv);

    return exit_code;
}
