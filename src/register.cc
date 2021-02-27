/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
register two source images.
**/

#include "image.h"
#include "log.h"
#include "register.h"

namespace {

const char *kInputFile1 = "stars1.CR2";
const char *kInputFile2 = "stars2.CR2";

class RegisterTwoImages {
public:
    RegisterTwoImages() = default;
    ~RegisterTwoImages() = default;

    Image image1_;
    Image image2_;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("register two images - work in progress.");

        load_image(kInputFile1, image1_);
        if (image1_.is_loaded_ == false) {
            return 1;
        }
        load_image(kInputFile2, image2_);
        if (image2_.is_loaded_ == false) {
            return 1;
        }

        LOG("success.");
        return 0;
    }

    void load_image(
        const char *filename,
        Image &image
    ) {
        image.load_raw(filename);
    }
};
}

int register_two_images(
    int argc,
    clo_argv_t argv
) {

    RegisterTwoImages reg;
    int exit_code = reg.run(argc, argv);

    return exit_code;
}
