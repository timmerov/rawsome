/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
register two source images.
**/

#include "canon.h"
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
    Plane luminance1_;
    Plane luminance2_;

    int run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;

        LOG("register two images - work in progress.");

        image1_.load_raw(kInputFile1);
        if (image1_.is_loaded_ == false) {
            return 1;
        }
        image2_.load_raw(kInputFile2);
        if (image2_.is_loaded_ == false) {
            return 1;
        }

        compute_luminance(image1_.planes_, luminance1_);
        compute_luminance(image2_.planes_, luminance2_);

        LOG("success.");
        return 0;
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
