/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
run the experiment of the day.
**/

#include "experiment.h"
#include "image.h"
#include "log.h"

#include <vector>

namespace {

class Experiment {
public:
    Experiment() = default;
    ~Experiment() = default;

    Image image_;

    void run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;
        LOG("experiment of the day.");

        image_.load_raw("vulpeculae.CR2");
        if (image_.is_loaded_ == false) {
            return;
        }

        image_.save_png("vulpeculae.png");
    }
};

}

int experiment_of_the_day(
    int argc,
    clo_argv_t argv
) {
    Experiment x;
    x.run(argc, argv);

    return 0;
}
