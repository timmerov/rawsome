/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
run the experiment of the day.
**/

#include "experiment.h"
#include "log.h"
#include "rs_png.h"

#include <vector>

namespace {

class Experiment {
public:
    Experiment() = default;
    ~Experiment() = default;

    void run(
        int argc,
        clo_argv_t argv
    ) {
        (void) argc;
        (void) argv;
        LOG("experiment of the day.");
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
