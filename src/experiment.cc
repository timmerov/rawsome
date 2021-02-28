/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
run the experiment of the day.
**/

#include "common.h"
#include "experiment.h"
#include "image.h"
#include "log.h"

#include <sstream>

namespace {

const char *kInputFile = "/home/timmer/Pictures/2021-02-27/mars800/IMG_1140.CR2";
const char *kOutputFile = "/home/timmer/Pictures/2021-02-27/mars800/experiment.png";

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

        image_.load_raw(kInputFile);

        determine_saturation();

        image_.planes_.scale_to_8bits();
        image_.save_png(kOutputFile);
    }

    void determine_saturation() {
        determine_saturation(image_.planes_.r_);
        determine_saturation(image_.planes_.g1_);
        determine_saturation(image_.planes_.g2_);
        determine_saturation(image_.planes_.b_);
    }

    void determine_saturation(
        Plane &plane
    ) {
        std::vector<int> histogram;
        histogram.resize(65536, 0);

        int minc = 65536;
        int maxc = 0;

        int sz = plane.width_ * plane.height_;
        for (int i = 0; i < sz; ++i) {
            int c = plane.samples_[i];
            c = pin_to_16bits(c);
            ++histogram[c];
            minc = std::min(minc, c);
            maxc = std::max(maxc, c);
        }

        std::stringstream ss;
        ss<<"min="<<minc<<" [";
        for (int i = 0; i < 16; ++i) {
            ss<<" "<<histogram[minc+i];
        }
        ss<<" ...]";
        LOG(ss.str());

        ss.str("");
        ss.clear();
        ss<<"max="<<maxc<<" [...";
        for (int i = -15; i <= 0; ++i) {
            ss<<" "<<histogram[maxc+i];
        }
        ss<<" ]";
        LOG(ss.str());
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
