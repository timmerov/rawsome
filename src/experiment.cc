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

        image_.load_raw("milkyway.CR2");
        if (image_.is_loaded_ == false) {
            return;
        }

        //crop_area_of_interest();
        to_rggb();

        image_.save_png("milkyway.png");
    }

    void crop_area_of_interest() {
        int x = 568;
        int y = 778;
        int sz = 150;
        int lf = x - sz/2;
        int tp = y - sz/2;
        lf &= ~1;
        tp &= ~1;
        sz &= ~1;
        int rt = lf + sz;
        int bt = tp + sz;
        image_.planes_.crop(lf, tp, rt, bt);
    }

    void to_rggb() {
        int wd = image_.planes_.r_.width_;
        int ht = image_.planes_.r_.height_;

        Planes rggb;
        rggb.init(2*wd, 2*ht);

        int sz = wd * ht * 4;
        for (int i = 0; i < sz; ++i) {
            rggb.r_.samples_[i] = 0;
            rggb.g1_.samples_[i] = 0;
            rggb.g2_.samples_[i] = 0;
            rggb.b_.samples_[i] = 0;
        }

        for (int y = 0; y < ht; ++y) {
            for (int x = 0; x < wd; ++x) {
                int r = image_.planes_.r_.get(x, y);
                int g1 = image_.planes_.g1_.get(x, y);
                int g2 = image_.planes_.g2_.get(x, y);
                int b = image_.planes_.b_.get(x, y);

                /**
                set the rggb bayer pattern.
                g2 is no longer used.
                **/
                //rggb.r_.set(2*x, 2*y, r);
                //rggb.g1_.set(2*x+1, 2*y, g1);
                //rggb.g1_.set(2*x, 2*y+1, g2);
                //rggb.b_.set(2*x+1, 2*y+1, b);
                rggb.g1_.set(2*x, 2*y, g2);
                rggb.r_.set(2*x+1, 2*y, b);
                rggb.b_.set(2*x, 2*y+1, r);
                rggb.g1_.set(2*x+1, 2*y+1, g1);
            }
        }

        image_.planes_ = std::move(rggb);
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
