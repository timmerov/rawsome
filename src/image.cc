/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
data structures for image processing.
**/

#include "canon.h"
#include "dump.h"
#include "image.h"
#include "log.h"
#include "planes.h"

#include <libraw/libraw.h>

#include <iomanip>

namespace {

void copy_raw_to_planes(
    LibRaw &raw_image,
    Planes &planes
) {
    LOG("copying raw to image...");
    /**
    raw_image_.imgdata.rawdata.raw_image is the raw samples in bayer format.
    evem rows: G B G B G B ...
    odd  rows: R G R G R G ...
    we want this pattern:
    every row: R G G B R G G B ...

    i don't know why the bayer pattern is GB/RG.
    i was expecting it to be RG/GB.
    but experimentation shows otherwise.
    maybe it has to do with the odd number of rows?
    and the pattern is bottom justified?
    dunno.
    weird.
    **/
    int raw_wd = raw_image.imgdata.sizes.width;
    int raw_ht = raw_image.imgdata.sizes.height;
    int wd = raw_wd / 2;
    int ht = raw_ht / 2;
    planes.init(wd, ht);

    for (int y = 0; y < ht; ++y) {
        for (int x = 0; x < wd; ++x) {
            int raw_idx = 2*x + 2*y*raw_wd;

            /** extract rggb from the GB/RG bayer pattern **/
            int g2 = raw_image.imgdata.rawdata.raw_image[raw_idx];
            int b  = raw_image.imgdata.rawdata.raw_image[raw_idx+1];
            int r  = raw_image.imgdata.rawdata.raw_image[raw_idx+raw_wd];
            int g1 = raw_image.imgdata.rawdata.raw_image[raw_idx+raw_wd+1];

            planes.r_.set(x, y, r);
            planes.g1_.set(x, y, g1);
            planes.g2_.set(x, y, g2);
            planes.b_.set(x, y, b);
        }
    }
}

} // anonymous namespace

void CameraParams::print() {
    auto tm = std::localtime(&timestamp_);

    LOG("make         : "<<make_);
    LOG("model        : "<<model_);
    LOG("lens         : "<<lens_);
    LOG("gamma        : "<<gamma0_<<" "<<gamma1_);
    LOG("white balance: R="<<wb_r_<<" B="<<wb_b_);
    LOG("iso          : "<<iso_);
    LOG("shutter      : "<<shutter_<<"s");
    LOG("aperture     : f/"<<aperture_);
    LOG("focal_len    : "<<focal_length_);
    LOG("timestamp    : "<<std::put_time(tm, "%c %Z"));
    LOG("temperature  : "<<temperature_);
}

void Image::load_raw(
    const char *filename
) {
    LOG("loading raw image from: \""<<filename<<"\"...");
    is_loaded_ = false;

    LibRaw raw_image;
    raw_image.open_file(filename);
    int raw_wd = raw_image.imgdata.sizes.width;
    int raw_ht = raw_image.imgdata.sizes.height;
    LOG("raw_wd="<<raw_wd);
    LOG("raw_ht="<<raw_ht);
    if (raw_wd <= 0 || raw_ht <= 0) {
        LOG("File not opened: "<<filename);
        return;
    }
    raw_image.unpack();
    raw_image.raw2image();
    //dump(raw_image);

    double gamma0 = raw_image.imgdata.params.gamm[0];
    double gamma1 = raw_image.imgdata.params.gamm[1];
    gamma0 = 1.0 / gamma0;

    double white_balance_r = raw_image.imgdata.color.cam_mul[0];
    double white_balance_g = raw_image.imgdata.color.cam_mul[1];
    double white_balance_b = raw_image.imgdata.color.cam_mul[2];

    camera_.make_ = raw_image.imgdata.idata.make;
    camera_.model_ = raw_image.imgdata.idata.model;
    camera_.lens_ = raw_image.imgdata.lens.Lens;
    camera_.make_ = raw_image.imgdata.idata.make;
    camera_.gamma0_ = gamma0;
    camera_.gamma1_ = gamma1;
    camera_.wb_r_ = white_balance_r / white_balance_g;
    camera_.wb_b_ = white_balance_b / white_balance_g;
    camera_.iso_ = raw_image.imgdata.other.iso_speed;
    camera_.shutter_ = raw_image.imgdata.other.shutter;
    camera_.aperture_ = raw_image.imgdata.other.aperture;
    camera_.focal_length_ = raw_image.imgdata.other.focal_len;
    camera_.timestamp_ = raw_image.imgdata.other.timestamp;
    camera_.temperature_ = raw_image.imgdata.other.CameraTemperature;

    copy_raw_to_planes(raw_image, planes_);

    RggbPixel black;
    determine_black(planes_, black, noise_);
    planes_.subtract(black);

    raw_image.recycle();
    is_loaded_ = true;
}
