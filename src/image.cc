/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
data structures for image processing.
**/

#include "canon.h"
#include "common.h"
#include "dump.h"
#include "image.h"
#include "log.h"
#include "planes.h"
#include "rs_png.h"

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

void write_plane(
    std::ofstream &out,
    const char *which,
    Plane &plane
) {
    int sz = plane.width_ * plane.height_ * sizeof(int);
    auto data = (const char *) plane.samples_.data();

    out<<which<<":"<<std::endl;
    out.write(data, sz);
    out<<std::endl;
}

void read_plane(
    std::ifstream &in,
    int wd,
    int ht,
    Plane &plane
) {
    plane.init(wd, ht);

    int sz = wd * ht * sizeof(int);
    auto data = (char *) plane.samples_.data();

    in.read(data, sz);

    std::string line;
    std::getline(in, line);
}

std::string trim(
    const std::string &s
) {
    auto len = s.size();
    auto first = s.find_first_not_of(" \t");
    if (first >= len) {
        /** it's all whitespace. **/
        return "";
    }

    auto last = s.find_last_not_of(" \t");
    len = last - first + 1;

    return s.substr(first, len);
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

void Image::save_png(
    const char *filename
) {
    LOG("saving to png: \""<<filename<<"\"...");

    int wd = planes_.r_.width_;
    int ht = planes_.r_.height_;
    Png png;
    png.init(wd, ht);
    for (int y = 0; y < ht; ++y) {
        for (int x = 0; x < wd; ++x) {
            int r = planes_.r_.get(x, y);
            int g = planes_.g1_.get(x, y);
            int b = planes_.b_.get(x, y);
            int idx = 3*x + y*png.stride_;
            png.data_[idx] = pin_to_8bits(r);
            png.data_[idx+1] = pin_to_8bits(g);
            png.data_[idx+2] = pin_to_8bits(b);
        }
    }
    png.write(filename);
}

void Image::load_rawsome(
    const char *filename
) {
    LOG("loading rawsome from: \""<<filename<<"\"...");
    is_loaded_ = false;
    planes_ = Planes();

    std::ifstream in(filename);
    if (in.is_open() == false) {
        return;
    }

    int wd = 0;
    int ht = 0;

    for(;;) {
        /** get data line by line. **/
        std::string line;
        std::getline(in, line);
        if (line.empty()) {
            break;
        }

        /** lines start with a command followed by a colon. **/
        auto len = line.size();
        auto found = line.find(":");
        if (found >= len) {
            /** a line with no colon? hrm... **/
            continue;
        }
        std::string cmd = line.substr(0, found);
        std::string arg = line.substr(found+1);

        /** commands and args may have leading and trailing whitespace. **/
        cmd = trim(cmd);
        arg = trim(arg);

        /** dispatch the command. **/
        if (cmd == "width") {
            wd = std::atoi(arg.c_str());
            LOG("width: "<<wd);
        } else if (cmd == "height") {
            ht = std::atoi(arg.c_str());
            LOG("height: "<<ht);
        } else if (cmd == "r") {
            read_plane(in, wd, ht, planes_.r_);
        } else if (cmd == "g1") {
            read_plane(in, wd, ht, planes_.g1_);
        } else if (cmd == "g2") {
            read_plane(in, wd, ht, planes_.g2_);
        } else if (cmd == "b") {
            read_plane(in, wd, ht, planes_.b_);
        }
    }

    is_loaded_ = true;
}

void Image::save_rawsome(
    const char *filename
) {
    LOG("saving to rawsome: \""<<filename<<"\"...");

    std::ofstream out(filename);
    if (out.is_open() == false) {
        return;
    }

    int wd = planes_.r_.width_;
    int ht = planes_.r_.height_;

    out<<"rawsome:"<<std::endl;
    out<<"version: 0"<<std::endl;
    out<<"width: "<<wd<<std::endl;
    out<<"height: "<<ht<<std::endl;
    out<<"format: rggb"<<std::endl;
    write_plane(out, "r", planes_.r_);
    write_plane(out, "g1", planes_.g1_);
    write_plane(out, "g2", planes_.g2_);
    write_plane(out, "b", planes_.b_);
    out<<"end:"<<std::endl;
}
