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

void rotate_90(
    Plane &src
) {
    int wd = src.width_;
    int ht = src.height_;
    Plane dst;
    dst.init(ht, wd);
    for (int y = 0; y < ht; ++y) {
        for (int x = 0; x < wd; ++x) {
            int p = src.get(x, y);
            dst.set(ht - y - 1, x, p);
        }
    }
    src = std::move(dst);
}

void rotate_180(
    Plane &src
) {
    int wd = src.width_;
    int ht = src.height_;
    Plane dst;
    dst.init(wd, ht);
    for (int y = 0; y < ht; ++y) {
        for (int x = 0; x < wd; ++x) {
            int p = src.get(x, y);
            dst.set(wd - x - 1, ht - y - 1, p);
        }
    }
    src = std::move(dst);
}

void rotate_270(
    Plane &src
) {
    int wd = src.width_;
    int ht = src.height_;
    Plane dst;
    dst.init(ht, wd);
    for (int y = 0; y < ht; ++y) {
        for (int x = 0; x < wd; ++x) {
            int p = src.get(x, y);
            dst.set(y, wd - x - 1, p);
        }
    }
    src = std::move(dst);
}

void rotate(
    Plane &plane,
    int rotation
) {
    switch (rotation) {
        default:
            break;
        case 90:
            rotate_90(plane);
            break;
        case 180:
            rotate_180(plane);
            break;
        case 270:
            rotate_270(plane);
            break;
    }
}

void rotate(
    Planes &planes,
    int rotation
) {
    rotate(planes.r_, rotation);
    rotate(planes.g1_, rotation);
    rotate(planes.g2_, rotation);
    rotate(planes.b_, rotation);
}

void fix_bad_pixels(
    Planes &planes,
    int noise
) {
    LOG("fixing bad pixels...");

    /**
    sometimes pixels go bloop.
    they randomly dump charge into the sensor.
    no idea why.
    this often shows up with negative values when we convert to srgb.
    we focus on fixing those bad pixels.
    **/

    /** white balance **/
    double balance_r = 2.27051;
    double balance_b = 1.54199;

    /** srgb conversion **/
    double mat[3][3] = {
        {+1.901824, -0.972035, +0.070211},
        {-0.229410, +1.659384, -0.429974},
        {+0.042001, -0.519143, +1.477141}
    };

    /** combine white balance and srgb **/
    for (int i = 0; i < 3; ++i) {
        mat[i][0] *= balance_r;
        mat[i][2] *= balance_b;
    }

    /** small negative values are allowed. **/
    double threshold = -3.0 * noise;

    /** stats **/
    int nfixed = 0;
    int nnotfixed = 0;

    /**
    assume we're going to crop the image.
    so garbage pixels on the border don't matter.
    **/
    int wd = planes.r_.width_;
    int ht = planes.r_.height_;
    for (int y = 1; y < ht-1; ++y) {
        for (int x = 1; x < wd-1; ++x) {
            /** get the source pixel **/
            double in_r = planes.r_.get(x, y);
            double in_g1 = planes.g1_.get(x, y);
            double in_g2 = planes.g2_.get(x, y);
            double in_b = planes.b_.get(x, y);

            /** pin to 0 **/
            in_r = std::max(0.0, in_r);
            in_g1 = std::max(0.0, in_g1);
            in_g2 = std::max(0.0, in_g2);
            in_b = std::max(0.0, in_b);

            /** compute srgb from g1 **/
            double out_r1 = mat[0][0]*in_r + mat[0][1]*in_g1 + mat[0][2]*in_b;
            double out_g1 = mat[1][0]*in_r + mat[1][1]*in_g1 + mat[1][2]*in_b;
            double out_b1 = mat[2][0]*in_r + mat[2][1]*in_g1 + mat[2][2]*in_b;
            bool bad1 = (out_r1 < threshold || out_g1 < threshold || out_b1 < threshold);

            /** compute srgb from g2 **/
            double out_r2 = mat[0][0]*in_r + mat[0][1]*in_g2 + mat[0][2]*in_b;
            double out_g2 = mat[1][0]*in_r + mat[1][1]*in_g2 + mat[1][2]*in_b;
            double out_b2 = mat[2][0]*in_r + mat[2][1]*in_g2 + mat[2][2]*in_b;
            bool bad2 = (out_r2 < threshold || out_g2 < threshold || out_b2 < threshold);

            if (bad1 && bad2 == false) {
                /** fix by replacing a bad green with a good green. **/
                planes.g1_.set(x, y, (int) in_g2);
                bad1 = false;
                ++nfixed;
            }
            else if (bad2 && bad1 == false) {
                /** fix by replacing a bad green with a good green. **/
                planes.g2_.set(x, y, (int) in_g1);
                bad2 = false;
                ++nfixed;
            } else if (in_r >= in_b && in_r >= in_g1 && in_r >= in_g2) {
                /** fix by replacing a bad red with the average of its neighbors. **/
                int r0 = planes.r_.get(x, y-1);
                int r1 = planes.r_.get(x-1, y);
                int r2 = planes.r_.get(x+1, y);
                int r3 = planes.r_.get(x, y+1);
                int r = (r0 + r1 + r2 + r3 + 2) / 4;
                planes.r_.set(x, y, r);
                bad1 = false;
                bad2 = false;
                ++nfixed;
            } else if (in_b >= in_r && in_b >= in_g1 && in_b >= in_g2) {
                /** fix by replacing a bad blue with the average of its neighbors. **/
                int b0 = planes.b_.get(x, y-1);
                int b1 = planes.b_.get(x-1, y);
                int b2 = planes.b_.get(x+1, y);
                int b3 = planes.b_.get(x, y+1);
                int b = (b0 + b1 + b2 + b3 + 2) / 4;
                planes.b_.set(x, y, b);
                bad1 = false;
                bad2 = false;
                ++nfixed;
            } else {
                /** unfixable pixel. **/
            }

            if (bad1 || bad2) {
                /*LOG("suspect pixel at "<<x<<","<<y<<" in: "<<in_r<<","<<in_g1<<","<<in_g2<<","<<in_b
                    <<" out1: "<<out_r1<<","<<out_g1<<","<<out_b1
                    <<" out2: "<<out_r2<<","<<out_g2<<","<<out_b2);*/
                ++nnotfixed;
            }
        }
    }
    LOG("fixed "<<nfixed<<" suspect pixels.");
    LOG("did not fix "<<nnotfixed<<" suspect pixels.");
}

} // anonymous namespace

void CameraParams::print() {
    auto tm = std::localtime(&timestamp_);

    LOG("make        : "<<make_);
    LOG("model       : "<<model_);
    LOG("lens        : "<<lens_);
    LOG("gamma       : "<<gamma0_<<" "<<gamma1_);
    LOG("iso         : "<<iso_);
    LOG("shutter     : "<<shutter_<<"s");
    LOG("aperture    : f/"<<aperture_);
    LOG("focal_len   : "<<focal_length_);
    LOG("timestamp   : "<<std::put_time(tm, "%c %Z"));
    LOG("balance     : "<<balance_[0]<<" "<<balance_[1]<<" "<<balance_[2]<<" "<<balance_[3]);
}

void Image::load_raw(
    const char *filename
) {
    LOG("loading raw image from: \""<<filename<<"\"...");
    is_loaded_ = false;

    /**
    use LibRaw to load the raw image.
    **/
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
    copy_raw_to_planes(raw_image, planes_);

    /**
    extract useful information from libraw.
    **/
    double white_balance_r = raw_image.imgdata.color.cam_mul[0];
    double white_balance_g = raw_image.imgdata.color.cam_mul[1];
    double white_balance_b = raw_image.imgdata.color.cam_mul[2];
    double gamma0 = raw_image.imgdata.params.gamm[0];
    double gamma1 = raw_image.imgdata.params.gamm[1];
    gamma0 = 1.0 / gamma0;

    /**
    from libraw documentation:
    int flip; Image orientation (
        0 if does not require rotation;
        3 if requires 180-deg rotation;
        5 if 90 deg counterclockwise,
        6 if 90 deg clockwise).
    **/

    int rotation = 0;
    switch (raw_image.imgdata.sizes.flip) {
        default:
            rotation = 0;
            break;
        case 3:
            rotation = 180;
            break;
        case 5:
            rotation = 270;
            break;
        case 6:
            rotation = 90;
            break;
    }

    camera_.make_ = raw_image.imgdata.idata.make;
    camera_.model_ = raw_image.imgdata.idata.model;
    camera_.lens_ = raw_image.imgdata.lens.Lens;
    camera_.make_ = raw_image.imgdata.idata.make;
    camera_.gamma0_ = gamma0;
    camera_.gamma1_ = gamma1;
    camera_.iso_ = raw_image.imgdata.other.iso_speed;
    camera_.shutter_ = raw_image.imgdata.other.shutter;
    camera_.aperture_ = raw_image.imgdata.other.aperture;
    camera_.focal_length_ = raw_image.imgdata.other.focal_len;
    camera_.timestamp_ = raw_image.imgdata.other.timestamp;
    for (int i = 0; i < 4; ++i) {
        camera_.balance_[i] = raw_image.imgdata.other.analogbalance[i];
    }

    LOG("flip="<<raw_image.imgdata.sizes.flip);

    /** done with libraw **/
    raw_image.recycle();

    /**
    the camera has quite a few black pixels.
    why? dunno.
    use them to determine the black level.
    subtract this from all pixels.
    while we're here, also determine noise.
    and remove the black pixels from the image.
    **/
    RggbPixel black;
    determine_black(planes_, black, noise_);
    crop_black(planes_);
    planes_.subtract(black);

    /**
    remove bad pixels from the image.
    camera specific.
    **/
    fix_bad_pixels(planes_, noise_);

    /** white balance **/
    double wb_r = white_balance_r / white_balance_g;
    double wb_b = white_balance_b / white_balance_g;
    LOG("white balance r:"<<wb_r<<" b:"<<wb_b);
    planes_.r_.multiply(wb_r);
    planes_.b_.multiply(wb_b);

    rotate(planes_, rotation);

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
