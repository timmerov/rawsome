/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
read and write png files.
wrapper for libpng.
**/

#pragma GCC diagnostic ignored "-Wclobbered"

#include "log.h"
#include "rs_png.h"

#include <png.h>

#include <fstream>


namespace {
const int kPngSigSize = 8;

static void readData(
    png_structp png,
    png_bytep data,
    png_size_t len
) {
    auto f = (std::fstream *) png_get_io_ptr(png);
    f->read((char *) data, len);
}

static void writeData(
    png_structp png,
    png_bytep data,
    png_size_t len
) {
    auto f = (std::fstream *) png_get_io_ptr(png);
    f->write((char *) data, len);
}
} // anonymous namespace

Png::Png() :
    wd_(0),
    ht_(0),
    stride_(0),
    data_(nullptr)
{
}

Png::~Png() {
    destruct();
}

void Png::destruct() {
    delete[] data_;
    data_ = nullptr;
}

bool Png::read(
    const char *filename
) {
    bool result = true;
    std::fstream f;
    png_byte sig[kPngSigSize];
    png_structp png = nullptr;
    png_infop info = nullptr;
    int depth = 0;
    int channels = 0;
    int color_type = 0;
    int size = 0;
    png_byte **rows = nullptr;

    destruct();
    wd_ = 0;
    ht_ = 0;
    stride_ = 0;

    /*if (result)*/ {
        f.open(filename, std::fstream::in);
        if (f.is_open() == false) {
            LOG("Failed to open file \"" << filename << "\"");
            result = false;
        }
    } if (result) {
        f.read((char *) sig, sizeof(sig));
        if (f.good() == false) {
            LOG("Failed to read file signature.");
            result = false;
        }
        //agm::log::bytes(sig, kPngSigSize);
    } if (result) {
        if (png_sig_cmp(sig, 0, kPngSigSize) != 0) {
            LOG("File signature not PNG.");
            result = false;
        }
    } if (result) {
        png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (png == nullptr) {
            LOG("Failed to create PNG read struct.");
                result = false;
        }
    } if (result) {
        info = png_create_info_struct(png);
        if (info == nullptr) {
            LOG("Failed to create PNG info struct.");
            result = false;
        }
    } if (result) {
        // in case something goes wrong in the subsequent reads.
        if (setjmp(png_jmpbuf(png))) {
            LOG("Failed reading PNG file.");
            result = false;
        }
    } if (result) {
        // read from the file.
        png_set_read_fn(png, (png_voidp) &f, readData);

        png_set_sig_bytes(png, sizeof(sig));
        png_read_info(png, info);

        wd_ = png_get_image_width(png, info);
        ht_ = png_get_image_height(png, info);
        depth = png_get_bit_depth(png, info);
        channels = png_get_channels(png, info);
        color_type = png_get_color_type(png, info);

        // tell png lib we want 8 bits of r g b.
        switch (color_type) {
        case PNG_COLOR_TYPE_PALETTE:
            png_set_palette_to_rgb(png);
            channels = 3;
            break;

        case PNG_COLOR_TYPE_GRAY:
            if (depth < 8) {
                png_set_expand_gray_1_2_4_to_8(png);
                depth = 8;
            }
            break;
        }
        if (depth == 16) {
            png_set_strip_16(png);
        }

        // in case we made any changes.
        png_read_update_info(png, info);

        stride_ = wd_ * depth * channels / 8;
        size = stride_ * ht_;

        rows = new(std::nothrow) png_byte*[ht_];
        data_ = new(std::nothrow) png_byte[size];

        auto row = data_;
        for (int i = 0; i < ht_; ++i, row += stride_) {
            rows[i] = row;
        }

        png_read_image(png, rows);
    }

    LOG("resolution=" << wd_ << "x" << ht_);
    LOG("channels=" << channels);
    LOG("bits/channel=" << depth);
    LOG("stride=" << stride_);

    /*if (rows && data_) {
        agm::log::bytes(rows[ht_/2], std::min(stride_,256));
    }*/

    // clean up
    delete[] rows;
    if (png) {
        png_destroy_read_struct(&png, &info, nullptr);
    }

    return result;
}

void Png::init(
    int width,
    int height
) {
    destruct();
    wd_ = width;
    ht_ = height;
    stride_ = (3 * wd_ + 31) / 32 * 32;
    int size = stride_ * ht_;
    data_ = new(std::nothrow) png_byte[size];
}

bool Png::write(
    const char *filename
) {
    bool result = true;
    std::fstream f;
    //png_byte sig[kPngSigSize];
    png_structp png = nullptr;
    png_infop info = nullptr;

    LOG("resolution=" << wd_ << "x" << ht_);
    LOG("stride=" << stride_);

    /*if (result)*/ {
        f.open(filename, std::fstream::out);
        if (f.is_open() == false) {
            LOG("Failed to open file \"" << filename << "\"");
            result = false;
        }
    } if (result) {
        png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (png == nullptr) {
            LOG("Failed to create PNG write struct.");
            result = false;
        }
    } if (result) {
        info = png_create_info_struct(png);
        if (info == nullptr) {
            LOG("Failed to create PNG info struct.");
            result = false;
        }
    } if (result) {
        if (setjmp(png_jmpbuf(png))) {
            LOG("Failed writing PNG file.");
            result = false;
        }
    } if (result) {
        png_set_write_fn(png, (png_voidp) &f, writeData, nullptr);

        png_set_IHDR(png, info, wd_, ht_, 8,
            PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png, info);

        auto row = data_;
        for (int i = 0; i < ht_; ++i, row += stride_) {
            png_write_row(png, row);
        }
        png_write_end(png, nullptr);
    }

    // clean up
    if (png) {
        png_destroy_write_struct(&png, &info);
    }

    return result;
}
