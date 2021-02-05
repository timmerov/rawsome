/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
read and write png files.
wrapper for libpng.
**/

#include <png.h>

class Png {
public:
    Png();
    Png(const Png &) = delete;
    ~Png();

    void destruct();

    int wd_;
    int ht_;
    int stride_;
    png_byte *data_;

    bool read(const char *filename);

    // prep for writing.
    void init(int width, int height);

    bool write(const char *filename);
};
