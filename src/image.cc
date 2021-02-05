/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
data structures for image processing.
**/

#include "image.h"

namespace {

}

Plane::Plane() {
}

Plane::~Plane() {
}

void Plane::init(
    int wd,
    int ht
) {
    width_ = wd;
    height_ = ht;
    int sz = wd * ht;
    samples_.resize(sz);
}

void Plane::set(
    int x,
    int y,
    int value
) {
    int idx = x + y * width_;
    samples_[idx] = value;
}

int Plane::get(
    int x,
    int y
) {
    int idx = x + y * width_;
    int value = samples_[idx];
    return value;
}

void Plane::scale(
    int low,
    double factor
) {
    int sz = width_ * height_;
    for (int i = 0; i < sz; ++i) {
        int x = samples_[i];
        x -= low;
        if (x < 0) {
            x = 0;
        }
        x = x * factor;
        samples_[i] = x;
    }
}

void Plane::crop(
    int left,
    int top,
    int right,
    int bottom
) {
    int new_wd = right - left;
    int new_ht = bottom - top;
    Plane dst;
    dst.init(new_wd, new_ht);
    for (int y = 0; y < new_ht; ++y) {
        for (int x = 0; x < new_wd; ++x) {
            int c = get(x + left, y + top);
            dst.set(x, y, c);
        }
    }
    *this = std::move(dst);
}

void Plane::transpose() {
    Plane dst;
    dst.init(height_, width_);
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            int c = get(x,y);
            dst.set(y, x, c);
        }
    }
    *this = std::move(dst);
}

void Plane::interpolate_horz_1331() {
    /**
    given N pixels: a b c d ... w x y z
    we are going to output: b m c m d ... w m x m y
    where n are newly interpolated pixels.
    we are going to discard the left-most and right-most pixels.
    we will have N-2 (b-y) + N-3 = 2N-5 pixels.
    **/
    int new_wd = 2*width_ - 5;
    Plane dst;
    dst.init(new_wd, height_);
    for (int y = 0; y < height_; ++y) {
        int c1 = get(0, y);
        int c2 = get(1, y);
        int c3 = get(2, y);
        int dstx = 0;
        for (int x = 3; x < width_; ++x) {
            int c0 = c1;
            c1 = c2;
            c2 = c3;
            c3 = get(x, y);
            int new_c = (c0 + 3*c1 + 3*c2 + c3) / 8;
            dst.set(dstx, y, c1);
            dst.set(dstx+1, y, new_c);
            dstx += 2;
        }
        dst.set(dstx, y, c2);
    }
    *this = std::move(dst);
}

void Plane::interpolate_horz_1331(
    int sat
) {
    /**
    given N pixels: a b c d ... w x y z
    we are going to output: b m c m d ... w m x m y
    where n are newly interpolated pixels.
    we are going to discard the left-most and right-most pixels.
    we will have N-2 (b-y) + N-3 = 2N-5 pixels.
    **/
    int new_wd = 2*width_ - 5;
    Plane dst;
    dst.init(new_wd, height_);
    for (int y = 0; y < height_; ++y) {
        int c1 = get(0, y);
        int c2 = get(1, y);
        int c3 = get(2, y);
        int bits = 0;
        if (c1 >= sat) {
            bits |= 2;
        }
        if (c2 >= sat) {
            bits |= 4;
        }
        if (c3 >= sat) {
            bits |= 8;
        }
        int dstx = 0;
        for (int x = 3; x < width_; ++x) {
            int c0 = c1;
            c1 = c2;
            c2 = c3;
            c3 = get(x, y);
            bits >>= 1;
            if (c3 >= sat) {
                bits |= 8;
            }
            int new_c;
            if (bits == 0) {
                new_c = (c0 + 3*c1 + 3*c2 + c3 + 4) / 8;
            } else {
                int mid_bits = bits & (2|4);
                switch (mid_bits) {
                case 0:
                    new_c = (c1 + c2 + 1) / 2;
                    break;
                case 2:
                    new_c = c2;
                    break;
                case 4:
                    new_c = c1;
                    break;
                case 2+4:
                    new_c = sat;
                    break;
                }
            }
            dst.set(dstx, y, c1);
            dst.set(dstx+1, y, new_c);
            dstx += 2;
        }
        dst.set(dstx, y, c2);
    }
    *this = std::move(dst);
}

void Plane::downsample(
    Plane &src
) {
    int wd = (src.width_ + 1) / 2;
    int ht = (src.height_ + 1) / 2;
    init(wd, ht);
    for (int y = 0; y < ht; ++y) {
        int src_y0 = 2*y;
        int src_y1 = std::min(src_y0 + 1, src.height_);
        for (int x = 0; x < wd; ++x) {
            int src_x0 = 2*x;
            int src_x1 = std::min(src_x0 + 1, src.width_);
            int c00 = src.get(src_x0, src_y0);
            int c10 = src.get(src_x1, src_y0);
            int c01 = src.get(src_x0, src_y1);
            int c11 = src.get(src_x1, src_y1);
            int c = (c00 + c01 + c10 + c11 + 2) / 4;
            set(x, y, c);
        }
    }
}

void Plane::gaussian(
    int n
) {
    gaussian_horz(n);
    transpose();
    gaussian_horz(n);
    transpose();
}

void Plane::gaussian_horz(
    int n
) {
    for (int y = 0; y < height_; ++y) {
        int div = 1;
        for (int i = 0; i < n; ++i) {
            int c1 = get(0, y);
            int c2 = c1;
            for (int x = 1; x < width_; ++x) {
                int c0 = c1;
                c1 = c2;
                c2 = get(x, y);
                int c = c0 + 2*c1 + c2;
                set(x-1, y, c);
            }
            int c = c1 + 3*c2;
            set(width_-1, y, c);
            div *= 4;
            if (div >= 8196) {
                int round = div / 2;
                for (int x = 0; x < width_; ++x) {
                    int c = get(x, y);
                    c = (c + round) / div;
                    set(x, y, c);
                }
                div = 1;
            }
        }
        if (div > 1) {
            int round = div / 2;
            for (int x = 0; x < width_; ++x) {
                int c = get(x, y);
                c = (c + round) / div;
                set(x, y, c);
            }
        }
    }
}

Image::Image() {
}

Image::~Image() {
}

void Image::init(
    int wd,
    int ht
) {
    r_.init(wd, ht);
    g1_.init(wd, ht);
    g2_.init(wd, ht);
    b_.init(wd, ht);
}

void Image::crop(
    int left,
    int top,
    int right,
    int bottom
) {
    r_.crop(left, top, right, bottom);
    g1_.crop(left, top, right, bottom);
    g2_.crop(left, top, right, bottom);
    b_.crop(left, top, right, bottom);
}

void Image::transpose() {
    r_.transpose();
    g1_.transpose();
    g2_.transpose();
    b_.transpose();
}

void Image::interpolate_horz_1331() {
    r_.interpolate_horz_1331();
    g1_.interpolate_horz_1331();
    g2_.interpolate_horz_1331();
    b_.interpolate_horz_1331();
}
