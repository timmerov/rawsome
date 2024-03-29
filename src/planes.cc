/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
data structures for image processing.
**/

#include "common.h"
#include "planes.h"

#include <cmath>
#include <thread>

namespace {

/** derived from dcraw. **/
void display_gamma_curve(
    std::vector<int> &curve,
    double pwr,
    double ts,
    int imax
){
    curve.resize(0x10000);

    #define SQR(x) ((x)*(x))
    double g2 = 0.0;
    double g3 = 0.0;
    double g4 = 0.0;

    double bnd[2] = {0.0, 0.0};
    double r;

    pwr = pwr;
    ts = ts;
    g2 = g3 = g4 = 0;
    bnd[ts >= 1] = 1;
    if ((ts-1)*(pwr-1) <= 0) {
        for (int i = 0; i < 48; ++i) {
            g2 = (bnd[0] + bnd[1])/2;
            bnd[(std::pow(g2/ts,-pwr) - 1)/pwr - 1/g2 > -1] = g2;
        }
        g3 = g2 / ts;
        g4 = g2 * (1/pwr - 1);
    }
    for (int i = 0; i < 0x10000; ++i) {
        curve[i] = 0xffff;
        r = (double) i / imax;
        if (r < 1) {
            curve[i] = 0x10000 * (r < g3 ? r*ts : std::pow(r,pwr)*(1+g4)-g4);
        }
    }
}

void user_gamma_curve(
    std::vector<int> &curve,
    double pwr
){
    curve.resize(0x10000);

    for (int i = 0; i < 0x10000; ++i) {
        double in = double(i) / double(0x0FFFF);
        double out = std::pow(in, pwr);
        int x = int(out * double(0x0FFFF));
        if (x < 0) {
            x = 0;
        } else if (x > 0x0FFFF) {
            x = 0x0FFFF;
        }
        curve[i] = x;
    }
}

void interpolate_horz_1331_thread(
    Plane *src,
    Plane *dst,
    int y0,
    int dy
) {
    int ht = src->height_;
    int wd = src->width_;
    for (int y = y0; y < ht; y += dy) {
        int c1 = src->get(0, y);
        int c2 = src->get(1, y);
        int c3 = src->get(2, y);
        int dstx = 0;
        for (int x = 3; x < wd; ++x) {
            int c0 = c1;
            c1 = c2;
            c2 = c3;
            c3 = src->get(x, y);
            /**
            r is the midpoint of segment c1.c2.
            s is the extrapolated value of segment c0.c1.
            t is the extrapolated value of segment c2.c3.
            if s,t are both greater than r...
            or s,t are both smaller than r...
            then average r with the closer of s,t.
            otherwise use r.
            **/
            int r = c1 + c2;
            int s = 3*c1 - c0;
            int t = 3*c2 - c3;
            int new_c;
            if (s >= r && t >= r) {
                s = std::min(s, t);
                new_c = (r + s + 2) / 4;
            } else if (s <= r && t <= r) {
                s = std::max(s, t);
                new_c = (r + s + 2) / 4;
            } else {
                new_c = (r + 1) / 2;
            }
            dst->set(dstx, y, c1);
            dst->set(dstx+1, y, new_c);
            dstx += 2;
        }
        dst->set(dstx, y, c2);
    }
}

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

void Plane::subtract(
    int delta
) {
    int sz = width_ * height_;
    for (int i = 0; i < sz; ++i) {
        samples_[i] -= delta;
    }
}

void Plane::multiply(
    double factor
) {
    int sz = width_ * height_;
    for (int i = 0; i < sz; ++i) {
        int x = samples_[i];
        x = std::round(x * factor);
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

namespace {
void transpose_thread(
    Plane *src,
    Plane *dst,
    int y0,
    int dy
) {
    int ht = src->height_;
    int wd = src->width_;
    for (int y = y0; y < ht; y += dy) {
        for (int x = 0; x < wd; ++x) {
            int c = src->get(x,y);
            dst->set(y, x, c);
        }
    }
}
}

void Plane::transpose_mt() {
    Plane dst;
    dst.init(height_, width_);
    std::thread th0(transpose_thread, this, &dst, 0, 8);
    std::thread th1(transpose_thread, this, &dst, 1, 8);
    std::thread th2(transpose_thread, this, &dst, 2, 8);
    std::thread th3(transpose_thread, this, &dst, 3, 8);
    std::thread th4(transpose_thread, this, &dst, 4, 8);
    std::thread th5(transpose_thread, this, &dst, 5, 8);
    std::thread th6(transpose_thread, this, &dst, 6, 8);
    std::thread th7(transpose_thread, this, &dst, 7, 8);
    th0.join();
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
    *this = std::move(dst);
}

void Plane::interpolate_1331() {
    transpose();
    interpolate_horz_1331();
    transpose();
    interpolate_horz_1331();
}

void Plane::interpolate_1331_mt() {
    transpose_mt();
    interpolate_horz_1331_mt();
    transpose_mt();
    interpolate_horz_1331_mt();
}

void Plane::interpolate_horz_1331() {
    int new_wd = 2*width_ - 5;
    Plane dst;
    dst.init(new_wd, height_);
    interpolate_horz_1331_thread(this, &dst, 0, 1);
    *this = std::move(dst);
}

void Plane::interpolate_horz_1331_mt() {
    int new_wd = 2*width_ - 5;
    Plane dst;
    dst.init(new_wd, height_);
    std::thread th0(interpolate_horz_1331_thread, this, &dst, 0, 8);
    std::thread th1(interpolate_horz_1331_thread, this, &dst, 1, 8);
    std::thread th2(interpolate_horz_1331_thread, this, &dst, 2, 8);
    std::thread th3(interpolate_horz_1331_thread, this, &dst, 3, 8);
    std::thread th4(interpolate_horz_1331_thread, this, &dst, 4, 8);
    std::thread th5(interpolate_horz_1331_thread, this, &dst, 5, 8);
    std::thread th6(interpolate_horz_1331_thread, this, &dst, 6, 8);
    std::thread th7(interpolate_horz_1331_thread, this, &dst, 7, 8);
    th0.join();
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
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
    gaussian_horz_mt(n);
    transpose();
    gaussian_horz_mt(n);
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

namespace {
void gaussian_horz_thread(
    Plane *plane,
    int n,
    int y0,
    int dy
) {
    int ht = plane->height_;
    int wd = plane->width_;
    for (int y = y0; y < ht; y += dy) {
        int div = 1;
        for (int i = 0; i < n; ++i) {
            int c1 = plane->get(0, y);
            int c2 = c1;
            for (int x = 1; x < wd; ++x) {
                int c0 = c1;
                c1 = c2;
                c2 = plane->get(x, y);
                int c = c0 + 2*c1 + c2;
                plane->set(x-1, y, c);
            }
            int c = c1 + 3*c2;
            plane->set(wd-1, y, c);
            div *= 4;
            if (div >= 8196) {
                int round = div / 2;
                for (int x = 0; x < wd; ++x) {
                    int c = plane->get(x, y);
                    c = (c + round) / div;
                    plane->set(x, y, c);
                }
                div = 1;
            }
        }
        if (div > 1) {
            int round = div / 2;
            for (int x = 0; x < wd; ++x) {
                int c = plane->get(x, y);
                c = (c + round) / div;
                plane->set(x, y, c);
            }
        }
    }
}
}

void Plane::gaussian_horz_mt(
    int n
) {
    std::thread th0(gaussian_horz_thread, this, n, 0, 8);
    std::thread th1(gaussian_horz_thread, this, n, 1, 8);
    std::thread th2(gaussian_horz_thread, this, n, 2, 8);
    std::thread th3(gaussian_horz_thread, this, n, 3, 8);
    std::thread th4(gaussian_horz_thread, this, n, 4, 8);
    std::thread th5(gaussian_horz_thread, this, n, 5, 8);
    std::thread th6(gaussian_horz_thread, this, n, 6, 8);
    std::thread th7(gaussian_horz_thread, this, n, 7, 8);
    th0.join();
    th1.join();
    th2.join();
    th3.join();
    th4.join();
    th5.join();
    th6.join();
    th7.join();
}

void Plane::apply_display_gamma(
    double pwr,
    double ts,
    int white
) {
    std::vector<int> curve;
    display_gamma_curve(curve, pwr, ts, white);
    apply_display_gamma(curve);
}

void Plane::apply_display_gamma(
    std::vector<int> &curve
) {
    int sz = width_ * height_;
    for (int i = 0; i < sz; ++i) {
        int c = samples_[i];
        c = pin_to_16bits(c);
        c = curve[c];
        samples_[i] = c;
    }
}

void Plane::scale_to_8bits() {
    double factor = 255.0/65535.0;
    multiply(factor);
}

Planes::Planes() {
}

Planes::~Planes() {
}

void Planes::init(
    int wd,
    int ht
) {
    r_.init(wd, ht);
    g1_.init(wd, ht);
    g2_.init(wd, ht);
    b_.init(wd, ht);
}

void Planes::set(
    int x,
    int y,
    RggbPixel p
) {
    r_.set(x, y, p.r_);
    g1_.set(x, y, p.g1_);
    g2_.set(x, y, p.g2_);
    b_.set(x, y, p.b_);
}

void Planes::set(
    int x,
    int y,
    int value
) {
    r_.set(x, y, value);
    g1_.set(x, y, value);
    g2_.set(x, y, value);
    b_.set(x, y, value);
}

void Planes::subtract(
    RggbPixel &delta
) {
    r_.subtract(delta.r_);
    g1_.subtract(delta.g1_);
    g2_.subtract(delta.g2_);
    b_.subtract(delta.b_);
}

void Planes::multiply(
    RggbDouble &factor
) {
    r_.multiply(factor.r_);
    g1_.multiply(factor.g1_);
    g2_.multiply(factor.g2_);
    b_.multiply(factor.b_);
}

void Planes::multiply4(
    double factor
) {
    r_.multiply(factor);
    g1_.multiply(factor);
    g2_.multiply(factor);
    b_.multiply(factor);
}

void Planes::multiply3(
    double factor
) {
    r_.multiply(factor);
    g1_.multiply(factor);
    b_.multiply(factor);
}

void Planes::crop(
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

void Planes::transpose() {
    r_.transpose_mt();
    g1_.transpose_mt();
    g2_.transpose_mt();
    b_.transpose_mt();
}

void Planes::interpolate_1331() {
    /**
    applying the 1331 filter is really fast for horizontal pixels.
    so we transpose while small.
    add pixels horizontally.
    tanspose again while medium sized.
    add pixels horizontally again.
    **/
    r_.interpolate_1331();
    g1_.interpolate_1331();
    g2_.interpolate_1331();
    b_.interpolate_1331();
}

void Planes::apply_user_gamma(
    double pwr
) {
    /**
    luminance_out = luminance_in ^ user_gamma
    rgb_out = rgb_in * luminance_out / luminance_in
    **/

    /** precompute the power curve. **/
    std::vector<int> curve;
    user_gamma_curve(curve, pwr);

    /** preserve luminance. **/
    int sz = r_.width_ * r_.height_;
    for (int i = 0; i < sz; ++i) {
        double r = r_.samples_[i];
        double g = g1_.samples_[i];
        double b = b_.samples_[i];
        /** convert srgb to luminance. **/
        double lum_in = 0.2125*r + 0.7154*g + 0.0721*b;
        int idx = int(lum_in);
        if (idx < 0 || idx >= 0x10000) {
            continue;
        }
        double lum_out = curve[idx];
        double factor = lum_out / lum_in;
        r *= factor;
        g *= factor;
        b *= factor;
        r_.samples_[i] = r;
        g1_.samples_[i] = g;
        b_.samples_[i] = b;
    }
}

void Planes::apply_display_gamma(
    double pwr,
    double ts,
    int white
) {
    std::vector<int> curve;
    display_gamma_curve(curve, pwr, ts, white);
    r_.apply_display_gamma(curve);
    g1_.apply_display_gamma(curve);
    b_.apply_display_gamma(curve);
}

void Planes::scale_to_8bits() {
    double factor = 255.0/65535.0;
    multiply3(factor);
}
