/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
common classes.
**/

const int kSaturated = 0x40000000;

class RggbPixel {
public:
    int r_ = 0;
    int g1_ = 0;
    int g2_ = 0;
    int b_ = 0;
};

class RggbDouble {
public:
    double r_ = 1.0;
    double g1_ = 1.0;
    double g2_ = 1.0;
    double b_ = 1.0;
};

int pin_to_8bits(int x);
int pin_to_16bits(int x);
