/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

/**
common utility functions.
**/

#include "common.h"

int pin_to_8bits(
    int x
) {
    if (x < 0) {
        return 0;
    }
    if (x > 255) {
        return 255;
    }
    return x;
}

int pin_to_16bits(
    int x
) {
    if (x < 0) {
        return 0;
    }
    if (x > 65535) {
        return 65535;
    }
    return x;
}
