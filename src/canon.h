/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
functions specific to canon.
**/

#include "common.h"
#include "image.h"

void determine_black(Image &image, RggbPixel &black, int &noise);
