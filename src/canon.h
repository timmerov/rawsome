/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

/**
functions specific to canon.
**/

#include "common.h"
#include "planes.h"

void determine_black(Planes &planes, RggbPixel &black, int &noise);

void crop_black(Planes &planes);

