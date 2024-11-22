//
// IIR Coefficients
// Copyright (C) 2024 Jeff Plaisance
//
// This software is distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE or copy at http://boost.org/LICENSE_1_0.txt
//

#pragma once

struct BiquadFilterCoefficients {
  float b0;
  float b1;
  float b2;
  float a1;
  float a2;
};

BiquadFilterCoefficients firstOrderLowPassFilterCoefficients(float samplePeriod, float freq);

BiquadFilterCoefficients firstOrderHighPassFilterCoefficients(float samplePeriod, float freq);

BiquadFilterCoefficients secondOrderLowPassFilterCoefficients(float samplePeriod, float freq, float q);

BiquadFilterCoefficients secondOrderHighPassFilterCoefficients(float samplePeriod, float freq, float q);

BiquadFilterCoefficients firstOrderLowShelfFilterCoefficients(float samplePeriod, float freq, float gain);

BiquadFilterCoefficients firstOrderHighShelfFilterCoefficients(float samplePeriod, float freq, float gain);

BiquadFilterCoefficients secondOrderLowShelfFilterCoefficients(float samplePeriod, float freq, float q, float gain);

BiquadFilterCoefficients secondOrderHighShelfFilterCoefficients(float samplePeriod, float freq, float q, float gain);
