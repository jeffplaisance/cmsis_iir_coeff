//
// IIR Coefficients
// Copyright (C) 2024 Jeff Plaisance
//
// This software is distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE or copy at http://boost.org/LICENSE_1_0.txt
//

#include "iir_coefficients.h"
#include <cmath>

static constexpr float pi = M_PI;

BiquadFilterCoefficients firstOrderLowPassFilterCoefficients(float samplePeriod, float freq) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float cosOmega0p1 = cosOmega0 + 1.0f;

  float a0inverse = 1.0f / (sinOmega0 + cosOmega0p1);
  ret.a1 = -(sinOmega0 - cosOmega0p1) * a0inverse;
  ret.a2 = 0;
  ret.b0 = sinOmega0 * a0inverse;
  ret.b1 = ret.b0;
  ret.b2 = 0;
  return ret;
}

BiquadFilterCoefficients firstOrderHighPassFilterCoefficients(float samplePeriod, float freq) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float cosOmega0p1 = cosOmega0 + 1.0f;

  float a0inverse = 1.0f / (sinOmega0 + cosOmega0p1);
  ret.a1 = -(sinOmega0 - cosOmega0p1) * a0inverse;
  ret.a2 = 0;
  ret.b0 = cosOmega0p1 * a0inverse;
  ret.b1 = -ret.b0;
  ret.b2 = 0;
  return ret;
}

BiquadFilterCoefficients secondOrderLowPassFilterCoefficients(float samplePeriod, float freq, float q) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float alpha = sinOmega0 * 0.5f / q;

  float a0inverse = 1.0f / (1.0f + alpha);
  ret.a1 = 2.0f * cosOmega0 * a0inverse;
  ret.a2 = -(1.0f - alpha) * a0inverse;
  ret.b1 = (1.0f - cosOmega0) * a0inverse;
  ret.b0 = ret.b1 * 0.5f;
  ret.b2 = ret.b0;
  return ret;
}

BiquadFilterCoefficients secondOrderHighPassFilterCoefficients(float samplePeriod, float freq, float q) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float alpha = sinOmega0 * 0.5f / q;

  float a0inverse = 1.0f / (1.0f + alpha);
  ret.a1 = 2.0f * cosOmega0 * a0inverse;
  ret.a2 = -(1.0f - alpha) * a0inverse;
  ret.b1 = -(1.0f + cosOmega0) * a0inverse;
  ret.b0 = -ret.b1 * 0.5f;
  ret.b2 = ret.b0;
  return ret;
}

BiquadFilterCoefficients firstOrderLowShelfFilterCoefficients(float samplePeriod, float freq, float gain) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float cosOmega0p1 = cosOmega0 + 1.0f;
  float A = ::sqrtf(gain);
  float AInverse = 1.0f / A;
  float sinOmega0A = sinOmega0 * A;
  float sinOmega0AInverse = sinOmega0 * AInverse;

  float a0inverse = 1.0f / (sinOmega0AInverse + cosOmega0p1);
  ret.a1 = -(sinOmega0AInverse - cosOmega0p1) * a0inverse;
  ret.a2 = 0;
  ret.b0 = (sinOmega0A + cosOmega0p1) * a0inverse;
  ret.b1 = (sinOmega0A - cosOmega0p1) * a0inverse;
  ret.b2 = 0;

  return ret;
}

BiquadFilterCoefficients firstOrderHighShelfFilterCoefficients(float samplePeriod, float freq, float gain) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);
  float cosOmega0p1 = cosOmega0 + 1.0f;
  float A = ::sqrtf(gain);
  float AInverse = 1.0f / A;
  float cosOmega0p1A = cosOmega0p1 * A;
  float cosOmega0p1AInverse = cosOmega0p1 * AInverse;

  float a0inverse = 1.0f / (sinOmega0 + cosOmega0p1AInverse);
  ret.a1 = -(sinOmega0 - cosOmega0p1AInverse) * a0inverse;
  ret.a2 = 0;
  ret.b0 = (sinOmega0 + cosOmega0p1A) * a0inverse;
  ret.b1 = (sinOmega0 - cosOmega0p1A) * a0inverse;
  ret.b2 = 0;

  return ret;
}

BiquadFilterCoefficients secondOrderLowShelfFilterCoefficients(float samplePeriod, float freq, float q, float gain) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);

  float A = ::sqrtf(gain);
  float Ap1 = A + 1.0f;
  float Am1 = A - 1.0f;

  float cosOmega0Ap1 = cosOmega0 * Ap1;
  float cosOmega0Am1 = cosOmega0 * Am1;

  float alpha = sinOmega0 * 0.5f / q;

  float twoSqrtAAlpha = 2.0f * ::sqrtf(A) * alpha;

  float Ap1pCosOmega0Am1 = Ap1 + cosOmega0Am1;
  float Ap1mCosOmega0Am1 = Ap1 - cosOmega0Am1;

  float a0inverse = 1.0f / (Ap1pCosOmega0Am1 + twoSqrtAAlpha);
  ret.a1 = 2.0f * (Am1 + cosOmega0Ap1) * a0inverse;
  ret.a2 = -(Ap1pCosOmega0Am1 - twoSqrtAAlpha) * a0inverse;
  ret.b0 = A * (Ap1mCosOmega0Am1 + twoSqrtAAlpha) * a0inverse;
  ret.b1 = 2.0f * A * (Am1 - cosOmega0Ap1) * a0inverse;
  ret.b2 = A * (Ap1mCosOmega0Am1 - twoSqrtAAlpha) * a0inverse;

  return ret;
}

BiquadFilterCoefficients secondOrderHighShelfFilterCoefficients(float samplePeriod, float freq, float q, float gain) {
  BiquadFilterCoefficients ret{};
  float omega0 = 2.0f * pi * freq * samplePeriod;
  float cosOmega0 = ::cosf(omega0);
  float sinOmega0 = ::sinf(omega0);

  float A = ::sqrtf(gain);
  float Ap1 = A + 1.0f;
  float Am1 = A - 1.0f;

  float cosOmega0Ap1 = cosOmega0 * Ap1;
  float cosOmega0Am1 = cosOmega0 * Am1;

  float alpha = sinOmega0 * 0.5f / q;

  float twoSqrtAAlpha = 2.0f * ::sqrtf(A) * alpha;

  float Ap1pCosOmega0Am1 = Ap1 + cosOmega0Am1;
  float Ap1mCosOmega0Am1 = Ap1 - cosOmega0Am1;

  float a0inverse = 1.0f / (Ap1mCosOmega0Am1 + twoSqrtAAlpha);
  ret.a1 = -2.0f * (Am1 - cosOmega0Ap1) * a0inverse;
  ret.a2 = -(Ap1mCosOmega0Am1 - twoSqrtAAlpha) * a0inverse;
  ret.b0 = A * (Ap1pCosOmega0Am1 + twoSqrtAAlpha) * a0inverse;
  ret.b1 = -2.0f * A * (Am1 + cosOmega0Ap1) * a0inverse;
  ret.b2 = A * (Ap1pCosOmega0Am1 - twoSqrtAAlpha) * a0inverse;

  return ret;
}
