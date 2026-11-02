// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once
#if defined(__wasm_simd128__) && !defined(OCCT_GEOMALGOEXT_SCALAR)
  #include <wasm_simd128.h>
#else
  #include <cstring>

namespace GeomAlgoExt
{
struct v128_t
{
  double lane[2];
};

inline v128_t wasm_f64x2_splat(double a)
{
  return {{a, a}};
}

inline v128_t wasm_v128_load(const void* p)
{
  v128_t a;
  std::memcpy(&a, p, sizeof(a));
  return a;
}

inline void wasm_v128_store(void* p, v128_t a)
{
  std::memcpy(p, &a, sizeof(a));
}

inline v128_t wasm_f64x2_add(v128_t a, v128_t b)
{
  volatile double x = a.lane[0] + b.lane[0], y = a.lane[1] + b.lane[1];
  return {{x, y}};
}

inline v128_t wasm_f64x2_sub(v128_t a, v128_t b)
{
  volatile double x = a.lane[0] - b.lane[0], y = a.lane[1] - b.lane[1];
  return {{x, y}};
}

inline v128_t wasm_f64x2_mul(v128_t a, v128_t b)
{
  volatile double x = a.lane[0] * b.lane[0], y = a.lane[1] * b.lane[1];
  return {{x, y}};
}

inline double wasm_f64x2_extract_lane(v128_t a, int lane)
{
  return a.lane[lane];
}
} // namespace GeomAlgoExt
#endif
