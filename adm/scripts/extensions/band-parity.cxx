// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

// Parity harness of the banded solver GeomAlgoExt::Band: factors and solves seeded random banded
// systems, ordinary and paired, and prints one FNV-1a digest of every bit of the factors and
// solutions. The digest must be the same for a native build, the WebAssembly SIMD build and the
// forced scalar build (see dox/extensions/utilities.md, "Portable SIMD").

#include <GeomAlgoExt_Band.hxx>
#include <cstdint>
#include <cstdio>
#include <cstring>

int main()
{
  std::uint64_t digest = 14695981039346656037ull, seed = 0x1938ad37;
  auto          sample = [&]() {
    seed = seed * 6364136223846793005ull + 1;
    return double(seed >> 12) / double(1ull << 52) - .5;
  };
  auto record = [&](double x) {
    std::uint64_t b;
    std::memcpy(&b, &x, 8);
    for (int k = 0; k < 8; ++k)
    {
      digest ^= (b >> (8 * k)) & 255;
      digest *= 1099511628211ull;
    }
  };
  for (bool paired : {false, true})
    for (int n : {7, 16, 37, 64})
      for (int half : {1, 3, 6})
      {
        GeomAlgoExt::Band   a(n, half, paired);
        std::vector<double> x(n);
        for (int i = 0; i < n; ++i)
        {
          a.add(i, i, double(half + 1));
          x[i] = sample();
          for (int j = std::max(0, i - half); j < i; ++j)
            a.add(i, j, sample());
        }
        if (!a.factor())
          return 1;
        for (double v : a.data)
          record(v);
        a.solve(x);
        for (double v : x)
          record(v);
      }
  std::printf("%016llx\n", static_cast<unsigned long long>(digest));
}
