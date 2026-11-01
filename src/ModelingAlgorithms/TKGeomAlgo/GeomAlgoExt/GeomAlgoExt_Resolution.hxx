// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#pragma once
namespace loft::tolerance {
// Inherited numerical allowances; empirical, without a proved error bound.
inline constexpr double fineAngle = 1e-12;
inline constexpr double parameterMerge = 1e-9;
}
