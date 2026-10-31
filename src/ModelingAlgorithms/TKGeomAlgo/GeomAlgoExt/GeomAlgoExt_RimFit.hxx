// Copyright (c) 2026 Timon Gentzsch
//
// This file is part of Open CASCADE Technology software library.
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult LICENSE_LGPL_21.txt for the complete license.

#pragma once

// The last of the way from a fair surface to the curves it is to pass through. A fit by least
// squares over a net weighs what it is asked against how the surface bends, and what it leaves
// of a curve is what that net cannot follow fairly: a curve from a Boolean wavers about its
// faces by a share of its size, and a tolerance that is a length and not a share asks for more
// of that the larger the model is. Here the surface is moved onto the curves by the least
// movement of its poles that does it, with knots added where the net as it is cannot: a small
// change beside the curve and none away from it.

#include <GeomAlgoExt_SplineFit.hxx>

#include <gp_Pnt.hxx>
#include <gp_XY.hxx>

#include <functional>

namespace GeomAlgoExt {

/// A place that the surface is to pass through. `uv` is where to look for the surface's point
/// nearest to it, and is where that was found afterwards. With `across`, a direction of the
/// parameters that leaves the curve the place is on, the surface's tangent plane is not to
/// tilt that way; zero leaves it free. With `bend` as well, the surface's bend that way is not
/// to change either.
struct RimPlace {
    gp_Pnt point;
    gp_XY uv;
    gp_XY across;
    bool bend;
};

/// The parameters of the point of a surface nearest to `point`, looked for from `uv`: steps of
/// Gauss and Newton, which for a point that is nearly on the surface are as good as exact.
Standard_EXPORT gp_XY nearest(const occ::handle<Geom_Surface> &surface, gp_XY uv, const gp_Pnt &point);

/// The net of a surface: the reverse of `surfaceOf`.
Standard_EXPORT SplineNet netOf(const occ::handle<Geom_BSplineSurface> &surface);

/// `placesFor` gives the places for a net whose shortest span is `span`: some to every span
/// that a curve passes through, and is asked again whenever the net is made finer.
/// Move the poles of `net` along `up`, a direction of unit length that its surface faces,
/// until the surface is within `within` of every place, by the least movement that does it;
/// where the net cannot, each span that a place still misses in is halved, up to `mostSpans`
/// either way. A surface within `enough` of every place to begin with is left as it is.
/// Returns the largest distance that is left.
Standard_EXPORT double refineRim(const char *operation, SplineNet &net,
                  const std::function<std::vector<RimPlace>(double span)> &placesFor, const gp_XYZ &up,
                  double within, double enough, int mostSpans, const Message_ProgressRange& range = {}, StageObserver observer = nullptr);

} // namespace GeomAlgoExt
