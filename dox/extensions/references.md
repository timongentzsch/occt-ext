# References

Numbered entries are cited by number in [utilities.md](utilities.md) and [corrections.md](corrections.md).
Entries marked *read* were read in full or in the sections named.

[1] Rida T. Farouki (2012), “The Bernstein polynomial basis: A centennial
retrospective”, Computer Aided Geometric Design 29(6), 379–419.
[Author manuscript](https://www.researchgate.net/profile/James-Peters-3/post/Since_the_Bernstein_basis_is_more_stable_than_the_1_x_x2_x3_basis_why_do_we_use_the_second_one/attachment/59d62e4cc49f478072e9ef35/AS:273573887840257@1442236505726/download/Farouki-on-Berstein-basis-2012.pdf),
degree elevation/product formulas and convex hull sections. The floating-point
implementation does not inherit an exact-arithmetic theorem without rounding analysis.

[2] Siegfried M. Rump (2010), “Verification methods: Rigorous results using
floating-point arithmetic”, Acta Numerica 19, 287–449.
[Author manuscript](https://www.tuhh.de/ti3/rump/intlab/ActaNumerica2010.pdf),
sections 2, 4 and 5.5 on IEEE arithmetic, directed rounding and inclusion.
The finite stored-double helper is narrower than general interval arithmetic.

[3] Carl de Boor (1972), “On calculating with B-splines”, Journal of Approximation
Theory 6(1), 50–62. DOI [10.1016/0021-9045(72)90080-9](https://doi.org/10.1016/0021-9045(72)90080-9).
[Read paper](https://ftp.cs.wisc.edu/Approx/bsplcalc.pdf), recurrence and convex combinations.

[4] John Skilling (2004), “Programming the Hilbert curve”, AIP Conference
Proceedings 707, 381. DOI [10.1063/1.1751381](https://doi.org/10.1063/1.1751381).
[Read transcription](https://doczz.net/doc/7960692/programming-the-hilbert-curve).
This supports the locality construction concept. The implemented two-dimensional
quadrant rotation is not a claim to reproduce Skilling's global Gray-code algorithm;
no speed theorem follows from the citation.

[5] Michael E. Hohmeyer and Brian A. Barsky (1989), “Rational continuity: parametric, geometric, and Frenet
frame continuity of rational curves”, ACM Transactions on Graphics 8(4), 335–359. The beta-constraints for G^n
(equations 105 to 109) are the end conditions of `BRepCurveExt_Bridge`; the three free scalars at an end are its
tensions.

[6] Standard constructions without a single source. A circle tangent to two curves has its centre on the
intersection of their offsets by the radius; the offset point is C + r n, solved by Newton's method with the exact
derivative dO/du = C' + r (n x dT/du), started along the pieces as in the documentation of `ChFi2d_FilletAlgo`.
Chaining unordered curves grows greedily from nearest ends. A gap in a chart closes with a straight segment
between the surface parameters of its two ends.

## Cited by name in operations.md

* Jaiswal, Ghosal and Gurumoorthy, “Circular blending”, ASME DETC 2001 (theory, algorithm and results read;
  Hermann 1995 and Lukacs 1998 are cited second-hand there).
* Salvi and Várady (2014), “Multi-sided surfaces with curvature continuity” (statement of the Linkage Curve
  Theorem read, appendix proof not).
* Zang, Liu and Lai, GMP 2008, sections 1 to 3, on Hu's surface extension (Hu 2002 not obtained).
* Floater (2003), “One-to-one piecewise linear mappings over triangulations”, Mathematics of Computation 72,
  Theorem 4.1 (statement read, proof not). Karni, Gotsman and Gortler, and Weber and Zorin, on non-convex
  boundaries.
* Hermann, Peters and Strotman (2009), on the vertex-enclosure condition (all five pages read).
* Zou and Feng (2019), push-pull modelling of solids.
* Wang, Jüttler, Zheng and Liu (2008), “Computation of rotation minimizing frames”, ACM Transactions on
  Graphics 27(1), double reflection method.
* Siggel et al., Gordon surfaces, section 3.2.
* Dogan, Bernal and Hagwood (2025), seam alignment of closed sections by cyclic cross-correlation (the
  rotation proof by Kabsch / SVD read; this implementation fixes rotation).
* Sawhney and Crane (2017), “Boundary First Flattening”, sections 3 to 5.
* Piessens and de Doncker (1980, rev. 1983), QUADPACK `dqk15`, <https://netlib.org/quadpack/dqk15.f>.
