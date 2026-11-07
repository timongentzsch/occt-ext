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

[5] Standard constructions without a single source. Chaining unordered curves grows greedily from nearest ends. A gap in a chart closes with a straight segment
between the surface parameters of its two ends.

## Keys used in operations.md

Entries are cited as `[Key section, equation]`. "Read" states what was read.

**[PMC09]** Nicholas M. Patrikalakis, Takashi Maekawa and Wonjoon Cho (2009), *Shape Interrogation for Computer
Aided Design and Manufacturing*, MIT hyperbook edition,
<https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/>. Sections are pages `nodeN.html` of that site:
1.3.6 node14, 1.4.1 node16, 1.4.2 node17, 1.4.3 node18, 1.5 node20, 3.3 node29, 3.4 node30, 4.9.1 node53,
7.2.1 node135, 8.1.3.5 node156 (Linkage Curve Theorem), 9.7.1 node190, 11.1.4 node217, 11.3.1 node224,
11.3.2 node225, 11.3.7 node230, 11.6.2 node237, 11.6.3.4 node242. Used for the Bernstein and B-spline identities
(convex hull property, partition of unity, knot insertion, derivative nets), fundamental forms and normal curvature,
closest-point stationarity, ruled surfaces, normal offsets, focal distances and offset pipes.

**[Floater03]** Michael S. Floater (2003), "One-to-one piecewise linear mappings over triangulations", Mathematics of
Computation 72, Theorem 4.1 (statement read, proof not). Non-convex boundaries: Karni, Gotsman and Gortler; Weber
and Zorin.
**[HermannPetersStrotman09]** Hermann, Peters and Strotman (2009), vertex-enclosure condition (all five pages read).
**[SalviVarady14]** Salvi and Várady (2014), "Multi-sided surfaces with curvature continuity",
<https://salvi.chaosnet.org/texts/ribbons.pdf> (§3.3.2 and the statement of the Linkage Curve Theorem read, appendix
proof not).
**[Siggel18]** Siggel, Kleinert, Stollenwerk and Maierl (2018/2019), "TiGL: An Open Source Computational Geometry
Library for Parametric Aircraft Design", <https://arxiv.org/html/1810.10795v1> (also
<https://elib.dlr.de/124524/1/1810.10795.pdf>); sections 3.2.1-3.2.3 and Algorithm 1; implementation
`CTiglGordonSurfaceBuilder::CreateGordonSurface`,
<https://github.com/DLR-SC/tigl/blob/main/src/geometry/CTiglGordonSurfaceBuilder.cpp>, lines 139-188.
**[Gordon71]** W. J. Gordon (1971), "Blending-function methods of bivariate and multivariate interpolation and
approximation", SIAM Journal on Numerical Analysis 8, 158-177.
**[PieglTiller97]** Piegl and Tiller (1997), *The NURBS Book*, second edition, chapter 10 (§10.3, (10.7)),
<https://home.zcu.cz/~bastl/GM1/the-nurbs-book.pdf>.
**[EngleitnerJuttler18]** Engleitner and Jüttler (2018), G+S Report 73, §2-4, <https://www.gs.jku.at/pubs/NFNreport73.pdf>.
**[Floater97]** Floater (1997), "An O(h^2n) Hermite approximation for conic sections", §2, Theorem 3.1,
<https://numerika.fmf.uni-lj.si/bib2html/e-documents/clanki/FloaterAnO(h2n)Hermiteapproximationforconicsections.pdf>.
**[HuWangWang18]** Hu, Wang and Wang (2018), "Piecewise Möbius reparameterization of rational Bézier curves", §2.1-2.2,
<https://www.jcad.cn/cn/article/pdf/preview/10.3724/SP.J.1089.2018.16764.pdf>.
**[PegnaWolter92]** Pegna and Wolter (1992), Journal of Mechanical Design 114, linkage curve theorem.
**[WangJuttler08]** Wang, Jüttler, Zheng and Liu (2008), "Computation of rotation minimizing frames", ACM Transactions
on Graphics 27(1), double reflection method.
**[Dogan25]** Dogan, Bernal and Hagwood (2025), seam alignment of closed sections by cyclic cross-correlation (rotation
proof by Kabsch / SVD read; this implementation fixes rotation).
**[ZouFeng19]** Zou and Feng (2019), push-pull modelling of solids; preprint "Push-Pull Direct Modeling of Solid CAD
Models" (2018), §5.2, Algorithm 3, <https://arxiv.org/pdf/1811.04036>.
**[Federer59]** Federer (1959), "Curvature measures", Theorem 4.8 (reach).
**[SawhneyCrane17]** Sawhney and Crane (2017), "Boundary First Flattening", §4-7,
<https://www.cs.cmu.edu/~kmcrane/Projects/BoundaryFirstFlattening/paper.pdf>.
**[SpringbornSchroederPinkall08]** Springborn, Schröder and Pinkall (2008), conformal equivalence of triangle meshes,
appendix E. **[Pinkall93]** Pinkall and Polthier (1993), cotangent Laplacian.
**[Kobbelt97]** Kobbelt (1997), "Discrete fairing". **[BotschKobbelt04]** Botsch and Kobbelt (2004), "An intuitive
framework for real-time freeform modeling", SIGGRAPH.
**[Farin02]** Farin (2002), *Curves and Surfaces for CAGD*, §5.8. **[Gravesen14]** Gravesen, Evgrafov, Nguyen and
Nortoft (2014), Bernstein test of a Jacobian's sign in planar parametrisation.
**[Boyd04]** Boyd and Vandenberghe (2004), *Convex Optimization*, §8.1.1 (support function).
**[Gilbert66]** Gilbert (1966), IEEE Journal of Robotics and Automation 2, nearest point of a convex hull (as cited; not re-read).
**[Shewchuk97]** Shewchuk (1997), adaptive precision floating-point arithmetic, `ccwerrboundA` of `orient2d`.
**[SandeManniSpeleers20]** Sande, Manni and Speleers (2020), "Explicit error estimates for spline approximation of
arbitrary smoothness in isogeometric analysis", §1.2, (2)-(3), §6-7, <https://arxiv.org/html/1909.03559v3>.
**[NocedalWright06]** Nocedal and Wright (2006), *Numerical Optimization*, second edition, §10.2-10.3, (10.20)-(10.26)
(pp. 252-255) and §18.4, Algorithm 18.3, <https://mcube.lab.nycu.edu.tw/~cfung/docs/books/nocedal2006.pdf>.
**[HormannAgathos01]** Hormann and Agathos (2001), "The point in polygon problem for arbitrary polygons", §1-2, (1)-(3),
<https://www.inf.usi.ch/hormann/papers/Hormann.2001.TPI.pdf>.
**[DymLawrenceSiegel24]** Dym, Lawrence and Siegel (2024), "Equivariant frames and the impossibility of continuous
canonicalization", §2.2-3, <https://proceedings.mlr.press/v235/dym24a.html>.
**[GreengardSerkh18]** Greengard and Serkh (2018), "Zernike polynomials: evaluation, quadrature, and interpolation",
§2.3, 4, <https://arxiv.org/abs/1811.02720>.
**[Narayan22]** Akil Narayan (2022), "Regularization and relaxation", University of Utah lecture slides 6-7 (Tikhonov
regularization), <https://www.sci.utah.edu/~akil/docs/courses/2022spring/math6880/slides09.pdf>.
**[Mukundan04]** Mukundan, Ko, Maekawa, Sakkalis and Patrikalakis (2004), "Tracing surface intersections with validated
ODE system solver", §1, 2.2 (4)-(6), 3-4, Theorem 1, <https://web.mit.edu/deslab/pubs/files/SM04.pdf>.
**[Sage10]** Malakhaltsev and Vankerschaver (2010), *Sage Reference Manual: Differential Geometry of Parametrized
Surfaces* (Sage 9.3 edition), `principal_directions`, torus example,
<https://sporadic.stanford.edu/reference/riemannian_geometry/sage/geometry/riemannian_manifolds/parametrized_surface3d.html>.
**[OCCT8]** OPEN CASCADE, *OCCT Reference Manual* 8.0.1 (2026), class pages on <https://occt3d.com/dev/doc/refman/html/>: `GeomConvert`,
`GeomConvert_BSplineSurfaceToBezierSurface`, `Geom_BSplineSurface` (also
<https://dev.opencascade.org/doc/refman/html/class_geom___b_spline_surface.html>), `BRepOffsetAPI_MakeThickSolid`
(`MakeThickSolidByJoin`), `BRepOffsetAPI_MakeOffsetShape` (`PerformByJoin`), `BRep_Tool` (`IsClosed`,
`CurveOnSurface`), `BRepExtrema_DistShapeShape` (`Perform`), `BRepCheck_Shell` (`Closed`), `BRepLib`
(`OrientClosedSolid`), `BRepTools_History`; *Modeling Data* ("Orientation"); *Boolean Operations* ("Interferences");
`BRepTools_Modifier::Rebuild` at
<https://github.com/Open-Cascade-SAS/OCCT/blob/b8f597c6/src/ModelingData/TKBRep/BRepTools/BRepTools_Modifier.cxx>.

## Cited by name in operations.md (blend chapters)

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
