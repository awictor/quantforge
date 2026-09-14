# Changelog

All notable changes to QuantForge are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.576.0] - 2026-09-14

### Documentation
- README: documented `weighted_reservoir_sample` and `weighted_sample_with_replacement`
  alongside the uniform `reservoir_sample` in the streaming-statistics section -- the
  A-Res key rule, the empirical `k=1` inclusion probabilities matching `w_i / sum(w)`,
  zero-weight exclusion, and with-replacement frequency convergence. All snippet values
  verified live.

## [1.575.0] - 2026-09-14

### Added
- `weighted_reservoir.py`: `weighted_reservoir_sample` selects `k` distinct items with
  probability proportional to weight in a single streaming pass -- the
  Efraimidis-Spirakis A-Res algorithm, which assigns each item the key ``u^(1/w)`` and
  keeps the ``k`` largest via a size-``k`` min-heap (``O(n log k)`` time, ``O(k)``
  memory) -- plus `weighted_sample_with_replacement` for independent weight-proportional
  draws via cumulative-weight bisection. Cross-checked: the ``k=1`` inclusion probability
  matches the closed form ``w_i / sum(w)`` empirically, ``k`` results are distinct,
  ``k`` past the positive-weight count returns all such items, zero-weight items are
  never selected, and the with-replacement frequencies converge to the weight fractions.

## [1.574.0] - 2026-09-12

### Added
- `entropy_pooling.py`: `entropy_pooling_mean` reweights a scenario set to satisfy
  a target-mean view with minimal Kullback-Leibler divergence from the prior
  (Meucci 2008) -- an exponential tilt ``p_i ∝ q_i exp(lambda x_i)`` with
  ``lambda`` solved by bisection -- plus `relative_entropy` (the KL divergence).
  Cross-checked: the posterior sums to 1 and hits the target mean exactly, a view
  equal to the prior mean returns the prior, the relative entropy is positive under
  a view and zero without, a stronger view costs more entropy, and an infeasible
  target (outside the scenario range) is rejected.

## [1.721.0] - 2026-09-12

### Documentation
- README distributional-tests section now documents the classical parametric tests
  (`chi_square_gof_test`, `chi_square_independence_test`, `one_way_anova`,
  `two_sample_t_test`, `binomial_test`) with a worked example, noting the two-group
  ANOVA reproduces the pooled t-test (`F = t^2`) and the binomial test is exact.

## [1.787.0] - 2026-09-12

### Documentation
- README gains a "Survival analysis (Kaplan-Meier / Nelson-Aalen)" section
  documenting `kaplan_meier`, `nelson_aalen` and `survival_at` with a worked
  right-censored example. TOC regenerated.

## [3.60.0] - 2026-09-14

### Documentation
- README: documented the strong PRNGs in the Quasi-Monte Carlo section, next to the
  variate samplers — `PCG32` and `Xorshift128Plus` with the reference-vector and
  bias-free-randint notes. All snippet values verified live.

## [3.59.0] - 2026-09-14

### Added
- `pcg.py`: statistically strong PRNGs — `PCG32` (O'Neill's PCG-XSH-RR, with unbiased
  `randint`) and `Xorshift128Plus`. Cross-checked: PCG32 reproduces O'Neill's canonical
  reference output vector exactly, both are uniform (mean ~0.5, variance ~1/12), `randint`
  is bias-free over a small range, and both streams are reproducible per seed.

## [3.58.0] - 2026-09-14

### Documentation
- README: documented the variate samplers in the Quasi-Monte Carlo section, next to the
  alias sampler — `sample_normal`, `sample_exponential`, `sample_gamma`, and
  `sample_poisson` with worked moments. All snippet values verified live.

## [3.57.0] - 2026-09-14

### Added
- `samplers.py`: random variate samplers — `sample_normal` (Box-Muller),
  `sample_exponential` (inverse CDF), `sample_gamma` (Marsaglia-Tsang with small-shape
  boost), and `sample_poisson` (Knuth), all on a deterministic seeded stream.
  Cross-checked: the sample mean and variance match each distribution's over 100k draws
  (including gamma shapes below and above 1), and every sampler is reproducible per seed.

## [3.56.0] - 2026-09-14

### Documentation
- README: documented `AliasSampler` in the Quasi-Monte Carlo section, next to Latin
  hypercube sampling — the O(1)-draw alias method with worked frequencies. All snippet
  values verified live.

## [3.55.0] - 2026-09-14

### Added
- `alias_sampler.py`: `AliasSampler` — Walker's alias method for O(1)-per-draw categorical
  sampling after O(n) setup, with a deterministic seeded stream. Cross-checked: empirical
  draw frequencies match the (normalized) weights over 200k samples for skewed, uniform,
  and degenerate distributions, and results are reproducible per seed.

## [3.54.0] - 2026-09-14

### Documentation
- README: documented the Diophantine tools in the Number theory section, after
  combinatorial ranking — `linear_diophantine`, `sqrt_continued_fraction`, and
  `pell_fundamental` with worked values. All snippet values verified live.

## [3.53.0] - 2026-09-14

### Added
- `diophantine.py`: integer-equation tools — `linear_diophantine` (all solutions of
  `ax + by = c`), `sqrt_continued_fraction` (periodic CF of an irrational square root),
  and `pell_fundamental` (smallest solution of `x^2 - n y^2 = 1`). Cross-checked: linear
  solutions satisfy the equation and the homogeneous step preserves it, the sqrt CFs match
  known expansions (and their convergents converge to the root), and the Pell fundamentals
  match published values including the famous n=61 case.

## [3.52.0] - 2026-09-14

### Documentation
- README: documented the prime sieves in the Number theory section, after
  factorization — `primes_up_to`, `prime_count`, `nth_prime`, and
  `smallest_prime_factors` with worked values. All snippet values verified live.

## [3.51.0] - 2026-09-14

### Added
- `sieve.py`: prime sieves — `primes_up_to` (Sieve of Eratosthenes), `prime_count`,
  `nth_prime`, and `smallest_prime_factors`. Cross-checked: the sieve agrees with
  `is_prime`, the counts hit pi(10)=4 / pi(100)=25 / pi(1000)=168, the n-th prime matches
  known values (10000th = 104729), and the smallest-prime-factor table factors 2000 random
  numbers identically to `factorize`.

## [3.50.0] - 2026-09-14

### Documentation
- README: documented the quadratic-residue / multiplicative-group tools in the Number
  theory modular subsection — `legendre_symbol`, `tonelli_shanks`,
  `multiplicative_order`, and `primitive_root` with worked values. All snippet values
  verified live.

## [3.49.0] - 2026-09-14

### Added
- `modular2.py`: quadratic residues and multiplicative structure mod a prime —
  `legendre_symbol`, `jacobi_symbol`, `tonelli_shanks` (modular square root),
  `multiplicative_order`, and `primitive_root`. Cross-checked: the Legendre symbol matches
  a brute residue set, Jacobi agrees with Legendre and is multiplicative, Tonelli-Shanks
  returns a valid square root (and rejects non-residues), orders divide p-1, and a
  primitive root generates the whole multiplicative group.

## [3.48.0] - 2026-09-14

### Documentation
- README: documented combinatorial ranking in the Number theory section, after numerals —
  `gray_code`, `permutation_rank`/`unrank`, and `combination_rank`/`unrank` with worked
  values. All snippet values verified live.

## [3.47.0] - 2026-09-14

### Added
- `combinatorics_rank.py`: combinatorial ranking bijections — `gray_code`/`gray_decode`,
  `permutation_rank`/`permutation_unrank` (Lehmer code), and
  `combination_rank`/`combination_unrank` (combinatorial number system). Cross-checked:
  consecutive Gray codes differ in one bit and round-trip, and the permutation and
  combination rank/unrank pairs match `itertools`' lexicographic order exactly.

## [3.46.0] - 2026-09-14

### Documentation
- README: documented the circle-geometry routines in the Computational geometry section —
  `circle_from_3points`, `circle_line_intersection`, `circle_circle_intersection`, and
  `point_in_circle` with worked values. All snippet values verified live.

## [3.45.0] - 2026-09-14

### Added
- `circle.py`: circle geometry — `circle_from_3points` (circumcircle), `point_in_circle`,
  `circle_line_intersection`, and `circle_circle_intersection`. Cross-checked: the
  circumcircle passes through all three points (collinear input raises), every line/circle
  intersection lies on both the circle and the line, and the tangent/miss/separate cases
  return the right point counts.

## [3.44.0] - 2026-09-14

### Documentation
- README: documented `douglas_peucker` polyline simplification in the Computational
  geometry section — the straight-run and keep-the-peak examples. All snippet values
  verified live.

## [3.43.0] - 2026-09-14

### Added
- `simplify.py`: `douglas_peucker` polyline simplification within a perpendicular
  tolerance. Cross-checked: a straight line collapses to its two endpoints, the endpoints
  are always kept and the result is a subsequence of the input, every dropped point stays
  within epsilon of its retained chord over 500 random polylines, and the point count
  decreases monotonically as epsilon grows.

## [3.42.0] - 2026-09-14

### Documentation
- README: documented the point-proximity distances in the Computational geometry section
  — `point_to_line_distance`, `point_segment_distance`, `closest_point_on_segment`, and
  `point_polyline_distance` with worked values. All snippet values verified live.

## [3.41.0] - 2026-09-14

### Added
- `geometry_dist.py`: point distances — `point_to_line_distance` (perpendicular to an
  infinite line), `closest_point_on_segment`, `point_segment_distance` (clamped), and
  `point_polyline_distance`. Cross-checked: perpendicular distances match known cases,
  the segment distance clamps to endpoints, and the closest point beats a 1000-step
  brute-force scan of the segment over 2000 random cases.

## [3.40.0] - 2026-09-14

### Documentation
- README: documented polygon triangulation in the Computational geometry section —
  `ear_clipping_triangulate`, `signed_area`, `is_clockwise`, and `is_convex_polygon` with
  worked values. All snippet values verified live.

## [3.39.0] - 2026-09-14

### Added
- `triangulate.py`: polygon triangulation and orientation tests — `ear_clipping_triangulate`,
  `signed_area`, `is_clockwise`, and `is_convex_polygon`. Cross-checked: a polygon triangulates
  into exactly `n - 2` triangles whose areas sum to the polygon's own area (square, concave
  L-shape, and 300 random convex polygons), and the orientation/convexity tests agree with the
  shoelace sign.

## [3.38.0] - 2026-09-14

### Documentation
- README: documented `RunningCovariance` in the Numerical utilities section, next to the
  streaming moments — the one-pass covariance/correlation and the exact merge. All snippet
  values verified live.

## [3.37.0] - 2026-09-14

### Added
- `online_cov.py`: `RunningCovariance` — single-pass covariance and Pearson correlation
  via the Welford co-moment update, with exact `+` merge for chunked/parallel data.
  Cross-checked: it matches batch covariance/correlation/variance over 500 random paired
  samples, returns +/-1 for perfectly linear data, and two merged accumulators reproduce
  the full-sample result.

## [3.36.0] - 2026-09-14

### Documentation
- README: documented `smith_waterman` local alignment in the String algorithms section,
  after Needleman-Wunsch — the best-substring-pair framing with the HELLO example. All
  snippet values verified live.

## [3.35.0] - 2026-09-14

### Added
- `alignment.py`: `smith_waterman` local alignment — the best-scoring pair of substrings
  (scores floored at zero, traceback from the peak cell). Cross-checked: it recovers an
  embedded common substring, its score always equals the returned alignment, the aligned
  regions are substrings of the inputs, and no positive common region yields score 0.

## [3.34.0] - 2026-09-14

### Documentation
- README: documented Damerau-Levenshtein and Needleman-Wunsch in the String algorithms
  section, next to the edit-distance functions — with the transposition and alignment
  examples. All snippet values verified live.

## [3.33.0] - 2026-09-14

### Added
- `alignment.py`: `damerau_levenshtein` (true unrestricted distance, adjacent
  transpositions count as one edit) and `needleman_wunsch` (global alignment score with
  traceback). Cross-checked: a transposition is a single edit (Levenshtein would charge
  two), the DL distance never exceeds Levenshtein and matches a brute-force BFS over edit
  operations, and the Needleman-Wunsch score equals the scored alignment while the gapped
  strings reconstruct the inputs.

## [3.32.0] - 2026-09-14

### Documentation
- README: documented the classic sequence algorithms in the String algorithms section —
  `longest_increasing_subsequence`, `maximum_subarray`, and `longest_run` with worked
  values. All snippet values verified live.

## [3.31.0] - 2026-09-14

### Added
- `sequences.py`: classic sequence algorithms — `longest_increasing_subsequence`
  (patience sorting, O(n log n), with reconstruction), `maximum_subarray` (Kadane, with
  bounds), and `longest_run`. Cross-checked: LIS length matches a brute-force DP over 1000
  random arrays (both strict and non-strict, result is a valid increasing subsequence),
  Kadane matches a brute all-slices scan with correct bounds, and the longest run is
  maximal.

## [3.30.0] - 2026-09-14

### Documentation
- README: documented the interval-set operations in the Range-query structures section —
  `merge_intervals`, `total_covered_length`, `intervals_intersection`, and `max_overlap`
  with worked values. All snippet values verified live.

## [3.29.0] - 2026-09-14

### Added
- `interval_set.py`: operations on sets of 1-D intervals — `merge_intervals`,
  `total_covered_length`, `intervals_intersection`, `intervals_union`, and `max_overlap`.
  Cross-checked: the union and intersection match an integer point-set computation over
  500 random collections, merging collapses overlapping and touching ranges, and the
  maximum-overlap sweep matches a brute-force point scan.

## [3.28.0] - 2026-09-14

### Documentation
- README: documented the linear-time selection routines in the Range-query structures
  section — `kth_smallest`, `median`, and `top_k` with worked values. All snippet values
  verified live.

## [3.27.0] - 2026-09-14

### Added
- `selection.py`: linear-time selection — `kth_smallest` (quickselect with the
  median-of-medians pivot), `median`, and `top_k`. Cross-checked: the k-th order statistic
  matches a full sort for every k across hundreds of random arrays (with and without
  duplicates), the median matches `statistics.median`, and top-k matches the sorted
  slices.

## [3.26.0] - 2026-09-14

### Documentation
- README: documented the numeral-system conversions in the Number theory section —
  `to_base`, `from_base`, `to_roman`, and `from_roman` with worked values. All snippet
  values verified live.

## [3.25.0] - 2026-09-14

### Added
- `numeral.py`: numeral-system conversions — `to_base`/`from_base` (any base 2..36) and
  `to_roman`/`from_roman`. Cross-checked: base conversion round-trips and matches Python's
  `int(s, base)` over 5000 random values, and Roman numerals round-trip for every integer
  1..3999 plus the classic subtractive cases.

## [3.24.0] - 2026-09-14

### Documentation
- README: documented the Bernoulli numbers and Faulhaber's formula in the Number theory
  section, after combinatorics — `bernoulli_number`, `faulhaber`, and `bernoulli_sequence`
  with worked values. All snippet values verified live.

## [3.23.0] - 2026-09-14

### Added
- `bernoulli.py`: exact Bernoulli numbers and Faulhaber's sum-of-powers formula —
  `bernoulli_number`, `bernoulli_sequence`, and `faulhaber` (all via `fractions`).
  Cross-checked: the Bernoulli numbers match known values (B_1 = +1/2, odd > 1 vanish),
  Faulhaber's closed form equals the direct power sum for every degree up to 7, and it
  reproduces the classic sum-of-k / k^2 / k^3 identities.

## [3.22.0] - 2026-09-14

### Documentation
- README: added an "Interval arithmetic" section documenting the `Interval` class — the
  enclosure guarantee, arithmetic, and helper methods with worked values. TOC
  regenerated. All snippet values verified live.

## [3.21.0] - 2026-09-14

### Added
- `interval.py`: an `Interval` class for interval arithmetic — `+ - * /`, integer powers,
  `exp`/`log`/`sqrt`, and `width`/`midpoint`/`contains`/`intersect`. Cross-checked: the
  enclosure property holds for thousands of random samples across every operation (the
  true value always lands inside the computed interval), the even-power-straddling-zero
  and sign cases are correct, and division by an interval containing zero raises.

## [3.20.0] - 2026-09-14

### Documentation
- README: documented the signal feature descriptors in the Spectral analysis section —
  `zero_crossing_rate`, `rms`, `crest_factor`, `spectral_centroid`, `spectral_bandwidth`,
  and `spectral_flatness` with worked values. All snippet values verified live.

## [3.19.0] - 2026-09-14

### Added
- `signal_features.py`: time- and frequency-domain signal descriptors —
  `zero_crossing_rate`, `rms`, `crest_factor`, `spectral_centroid`, `spectral_bandwidth`,
  and `spectral_flatness`. Cross-checked: a sine's RMS is 1/sqrt(2) and its crest factor
  sqrt(2), a pure tone's spectral centroid lands at its frequency with a far narrower
  bandwidth than noise, and spectral flatness is near 0 for a tone and high for
  white noise.

## [3.18.0] - 2026-09-14

### Documentation
- README: documented the cross-spectrum and coherence in the Spectral analysis section,
  next to Welch — `coherence` and `cross_spectral_density` with the shared-band and
  independent-noise examples. All snippet values verified live.

## [3.17.0] - 2026-09-14

### Added
- `coherence.py`: Welch cross-spectral density and magnitude-squared coherence —
  `cross_spectral_density` and `coherence`. Cross-checked: identical signals have
  coherence 1 at every frequency, independent noise averages near zero, two linearly
  related signals reach ~1 in their shared band with the cross-spectrum peaking there, and
  all coherence values stay in [0, 1].

## [3.16.0] - 2026-09-14

### Documentation
- README: documented partial and semi-partial correlation in the rank-dependence section,
  next to Kendall's tau — `partial_correlation` and `semipartial_correlation` with the
  spurious-correlation example. All snippet values verified live.

## [3.15.0] - 2026-09-14

### Added
- `partial_corr.py`: `partial_correlation` and `semipartial_correlation` — the
  correlation of two variables after regressing out one or more controls (OLS residuals).
  Cross-checked: a spurious correlation driven by a common cause drops from 0.92 to ~0
  when that cause is controlled, the result matches the textbook three-variable partial
  formula, an uncorrelated control leaves the correlation unchanged, and multiple controls
  work.

## [3.14.0] - 2026-09-14

### Documentation
- README: documented the categorical association measures in the hypothesis-testing
  section, next to the chi-square independence test — `cramers_v`, `phi_coefficient`,
  `tschuprow_t`, and `contingency_coefficient` with worked values. All snippet values
  verified live.

## [3.13.0] - 2026-09-14

### Added
- `association.py`: association strength for categorical contingency tables — `cramers_v`,
  `phi_coefficient`, `tschuprow_t`, and `contingency_coefficient`, all built from the
  Pearson chi-square. Cross-checked: independence gives 0 and perfect association 1, phi
  equals Cramer's V and the closed-form 2x2 phi, the chi-square matches the existing
  independence test, and V equals Tschuprow's T on a square table.

## [3.12.0] - 2026-09-14

### Documentation
- README: documented the Qn/Sn/biweight robust scale estimators in the Robust scale and
  location section — the location-free framing, contamination example, and finite-sample
  caveat. All snippet values verified live.

## [3.11.0] - 2026-09-14

### Added
- `robust_scale.py`: high-breakdown scale estimators — `qn_scale` and `sn_scale`
  (Rousseeuw-Croux, 50% breakdown) and `biweight_midvariance`. Cross-checked: all are
  consistent with the true sigma on clean Gaussian data, stay near the clean scale when
  20% gross outliers destroy the standard deviation, and Sn matches robustbase's
  asymptotic value on 1..10. (The heaviest O(n^2) Gaussian-consistency checks are marked
  slow.)

## [3.10.0] - 2026-09-14

### Documentation
- README: documented the 2-D grid interpolators in the Numerical utilities section, next
  to the 1-D interpolation methods — `bilinear_interp` and `nearest_interp` with worked
  values. All snippet values verified live.

## [3.9.0] - 2026-09-14

### Added
- `grid_interp.py`: 2-D regular-grid interpolation — `bilinear_interp` and
  `nearest_interp`, both clamping out-of-range queries to the grid edge. Cross-checked:
  bilinear is exact at grid nodes, reproduces a planar field everywhere, gives the
  four-corner mean at a cell center, and nearest-neighbour snaps to the closest node.

## [3.8.0] - 2026-09-14

### Documentation
- README: documented the SDE integrators in the Monte Carlo section — `euler_maruyama`,
  `milstein`, and `gbm_paths` with the general-diffusion framing and the strong-order
  contrast. All snippet values verified live.

## [3.7.0] - 2026-09-14

### Added
- `sde.py`: one-dimensional SDE integrators — `euler_maruyama` (strong order 0.5),
  `milstein` (strong order 1.0, with the `b b'` correction), and a `gbm_paths` convenience
  wrapper. Cross-checked: zero diffusion reduces to the deterministic ODE, the GBM
  terminal sample mean approaches the analytic `x0 exp(mu t)`, Milstein's strong error is
  smaller than Euler's against the exact GBM path on shared increments, and paths are
  reproducible per seed.

## [3.6.0] - 2026-09-14

### Documentation
- README: added a "Range-query structures" section documenting `FenwickTree` and
  `SegmentTree` with worked values. TOC regenerated. All snippet values verified live.

## [3.5.0] - 2026-09-14

### Added
- `fenwick.py`: `FenwickTree` (binary indexed tree — O(log n) point update and
  prefix/range sums) and `SegmentTree` (generic associative range query with point
  update; sum by default, min/max via `combine`). Cross-checked: both match a brute-force
  recompute over 200 random arrays with interleaved updates, for sum, minimum, and maximum
  queries.

## [3.4.0] - 2026-09-14

### Documentation
- README: documented the weighted descriptive statistics in the Numerical utilities
  section, next to the streaming moments — `weighted_mean`, `weighted_std`,
  `weighted_median`, and the rest with worked values. All snippet values verified live.

## [3.3.0] - 2026-09-14

### Added
- `weighted_stats.py`: weighted descriptive statistics — `weighted_mean`,
  `weighted_variance`, `weighted_std` (reliability-weight bias correction),
  `weighted_quantile`, and `weighted_median`. Cross-checked: equal weights reproduce the
  ordinary mean/variance/std, integer weights equal the statistics of the replicated
  sample, and a known weighted mean and skewed median come out right.

## [3.2.0] - 2026-09-13

### Documentation
- README: added a "Color spaces" section documenting `rgb_to_hsv`, `rgb_to_hsl`,
  `rgb_to_hex`, `hex_to_rgb`, and their inverses with worked values. TOC regenerated. All
  snippet values verified live.

## [3.1.0] - 2026-09-13

### Added
- `color.py`: color-space conversions — `rgb_to_hsv`, `hsv_to_rgb`, `rgb_to_hsl`,
  `hsl_to_rgb`, `rgb_to_hex`, and `hex_to_rgb`. Cross-checked: the HSV and HSL conversions
  match Python's `colorsys` over 2000 random colors each, all conversions round-trip, the
  primary colors give their known values, and hex round-trips exactly over 500 codes.

## [3.0.0] - 2026-09-13

### Documentation
- README: documented the vector-algebra helpers in the Quaternions section —
  `cross`, `angle_between`, `vector_project`, `reflect`, and the rest with worked values
  and the naming-collision note. All snippet values verified live. Version rolls to 3.0.0.

## [2.99.0] - 2026-09-13

### Added
- `vector3.py`: vector algebra — `dot`, `cross`, `norm`, `normalize`, `angle_between`,
  `vector_project`, `vector_reject`, and `reflect` (the projection helpers are prefixed to
  avoid clashing with the existing PCA `project`). Cross-checked: the cross product is anticommutative
  and perpendicular to both inputs, angles come out at 90/0/180 degrees for the obvious
  cases, projection plus rejection reconstruct the original vector (rejection perpendicular
  to the target), and reflection flips the normal component while preserving length.

## [2.98.0] - 2026-09-13

### Documentation
- README: added a "Quaternions" section documenting `axis_angle_to_quat`,
  `rotate_vector`, `slerp`, `quat_multiply`, and the rest with worked values. TOC
  regenerated. All snippet values verified live.

## [2.97.0] - 2026-09-13

### Added
- `quaternion.py`: unit quaternions for 3-D rotation — `quat_multiply`, `quat_normalize`,
  `quat_conjugate`, `axis_angle_to_quat`, `quat_to_axis_angle`, `rotate_vector`, and
  `slerp`. Cross-checked: known rotations (90 deg about z maps x to y), rotation preserves
  vector length, axis-angle round-trips, composition equals the Hamilton product, and
  slerp hits its endpoints, the 45-degree midpoint, and stays unit-norm throughout.

## [2.96.0] - 2026-09-13

### Documentation
- README: added a "Geodesy" section documenting `haversine_distance`, `initial_bearing`,
  `destination_point`, and `cross_track_distance` with worked values. TOC regenerated. All
  snippet values verified live.

## [2.95.0] - 2026-09-13

### Added
- `geo.py`: great-circle geodesy on a spherical Earth — `haversine_distance`,
  `initial_bearing`, `destination_point`, and `cross_track_distance`. Cross-checked:
  known city-pair distances (London-Paris ~344 km, NYC-LA ~3936 km), distance symmetry,
  a destination/distance/bearing round-trip over 500 random points, cardinal bearings
  (due north 0, due east 90), and a cross-track distance of ~111 km per degree off a
  meridian.

## [2.94.0] - 2026-09-13

### Documentation
- README: documented the von Mises distribution in the Directional statistics section —
  `von_mises_pdf`, `von_mises_fit`, and the Bessel functions with worked values. All
  snippet values verified live.

## [2.93.0] - 2026-09-13

### Added
- `von_mises.py`: the von Mises distribution (circular normal) — `bessel_i0`, `bessel_i1`,
  `von_mises_pdf`, and `von_mises_fit` (maximum-likelihood mean and concentration).
  Cross-checked: the Bessel functions match published values, the density integrates to 1
  for every concentration and reduces to uniform at kappa=0 with its mode at the mean, and
  the fit recovers the parameters of a simulated sample (its kappa satisfies I1/I0 = R).

## [2.92.0] - 2026-09-13

### Documentation
- README: added a "Directional statistics" section documenting `circular_mean`,
  `resultant_length`, `circular_variance`, `circular_std`, and `rayleigh_test` with worked
  values. TOC regenerated. All snippet values verified live.

## [2.91.0] - 2026-09-13

### Added
- `circular_stats.py`: directional statistics for angles — `circular_mean`,
  `resultant_length`, `circular_variance`, `circular_std`, and `rayleigh_test`.
  Cross-checked: the mean wraps correctly (350 deg & 10 deg -> 0), the resultant length is
  1 for identical angles and ~0 for a uniform spread, cancelling vectors raise, and the
  Rayleigh test rejects a concentrated sample (p ~ 0) while accepting a uniform one.

## [2.90.0] - 2026-09-13

### Documentation
- README: documented the Hamming(7,4) code and Luhn checksum in the Data compression
  section — the single-bit correction and check-digit examples, next to the checksums.
  All snippet values verified live.

## [2.89.0] - 2026-09-13

### Added
- `hamming_code.py`: error-detecting/correcting codes — `hamming74_encode`,
  `hamming74_decode` (single-bit correction), `luhn_checksum`, and `luhn_check_digit`.
  Cross-checked: all 16 Hamming(7,4) words round-trip and a flip at every one of the 7
  positions is both corrected and located; the Luhn checksum validates known card numbers
  and the computed check digit makes any payload valid over 1000 random cases.

## [2.88.0] - 2026-09-13

### Documentation
- README: documented `find_all_roots` and `count_sign_changes` in the Numerical utilities
  section, next to the nonlinear-system solvers — the grid-scan approach, the cubic/sine
  examples, and the even-multiplicity caveat. All snippet values verified live.

## [2.87.0] - 2026-09-13

### Added
- `root_scan.py`: find all sign-changing roots of a function on an interval —
  `find_all_roots` (grid scan for sign changes, then Brent refinement) and
  `count_sign_changes`. Cross-checked: recovers a cubic's three roots, the multiples of
  pi where sine vanishes, cosine's six roots on [0, 20], and two transcendental roots of
  e^x = 3x, all to residuals near machine precision.

## [2.86.0] - 2026-09-13

### Documentation
- README: documented the `hungarian` assignment and `knapsack_01` solvers in the
  optimization section with worked values. All snippet values verified live.

## [2.85.0] - 2026-09-13

### Added
- `assignment.py`: combinatorial optimization — `hungarian` (minimum-cost assignment,
  O(n^3) Kuhn-Munkres) and `knapsack_01` (0/1 knapsack DP with a chosen-item trace).
  Cross-checked: Hungarian matches a brute-force permutation search over 300 random
  matrices and returns a valid permutation, and the knapsack value and item trace match a
  brute-force subset search over 500 random instances.

## [2.84.0] - 2026-09-13

### Documentation
- README: documented the `linprog` two-phase simplex solver in the optimization section,
  next to `nnls` — the constraint format, the max/min examples, and the phase-1 / Bland's
  rule notes. All snippet values verified live.

## [2.83.0] - 2026-09-13

### Added
- `linprog.py`: linear programming by two-phase simplex — `linprog` handles `<=`, `>=`,
  and `=` constraints (with a phase-1 artificial-variable stage) and Bland's anti-cycling
  rule. Cross-checked against known optima and a brute-force vertex enumeration; a diet-
  style minimum satisfies every constraint; and unbounded / infeasible programs raise.

## [2.82.0] - 2026-09-13

### Documentation
- README: documented the general eigenvalue solver in the Matrix utilities section, next
  to power iteration — `characteristic_polynomial`, `eigenvalues_general`, and
  `determinant_from_charpoly` with worked values including a complex pair. All snippet
  values verified live.

## [2.81.0] - 2026-09-13

### Added
- `eigen_general.py`: general (non-symmetric) eigenvalues via Faddeev-LeVerrier —
  `characteristic_polynomial`, `eigenvalues_general` (roots the char-poly, so complex
  conjugate pairs come back correctly), and `determinant_from_charpoly`. Cross-checked:
  the char-poly matches known cases, a rotation matrix yields the pair ±i, eigenvalues
  agree with `jacobi_eigen` on symmetric matrices, their sum equals the trace and their
  product equals the determinant, and the char-poly determinant matches the LU one.

## [2.80.0] - 2026-09-13

### Documentation
- README: documented the nonlinear-system solvers in the Numerical utilities section —
  `newton_system` and `broyden` with the circle-line example, next to the scalar root
  finders. All snippet values verified live.

## [2.79.0] - 2026-09-13

### Added
- `newton_system.py`: root finding for nonlinear systems — `newton_system`
  (finite-difference Jacobian + linear solve each step) and `broyden` (rank-1
  quasi-Newton inverse-Jacobian updates). Cross-checked: both solve a circle-line
  intersection, a transcendental system, and a 3-variable polynomial system to residuals
  near machine precision, Newton and Broyden agree, an at-root input returns in zero
  iterations, and a singular Jacobian raises.

## [2.78.0] - 2026-09-13

### Documentation
- README: documented the checksums and hashes in the Data compression section — `crc32`,
  `adler32`, and `fnv1a_32` with worked values and the zlib-match note. All snippet
  values verified live.

## [2.77.0] - 2026-09-13

### Added
- `checksums.py`: checksums and non-cryptographic hashes — `crc32` (IEEE 802.3),
  `adler32`, and `fnv1a_32`. Cross-checked: CRC-32 and Adler-32 match Python's `zlib`
  over 1000 random inputs each, FNV-1a reproduces its published test vectors, and str and
  bytes inputs hash identically.

## [2.76.0] - 2026-09-13

### Documentation
- README: documented LZW compression and delta coding in the Data compression section —
  the on-the-fly dictionary, the repetitive-text ratio, and the slowly-varying delta
  example. All snippet values verified live.

## [2.75.0] - 2026-09-13

### Added
- `lzw.py`: LZW dictionary compression and integer delta coding — `lzw_compress`,
  `lzw_decompress`, `delta_encode`, and `delta_decode`. Cross-checked: LZW round-trips
  exactly (including the code-equals-next-code special case and 500 random strings) and
  emits fewer codes than characters on repetitive text; delta coding round-trips and
  turns a slowly-varying series into small differences.

## [2.74.0] - 2026-09-13

### Documentation
- README: documented the Burrows-Wheeler transform and move-to-front coding in the Data
  compression section — the bzip2-pipeline framing, the `banana` example, and the
  clustering-to-zeros behavior. All snippet values verified live.

## [2.73.0] - 2026-09-13

### Added
- `bwt.py`: Burrows-Wheeler transform and move-to-front coding (bzip2 building blocks) —
  `bwt_transform`, `bwt_inverse`, `move_to_front_encode`, and `move_to_front_decode`.
  Cross-checked: the BWT round-trips exactly (known and 500 random strings, `banana` ->
  `nnbaaa`), move-to-front round-trips and turns clustered input into mostly zeros, and
  the full BWT+MTF pipeline inverts back to the original.

## [2.72.0] - 2026-09-13

### Documentation
- README: added a "Data compression" section documenting `huffman_encode`,
  `huffman_decode`, `run_length_encode`, and `run_length_decode` with worked values. TOC
  regenerated. All snippet values verified live.

## [2.71.0] - 2026-09-13

### Added
- `compression.py`: lossless compression primitives — `huffman_codebook`,
  `huffman_encode`, `huffman_decode`, `run_length_encode`, and `run_length_decode`.
  Cross-checked: Huffman round-trips exactly, the codes are prefix-free, the expected
  code length lands within `[H, H+1)` of the entropy and beats fixed-length coding,
  frequent symbols get the shortest codes, and RLE round-trips.

## [2.70.0] - 2026-09-13

### Documentation
- README: documented the fuzzy string-similarity scores in the String algorithms
  section — `jaro`, `jaro_winkler`, `dice_coefficient`, and `jaccard_similarity` with
  worked values. All snippet values verified live.

## [2.69.0] - 2026-09-13

### Added
- `fuzzy_match.py`: fuzzy string similarity — `jaro`, `jaro_winkler`,
  `dice_coefficient` (character bigrams), and `jaccard_similarity` (token/character
  sets). Cross-checked: Jaro and Jaro-Winkler reproduce Winkler's published reference
  values (MARTHA/MARHTA, DWAYNE/DUANE, DIXON/DICKSONX), Jaro is symmetric and bounded in
  [0, 1] over random pairs, the Winkler prefix boost never lowers the score, and Dice /
  Jaccard match known values.

## [2.68.0] - 2026-09-13

### Documentation
- README: added a "String algorithms" section documenting `levenshtein`,
  `hamming_distance`, `longest_common_subsequence`, `longest_common_substring`, and
  `kmp_search` with worked values. TOC regenerated. All snippet values verified live.

## [2.67.0] - 2026-09-13

### Added
- `strings.py`: string algorithms — `levenshtein`, `hamming_distance`,
  `longest_common_subsequence`, `longest_common_substring`, and `kmp_search`
  (Knuth-Morris-Pratt). Cross-checked: the edit distance matches a brute-force recursion
  and is symmetric, the LCS is a subsequence of both inputs, and KMP matches a naive scan
  (including overlapping matches) across hundreds of random cases.

## [2.66.0] - 2026-09-13

### Documentation
- README: documented Bellman-Ford, Floyd-Warshall, and A* in the Graph algorithms
  section — negative edges, all-pairs distances, and heuristic search with worked values.
  All snippet values verified live.

## [2.65.0] - 2026-09-13

### Added
- `graph4.py`: shortest paths beyond Dijkstra — `bellman_ford` (negative edges + negative
  cycle detection), `floyd_warshall` (all-pairs), and `a_star` (heuristic point-to-point
  search). Cross-checked: Bellman-Ford handles a negative edge, detects a negative cycle,
  and agrees with Dijkstra on non-negative graphs; Floyd-Warshall matches per-source
  Dijkstra over all pairs; and A* finds the optimal grid path (matching Dijkstra) and
  reduces to Dijkstra with a zero heuristic.

## [2.64.0] - 2026-09-13

### Documentation
- README: documented the centrality measures in the Graph algorithms section —
  `pagerank`, `degree_centrality`, `closeness_centrality`, and `betweenness_centrality`
  with worked star/path values. All snippet values verified live.

## [2.63.0] - 2026-09-13

### Added
- `graph3.py`: graph centrality measures — `pagerank` (power iteration with teleportation
  and dangling-node handling), `degree_centrality`, `closeness_centrality`, and
  `betweenness_centrality` (Brandes). Cross-checked: PageRank sums to 1 and is uniform on
  a symmetric ring, and the centralities match the known star/path structures (the star
  center dominates degree/closeness/betweenness; the path midpoint has the highest
  betweenness).

## [2.62.0] - 2026-09-13

### Documentation
- README: documented union-find, MST, and max flow in the Graph algorithms section —
  `UnionFind`, `minimum_spanning_tree`, and `max_flow` with worked values. All snippet
  values verified live.

## [2.61.0] - 2026-09-13

### Added
- `graph2.py`: weighted-graph structure algorithms — `UnionFind` (path compression +
  union by rank), `minimum_spanning_tree` (Kruskal), and `max_flow` (Edmonds-Karp).
  Cross-checked: the MST weight matches a brute-force minimum over all spanning trees and
  the tree is spanning and acyclic, and the maximum flow reproduces the CLRS classic
  network (23) plus series-bottleneck and parallel-path cases.

## [2.60.0] - 2026-09-13

### Documentation
- README: added a "Graph algorithms" section documenting `dijkstra`, `shortest_path`,
  `bfs`, `connected_components`, and `topological_sort` with worked values and the
  adjacency-dict conventions. TOC regenerated. All snippet values verified live.

## [2.59.0] - 2026-09-13

### Added
- `graph.py`: graph algorithms — `dijkstra` and `shortest_path` (non-negative weights),
  `bfs` (unweighted hop distances), `connected_components` (undirected), and
  `topological_sort` (Kahn's algorithm with cycle detection). Cross-checked: Dijkstra
  matches Bellman-Ford on 100 random weighted graphs, path reconstruction and BFS hops
  are correct, components partition an undirected graph, the topological order respects
  every edge, and cycles / negative weights / bad sources raise.

## [2.58.0] - 2026-09-13

### Documentation
- README: documented the point-set extremal measures in the Computational geometry
  section — `bounding_box`, `polygon_diameter`, and `min_enclosing_circle` with worked
  values. All snippet values verified live.

## [2.57.0] - 2026-09-13

### Added
- `geometry3.py`: extremal measures of a point set — `bounding_box`, `polygon_diameter`
  (farthest pair, via the convex hull), and `min_enclosing_circle` (Welzl's algorithm).
  Cross-checked: the diameter matches a brute-force farthest-pair scan, the minimum
  enclosing circle contains every point yet cannot be shrunk 1% without excluding one,
  and it reproduces the known circle for two points and the square's circumcircle.

## [2.56.0] - 2026-09-13

### Documentation
- README: documented segment intersection and polygon clipping in the Computational
  geometry section — the crossing/boolean tests, perimeter, and Sutherland-Hodgman
  window clipping. All snippet values verified live.

## [2.55.0] - 2026-09-13

### Added
- `geometry2.py`: segment intersection and convex polygon clipping — `segments_intersect`,
  `segment_intersection`, `polygon_perimeter`, and `clip_polygon` (Sutherland-Hodgman).
  Cross-checked: crossing segments meet at the expected point, parallel/disjoint return
  no intersection, shared-endpoint and T-junctions are detected, the perimeter matches
  known polygons, and clipping a square by a smaller window (or by a diagonal
  half-plane, or a partial overlap) yields the correct clipped area.

## [2.54.0] - 2026-09-13

### Documentation
- README: added a "Computational geometry" section documenting `convex_hull`,
  `polygon_area`, `polygon_centroid`, `point_in_polygon`, and `closest_pair` with worked
  values. TOC regenerated. All snippet values verified live.

## [2.53.0] - 2026-09-13

### Added
- `geometry.py`: planar computational geometry — `convex_hull` (Andrew's monotone
  chain), `polygon_area` and `polygon_centroid` (shoelace), `point_in_polygon`
  (ray casting), and `closest_pair` (divide-and-conquer). Cross-checked: the hull of a
  square plus interior points is its four corners and every point of a random cloud lies
  inside its own hull; areas and centroids match known values orientation-independently;
  point-in-polygon handles a concave L-shape and on-edge points; and closest_pair matches
  a brute-force search.

## [2.52.0] - 2026-09-13

### Documentation
- README: documented the exact integer/rational linear algebra in the Matrix utilities
  section, next to LU — the fraction-free determinant, exact solve/inverse, and the
  Hilbert-matrix contrast with floating solvers. All snippet values verified live.

## [2.51.0] - 2026-09-13

### Added
- `bareiss.py`: exact integer/rational linear algebra — `bareiss_determinant`
  (fraction-free elimination, exact integer determinant), `rational_solve`, and
  `rational_inverse` (exact `Fraction` results). Cross-checked: the determinant matches a
  cofactor expansion over 200 random matrices, the solver's residual `A x - b` is exactly
  zero, an inverse times its matrix is the exact identity, a Hilbert system (which floats
  botch) solves to exact ones, and singular matrices raise.

## [2.50.0] - 2026-09-13

### Documentation
- README: documented Latin hypercube sampling in the Quasi-Monte Carlo section — the
  stratification, the discrepancy comparison against random sampling, and the maximin
  variant. All snippet values verified live.

## [2.49.0] - 2026-09-13

### Added
- `lhs.py`: Latin hypercube sampling and uniformity diagnostics — `latin_hypercube`
  (stratified, with a centered variant), `maximin_lhs` (most-spread over several tries),
  and `l2_star_discrepancy` (Warnock's formula). Cross-checked: every axis has exactly
  one sample per bin, centered points sit at bin centers, an LHS design has lower L2 star
  discrepancy than plain random sampling on average, and the maximin search improves the
  minimum inter-point distance.

## [2.48.0] - 2026-09-13

### Documentation
- README: documented the combinatorics functions in the Number theory section —
  `binomial`, `multinomial`, `stirling_second`, `bell`, `catalan`, `partition_count`,
  and `derangements` with worked values. All snippet values verified live.

## [2.47.0] - 2026-09-13

### Added
- `combinatorics.py`: exact integer combinatorics — `binomial`, `multinomial`,
  `stirling_second`, `bell`, `catalan`, `partition_count`, and `derangements`.
  Cross-checked: the binomial row sums to 2^n, the multinomial gives 34650 for
  "mississippi", Bell numbers equal the summed Stirling row, the Catalan and
  integer-partition sequences match their known values, and derangements equal
  round(n!/e).

## [2.46.0] - 2026-09-13

### Documentation
- README: documented the modular-arithmetic routines in the Number theory section —
  `extended_gcd`, `mod_inverse`, `chinese_remainder`, `mod_pow`, and `discrete_log` with
  worked values. All snippet values verified live.

## [2.45.0] - 2026-09-13

### Added
- `modular.py`: modular arithmetic — `extended_gcd`, `mod_inverse`,
  `chinese_remainder`, `mod_pow` (negative exponents via inverse), and `discrete_log`
  (baby-step giant-step). Cross-checked: the extended GCD satisfies Bezout's identity and
  matches `math.gcd`, the inverse gives `a*inv == 1 (mod m)`, CRT reconstructs a known
  value (and the classic (2,3,2) mod (3,5,7) = 23), and the discrete log inverts
  exponentiation on a prime field (3^x = 13 mod 17 -> 4).

## [2.44.0] - 2026-09-13

### Documentation
- README: added a "Number theory" section documenting `is_prime`, `factorize`,
  `euler_totient`, `divisors`, `gcd`, and `lcm` — the deterministic-primality guarantee,
  the Mersenne/Carmichael/semiprime examples, and the derived identities. TOC regenerated.
  All snippet values verified live.

## [2.43.0] - 2026-09-13

### Added
- `number_theory.py`: integer number theory — `gcd`, `lcm`, `is_prime` (deterministic
  Miller-Rabin, exact past 64-bit), `factorize` (trial division + Pollard's rho),
  `divisors`, and `euler_totient`. Cross-checked: primality matches brute force for all
  n < 2001 and handles the Mersenne prime 2^61-1 and the Carmichael number 561;
  factorization reconstructs the input (including a ~10^18 semiprime) with certified-prime
  factors; the totient matches a brute-force coprime count; and gcd*lcm == a*b over
  1000 random pairs.

## [2.42.0] - 2026-09-13

### Documentation
- README: documented continued-fraction expansion and best rational approximation next
  to Padé — the number-vs-function contrast, pi/sqrt(2) expansions and convergents, and
  the `best_rational` bounds. All snippet values verified live.

## [2.41.0] - 2026-09-13

### Added
- `continued_fraction.py`: continued-fraction expansion and best rational approximation
  — `cf_expansion`, `convergents`, and `best_rational`. Cross-checked: the expansions of
  415/93, pi, and sqrt(2) match their known coefficients; pi's convergents are
  22/7, 333/106, 355/113; each convergent beats every fraction with a smaller
  denominator; and `best_rational` matches Python's `Fraction.limit_denominator` on 3000
  random values across denominator bounds.

## [2.40.0] - 2026-09-13

### Documentation
- README: documented matched filtering and peak detection in the spectral section, after
  cross-correlation — the optimal-detector framing, the peak-finding filters, the
  multi-detection example, and the short-template caveat. All snippet values verified
  live.

## [2.39.0] - 2026-09-13

### Added
- `matched_filter.py`: matched filtering and peak detection — `matched_filter` (max-SNR
  detection statistic for a known template), `normalized_matched_filter` (correlation
  coefficient in [-1, 1]), `find_peaks` (local maxima with height/distance filtering),
  and `detect_template` (template occurrences above a correlation threshold).
  Cross-checked: the matched filter peaks at the planted template location in noise, the
  normalized response is exactly 1 at a scaled match, peak finding respects height and
  minimum-distance constraints, and a distinctive template is detected at all its planted
  offsets.

## [2.38.0] - 2026-09-13

### Documentation
- README: documented total-variation denoising in the spectral section, after the
  order-statistic filters — the objective, the `lam` behavior, the noisy-step recovery
  numbers, and the exact-minimizer note. All snippet values verified live.

## [2.37.0] - 2026-09-13

### Added
- `tv_denoise.py`: total-variation denoising via Condat's exact O(n) 1-D algorithm —
  `tv_denoise` and `tv_total_variation`. Cross-checked: `lam=0` returns the input, a
  large `lam` collapses to the mean, a noisy step is recovered with sharp edges and much
  lower total variation, the output is the true minimizer (no coordinate perturbation
  lowers the objective across thousands of trials), and it matches an independent
  subgradient-descent solver.

## [2.36.0] - 2026-09-13

### Documentation
- README: documented the nonlinear order-statistic filters in the spectral section,
  after the sample-rate subsection — the spike-removal / edge-preservation contrast with
  linear filters, the rank generalization, and the surgical Hampel behavior. All snippet
  values verified live.

## [2.35.0] - 2026-09-13

### Added
- `median_filter.py`: nonlinear order-statistic filters — `median_filter` (sliding
  median), `rank_filter` (any percentile, incl. min/max), and `hampel_filter` (MAD-based
  outlier replacement). Cross-checked: the median filter removes an impulsive spike,
  preserves a step edge that a moving average would smear, and matches a brute-force
  centered median in the interior; the rank filter gives correct min/median/max; and the
  Hampel filter flags only the injected outlier while leaving clean data untouched.

## [2.34.0] - 2026-09-13

### Documentation
- README: documented the Butterworth IIR filters in the spectral section, between the
  FIR and sample-rate subsections — the FIR/IIR trade-off, the biquad-cascade design,
  the -3 dB / stopband landmarks, and the stability note. All snippet values verified
  live.

## [2.33.0] - 2026-09-13

### Added
- `iir_filter.py`: Butterworth IIR filter design via the bilinear transform as a biquad
  cascade — `butter_lowpass`, `butter_highpass`, `sosfilt` (Direct Form II transposed
  cascade), and `iir_frequency_response`. Cross-checked: the low-pass has unit DC gain,
  exactly -3 dB (0.7071) at the cutoff, a monotone maximally-flat rolloff, and a deep
  stopband; the high-pass is its mirror; higher order rolls off faster; every section is
  stable (poles inside the unit circle) for both even and odd orders; and in the time
  domain the filter cleanly separates a two-tone signal.

## [2.32.0] - 2026-09-13

### Documentation
- README: documented sample-rate conversion in the spectral section, next to the FIR
  filters — the aliasing pitfall, sinc reconstruction, up/downsampling, the anti-alias
  suppression example, and the group-delay note. All snippet values verified live.

## [2.31.0] - 2026-09-13

### Added
- `sample_rate.py`: sample-rate conversion — `sinc_interp` (Whittaker-Shannon
  reconstruction at arbitrary points), `upsample` and `downsample` (integer factors with
  anti-imaging / anti-alias FIR filters), and `resample_rational` (rational `up/down`
  ratio). Cross-checked: sinc interpolation is exact at integer positions and
  reconstructs a band-limited tone between them, up/downsampling preserve a low tone
  (accounting for the filter's group delay), and downsampling a tone above the new
  Nyquist suppresses it (rms 0.001) instead of aliasing it back.

## [2.30.0] - 2026-09-13

### Documentation
- README: documented the parametric AR spectral estimator in the spectral section,
  next to Welch — the contrast with the periodogram, the AR(2) recovery and short-record
  two-sinusoid resolution, and the Burg/Yule-Walker options. All snippet values verified
  live.

## [2.29.0] - 2026-09-13

### Added
- `ar_spectrum.py`: parametric (autoregressive) power-spectral-density estimation —
  `burg` (Burg's forward-backward method, best for short records), `ar_psd` (PSD from
  AR coefficients), and `ar_spectrum` (fit + evaluate, Burg or Yule-Walker).
  Cross-checked: Burg recovers a known AR(2) process's coefficients, the reflection
  coefficients stay stable (|k| < 1), the spectrum of two sinusoids in noise peaks at
  both frequencies, and a positive/negative AR(1) coefficient gives a low-pass/high-pass
  shape.

## [2.28.0] - 2026-09-13

### Documentation
- README: documented the STFT and spectrogram in the spectral section, next to Welch —
  the time-frequency contrast, the frame-count/hop relationship, the tone and chirp
  behavior, and the overlap-add reconstruction. All snippet values verified live.

## [2.27.0] - 2026-09-13

### Added
- `stft.py`: short-time Fourier transform, spectrogram, and inverse STFT for
  time-frequency analysis — `stft` (per-frame windowed FFTs), `spectrogram`
  (`|STFT|^2` over the non-redundant half), and `istft` (weighted overlap-add
  reconstruction). Cross-checked: `istft(stft(x))` reconstructs the signal to ~1e-14 in
  the interior, a pure tone concentrates in its own frequency bin, a rising chirp's
  spectrogram peak climbs from early to late frames, and the frame count follows the hop.

## [2.26.0] - 2026-09-13

### Documentation
- README: documented the cepstrum and quefrency analysis in the spectral section —
  the log-spectrum-inverse interpretation, the echo-delay and impulse-train pitch
  examples, and the role of `min_quefrency`. All snippet values verified live.

## [2.25.0] - 2026-09-13

### Added
- `cepstrum.py`: real cepstrum, power cepstrum, and quefrency analysis for echo and
  pitch detection — `real_cepstrum`, `power_cepstrum`, and `fundamental_quefrency`.
  Cross-checked: a signal with an echo delayed by D samples produces a cepstral peak
  exactly at quefrency D, a periodic impulse train (the voiced-speech model) peaks at
  its period, the power cepstrum is non-negative, and the quefrency search band is
  respected.

## [2.24.0] - 2026-09-13

### Documentation
- README: documented the cross-correlation / lead-lag tools in the spectral section,
  next to convolution and autocorrelation — the raw and normalized forms, the
  peak-lag detector, and the known-delay recovery example. All snippet values verified
  against a live run.

## [2.23.0] - 2026-09-13

### Added
- `cross_correlation.py`: FFT cross-correlation and lead-lag detection between two
  series — `cross_correlation` (raw, over positive and negative lags),
  `normalized_cross_correlation` (correlation coefficient in [-1, 1]), and
  `lag_at_max_correlation` (the lag of best alignment). Cross-checked: the FFT result
  matches a direct brute-force cross-correlation sum to ~1e-14, a known delay between
  two series is recovered exactly, identical series peak at lag 0 with coefficient 1,
  and the normalized values stay within [-1, 1].

## [2.22.0] - 2026-09-13

### Documentation
- README: documented the Hilbert transform and analytic signal in the spectral
  section — the phasor interpretation, the FFT one-sided construction, the AM-envelope
  recovery and flat instantaneous-frequency examples, and the power-of-two constraint.
  All snippet values verified against a live run.

## [2.21.0] - 2026-09-13

### Added
- `hilbert.py`: Hilbert transform and the analytic signal via the FFT one-sided
  construction — `analytic_signal`, `hilbert_transform`, `envelope`,
  `instantaneous_phase` (unwrapped), and `instantaneous_frequency`. Cross-checked: the
  analytic real part reproduces the input, the Hilbert transform of a cosine equals the
  sine to ~1e-14, the envelope of a pure tone is its amplitude and the envelope of an
  AM signal recovers the modulating waveform exactly, and the instantaneous frequency
  of a single tone is flat at that tone's frequency.

## [2.20.0] - 2026-09-13

### Documentation
- README: documented the windowed-sinc FIR filters in the spectral section —
  Nyquist-normalized cutoffs, linear-phase symmetry, unit DC gain, the low/high/band
  design, and the two-tone separation example. All snippet values verified against a
  live run.

## [2.19.0] - 2026-09-13

### Added
- `fir_filter.py`: windowed-sinc FIR filter design — `fir_lowpass`, `fir_highpass`,
  `fir_bandpass` (Nyquist-normalized cutoffs, windowed) and `fir_apply` (convolution).
  Cross-checked: the lowpass has unit DC gain and linear-phase symmetry, passbands pass
  and stopbands block by orders of magnitude, the highpass has zero DC gain, the
  bandpass passes only its middle band, and a two-tone signal keeps the low tone while
  the high one is removed.

## [2.18.0] - 2026-09-13

### Documentation
- README: documented the window functions (`hann`/`hamming`/`blackman`/`bartlett`/
  `rectangular`) and `apply_window` in the spectral-analysis section.

## [2.17.0] - 2026-09-13

### Added
- `windows.py`: spectral window functions `hann`, `hamming`, `blackman`, `bartlett`,
  `rectangular`, and the `apply_window` helper. Cross-checked: correct endpoint/center
  values, all symmetric, and applying a Hann window before the FFT suppresses the
  spectral leakage of an off-bin sinusoid by orders of magnitude.

## [2.16.0] - 2026-09-13

### Documentation
- README: documented `goertzel` and `goertzel_power` in the spectral-analysis section.
  Every snippet value verified.

## [2.15.0] - 2026-09-13

### Added
- `goertzel.py`: the Goertzel single-frequency DFT — `goertzel` (the complex `X[k]`
  coefficient in `O(n)`) and `goertzel_power` (`|X[k]|^2` for tone detection).
  Cross-checked: matches the FFT and direct DFT at every bin, power equals the
  coefficient magnitude squared, a pure tone peaks at its bin, and the DC bin equals the
  sum. Caught a phase error in the termination formula (magnitude was right but the
  complex value was off by `e^{jw}`); fixed to match the FFT coefficient exactly.

## [2.14.0] - 2026-09-13

### Documentation
- README: documented `dct` and `idct` in the spectral-analysis section, alongside the
  FFT.

## [2.13.0] - 2026-09-13

### Added
- `dct.py`: the orthonormal discrete cosine transform `dct` (DCT-II) and its inverse
  `idct` (DCT-III). Cross-checked: `idct(dct(x)) == x` to machine precision, energy is
  preserved (Parseval), a constant maps to a pure DC coefficient, a smooth signal packs
  >95% of its energy into the first four coefficients (compaction), a single cosine mode
  isolates one coefficient, and the transform is linear.

## [2.12.0] - 2026-09-13

### Documentation
- README: documented `poisson2d` in the numerical-utilities section, completing the
  parabolic/hyperbolic/elliptic PDE trio.

## [2.11.0] - 2026-09-13

### Added
- `poisson2d.py`: `poisson2d`, a 2-D Poisson/Laplace solver (`u_xx + u_yy = f`) on a
  rectangle with Dirichlet boundaries by successive over-relaxation. Cross-checked: a
  linear (harmonic) boundary is recovered exactly, a `sin*sin` source matches its
  analytic solution to `O(h^2)`, the maximum principle holds (interior stays within the
  boundary range), and SOR converges about 7x faster than plain Gauss-Seidel.

## [2.10.0] - 2026-09-13

### Documentation
- README: documented `wave_equation` in the numerical-utilities section, after the heat
  equation.

## [2.9.0] - 2026-09-13

### Added
- `wave_equation.py`: `wave_equation`, a 1-D wave-equation solver (`u_tt = c^2 u_xx`) by
  the explicit central-difference (leapfrog) scheme with fixed ends, initial shape and
  velocity. Cross-checked: a standing wave matches the d'Alembert solution, returns to
  its start after a full period, a localized pulse splits into two half-height
  travelling waves, and the CFL condition `c dt/dx <= 1` is enforced.

## [2.8.0] - 2026-09-13

### Documentation
- README: documented `heat_equation_cn` in the numerical-utilities section, after the
  ODE/symplectic solvers.

## [2.7.0] - 2026-09-13

### Added
- `heat_equation.py`: `heat_equation_cn`, a general 1-D heat/diffusion solver
  (`u_t = alpha u_xx`) by Crank-Nicolson with Dirichlet boundaries — second-order in
  space and time and unconditionally stable. Cross-checked: a sine mode decays to its
  analytic `exp(-alpha (pi/L)^2 t)` profile, boundaries are held, the steady state is
  the exact linear profile, a Gaussian bump diffuses, and huge time steps stay stable.

## [2.6.0] - 2026-09-13

### Documentation
- README: documented `velocity_verlet` and `leapfrog` in the numerical-utilities
  section, after the ODE/BVP solvers.

## [2.5.0] - 2026-09-13

### Added
- `symplectic.py`: `velocity_verlet` and `leapfrog`, symplectic integrators for
  Hamiltonian systems that keep the total energy *bounded* over long runs (unlike a
  drifting general RK). Cross-checked: harmonic-oscillator energy oscillates within a
  tiny band rather than drifting, one period returns to the start, leapfrog equals
  velocity Verlet at unit mass, a circular orbit keeps its radius over 10k steps, and
  the scheme is time-reversible.

## [2.4.0] - 2026-09-13

### Documentation
- README: documented `shooting_bvp` in the numerical-utilities section, after the ODE
  solvers. Every snippet value verified.

## [2.3.0] - 2026-09-13

### Added
- `bvp.py`: `shooting_bvp`, the shooting method for two-point boundary-value problems
  (`y'' = f(t, y, y')` with fixed `y(a)` and `y(b)`) — guesses the initial slope,
  integrates with adaptive RK45, and Brent-root-finds the slope that hits the far
  boundary. Cross-checked: recovers the analytic slope and solution for `y''=y`
  (sinh), `y''=-y` (sine), `y''=6t` (cubic, every computed point on `t^3`), and a
  problem with a first-derivative term.

## [2.2.0] - 2026-09-13

### Documentation
- README: documented `rk4` and `rk45` in the numerical-utilities section. Every snippet
  value verified.

## [2.1.0] - 2026-09-13

### Added
- `ode.py`: initial-value ODE solvers — `rk4` (fixed-step fourth-order Runge-Kutta) and
  `rk45` (adaptive Dormand-Prince with an embedded error estimate). Handle scalar or
  vector systems. Cross-checked: recover `exp(t)` and `exp(-t)` to ~1e-8, integrate the
  harmonic-oscillator system back to `[1, 0]` at `2*pi`, match the logistic ODE's closed
  form, and the scalar and vector paths agree.

## [2.0.1] - 2026-09-13

### Documentation
- README: documented `richardson_extrapolate` and `richardson_table` in the
  numerical-utilities section. Every snippet value verified.

## [2.0.0] - 2026-09-13

### Added
- `richardson.py`: general Richardson extrapolation — `richardson_extrapolate` (accelerate
  a sequence of step-halved estimates to the `h -> 0` limit for any leading error order
  `p`) and `richardson_table` (the full convergence tableau). Cross-checked: a first-order
  forward difference and second-order central difference extrapolate to the exact
  derivative, the trapezoid sequence extrapolates to the Romberg integral, and it is
  exact on the assumed error model.

## [1.999.0] - 2026-09-13

### Documentation
- README: documented `laplace_inversion` in the numerical-utilities section. Every
  snippet value verified.

## [1.998.0] - 2026-09-13

### Added
- `laplace_inversion.py`: `laplace_inversion`, real-arithmetic numerical inversion of a
  Laplace transform by the Gaver-Stehfest algorithm. Cross-checked against known
  transform pairs (constant, ramp, `t^2`, decaying exponential, `sqrt(t)`, cosine),
  recovering the time function to ~1e-3 or better on smooth, non-oscillatory targets.

## [1.997.0] - 2026-09-13

### Documentation
- README: documented `bivariate_normal_cdf` and `trivariate_normal_cdf` in the
  special-functions section. Every snippet value verified.

## [1.996.0] - 2026-09-13

### Added
- `multivariate_normal_cdf.py`: public `bivariate_normal_cdf` (Drezner-Wesolowsky, the
  routine behind the American/compound models) and `trivariate_normal_cdf` (Genz
  reduction to a 1-D integral of the bivariate CDF). Cross-checked: the bivariate matches
  the orthant formula `1/4 + asin(rho)/(2pi)` to ~1e-9, both reduce to the product form at
  zero correlation, the trivariate collapses to the bivariate when the third bound is
  infinite, is monotone in correlation, and is permutation-symmetric.

## [1.995.0] - 2026-09-13

### Documentation
- README: documented `integrate2d_gauss` and `integrate2d_simpson` in the
  numerical-utilities quadrature section. Every snippet value verified.

## [1.994.0] - 2026-09-13

### Added
- `integrate2d.py`: two-dimensional integration over a rectangle — `integrate2d_gauss`
  (tensor Gauss-Legendre) and `integrate2d_simpson` (composite Simpson via Fubini).
  Cross-checked: exact on separable polynomials up to the Gauss degree, recovers the
  2-D Gaussian integral (pi) and sin*cos analytically, and the two rules agree to
  machine precision on a smooth polynomial.

## [1.993.0] - 2026-09-13

### Documentation
- README: documented `dagostino_k2` in the goodness-of-fit section.

## [1.992.0] - 2026-09-13

### Added
- `dagostino.py`: `dagostino_k2`, the D'Agostino-Pearson K^2 omnibus normality test
  (skewness Z + kurtosis Z, chi-square with 2 df). Cross-checked: null rejection rate
  ~0.05, rejects skewed (exponential) data via `z_skew` and heavy tails via `z_kurt`,
  gives a negative kurtosis Z for light-tailed uniform data, passes clean normals, and
  has high power against exponential samples.

## [1.991.0] - 2026-09-13

### Documentation
- README: documented `levene_test` and `bartlett_test` in the hypothesis-test section.

## [1.990.0] - 2026-09-13

### Added
- `variance_tests.py`: `levene_test` (Levene / median-centered Brown-Forsythe, F-test on
  absolute deviations) and `bartlett_test` (likelihood-ratio, chi-square) for equal
  variance across groups. Cross-checked: non-significant on equal-variance groups,
  strongly significant on 1x/3x/6x spreads, both `center` variants work, and the null
  rejection rate is calibrated.

## [1.989.0] - 2026-09-13

### Documentation
- README: documented `anderson_darling_ksample` in the goodness-of-fit section. Every
  snippet value verified.

## [1.988.0] - 2026-09-13

### Added
- `anderson_darling_ksample.py`: `anderson_darling_ksample`, the Scholz-Stephens
  k-sample Anderson-Darling test for a common distribution (tie-corrected, with the
  standardized statistic and an interpolated p-value). Cross-checked: reproduces the
  Scholz-Stephens worked example (A2akN ~ 8.36, p < 0.01), is non-significant under a
  shared distribution, flags differing means, and is calibrated under the null.

## [1.987.0] - 2026-09-13

### Documentation
- README: documented `permutation_test` and `paired_permutation_test` in the
  bootstrap/resampling section.

## [1.986.0] - 2026-09-13

### Added
- `permutation_test.py`: `permutation_test` (two-sample, any statistic, label shuffling)
  and `paired_permutation_test` (sign-flip on within-pair differences), both seeded and
  deterministic. Cross-checked: detects a mean or median difference, is non-significant
  under the null, calibrated rejection rate, custom statistics work, and the paired form
  flags a consistent within-pair shift.

## [1.985.0] - 2026-09-13

### Documentation
- README: documented `ecdf`, `quantile` and `qq_points` in the robust-statistics
  section. Every snippet value verified.

## [1.984.0] - 2026-09-13

### Added
- `ecdf.py`: the empirical CDF (`ecdf`), sample `quantile` (linear/lower/higher/nearest
  interpolation), and `qq_points` for quantile-quantile comparison. Cross-checked: the
  ECDF runs 0->1 monotonically, the linear quantile matches `statistics.quantiles`,
  quantiles recover the normal median and 97.5% point, and a Q-Q pairing sits on the
  diagonal for matched samples and shows the right slope under a scale difference.

## [1.983.0] - 2026-09-13

### Documentation
- README: documented `wasserstein_distance` / `wasserstein1_sorted` in the
  information-theory section. Every snippet value verified.

## [1.982.0] - 2026-09-13

### Added
- `wasserstein.py`: 1-D Wasserstein (earth-mover) distance — `wasserstein_distance`
  (any sizes, any `p`, via merged empirical CDFs) and `wasserstein1_sorted` (fast
  equal-length sorted-difference form). Cross-checked: a constant shift gives exactly
  that distance, identical samples give zero, point masses at 0 and 1 give 1, it is
  symmetric, `W2 >= W1`, and a uniform shift of 2 is recovered on large samples.

## [1.981.0] - 2026-09-13

### Documentation
- README: documented `kl_divergence`, `jensen_shannon_divergence`, `hellinger_distance`,
  `total_variation_distance` and `bhattacharyya_distance` in the information-theory
  section. Every snippet value verified.

## [1.980.0] - 2026-09-13

### Added
- `divergences.py`: distances between discrete distributions — `kl_divergence`,
  `jensen_shannon_divergence`, `hellinger_distance`, `total_variation_distance` and
  `bhattacharyya_distance` (all on auto-normalized weight vectors). Cross-checked:
  identical distributions give zero, KL is asymmetric while JS is symmetric and bounded
  by log 2, Hellinger/TV stay in [0, 1], disjoint supports give TV = 1 and infinite
  Bhattacharyya, and the Pinsker inequality (TV <= sqrt(KL/2)) holds.

## [1.979.0] - 2026-09-13

### Documentation
- README: documented `poisson_regression` / `poisson_predict` alongside logistic
  regression in the ML section.

## [1.978.0] - 2026-09-13

### Added
- `poisson_regression.py`: `poisson_regression` (Poisson GLM with a log link, fit by
  Fisher-scoring IRLS) and `poisson_predict`. Cross-checked: recovers a known
  log-linear rate from sampled counts, the intercept-only fit reproduces log of the
  mean count, predictions are always positive, and the full model's log-likelihood
  beats the null.

## [1.977.0] - 2026-09-13

### Documentation
- README: documented `bayesian_linear_regression` and `bayesian_predict` in the
  OLS-regression section, alongside ridge.

## [1.976.0] - 2026-09-13

### Added
- `bayesian_regression.py`: `bayesian_linear_regression` (conjugate Gaussian posterior
  over the coefficients) and `bayesian_predict` (predictive mean and variance).
  Cross-checked: a weak prior approaches OLS, the posterior mean equals ridge with a
  penalized intercept, posterior std shrinks with more data, a stronger prior shrinks
  the coefficients, and the predictive variance widens on extrapolation.

## [1.975.0] - 2026-09-13

### Documentation
- README: documented `lowess` alongside the Savitzky-Golay filter in the smoothing
  section.

## [1.974.0] - 2026-09-13

### Added
- `lowess.py`: `lowess`, locally-weighted scatterplot smoothing (tricube kernel, local
  linear fits, Cleveland robustifying iterations). Cross-checked: reproduces a straight
  line, smooths a noisy sine below the raw error, a larger `frac` is smoother, it keeps
  input order on unsorted data, and the robust iterations pin an outlier's neighbours to
  the true trend (caught a degenerate median-scale case that let a lone outlier drag the
  fit).

## [1.973.0] - 2026-09-13

### Documentation
- README: documented `ransac_line` in the robust-regression section. Every snippet value
  verified.

## [1.972.0] - 2026-09-13

### Added
- `ransac.py`: `ransac_line`, RANSAC robust line fitting — minimal 2-point samples,
  largest inlier consensus set, then a least-squares refit on the inliers. Cross-checked:
  recovers the true line with 60% of the data corrupted (beyond any median estimator's
  breakdown), is exact on clean data, deterministic for a fixed seed, and locks onto the
  majority structure amid gross contamination.

## [1.971.0] - 2026-09-13

### Documentation
- README: documented `huber_regression` in the robust-regression section. Every snippet
  value verified.

## [1.970.0] - 2026-09-13

### Added
- `huber_regression.py`: `huber_regression`, a robust M-estimator by IRLS — Huber loss
  (quadratic near zero, linear beyond `delta`) with a MAD-based residual scale.
  Cross-checked: matches OLS on clean data, resists vertical outliers that drag the OLS
  slope (1.997 vs 1.65 at the true 2), approaches OLS exactly as `delta` grows, and
  recovers an exact line.

## [1.969.0] - 2026-09-13

### Documentation
- README: documented `repeated_median_regression` in the robust-regression section.
  Every snippet value verified.

## [1.968.0] - 2026-09-13

### Added
- `siegel_regression.py`: Siegel's `repeated_median_regression`, a 50%-breakdown robust
  line fit (median of per-point median slopes). Cross-checked: exact on a clean line,
  recovers the true slope/intercept with 40% of the data corrupted (past Theil-Sen's
  ~29% limit), matches Theil-Sen on clean noisy data, and is unmoved by a single gross
  outlier.

## [1.967.0] - 2026-09-13

### Documentation
- README: documented `dbscan` alongside k-means in the clustering section.

## [1.966.0] - 2026-09-13

### Added
- `dbscan.py`: `dbscan` density-based clustering — no preset cluster count, finds
  arbitrarily-shaped clusters, and labels low-density points as noise (`-1`).
  Cross-checked: separates two blobs with zero noise, flags a far outlier as noise,
  clusters two non-convex concentric rings (which k-means cannot), collapses to all-noise
  at tiny `eps` and one cluster at huge `eps`, and is deterministic.

## [1.965.0] - 2026-09-13

### Documentation
- README: documented `fit_gradient_boost` / `predict_gradient_boost` alongside the
  regression tree in the ML section.

## [1.964.0] - 2026-09-13

### Added
- `gradient_boost.py`: gradient-boosted regression trees (`fit_gradient_boost`,
  `predict_gradient_boost`) — additive shallow trees fit to residuals with a shrinkage
  learning rate. Cross-checked: beats a single shallow tree, training error falls
  monotonically with the estimator count, a zero learning rate returns the mean, and it
  generalizes on held-out data.

## [1.963.0] - 2026-09-13

### Documentation
- README: documented `fit_regression_tree` / `predict_regression_tree` alongside the
  classification tree and random forest. Every snippet value verified.

## [1.962.0] - 2026-09-13

### Added
- `regression_tree.py`: CART regression tree (`fit_regression_tree`,
  `predict_regression_tree`) — recursive squared-error-reduction splits with mean-valued
  leaves. Cross-checked: recovers a step function exactly, collapses a constant target to
  a single leaf, lowers MSE monotonically with depth on a parabola, and picks the
  relevant feature in a two-feature problem.

## [1.961.0] - 2026-09-13

### Documentation
- README: documented `kde`, `kde_function` and the bandwidth rules in the
  robust-statistics section.

## [1.960.0] - 2026-09-13

### Added
- `kde.py`: Gaussian kernel density estimation — `kde` (pointwise), `kde_function`
  (callable estimator), and the `silverman_bandwidth` / `scott_bandwidth` rules.
  Cross-checked: the estimate integrates to one and is non-negative, recovers the
  N(0,1) and N(5,2) densities at their peaks, and resolves a bimodal sample with both
  modes above the valley.

## [1.959.0] - 2026-09-13

### Documentation
- README: documented `logsumexp`, `softmax` and `log_softmax` in the numerical-utilities
  section. Every snippet value verified.

## [1.958.0] - 2026-09-13

### Added
- `logsumexp.py`: numerically stable `logsumexp` (max-shifted, optional weights),
  `softmax` and `log_softmax`. Cross-checked: matches the naive formula on small
  inputs, never overflows on ~1000-scale values (where the naive version returns inf)
  or underflows on very negative ones, softmax sums to one and stays stable at large
  scale, and log-softmax is consistent with log(softmax).

## [1.957.0] - 2026-09-13

### Documentation
- README: documented `hyperdual_derivatives` and `second_derivative` in the
  numerical-utilities section. Every snippet value verified.

## [1.956.0] - 2026-09-13

### Added
- `hyperdual.py`: hyperdual numbers for exact first *and* second derivatives in one
  evaluation — a `HyperDual` type plus `hyperdual_derivatives` (returns
  `(f, f', f'')`) and `second_derivative`. Cross-checked: first and second derivatives
  of elementary and composite functions are exact, the second derivative agrees with
  `ridders_second_derivative`, and the first derivative matches the plain dual number.

## [1.955.0] - 2026-09-13

### Documentation
- README: documented `dual_gradient` and `dual_newton` in the numerical-utilities
  section. Every snippet value verified.

## [1.954.0] - 2026-09-13

### Added
- `dual_calculus.py`: `dual_gradient` (exact multivariate gradient via forward-mode
  autodiff) and `dual_newton` (Newton's method that gets `f'` from autodiff, so no
  hand-coded derivative is needed). Cross-checked: the gradient matches the analytic
  and numerical gradients including transcendentals, and `dual_newton` finds sqrt(2)
  and ln(2) to machine precision, agreeing with the library's `newton`.

## [1.953.0] - 2026-09-13

### Documentation
- README: documented `dual_derivative` and the `Dual` type in the numerical-utilities
  section. Every snippet value verified.

## [1.952.0] - 2026-09-13

### Added
- `dual.py`: forward-mode automatic differentiation with a `Dual` number type and the
  `dual_derivative` helper (plus dual-aware `exp`/`log`/`sqrt`/`sin`/`cos`/`tan`/`tanh`).
  Cross-checked: exact (zero-error) derivatives of elementary and composite functions,
  correct product/quotient/chain rules, a dual exponent (`x**x`), and agreement with the
  complex-step derivative to machine precision.

## [1.951.0] - 2026-09-13

### Documentation
- README: documented `CountMinSketch` and `BloomFilter` in the numerical-utilities
  section.

## [1.950.0] - 2026-09-13

### Added
- `count_min.py`: `CountMinSketch` (streaming frequency estimation, one-sided
  overestimate) and `BloomFilter` (approximate membership with no false negatives),
  both on a stable SHA-1 double hash. Cross-checked: Count-Min never underestimates and
  nails a heavy hitter, and the Bloom filter has zero false negatives with a
  false-positive rate near its target.

## [1.949.0] - 2026-09-13

### Documentation
- README: documented `HyperLogLog` in the numerical-utilities section.

## [1.948.0] - 2026-09-13

### Added
- `hyperloglog.py`: `HyperLogLog`, a streaming distinct-count (cardinality) estimator
  using `2^p` small registers and a stable SHA-1 hash, with mergeable sketches.
  Cross-checked: relative error under 5% from 100 to 100,000 distinct items, duplicates
  don't inflate the count, `merge` gives the union cardinality, small-range linear
  counting is near-exact, and it is deterministic across runs.

## [1.947.0] - 2026-09-13

### Documentation
- README: documented `EWMAStats` and `ewma` in the numerical-utilities section. Every
  snippet value verified.

## [1.946.0] - 2026-09-13

### Added
- `ewma.py`: `EWMAStats`, a streaming exponentially-weighted mean/variance/volatility
  accumulator (the single-series RiskMetrics recursion), plus the `ewma` batch helper.
  Cross-checked: constant series gives zero variance, the batch output matches the
  manual recursion, it tracks a step change with the expected lag, recovers the vol of
  N(0,2), and reacts sharply to a volatility regime change.

## [1.945.0] - 2026-09-13

### Documentation
- README: documented `P2Quantile` and `reservoir_sample` in the numerical-utilities
  section.

## [1.944.0] - 2026-09-13

### Added
- `streaming_quantile.py`: `P2Quantile` (Jain-Chlamtac P-square one-pass quantile
  estimator, O(1) memory) and `reservoir_sample` (Vitter uniform sampling from a stream
  of unknown length). Cross-checked: P-square tracks the uniform median, the true 95th
  percentile to 0.01, and the normal 90th percentile; the reservoir gives uniform
  coverage and is deterministic. Caught an LCG low-bit bias -- `state % (i+1)` gave a
  badly skewed reservoir (the LCG's lowest bit merely alternates); switched to a
  high-bit float index.

## [1.943.0] - 2026-09-13

### Documentation
- README: documented `RunningMoments` (streaming, mergeable moments) in the
  numerical-utilities section. Every snippet value verified.

## [1.942.0] - 2026-09-13

### Added
- `running_moments.py`: `RunningMoments`, a one-pass accumulator for mean, variance,
  skewness and excess kurtosis via Terriberry's stable central-moment recurrences, with
  a `+` operator that merges two accumulators exactly (parallel reduction).
  Cross-checked: all four moments match the batch statistics, splitting a sample and
  merging the parts reproduces the full accumulator (M2/M3/M4 to ~1e-10), a 10-way merge
  equals a single pass, symmetric data has zero skew and normal data ~0 excess kurtosis.

## [1.941.0] - 2026-09-13

### Documentation
- README: documented `neumaier_sum`, `kahan_sum`, `accurate_dot` and `welford` in the
  numerical-utilities section. Every snippet value verified.

## [1.940.0] - 2026-09-13

### Added
- `compensated.py`: numerically stable summation and statistics — `kahan_sum`,
  `neumaier_sum`, `accurate_dot` (compensated dot product) and `welford` (one-pass
  stable mean/variance). Cross-checked: Neumaier handles catastrophic cancellation
  (`[1, 1e100, 1, -1e100] -> 2`), matches `math.fsum` on scaled sequences, Kahan beats
  the naive running sum of a million `0.1`s, and Welford recovers the variance of
  large-mean data where the textbook `E[x^2]-E[x]^2` collapses to zero.

## [1.939.0] - 2026-09-13

### Documentation
- README: documented `power_iteration`, `inverse_iteration` and `rayleigh_quotient` in
  the matrix-utilities section. Every snippet value verified.

## [1.938.0] - 2026-09-13

### Added
- `power_iteration.py`: individual-eigenpair methods — `power_iteration` (dominant
  eigenvalue/vector), `inverse_iteration` (eigenvalue nearest a shift, e.g. the
  smallest), and `rayleigh_quotient`. Cross-checked against the full `jacobi_eigen`
  spectrum: power iteration recovers the largest eigenvalue and inverse iteration the
  smallest and any shift-targeted one, the eigenvector residual `||A v - lam v||` is
  ~1e-5, and the Rayleigh quotient is exact on a true eigenvector.

## [1.937.0] - 2026-09-13

### Documentation
- README: documented `conjugate_gradient`, `gauss_seidel` and `jacobi` in the
  matrix-utilities section. Every snippet value verified.

## [1.936.0] - 2026-09-13

### Added
- `conjugate_gradient.py`: iterative linear solvers — `conjugate_gradient` (for
  symmetric positive-definite systems), `gauss_seidel` and `jacobi` (stationary
  iterations for diagonally dominant systems). Cross-checked: CG matches the direct
  solve in at most n iterations, resolves an identity system in one step, Gauss-Seidel
  converges in about half the iterations of Jacobi, all three agree, and a 30-D SPD
  system is solved to ~1e-11.

## [1.935.0] - 2026-09-13

### Documentation
- README: documented `bfgs` alongside the other multivariate optimizers in the
  numerical-utilities section, framed as the smooth-objective / local-polish choice.
  Every snippet value verified.

## [1.934.0] - 2026-09-13

### Added
- `bfgs.py`: `bfgs`, a quasi-Newton minimizer building an inverse-Hessian approximation
  from successive numerical gradients, with an Armijo backtracking line search.
  Cross-checked: converges on a quadratic bowl in a couple of iterations (superlinear),
  reaches the Rosenbrock minimum, matches the exact linear solve on a 5-D quadratic to
  1e-6, and drives a smooth convex objective's gradient to zero.

## [1.933.0] - 2026-09-13

### Documentation
- README: documented `simulated_annealing` alongside the other multivariate optimizers
  in the numerical-utilities section, with guidance on choosing it versus differential
  evolution. Snippet value verified.

## [1.932.0] - 2026-09-13

### Added
- `simulated_annealing.py`: `simulated_annealing`, a single-point global optimizer that
  accepts uphill moves by the Metropolis criterion with geometric cooling, reproducible
  via a seeded LCG. Cross-checked: minimizes sphere and Rosenbrock, climbs out of the
  multimodal Rastrigin's local minima where Nelder-Mead stalls, is deterministic for a
  fixed seed, and respects the box bounds.

## [1.931.0] - 2026-09-13

### Documentation
- README: documented `nnls` (non-negative least squares) in the numerical-utilities
  section. Every snippet value verified.

## [1.930.0] - 2026-09-13

### Added
- `nnls.py`: `nnls`, non-negative least squares by the Lawson-Hanson active-set
  method (`min ||A x - b||^2` subject to `x >= 0`). Cross-checked: reproduces the exact
  fit when the unconstrained solution is already non-negative, clamps an otherwise-
  negative coefficient to zero, every solution is non-negative over 200 random cases,
  and the KKT optimality conditions hold (gradient ~0 on the active coefficients,
  <=0 on the zeroed ones) on 200 more.

## [1.929.0] - 2026-09-13

### Documentation
- README: documented `nelder_mead` and `differential_evolution` in the
  numerical-utilities section. Every snippet value verified.

## [1.928.0] - 2026-09-13

### Added
- `differential_evolution.py`: `differential_evolution`, a derivative-free global
  optimizer (Storn-Price), reproducible via a seeded LCG. Also exposes the previously
  internal `nelder_mead` local minimizer publicly. Cross-checked: finds the global
  minimum of sphere / Rosenbrock / Rastrigin / Ackley, beats Nelder-Mead on the
  multimodal Rastrigin (where the local method gets stuck), is deterministic for a
  fixed seed, and respects the box bounds.

## [1.927.0] - 2026-09-13

### Documentation
- README: documented `principal_components_regression` in the OLS-regression section.

## [1.926.0] - 2026-09-13

### Added
- `pcr.py`: `principal_components_regression` — regress on the top-``k`` principal
  components of the predictors and map the fit back, a stable alternative to OLS under
  collinearity. Cross-checked: retaining all components reproduces OLS coefficients and
  intercept exactly, the retained explained-variance fraction is monotone in ``k`` and
  reaches 1, and a single predictor matches OLS.

## [1.925.0] - 2026-09-13

### Documentation
- README: documented `RecursiveLeastSquares` and `recursive_least_squares` in the
  OLS-regression section. Every snippet value verified.

## [1.924.0] - 2026-09-13

### Added
- `rls.py`: `RecursiveLeastSquares` (online OLS updated one observation at a time via
  the Sherman-Morrison identity, with an optional forgetting factor) and the
  `recursive_least_squares` batch wrapper. Cross-checked: with `forgetting = 1` the
  estimate matches batch OLS, it recovers an exact line through noise-free points,
  `predict` works, and a forgetting factor < 1 tracks a mid-stream slope regime change.

## [1.923.0] - 2026-09-13

### Documentation
- README: documented `weighted_least_squares` and `generalized_least_squares` in the
  OLS-regression section.

## [1.922.0] - 2026-09-13

### Added
- `wls.py`: `weighted_least_squares` (minimize weighted squared residuals) and
  `generalized_least_squares` (known error covariance, solved by Cholesky whitening).
  Cross-checked: equal weights reproduce OLS coefficients and standard errors exactly,
  a diagonal GLS covariance equals WLS with `1/diag` weights, an identity covariance
  reduces to OLS, and GLS with an AR(1) covariance recovers the true slope.

## [1.921.0] - 2026-09-13

### Documentation
- README: documented `white_hc0` and `newey_west` robust OLS standard errors in the
  OLS-regression section.

## [1.920.0] - 2026-09-13

### Added
- `ols_hac.py`: robust OLS standard errors via the sandwich estimator — `white_hc0`
  (heteroskedasticity-consistent) and `newey_west` (Bartlett-weighted HAC, valid under
  autocorrelation too). Cross-checked: coefficients match plain OLS, `newey_west` with
  0 lags equals White exactly, the White SE grows under heteroskedastic errors where
  the OLS SE is wrong, and the Newey-West SE grows further under AR(1) errors.

## [1.919.0] - 2026-09-13

### Documentation
- README: documented `thiele_interpolate` / `thiele_coefficients` / `thiele_eval`
  (rational interpolation) in the numerical-utilities section. Every snippet value
  verified.

## [1.918.0] - 2026-09-13

### Added
- `thiele.py`: Thiele's continued-fraction rational interpolation
  (`thiele_coefficients`, `thiele_eval`, `thiele_interpolate`) via reciprocal
  differences. Cross-checked: passes through every node, recovers a known rational
  function `(2x+1)/(x^2+1)` to machine precision, and interpolates `tan` (which has
  poles a polynomial can't model). Caught an off-by-one index bug in the
  continued-fraction evaluation before shipping.

## [1.917.0] - 2026-09-13

### Documentation
- README: documented the `poly_*` dense-polynomial-algebra helpers alongside
  `polynomial_roots` in the numerical-utilities section. Every snippet value verified.

## [1.916.0] - 2026-09-13

### Added
- `polynomial.py`: dense polynomial arithmetic on coefficient lists — `poly_add`,
  `poly_sub`, `poly_mul`, `poly_divmod`, `poly_derivative`, `poly_integral`,
  `poly_eval` and `poly_gcd`. Cross-checked: the product matches the FFT convolution,
  the division identity `num = q*den + r` holds over 500 random pairs, derivative and
  integral invert, and `gcd(p, p')` recovers a repeated root. Caught a non-terminating
  loop in division by a constant (float round-off left a tiny leading residue that the
  exact-zero trim never dropped); forced the cancelled leading term to zero.

## [1.915.0] - 2026-09-13

### Documentation
- README: documented `polynomial_roots` alongside the scalar root finders in the
  numerical-utilities section. Every snippet value verified.

## [1.914.0] - 2026-09-13

### Added
- `polyroots.py`: `polynomial_roots`, an all-roots solver by the Durand-Kerner
  (Weierstrass) simultaneous iteration — finds every real and complex root at once, no
  deflation. Cross-checked: exact on quadratics/cubics, complex conjugate pairs, the
  5th roots of unity and a double root, recovers the roots of 200 random polynomials,
  and every returned root drives the polynomial to ~1e-15.

## [1.913.0] - 2026-09-13

### Documentation
- README: documented `barycentric_weights`, `barycentric_eval`, `chebyshev_nodes` and
  `chebyshev_barycentric_weights` in the numerical-utilities section, including the
  Runge-phenomenon comparison. Every snippet value verified.

## [1.912.0] - 2026-09-13

### Added
- `barycentric.py`: stable, reusable barycentric Lagrange interpolation
  (`barycentric_weights`, `barycentric_eval`) plus `chebyshev_nodes` and the
  closed-form `chebyshev_barycentric_weights`. Cross-checked: matches Neville, exact at
  nodes, the closed-form Chebyshev weights agree with the general ones, and on the
  Runge function Chebyshev nodes converge (max error ~0.02) where equispaced nodes blow
  up (~59) — with exp interpolated to machine precision on 25 Chebyshev points.

## [1.911.0] - 2026-09-13

### Documentation
- README: documented `neville`, `divided_differences` and `newton_polynomial` in the
  numerical-utilities section, including the Richardson-extrapolation use. Every
  snippet value verified.

## [1.910.0] - 2026-09-13

### Added
- `polyinterp.py`: `neville` (polynomial interpolation at a point with an error
  estimate), `divided_differences` and `newton_polynomial` (Newton form).
  Cross-checked: exact recovery of a cubic and 200 random polynomials, Neville equals
  the Newton form, the leading divided difference equals the leading coefficient, an
  honest error estimate, and Richardson extrapolation via Neville recovers a
  derivative to 1e-13.

## [1.909.0] - 2026-09-13

### Documentation
- README: documented `pade`, `pade_eval` and `lentz_continued_fraction` in the
  numerical-utilities section. Every snippet value verified.

## [1.908.0] - 2026-09-13

### Added
- `pade.py`: `pade` (Pade [m/n] approximant from Taylor coefficients), `pade_eval`,
  and `lentz_continued_fraction` (general modified-Lentz evaluator). Cross-checked:
  the exp [2/2] coefficients are exact and the approximant beats the degree-4 Taylor
  series, its own series matches the input coefficients through order m+n, a [4/4]
  approximant of ln(1+x) converges at x=2 where the Taylor series diverges, and Lentz
  reproduces tan(1) and the golden ratio.

## [1.907.0] - 2026-09-13

### Documentation
- README: documented `wynn_epsilon` and `euler_transform` in the numerical-utilities
  section. Every snippet value verified.

## [1.906.0] - 2026-09-13

### Added
- `series_transform.py`: `wynn_epsilon` (Wynn's epsilon algorithm — Shanks iterated to
  all orders) and `euler_transform` (alternating-series acceleration). Cross-checked:
  Wynn takes the Leibniz-pi partial sums to machine precision and is exact on a
  geometric series (picking the best even-column estimate to dodge post-convergence
  round-off), and the Euler transform recovers pi/4 and ln 2 to 1e-6 / 1e-9.

## [1.905.0] - 2026-09-13

### Documentation
- README: documented `aitken`, `shanks` and `steffensen` in the numerical-utilities
  section. Every snippet value verified.

## [1.904.0] - 2026-09-13

### Added
- `sequence_accel.py`: `aitken` (delta-squared acceleration), `shanks` (Shanks
  transform) and `steffensen` (derivative-free quadratic fixed-point solver).
  Cross-checked: Aitken cuts the Leibniz-pi partial-sum error by >100x and is exact on
  a geometric series, Steffensen finds the cos fixed point and sqrt(2) to machine
  precision in a handful of steps (vs ~69 for plain iteration), and Shanks matches
  Aitken.

## [1.903.0] - 2026-09-13

### Documentation
- README: documented `complex_step_derivative` / `complex_step_gradient` in the
  numerical-utilities section, including the `h = 1e-100` no-cancellation
  demonstration. Every snippet value verified.

## [1.902.0] - 2026-09-13

### Added
- `complex_step.py`: `complex_step_derivative` and `complex_step_gradient`, first
  derivatives by the complex-step method (`Im(f(x + i h)) / h`) with no subtractive
  cancellation. Cross-checked: derivatives of sin/exp/log/power/sqrt to machine
  precision, accuracy preserved even at `h = 1e-100` (where a finite difference would
  be pure round-off), agreement with the Ridders derivative, and an exact gradient.

## [1.901.0] - 2026-09-13

### Documentation
- README: documented `ridders_derivative` / `ridders_second_derivative` in the
  numerical-utilities section, including the numerical-vs-analytic Black-Scholes delta
  cross-check. Every snippet value verified.

## [1.900.0] - 2026-09-13

### Added
- `richardson_derivative.py`: `ridders_derivative` and `ridders_second_derivative`,
  high-accuracy numerical derivatives by Ridders' adaptive Richardson extrapolation
  with an error estimate. Cross-checked: first and second derivatives of sin/exp/log/
  power/reciprocal to near machine precision, an honest error estimate, and the
  numerical derivative of the Black-Scholes price matching the closed-form delta.

## [1.899.0] - 2026-09-13

### Documentation
- README: extended the agreement section with `cohen_kappa`, `weighted_kappa` and
  `fleiss_kappa` for categorical inter-rater agreement. Every snippet value verified.

## [1.898.0] - 2026-09-13

### Added
- `cohen_kappa.py`: chance-corrected categorical agreement — `cohen_kappa` (two
  raters, nominal), `weighted_kappa` (ordinal, linear/quadratic distance weights) and
  `fleiss_kappa` (m raters). Cross-checked: perfect agreement gives 1 and chance-level
  gives ~0, Cohen matches a hand computation (0.6154), a two-category linear weighted
  kappa equals Cohen's exactly, weighting rewards near-misses over wild disagreements,
  and Fleiss reproduces its textbook value (0.2519).

## [1.897.0] - 2026-09-13

### Documentation
- README: extended the robust-regression / agreement section with the `icc`
  intraclass correlation family for multi-rater reliability. Every snippet value
  verified.

## [1.896.0] - 2026-09-13

### Added
- `icc.py`: the `icc` intraclass correlation coefficients (Shrout-Fleiss ICC(1),
  ICC(2,1)/(2,k), ICC(3,1)/(3,k)) from a two-way ANOVA variance decomposition.
  Cross-checked: reproduces the canonical Shrout-Fleiss (1979) table values, gives 1
  under perfect agreement and ~0 under noise, satisfies the Spearman-Brown relation
  between single and average-rater forms, and consistency (ICC3) exceeds absolute
  agreement (ICC2) under a rater offset.

## [1.895.0] - 2026-09-13

### Documentation
- README: extended the robust-regression section with `bland_altman` and
  `concordance_correlation` for method-agreement analysis. Every snippet value
  verified.

## [1.894.0] - 2026-09-13

### Added
- `bland_altman.py`: `bland_altman` agreement analysis (bias, SD and 95% limits of
  agreement plus per-point means/diffs) and `concordance_correlation` (Lin's CCC).
  Cross-checked: identical methods give zero bias and CCC 1, a constant offset drops
  CCC below the Pearson value (0.5 vs 1), the CCC matches a reference formula on 300
  random cases and never exceeds |Pearson|, and about 95% of differences fall inside
  the limits of agreement.

## [1.893.0] - 2026-09-13

### Documentation
- README: extended the robust-regression section with `passing_bablok_regression`,
  the distribution-free method-comparison alternative to Deming. Every snippet value
  verified.

## [1.892.0] - 2026-09-13

### Added
- `passing_bablok.py`: `passing_bablok_regression`, the robust nonparametric
  method-comparison fit (shifted median of pairwise slopes with the K offset and the
  slope = -1 exclusion). Cross-checked: exact on a noiseless line, unmoved by a gross
  outlier, recovers the true slope on noisy data, exact on a decreasing line, and
  symmetric (reciprocal slope) fitting x on y; no index errors over 1000 random
  inputs.

## [1.891.0] - 2026-09-13

### Documentation
- README: extended the robust-regression section with `deming_regression` and
  `orthogonal_regression` for errors-in-variables fitting. Every snippet value
  verified.

## [1.890.0] - 2026-09-13

### Added
- `deming.py`: `deming_regression` (errors-in-variables line fitting with an error-
  variance ratio `lambda`) and `orthogonal_regression` (total least squares, the
  `lambda = 1` case). Cross-checked: exact on a noiseless line for every lambda,
  `lambda -> infinity` recovers the OLS slope, orthogonal regression equals the PCA
  first-eigenvector slope and is symmetric (`slope_xy * slope_yx = 1`), and the fit
  passes through the sample means.

## [1.889.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `fisher_exact_test` for small-count
  2x2 tables. Every snippet value verified.

## [1.888.0] - 2026-09-13

### Added
- `fisher_exact.py`: `fisher_exact_test` for a 2x2 table (two-sided / greater / less)
  via the hypergeometric distribution, returning the exact p-value and the sample odds
  ratio. Cross-checked against the classic tea-tasting table (OR = 9, two-sided
  p = 0.4857, one-sided 0.2429) and a known scipy case, plus the hypergeometric
  probabilities summing to one and p ~ 1 under independence.

## [1.887.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `mcnemar_test` and
  `cochran_q_test` for paired / repeated binary outcomes. Every snippet value
  verified.

## [1.886.0] - 2026-09-13

### Added
- `cochran_mcnemar.py`: `mcnemar_test` (paired 2x2, exact two-sided binomial plus the
  continuity-corrected chi-square) and `cochran_q_test` (k binary treatments over
  shared blocks). Cross-checked: McNemar's exact p-value matches the binomial CDF and
  the corrected chi-square matches hand computation, Cochran's Q on two treatments
  equals the uncorrected McNemar statistic exactly, and Q is calibrated under the
  null.

## [1.885.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `jonckheere_terpstra_test` for an
  ordered-trend alternative to Kruskal-Wallis. Every snippet value verified.

## [1.884.0] - 2026-09-13

### Added
- `jonckheere.py`: the `jonckheere_terpstra_test` for a monotone trend across ordered
  groups (a directional, more powerful alternative to Kruskal-Wallis when the groups
  have a natural order). Cross-checked: the null mean and tie-corrected variance match
  the closed-form formulas, a perfectly increasing arrangement gives the maximal
  statistic and a large positive z (decreasing gives the mirror-image negative z), the
  null z is calibrated (mean ~0, variance ~1 over 2000 replications), and a graded
  mean shift is detected.

## [1.883.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `dunn_test`, the post-hoc
  pairwise follow-up to Kruskal-Wallis. Every snippet value verified.

## [1.882.0] - 2026-09-13

### Added
- `dunn_test.py`: Dunn's post-hoc pairwise rank test after a Kruskal-Wallis rejection,
  using one pooled tie-corrected ranking with Holm/Bonferroni/none p-value adjustment.
  Cross-checked: identical groups give z = 0, well-separated groups are significant
  after adjustment, adjusted p-values never fall below the raw ones, and the pairwise
  z^2 equals the Kruskal-Wallis H exactly on two groups.

## [1.881.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `kruskal_wallis_test` and
  `friedman_test` for k-sample and repeated-measures rank comparisons. Every snippet
  value verified.

## [1.880.0] - 2026-09-13

### Added
- `kruskal_wallis.py`: the `kruskal_wallis_test` (tie-corrected k-sample rank test,
  nonparametric one-way ANOVA) and `friedman_test` (repeated-measures rank test).
  Cross-checked: Kruskal-Wallis reproduces the textbook H = 0.7714 (p = 0.68), equals
  the Mann-Whitney z^2 exactly on two groups, is calibrated under the null and detects
  a three-group shift; Friedman matches its textbook value and flags a consistently
  ordered treatment.

## [1.879.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `ansari_bradley_test` and
  `mood_test` for nonparametric equal-dispersion testing, contrasting them with a
  location rank test that misses a spread difference.

## [1.878.0] - 2026-09-13

### Added
- `scale_tests.py`: the nonparametric two-sample dispersion tests
  `ansari_bradley_test` (outside-in rank scores) and `mood_test` (squared rank
  deviations), each with a normal-approximation p-value. Cross-checked: the null z is
  calibrated (mean ~0, variance ~1 over 500 replications), both reject a 1x-vs-4x
  spread difference at equal location that a rank-sum location test would miss, the
  null rejection rate sits near alpha, and the Ansari statistic matches a hand
  computation.

## [1.877.0] - 2026-09-13

### Documentation
- README: extended the hypothesis-test section with `wilcoxon_signed_rank_test` and
  `sign_test` for paired nonparametric data. Every snippet value verified.

## [1.876.0] - 2026-09-13

### Added
- `wilcoxon.py`: the `wilcoxon_signed_rank_test` (paired / one-sample, tie- and
  zero-handled with a continuity-corrected normal approximation) and the `sign_test`
  (exact two-sided binomial plus a normal approximation). Cross-checked: the W+
  statistic matches hand computation and the classic all-positive textbook case
  (W+ = 45, p = 0.0092), paired mode equals the one-sample test on the differences,
  zeros are dropped, and the sign test's exact p-value matches the binomial CDF.

## [1.875.0] - 2026-09-13

### Documentation
- README: extended the robust scale/location section with `hodges_lehmann_location`
  and `hodges_lehmann_shift`. Every snippet value verified.

## [1.874.0] - 2026-09-13

### Added
- `hodges_lehmann.py`: the Hodges-Lehmann robust location (`hodges_lehmann_location`,
  median of Walsh averages) and two-sample shift (`hodges_lehmann_shift`, median of
  pairwise differences). Cross-checked: exact on symmetric data and a pure +5 shift,
  matches a brute-force reference on 500 random cases, barely moves under a gross
  outlier that drags the mean past 90, and the shift is antisymmetric and ~0 for two
  samples from the same distribution.

## [1.873.0] - 2026-09-13

### Documentation
- README: extended the goodness-of-fit section with `cramer_von_mises_2samp`,
  contrasting its whole-curve squared-gap statistic with KS's single largest gap.
  Every snippet value verified.

## [1.872.0] - 2026-09-13

### Added
- `cramer_von_mises.py`: the two-sample Cramer-von Mises test
  (`cramer_von_mises_2samp`, Anderson's 1962 rank statistic) with an asymptotic
  p-value from the limiting distribution (computed via its Bessel-K_{1/4} series).
  Cross-checked: the asymptotic tail reproduces the published critical values
  (0.461 -> 0.05, 0.743 -> 0.01, 0.347 -> 0.10), the p-value agrees with a
  permutation test, a location shift is rejected, and the null rejection rate sits
  near alpha.

## [1.871.0] - 2026-09-13

### Documentation
- README: extended the goodness-of-fit section with the `energy_distance` /
  `energy_test` distribution-free two-sample test. Every snippet value verified.

## [1.870.0] - 2026-09-13

### Added
- `energy_distance.py`: the Szekely-Rizzo `energy_distance` two-sample statistic
  (`2A - B - C`) and `energy_test`, its permutation test of equal distributions.
  Cross-checked: zero for identical samples, symmetric, matches a brute-force
  reference on 300 random cases, grows with mean separation, and the permutation test
  flags both a mean shift and a variance-only difference (same mean) that a t-test
  would miss.

## [1.869.0] - 2026-09-13

### Documentation
- README: extended the rank-dependence section with `chatterjee_xi` and
  `blomqvist_beta`, showing xi flag a V-shaped function where Kendall's tau is exactly
  zero. Every snippet value verified.

## [1.868.0] - 2026-09-13

### Added
- `chatterjee.py`: Chatterjee's (2020) `chatterjee_xi` rank correlation and
  Blomqvist's `blomqvist_beta` medial correlation. Cross-checked: xi approaches 1 for
  a noiseless monotone relation, is near zero under independence, reaches ~0.99 on a
  multi-period sine where Kendall's tau is nearly blind (nonmonotone functional
  dependence), and matches a reference implementation on 300 random cases; Blomqvist
  gives +/-1 at the comonotone/countermonotone extremes and ~0 under independence.

## [1.867.0] - 2026-09-13

### Documentation
- README: extended the rank-dependence section with distance correlation
  (`distance_correlation`, `distance_covariance`, `distance_variance`), showing it
  flag a symmetric parabola where Pearson is exactly zero. Every snippet value
  verified.

## [1.866.0] - 2026-09-13

### Added
- `distance_correlation.py`: Szekely-Rizzo distance correlation
  (`distance_correlation`, `distance_covariance`, `distance_variance`) via
  double-centered distance matrices — a dependence measure that is zero iff the
  variables are independent, detecting nonlinear structure Pearson misses.
  Cross-checked: 1 for a linear relation, near 0 under independence, ~0.5 on a
  parabola where Pearson is ~0, matches a brute-force reference on 200 random cases,
  and stays in [0, 1].

## [1.865.0] - 2026-09-13

### Documentation
- README: extended the rank-dependence section with `kendall_tau_b` (tie correction),
  `goodman_kruskal_gamma` and `kendall_tau_test`, contrasting tau-a and tau-b on an
  aligned-tie example. Every snippet value verified.

## [1.864.0] - 2026-09-13

### Added
- `kendall_test.py`: tie-corrected `kendall_tau_b`, Goodman-Kruskal
  `goodman_kruskal_gamma`, and `kendall_tau_test` (large-sample normal approximation
  for tau = 0). Cross-checked: tau-b reaches +/-1 for a perfect monotone relation
  even with aligned ties, matches the existing tau-a on 300 tie-free random cases,
  the p-value equals `erfc(|z|/sqrt2)`, and the S-variance uses `n(n-1)(2n+5)/18`.

## [1.863.0] - 2026-09-13

### Documentation
- README: extended the calibration section with Platt scaling (`platt_fit`,
  `platt_predict`, `platt_calibrate`) and guidance on choosing it versus
  `isotonic_fit`. Every snippet value verified.

## [1.862.0] - 2026-09-13

### Added
- `platt_scaling.py`: Platt sigmoid probability calibration (`platt_fit`,
  `platt_predict`, `platt_calibrate`) with Platt's smoothed targets and a Newton
  solver on the regularized logistic loss. Cross-checked: the fit is a genuine loss
  minimum (no local perturbation lowers it), the gradient is zero there, calibrated
  probabilities rise monotonically with the score and stay in [0, 1].

## [1.861.0] - 2026-09-13

### Documentation
- README: new "Forecast calibration (Brier decomposition)" section covering
  `brier_decomposition`, `reliability_curve` and `expected_calibration_error`, with a
  pointer to `isotonic_fit` for recalibration. Every snippet value verified.

## [1.860.0] - 2026-09-13

### Added
- `calibration.py`: Murphy's Brier-score decomposition (`brier_decomposition` into
  reliability / resolution / uncertainty), the `reliability_curve` calibration
  diagram, and `expected_calibration_error`. Cross-checked: the
  reliability-resolution+uncertainty identity reconstructs the raw Brier score
  exactly on 300 random cases, a perfect forecast has zero reliability, a
  constant-at-base-rate forecast has zero resolution, and the reliability curve of a
  well-calibrated forecaster sits on the diagonal (ECE ~0.005 vs ~0.42 miscalibrated).

## [1.859.0] - 2026-09-13

### Documentation
- README: new "Isotonic regression (monotone fit)" section covering
  `isotonic_regression` (PAVA) and `isotonic_fit`, framed around dose-response and
  probability calibration. Every snippet value verified.

## [1.858.0] - 2026-09-13

### Added
- `isotonic.py`: weighted isotonic (monotone) regression by pool-adjacent-violators
  (`isotonic_regression`) and `isotonic_fit`, which sorts by `x`, fits, realigns to
  the original order, and returns an interpolating predictor. Cross-checked against
  the independent max-min weighted-average formula on 300 random weighted cases, plus
  optimality (beats `sorted(y)` in weighted SSE) and decreasing fits.

## [1.857.0] - 2026-09-13

### Documentation
- README: extended the wavelet section with wavelet denoising — `wavelet_denoise`
  (VisuShrink soft/hard shrinkage) plus the exposed `mad_sigma`, `universal_threshold`
  and `soft_threshold` building blocks. Snippet values verified (MSE 0.164 → 0.102 on
  a noisy sinusoid, threshold ~1.53).

## [1.856.0] - 2026-09-13

### Added
- `wavelet_denoise.py`: Donoho-Johnstone wavelet shrinkage on Haar detail
  coefficients (`wavelet_denoise`) with `soft_threshold`/`hard_threshold` rules, the
  `mad_sigma` robust noise estimate (median absolute deviation), and the VisuShrink
  `universal_threshold` (`sigma * sqrt(2 log n)`). Cross-checked: MAD recovers the
  true noise scale of orthonormal Haar detail to within 8%, denoising cuts MSE
  against the clean signal, a zero threshold is the identity, and a huge threshold
  keeps only the coarse approximation.

## [1.855.0] - 2026-09-13

### Documentation
- README: new "Wavelet transform (Haar multiresolution)" section documenting
  `haar_dwt`/`haar_idwt` and `wavelet_energy`, contrasting the scale-localized
  wavelet view with the Fourier transform. Every snippet value verified to run
  (reconstruction error ~1e-16, single-level average/difference, and the
  Parseval energy split summing to one).

## [1.854.0] - 2026-09-13

### Added
- `wavelet.py`: the multilevel Haar discrete wavelet transform (`haar_dwt`,
  `haar_idwt`) and `wavelet_energy` (energy fraction per detail level and the coarse
  approximation). Cross-checked: the transform reconstructs exactly, a constant has
  no detail (all energy in the approximation), it preserves total energy (Parseval),
  the energy fractions sum to one, and a single level gives the scaled
  average/difference of each pair.

## [1.853.0] - 2026-09-13

### Documentation
- README spectral-analysis section now documents `welch_psd` with a worked example.

## [1.852.0] - 2026-09-13

### Added
- `spectral.py` gains `welch_psd`, Welch's averaged-periodogram power-spectral-density
  estimate: overlapping Hann-windowed segments averaged for a low-variance spectrum.
  Cross-checked: it peaks at a sinusoid's frequency (including a noisy one), and its
  variance is far below the raw periodogram's on white noise (the variance reduction
  Welch is for).

## [1.851.0] - 2026-09-13

### Documentation
- README spectral-analysis section now documents `convolve` and
  `fft_autocorrelation` with a worked example.

## [1.850.0] - 2026-09-13

### Added
- `convolution.py`: FFT-based `convolve` (full linear convolution / polynomial
  product in ``O(n log n)``) and `fft_autocorrelation` (Wiener-Khinchin
  autocorrelation, ``acf[0] = 1``). Cross-checked: the convolution matches a direct
  ``O(n m)`` sum, gives the right polynomial product and delta identity, and the FFT
  autocorrelation matches a direct autocovariance sum.

## [1.849.0] - 2026-09-13

### Documentation
- README spectral-analysis section now documents the public `fft` / `ifft` with a
  worked example.

## [1.848.0] - 2026-09-13

### Added
- `fft.py`: public radix-2 Cooley-Tukey `fft` and `ifft` (power-of-two lengths),
  exposing the ``O(n log n)`` transform used internally by the Carr-Madan pricer for
  general spectral work. Cross-checked: the FFT matches a direct DFT, ``ifft(fft(x))``
  round-trips, a unit impulse gives a flat spectrum and a constant a single spike,
  and the transform is linear.

## [1.847.0] - 2026-09-13

### Documentation
- README Hodrick-Prescott section now documents the Savitzky-Golay filter
  (`savgol_filter`, `savgol_coeffs`) with a worked example.

## [1.846.0] - 2026-09-13

### Added
- `savgol.py`: the Savitzky-Golay filter -- `savgol_coeffs` (convolution weights from
  the window's Vandermonde normal equations) and `savgol_filter` (smoothing or
  differentiation, with proper polynomial edge handling). Cross-checked: the
  smoothing weights sum to one and derivative weights to zero, it reproduces
  polynomials up to the fit degree exactly, the first derivative of a quadratic is
  recovered, and it cuts the variance of a noisy series by more than half.

## [1.845.0] - 2026-09-13

### Documentation
- README Kalman-filter section now documents the general `kalman_filter` and
  `kalman_smoother` with a worked constant-velocity example.

## [1.844.0] - 2026-09-13

### Added
- `kalman_filter.py`: the general multivariate linear-Gaussian `kalman_filter`
  (predict/update with the data log-likelihood) and `kalman_smoother` (RTS backward
  pass). Cross-checked: with scalar 1x1 matrices it reproduces the existing
  local-level filter to machine precision, the smoother's covariances never exceed
  the filter's, a 2-D constant-velocity model tracks a moving target and recovers its
  velocity, and the smoothed endpoint equals the filtered endpoint.

## [1.843.0] - 2026-09-13

### Documentation
- README Hidden Markov model section now documents `hmm_baum_welch` and
  `hmm_simulate` with a worked example.

## [1.842.0] - 2026-09-13

### Added
- `hmm.py` gains `hmm_baum_welch` (Baum-Welch EM parameter estimation from an
  observation sequence) and `hmm_simulate` (sample states and observations from an
  HMM). Cross-checked: on a long simulated sticky two-state chain it recovers the
  transition and emission matrices (up to a state permutation), the fitted
  log-likelihood beats the true model's on the data, and the log-likelihood is
  non-decreasing. Uses a strongly-asymmetric initialization to escape the symmetric
  saddle point that traps a uniform start.

## [1.841.0] - 2026-09-13

### Documentation
- README Hidden Markov model section now documents `hmm_posterior` with a worked
  example, contrasting the smoothed marginals with Viterbi's hard path.

## [1.840.0] - 2026-09-13

### Added
- `hmm.py` gains `hmm_posterior`, the scaled forward-backward smoothed posterior
  state probabilities ``P(state_t = i | obs)`` (one distribution per time step, each
  summing to one). Cross-checked: the marginals match a brute-force sum over all
  paths at every step, a deterministic HMM concentrates the posterior on the true
  state, and the rows normalize.

## [1.839.0] - 2026-09-13

### Documentation
- README gains a "Hidden Markov model" section documenting `hmm_forward` and
  `hmm_viterbi` with a worked example. TOC regenerated.

## [1.838.0] - 2026-09-13

### Added
- `hmm.py`: discrete hidden Markov model with `hmm_forward` (sequence log-likelihood
  by the scaled forward algorithm, underflow-safe) and `hmm_viterbi` (most-likely
  hidden-state path by log-space dynamic programming). Cross-checked: the forward
  log-likelihood matches a brute-force sum over all state paths, the Viterbi path and
  its log-probability match a brute-force argmax, a deterministic HMM recovers the
  exact path, and the scaled forward pass stays finite on a 1000-step sequence.

## [1.837.0] - 2026-09-13

### Documentation
- README gains a "Gaussian mixture model" section documenting `fit_gaussian_mixture`
  with a worked regime-clustering example. TOC regenerated.

## [1.836.0] - 2026-09-13

### Added
- `gmm.py`: `fit_gaussian_mixture`, a one-dimensional Gaussian-mixture fit by
  expectation-maximization (responsibility E-step, weighted-moment M-step). Returns
  the mixing weights, means, variances and log-likelihood. Cross-checked: it recovers
  well-separated components and their mixing proportions, the weights sum to one, a
  single component equals the sample mean and variance, and the log-likelihood is
  non-decreasing across iterations.

## [1.835.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `log_determinant`, `mvn_logpdf` and
  `mvn_pdf` with a worked example.

## [1.834.0] - 2026-09-13

### Added
- `mvn.py`: multivariate-normal density and log-determinant via Cholesky --
  `log_determinant` (``2 sum log L_ii``, stable), `mvn_logpdf` (full log-density with
  the Mahalanobis term as a triangular solve, no explicit inverse) and `mvn_pdf`.
  Cross-checked: the log-determinant matches ``ln 24`` for ``diag(2,3,4)`` and the LU
  determinant, the density reduces to the univariate normal, a 2-D density integrates
  to one, and it peaks at the mean.

## [1.833.0] - 2026-09-13

### Documentation
- README Markov-chains section now documents the continuous-time
  `generator_to_transition` and `generator_default_probability` with a worked rating
  generator.

## [1.832.0] - 2026-09-13

### Added
- `markov.py` gains continuous-time chain support: `generator_to_transition`
  (``P(t) = exp(Q t)`` from a rate generator, validated for zero row sums and
  non-negative off-diagonals) and `generator_default_probability` (cumulative
  default at each horizon from a rating generator). Cross-checked: ``P(0)`` is the
  identity, transition rows sum to one and are non-negative, the semigroup
  ``P(1) P(1) = P(2)`` holds, default probabilities are non-decreasing, and an
  absorbing state stays with probability one.

## [1.831.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `matrix_exp` with a worked
  Markov-generator example.

## [1.830.0] - 2026-09-13

### Added
- `matrix_exp.py`: `matrix_exp`, the matrix exponential by scaling-and-squaring with
  a Pade(6,6) approximant -- turns a continuous-time Markov generator ``Q`` into the
  transition matrix ``exp(Q t)``. Cross-checked: ``exp(0) = I``, a diagonal matrix
  exponentiates entrywise, a nilpotent matrix gives ``I + N``, ``exp(A) exp(-A) = I``,
  a rate-generator maps to a valid stochastic matrix (rows sum to one, non-negative),
  and small-norm results match the truncated series. Fixed a Pade denominator sign
  (``V - U``, not ``V - 2U``) caught by the diagonal check.

## [1.829.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `condition_number`, `matrix_rank`,
  `spectral_norm` and `frobenius_norm` with a worked example.

## [1.828.0] - 2026-09-13

### Added
- `svd.py` gains SVD-based matrix diagnostics: `condition_number`
  (``sigma_max/sigma_min``), `matrix_rank` (singular values above a relative
  tolerance), `spectral_norm` (largest singular value) and `frobenius_norm`. Cross-
  checked: the identity has condition 1, ``diag(1000, 1)`` gives 1000, a singular
  matrix is infinite, the rank counts the non-zero singular values, and the norms
  match their singular-value forms.

## [1.827.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `lu_decomposition`, `lu_solve` and
  `determinant` with a worked example.

## [1.826.0] - 2026-09-13

### Added
- `lu.py`: LU decomposition with partial pivoting (`lu_decomposition`, ``P A = L U``),
  the linear solve `lu_solve`, and `determinant` (signed product of the U pivots).
  Cross-checked: known 2x2/3x3 determinants (including a singular zero), ``P A = L U``
  to machine precision, the L/U triangular structure, an ``A x = b`` residual at
  machine precision, and a zero leading pivot handled by row swapping.

## [1.825.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `svd` and `pseudo_inverse` with a
  worked example, noting the numerical-rank and minimum-norm uses.

## [1.824.0] - 2026-09-13

### Added
- `svd.py`: one-sided Jacobi singular value decomposition (`svd`, ``A = U S V'``) and
  the Moore-Penrose `pseudo_inverse`. Cross-checked: ``U S V'`` reconstructs ``A`` and
  ``U`` has orthonormal columns to machine precision, the singular values are
  non-negative, descending, and equal the square roots of the eigenvalues of
  ``A' A``, and the pseudo-inverse gives the same least-squares solution as the QR
  solver.

## [1.823.0] - 2026-09-13

### Documentation
- README matrix-utilities section now documents `qr_decomposition` and `qr_solve`
  with a worked example, noting the stability advantage over the normal equations.

## [1.822.0] - 2026-09-13

### Added
- `qr.py`: Householder QR decomposition (`qr_decomposition`) and QR least squares
  (`qr_solve`). `A = Q R` with orthogonal ``Q`` and upper-triangular ``R``; solving
  ``R x = Q' b`` avoids the ill-conditioned normal-equations matrix and is the stabler
  regression solver. Cross-checked: ``Q R`` reconstructs ``A`` and ``Q' Q = I`` to
  machine precision, ``R`` is upper triangular, and `qr_solve` matches the OLS
  coefficients.

## [1.821.0] - 2026-09-13

### Documentation
- README numerics section now shows `gauss_kronrod` as the general-purpose adaptive
  quadrature, noting its embedded error estimate and sharp-peak resolution.

## [1.820.0] - 2026-09-13

### Added
- `gauss_kronrod.py`: `gauss_kronrod`, adaptive Gauss-Kronrod (G7-K15) quadrature
  with an embedded error estimate. The 15-point Kronrod rule reuses the 7-point
  Gauss nodes, so one evaluation set yields both the estimate and a local error
  bound; the interval with the largest error is bisected until the total falls below
  the tolerance (the scheme behind QUADPACK's QAG). Cross-checked to 1e-11 against
  closed forms, exact on a degree-7 polynomial, resolves a sharp peak that a fixed
  rule would miss, and the single-panel error estimate bounds the true error.

## [1.819.0] - 2026-09-13

### Documentation
- README numerical-utilities section now documents Chebyshev approximation
  (`chebyshev_fit`, `chebyshev_eval`, `chebyshev_derivative`) with a worked example.

## [1.818.0] - 2026-09-13

### Added
- `chebyshev.py`: Chebyshev polynomial approximation on an interval -- `chebyshev_fit`
  (coefficients from the Chebyshev-extrema DCT), `chebyshev_eval` (stable Clenshaw
  recurrence) and `chebyshev_derivative`. Cross-checked: it is exact for polynomials
  up to the truncation degree, reaches machine precision on ``exp`` and ``sin`` at
  modest degree, the coefficients decay geometrically for smooth functions, and the
  differentiated series matches ``d/dx exp = exp`` and ``d/dx sin = cos``. Fixed the
  Chebyshev-Lobatto top-coefficient half-weight that had aliased an exact
  degree-``n`` fit.

## [1.817.0] - 2026-09-13

### Documentation
- README covariance-shrinkage section now documents the RMT denoising routines
  (`marchenko_pastur_edge`, `clip_correlation_eigenvalues`) as an alternative to
  shrinkage.

## [1.816.0] - 2026-09-13

### Added
- `rmt.py`: random-matrix-theory denoising of a correlation matrix --
  `marchenko_pastur_edge` (the ``(1 + sqrt(N/T))^2`` noise-eigenvalue cutoff) and
  `clip_correlation_eigenvalues`, which replaces the sub-edge (noise) eigenvalues
  with their average while keeping the signal eigenvalues, preserving the trace and
  rescaling to unit diagonal. Cross-checked: the MP edge matches its formula, the
  denoised matrix has unit diagonal and is symmetric, the top signal eigenvalue is
  preserved while the noise bulk collapses, and the trace is unchanged.

## [1.815.0] - 2026-09-13

### Documentation
- README gains a "Market stress (turbulence / absorption ratio)" section
  documenting `turbulence_series` and `absorption_ratio` with a worked example. TOC
  regenerated.

## [1.814.0] - 2026-09-13

### Added
- `turbulence.py`: Kritzman-Li systemic-risk gauges -- `turbulence` /
  `turbulence_series` (the Mahalanobis distance of a return vector from its
  historical mean and covariance, spiking on unusual cross-asset moves) and
  `absorption_ratio` (the variance share of the top principal components). Cross-
  checked: in-sample multivariate-normal turbulence averages the number of assets,
  an outlier vector gives a far larger value, the absorption ratio is
  ``n_factors/n`` under equal variance and near one when a factor dominates, and
  stays in ``[0, 1]``.

## [1.813.0] - 2026-09-13

### Documentation
- README portfolio-optimization section now documents `risk_budget_weights` and
  `risk_contributions` next to risk parity.

## [1.812.0] - 2026-09-13

### Added
- `portopt.py` gains `risk_contributions` (each asset's ``w_i (C w)_i`` share of the
  portfolio variance) and `risk_budget_weights`, which generalizes risk parity to an
  arbitrary target risk-budget vector via the fixed point
  ``w_i <- sqrt(b_i w_i / (C w)_i)``. Cross-checked: an equal budget reproduces
  `risk_parity_weights`, the achieved percentage contributions match the target
  budget, the contributions sum to the portfolio variance, and the weights are
  positive and sum to one.

## [1.811.0] - 2026-09-13

### Documentation
- README portfolio-optimization section now documents the concentration measures
  (`herfindahl_index`, `effective_number_of_constituents`, `effective_number_of_bets`)
  with a worked example.

## [1.810.0] - 2026-09-13

### Added
- `concentration.py`: portfolio concentration measures -- `herfindahl_index`
  (``sum w_i^2``), `effective_number_of_constituents` (``1/HHI``) and Meucci's
  `effective_number_of_bets` (exp-entropy of the uncorrelated principal-component
  risk contributions). Cross-checked: equal weights give ENC = n, concentration
  drives ENC toward 1, ENB equals n for equal weights under an identity covariance
  (independent equal-risk factors) and falls toward 1 when one factor dominates, and
  ENB is bounded by n.

## [1.809.0] - 2026-09-13

### Documentation
- README factor-models section now documents `style_analysis` (Sharpe RBSA) next to
  the unconstrained factor regression, with a worked example.

## [1.808.0] - 2026-09-13

### Added
- `style_analysis.py`: Sharpe's returns-based style analysis (`style_analysis`) --
  explains a fund's returns as a long-only, fully-invested mix of index returns
  (weights non-negative, summing to one) by projected-gradient descent on the
  simplex, returning the implied style weights, the ``r_squared`` explained, and the
  selection ``tracking_error``. Cross-checked: it recovers a known blend, the weights
  lie on the simplex, an exact blend gives ``r_squared`` ~ 1, and a single matching
  index gets weight one.

## [1.807.0] - 2026-09-13

### Documentation
- README performance-metrics section now documents the CAPM measures
  (`market_beta`, `treynor_ratio`, `jensens_alpha`, `m_squared`) with a worked
  example.

## [1.806.0] - 2026-09-13

### Added
- `perfmetrics.py` gains CAPM-based performance measures: `market_beta`,
  `treynor_ratio` (excess return per unit of beta), `jensens_alpha` (CAPM
  risk-adjusted excess) and `m_squared` (Modigliani, the portfolio rescaled to the
  market's risk). Cross-checked: the market's beta against itself is 1 and a
  1.5x-levered portfolio recovers beta 1.5, Treynor matches its manual formula,
  Jensen's alpha is ~0 for a pure-beta portfolio and recovers an injected constant
  alpha, and M-squared equals the annualized Sharpe times the market volatility.

## [1.805.0] - 2026-09-13

### Documentation
- README performance-metrics section now documents `gain_to_pain_ratio`,
  `sterling_ratio` and `burke_ratio` with a worked example, noting the L1-vs-L2
  denominator distinction.

## [1.804.0] - 2026-09-13

### Added
- `perfmetrics.py` gains three drawdown-adjusted return ratios: `gain_to_pain_ratio`
  (Schwager, sum of returns over the summed absolute losses), `sterling_ratio`
  (annualized excess over the average drawdown plus a 10% margin), and `burke_ratio`
  (annualized excess over the root-sum-of-squared drawdowns). Cross-checked: each
  matches its component formula, the gain-to-pain ratio is infinite without losses,
  and the L2 Burke ratio penalizes a deep drawdown more than a shallow one.

## [1.803.0] - 2026-09-13

### Documentation
- README performance-metrics section now documents the drawdown risk measures
  (`average_drawdown`, `drawdown_at_risk`, `conditional_drawdown_at_risk`) with a
  worked example.

## [1.802.0] - 2026-09-13

### Added
- `cdar.py`: drawdown-based risk measures -- `average_drawdown`, `drawdown_at_risk`
  (DaR, the drawdown quantile) and `conditional_drawdown_at_risk` (CDaR, the mean of
  the worst-tail drawdowns, a coherent drawdown analogue of expected shortfall).
  Cross-checked: the ordering ``CDaR >= DaR >= average >= 0`` holds, all are bounded
  by the maximum drawdown, both are monotone in the confidence level, a monotonic
  decline gives CDaR equal to the maximum drawdown, and a monotonically rising path
  has zero drawdown.

## [1.801.0] - 2026-09-12

### Documentation
- README cointegration section now documents `kpss_test` as the confirmatory
  complement to ADF, with a worked example and the joint-interpretation note.

## [1.800.0] - 2026-09-12

### Added
- `kpss.py`: `kpss_test`, the Kwiatkowski-Phillips-Schmidt-Shin stationarity test
  (level ``"c"`` or trend ``"ct"``), the complement of ADF -- its null is
  stationarity, so a small p-value rejects it. The statistic uses the Newey-West
  long-run variance and is compared to the asymptotic critical values.
  Cross-checked: white noise is not rejected (and ADF rejects its unit root, the
  opposite conclusion), a random walk is rejected, and a linear-trend-plus-noise
  series is stationary under ``"ct"`` but rejected under ``"c"``.

## [1.799.0] - 2026-09-12

### Documentation
- README GARCH section now documents `arch_lm_test` with a worked example as the
  pre-fit check for volatility clustering.

## [1.798.0] - 2026-09-12

### Added
- `arch_test.py`: `arch_lm_test`, Engle's ARCH-LM test for conditional
  heteroskedasticity (volatility clustering). Regresses the squared series on its
  own lags and returns ``(n R^2, p_value)`` referenced to a chi-square with ``lags``
  degrees of freedom. Cross-checked: iid data holds the 5% size and is not rejected,
  while a simulated GARCH(1,1) series is strongly rejected.

## [1.797.0] - 2026-09-12

### Documentation
- README goodness-of-fit section now documents `anderson_darling_normal` with a
  worked example, noting its tail-sensitivity vs KS.

## [1.796.0] - 2026-09-12

### Added
- `anderson_darling.py`: `anderson_darling_normal`, the Anderson-Darling test of
  normality (mean/sd estimated), with the Stephens small-sample adjustment and the
  D'Agostino-Stephens p-value approximation. More tail-sensitive than
  Kolmogorov-Smirnov. Cross-checked: normal samples hold the 5% size and are not
  rejected, while exponential (skewed) and heavy-tailed mixtures are strongly
  rejected; the statistic is non-negative.

## [1.795.0] - 2026-09-12

### Documentation
- README serial-correlation section now documents the Wald-Wolfowitz `runs_test`
  and `runs_test_binary` with a worked example.

## [1.794.0] - 2026-09-12

### Added
- `runs_test.py`: the Wald-Wolfowitz runs test for sequence randomness --
  `runs_test` (dichotomizes a numeric series about its median) and
  `runs_test_binary` (two-symbol sequence), each returning a two-sided ``(z,
  p_value)``. Cross-checked: a random series holds its 5% size, sorted data gives
  too few runs (``z << 0``, clustering), an alternating series gives too many
  (``z >> 0``), and a known 3+3 binary case matches the closed-form z.

## [1.793.0] - 2026-09-12

### Documentation
- README gains a "Benford's law (first-digit anomaly detection)" section
  documenting `benford_chi_square`, `benford_mad`, `first_digit_distribution` and
  `benford_expected` with a worked example. TOC regenerated.

## [1.792.0] - 2026-09-12

### Added
- `benford.py`: Benford's-law first-digit analysis for anomaly detection --
  `benford_expected`, `first_digit`, `first_digit_distribution`,
  `benford_chi_square` (goodness-of-fit, 8 df) and `benford_mad` (Nigrini's mean
  absolute deviation). Cross-checked: the expected probabilities sum to one with
  ``P(1) = 0.301``, Fibonacci numbers and powers of two conform (not rejected, MAD
  below Nigrini's 0.006 threshold), and a uniform sample is strongly rejected.

## [1.791.0] - 2026-09-12

### Documentation
- README survival-analysis section now documents `median_survival_time` and
  `restricted_mean_survival_time` with a worked example.

## [1.790.0] - 2026-09-12

### Added
- `survival.py` gains `median_survival_time` (earliest time Kaplan-Meier drops to
  0.5, or ``None`` if it never does) and `restricted_mean_survival_time` (RMST, the
  area under the KM curve up to a horizon). Cross-checked: the no-censoring median
  matches the sorted median, RMST matches a by-hand step integral and the
  exponential closed form ``(1 - e^{-lambda tau})/lambda``, it grows with the horizon
  and is capped by it, and the median is ``None`` when every observation is censored.

## [1.789.0] - 2026-09-12

### Documentation
- README survival-analysis section now documents `log_rank_test` with a worked
  example.

## [1.788.0] - 2026-09-12

### Added
- `survival.py` gains `log_rank_test`, the Mantel-Cox log-rank test comparing two
  survival curves via the accumulated observed-minus-expected events and their
  hypergeometric variance, referenced to a chi-square(1). Cross-checked: two samples
  from the same hazard are not rejected, a threefold-faster hazard is strongly
  rejected, identical data gives a zero statistic, and the statistic is symmetric in
  the two groups.

## [1.786.0] - 2026-09-12

### Added
- `survival.py`: nonparametric survival estimators for right-censored data --
  `kaplan_meier` (product-limit survival), `nelson_aalen` (cumulative hazard), and
  `survival_at` to evaluate the step curve. Cross-checked: with no censoring the KM
  estimate equals ``1 - ECDF``, it is non-increasing in ``[0, 1]``, a censored
  worked example matches by hand, the Nelson-Aalen hazard is the running
  ``sum d_i/n_i`` and ``exp(-H)`` tracks KM while the risk set is large.

## [1.785.0] - 2026-09-12

### Documentation
- README cross-validation section now documents the generic information criteria
  (`gaussian_log_likelihood`, `aic`, `aicc`, `bic`, `hqic`) with a worked example.

## [1.784.0] - 2026-09-12

### Added
- `info_criteria.py`: generic model-selection criteria from a log-likelihood --
  `aic`, `aicc` (small-sample corrected), `bic`, `hqic`, plus
  `gaussian_log_likelihood` to turn a residual sum of squares into the log-likelihood
  they consume. Cross-checked: the Gaussian log-likelihood matches its direct MLE
  form, BIC exceeds AIC for ``n >= 8``, AICc approaches AIC as ``n`` grows, the
  penalties order ``AIC < HQIC < BIC``, more parameters raise every criterion, and a
  smaller RSS lowers the AIC.

## [1.783.0] - 2026-09-12

### Documentation
- README regression section now documents `elastic_net` next to LASSO, noting the
  L1/L2 blend and the grouping effect on correlated features.

## [1.782.0] - 2026-09-12

### Added
- `lasso.py` gains `elastic_net`, coordinate-descent regression with a mixed L1+L2
  penalty (`l1_ratio` blends the two). Cross-checked: ``l1_ratio = 1`` reproduces
  `lasso_regression` exactly, ``alpha = 0`` recovers OLS, it still zeros noise
  features at intermediate mixes, and on near-duplicate correlated predictors it
  shares weight between them (the grouping effect) where pure LASSO drops one.

## [1.781.0] - 2026-09-12

### Documentation
- README regression section now documents `lasso_regression` next to ridge, noting
  the exact-zero feature selection.

## [1.780.0] - 2026-09-12

### Added
- `lasso.py`: `lasso_regression`, L1-penalized least squares by cyclic coordinate
  descent with soft-thresholding on standardized features (unpenalized intercept,
  coefficients returned on the original scale). Unlike ridge it drives coefficients
  exactly to zero, so it selects features. Cross-checked: ``alpha = 0`` matches OLS,
  a moderate penalty zeros the noise features while keeping the real ones, a large
  penalty zeros all slopes leaving the intercept at ``mean(y)``, and the real-feature
  slope shrinks monotonically as the penalty grows.

## [1.779.0] - 2026-09-12

### Documentation
- README regression section now documents `quantile_regression` next to OLS and
  ridge, noting the median/tail fits and the residual-quantile property.

## [1.778.0] - 2026-09-12

### Added
- `quantile_regression.py`: linear `quantile_regression` fitting a conditional
  ``tau``-quantile by minimizing the pinball loss via iteratively-reweighted least
  squares. Cross-checked: the median fit recovers the true slope and intercept, the
  intercept rises with ``tau`` while the slope stays put under homoskedastic noise,
  the fitted line has the residual-quantile property (a fraction ``tau`` of points
  fall below it), and it beats OLS on the pinball loss at its target quantile.

## [1.777.0] - 2026-09-12

### Documentation
- README forecast-accuracy section now documents the probabilistic scoring rules
  (`pinball_loss`, `interval_score`, `interval_coverage`, `crps_ensemble`) with a
  worked example.

## [1.776.0] - 2026-09-12

### Added
- `prob_forecast.py`: proper scoring rules for probabilistic forecasts --
  `pinball_loss` (quantile loss), `interval_score` (Winkler score for prediction
  intervals), `interval_coverage`, and `crps_ensemble` (continuous ranked
  probability score of a sample forecast). Cross-checked: the pinball loss is
  minimized at the true quantile and is asymmetric in ``tau``; CRPS reduces to the
  absolute error for a single-member ensemble and is zero for a perfect forecast;
  coverage counts actuals inside their intervals; and the interval score rewards
  tight covering intervals while penalizing misses.

## [1.775.0] - 2026-09-12

### Documentation
- README forecast-accuracy section now documents `theil_u1` and `theil_u2` alongside
  the other error metrics.

## [1.774.0] - 2026-09-12

### Added
- `forecast_metrics.py` gains Theil's U statistics: `theil_u2` (forecast RMSE over
  the no-change naive RMSE) and `theil_u1` (the inequality coefficient bounded in
  ``[0, 1]``). Cross-checked: a perfect forecast gives zero, the persistence forecast
  gives ``U2 = 1`` exactly, a forecast better than naive gives ``U2 < 1`` and a worse
  one ``U2 > 1``, and ``U1`` stays in ``[0, 1]``.

## [1.773.0] - 2026-09-12

### Documentation
- README Newey-West HAC section now documents the forecast-combination routines
  (`simple_average_forecast`, `inverse_mse_weights`, `optimal_combination_weights`,
  `combine_forecasts`) with a worked example.

## [1.772.0] - 2026-09-12

### Added
- `forecast_combine.py`: forecast-combination weights -- `simple_average_forecast`
  (equal weights), `inverse_mse_weights` (proportional to ``1/MSE``),
  `optimal_combination_weights` (Bates-Granger minimum-variance from the error
  covariance) and a `combine_forecasts` helper. Cross-checked: all weight schemes
  sum to one, inverse-MSE favors the more accurate model and reduces to equal weights
  at equal variance, and the minimum-variance combination's error MSE never exceeds
  the best individual model's (dramatically better when the errors hedge).

## [1.771.0] - 2026-09-12

### Documentation
- README Newey-West HAC section now documents `diebold_mariano` with a worked
  example, noting it reuses the HAC variance and is antisymmetric in its arguments.

## [1.770.0] - 2026-09-12

### Added
- `forecast_test.py`: the Diebold-Mariano test of equal predictive accuracy
  (`diebold_mariano`) -- standardizes the mean loss differential of two forecasts by
  its Newey-West HAC standard error, with the Harvey-Leybourne-Newbold small-sample
  correction and a Student-t reference. Cross-checked: a lower-variance forecast
  gives a significantly negative statistic, the test is antisymmetric in its
  arguments (swapping flips the sign, same p-value), and two equally-accurate
  forecasts are not rejected.

## [1.769.0] - 2026-09-12

### Documentation
- README sample-risk-measures section now documents `expectile`, noting it is the
  only coherent-and-elicitable risk measure and is backtestable by a single scoring
  function.

## [1.768.0] - 2026-09-12

### Added
- `riskmeasures.py`: `expectile`, the ``tau``-expectile of a P&L sample -- the only
  risk measure that is both coherent (for ``tau >= 0.5``) and elicitable. Solved by
  bisection on the asymmetric-least-squares first-order condition. Cross-checked:
  the 0.5-expectile equals the mean loss, it is monotone in ``tau`` and above the
  mean for ``tau > 0.5``, and the first-order condition holds at the solution.

## [1.767.0] - 2026-09-12

### Documentation
- README portfolio-risk section now documents the VaR/ES backtests (`kupiec_pof`,
  `christoffersen_cc`, `christoffersen_independence`, `acerbi_szekely_es`) with a
  worked example.

## [1.766.0] - 2026-09-12

### Added
- `var_backtest.py`: Value-at-Risk / Expected-Shortfall backtests -- `kupiec_pof`
  (unconditional coverage), `christoffersen_independence` and `christoffersen_cc`
  (exception clustering / joint conditional coverage), and `acerbi_szekely_es` (the
  Acerbi-Szekely ES calibration statistic). Cross-checked on simulated normal losses:
  a correctly-specified 99% VaR/ES passes all three tests with the AS statistic near
  zero; too-low a VaR is rejected by Kupiec; and an understated (overstated) ES gives
  a positive (negative) AS statistic.

## [1.765.0] - 2026-09-12

### Documentation
- README entropy section now documents `mutual_information` and `transfer_entropy`
  with a worked example, noting the symmetry vs directionality and the lead-lag use.

## [1.764.0] - 2026-09-12

### Added
- `transfer_entropy.py`: histogram-based `mutual_information` (symmetric, zero iff
  independent) and Schreiber's directional `transfer_entropy` (lag-1). Cross-checked:
  mutual information is near zero for independent series, symmetric in its arguments,
  and large under dependence; transfer entropy is much larger in the driving
  direction of a coupled system (``Y_{t+1} = 0.6 X_t + noise`` gives
  ``TE_{X->Y} >> TE_{Y->X}``) and near zero for independent series.

## [1.763.0] - 2026-09-12

### Documentation
- README gains an "Entropy (time-series regularity)" section documenting
  `approximate_entropy`, `sample_entropy` and `permutation_entropy` with a worked
  example. TOC regenerated.

## [1.762.0] - 2026-09-12

### Added
- `entropy_ts.py`: time-series regularity measures -- `approximate_entropy`
  (Pincus), `sample_entropy` (Richman-Moorman) and `permutation_entropy`
  (Bandt-Pompe). Cross-checked: a regular sine has far lower approximate/sample
  entropy than white noise, permutation entropy is exactly zero for a monotone
  series and near one for noise, stays in ``[0, 1]``, and is invariant to a monotone
  transform of the series.

## [1.761.0] - 2026-09-12

### Documentation
- README long-memory section now documents `dfa_exponent` and `dfa_fluctuations`
  alongside the Hurst and GPH estimators.

## [1.760.0] - 2026-09-12

### Added
- `dfa.py`: detrended fluctuation analysis (`dfa_exponent`, `dfa_fluctuations`) --
  the scaling exponent from the RMS of the linearly-detrended integrated profile
  across window sizes, robust to slow trends where rescaled-range analysis is not.
  Cross-checked: white noise gives ``alpha ~ 0.5``, a random walk ``~ 1.5``,
  differenced noise ``< 0.5`` (anti-persistent), and integrating a series adds one
  to its exponent.

## [1.759.0] - 2026-09-12

### Documentation
- README long-memory section now documents `gph_estimate` and
  `fractional_integrate` alongside fractional differencing, closing the
  estimate-then-difference workflow.

## [1.758.0] - 2026-09-12

### Added
- `gph.py`: the Geweke-Porter-Hudak estimator of the long-memory parameter ``d``
  (`gph_estimate`, a log-periodogram regression with the asymptotic ``pi^2/6``
  standard error) and `fractional_integrate` (apply ``(1 - L)^{-d}``, the inverse of
  the fractional difference, for generating ARFIMA(0,d,0) series). Cross-checked:
  white noise gives ``d`` insignificantly different from zero, the estimator recovers
  a known memory parameter from a fractionally-integrated series, the standard error
  shrinks with the bandwidth, and integrate-then-difference round-trips to white
  noise.

## [1.757.0] - 2026-09-12

### Documentation
- README long-memory section now documents fractional differencing
  (`fractional_difference`, `fixed_width_fracdiff`, `fracdiff_weights`) next to the
  Hurst exponent.

## [1.756.0] - 2026-09-12

### Added
- `fracdiff.py`: fractional differencing for long-memory (ARFIMA) series --
  `fracdiff_weights` (binomial weights of ``(1 - L)^d``), `fractional_difference`
  (full expansion) and `fixed_width_fracdiff` (Lopez de Prado's fixed-window
  variant). Cross-checked: integer orders reproduce the standard differences
  (``d = 0`` identity, ``d = 1`` first difference, ``d = 2`` second difference), the
  weights follow the ``w_k = w_{k-1} (k-1-d)/k`` recursion and decay slowly for
  fractional ``d``, and the fixed-width filter trims to a constant memory window.

## [1.755.0] - 2026-09-12

### Documentation
- README Hawkes section now documents `hawkes_residuals` and `hawkes_gof_test`
  (time-rescaling diagnostics) with a worked example.

## [1.754.0] - 2026-09-12

### Added
- `hawkes.py` gains time-rescaling diagnostics: `hawkes_residuals` (the integrated
  intensity between events, i.i.d. unit exponentials under a correct model) and
  `hawkes_gof_test` (a one-sample Kolmogorov-Smirnov test of those residuals against
  the Exp(1) CDF). Cross-checked: the residual recursion matches a brute-force
  compensator to machine precision, the residuals of a well-fit path have mean one,
  and the GOF test accepts the true model (p ~ 0.6) while rejecting wrong parameters
  (p ~ 0).

## [1.753.0] - 2026-09-12

### Documentation
- README gains a "Hawkes self-exciting process" section documenting
  `hawkes_simulate`, `hawkes_fit`, `hawkes_intensity`, `hawkes_branching_ratio` and
  `hawkes_log_likelihood` with a worked example. TOC regenerated.

## [1.752.0] - 2026-09-12

### Added
- `hawkes.py`: exponential-kernel self-exciting (Hawkes) point process --
  `hawkes_intensity`, `hawkes_branching_ratio`, the exact recursive
  `hawkes_log_likelihood`, `hawkes_simulate` (Ogata thinning), and the
  maximum-likelihood `hawkes_fit`. Cross-checked: the recursive log-likelihood
  matches a brute-force sum term for term, the simulated event rate matches the
  stationary theory `mu / (1 - alpha/beta)` (and the Poisson limit at `alpha = 0`),
  and the MLE recovers the parameters of a long simulated path.

### Fixed
- The Hawkes Ogata simulation used a strict `t_i < t` upper-intensity bound, which
  dropped the just-fired event's own excitation and undersampled the process
  (event rate ~30% low at high branching). The bound now includes events at `t`, so
  the simulated rate matches the stationary mean.

## [1.751.0] - 2026-09-12

### Documentation
- README optimal-execution section now documents the spread decomposition
  (`quoted_spread`, `effective_spread`, `realized_spread`, `price_impact`) with a
  worked example and the `effective = realized + price_impact` identity.

## [1.750.0] - 2026-09-12

### Added
- `microstructure.py` gains the transaction-cost spread decomposition:
  `quoted_spread`, `effective_spread`, `realized_spread` and `price_impact`. Cross-
  checked: the decomposition identity `effective = realized + price_impact` holds
  term by term; with no midpoint move the realized spread equals the effective and
  the impact is zero; and a trade at the mid with a permanent move is pure price
  impact (zero effective, positive impact).

## [1.749.0] - 2026-09-12

### Documentation
- README optimal-execution section now documents the trade-sign classifiers
  (`tick_rule`, `quote_rule`, `lee_ready`) with a worked example, noting they
  produce the signed flow the order-flow measures consume.

## [1.748.0] - 2026-09-12

### Added
- `trade_sign.py`: trade-sign classifiers that infer the aggressor side from tape
  data -- `tick_rule` (sign by the last price change), `quote_rule` (side of the
  bid-ask midpoint), and `lee_ready` (the Lee-Ready hybrid: quote rule with a
  tick-rule tiebreak at the midpoint). These produce the signed order flow the
  microstructure measures consume. Cross-checked: the tick rule signs upticks and
  downticks and carries the prior sign on flat ticks, the quote rule is zero exactly
  at the midpoint, Lee-Ready resolves those to +/-1 via the tick rule, and the
  output feeds VPIN and the order-flow imbalance cleanly.

## [1.747.0] - 2026-09-12

### Documentation
- README optimal-execution section now documents the order-flow measures
  (`kyle_lambda_regression`, `order_flow_imbalance`, `vpin`) with a worked example.

## [1.746.0] - 2026-09-12

### Added
- `microstructure.py`: order-flow measures from trade data --
  `kyle_lambda_regression` (empirical Kyle's lambda, the OLS slope of price change
  on signed order flow), `order_flow_imbalance` (net signed volume over total), and
  `vpin` (volume-synchronized probability of informed trading). Cross-checked: the
  regression recovers a known price-impact coefficient from simulated data and is
  zero when the flow carries no impact; the imbalance and VPIN are correctly signed
  / bounded, with VPIN zero for balanced flow and one for one-sided flow.

## [1.745.0] - 2026-09-12

### Documentation
- README optimal-execution section now documents the empirical liquidity proxies
  (`roll_spread`, `amihud_illiquidity`, `corwin_schultz_spread`) alongside the
  theoretical impact models.

## [1.744.0] - 2026-09-12

### Added
- `liquidity.py`: low-frequency liquidity and transaction-cost proxies --
  `roll_spread` (Roll's effective spread from the serial covariance of price
  changes), `amihud_illiquidity` (absolute return per dollar of volume), and
  `corwin_schultz_spread` (high-low spread estimator). Cross-checked: Roll recovers
  a known bid-ask-bounce spread from a simulated tape (0.10 -> ~0.10), Amihud scales
  inversely with volume (10x less volume -> 10x the measure), and Corwin-Schultz is
  non-negative and rises monotonically with the true spread.

## [1.743.0] - 2026-09-12

### Documentation
- README realized-volatility section now documents `bns_jump_test` and
  `tripower_quarticity` with a worked example, noting the null size and jump power.

## [1.742.0] - 2026-09-12

### Added
- `jump_test.py`: the Barndorff-Nielsen-Shephard / Huang-Tauchen realized-volatility
  jump test (`bns_jump_test`) and its `tripower_quarticity` scale. The ratio
  statistic ``sqrt(n)(RV-BV)/RV`` standardized by the tripower quarticity is
  asymptotically standard normal under no jump; a large positive z rejects. Cross-
  checked on a simulated diffusion: under the null the statistic is centred at zero
  with ~5% rejection at the 5% level, the test has power ~1 against an added jump,
  and the tripower quarticity is consistent for ``sigma^4`` and jump-robust.

## [1.741.0] - 2026-09-12

### Documentation
- README realized-volatility section now documents `min_realized_variance`,
  `med_realized_variance` and `realized_quarticity` next to bipower variation.

## [1.740.0] - 2026-09-12

### Added
- `realized.py` gains `min_realized_variance` and `med_realized_variance` (the
  Andersen-Dobrev-Schaumburg MinRV/MedRV nearest-neighbour jump-robust
  integrated-variance estimators) and `realized_quarticity` (the `(n/3) sum r^4`
  estimator of integrated quarticity). Cross-checked on a simulated diffusion: with
  no jump MinRV, MedRV, bipower and RV all match the integrated variance; with a
  large jump the RV inflates ~26x while MinRV and MedRV stay near the truth; the
  quarticity is positive and consistent for `sigma^4`, and MedRV is unchanged by a
  single spiked return (median of three).

## [1.739.0] - 2026-09-12

### Documentation
- README realized-volatility section now documents `realized_kernel` next to the
  two-scale estimator, with a worked example and the noise-robustness note.

## [1.738.0] - 2026-09-12

### Added
- `realized_kernel.py`: `realized_kernel`, the Barndorff-Nielsen-Hansen-Lunde-
  Shephard realized-kernel estimator of integrated variance -- weighted intraday
  return autocovariances (flat-top Parzen kernel) that cancel the microstructure-
  noise bias in the naive realized variance. Cross-checked against a simulated
  Brownian-motion price process: it matches the realized variance and the true
  integrated variance without noise, and under i.i.d. noise its mean bias is far
  below the naive estimator's (which inflates ~10x), while staying non-negative.

## [1.737.0] - 2026-09-12

### Documentation
- README numerics section now shows `clenshaw_curtis` next to the other quadrature
  rules, noting its free order and Runge-function robustness.

## [1.736.0] - 2026-09-12

### Added
- `quadrature.py`: `clenshaw_curtis` quadrature -- samples the integrand at the
  Chebyshev extrema and combines them with the classic cosine-series weights.
  Spectrally accurate for smooth integrands like Gauss-Legendre, but with a free
  order ``n`` and nesting nodes, so it scales past the fixed 2-5 point Gauss rule.
  Cross-checked to 1e-12 against closed forms (``sin`` over a half period, ``exp``,
  a degree-7 polynomial exact at ``n = 8``, ``arctan'`` giving ``pi/4``, a Gaussian
  against ``erf``), against Romberg on a damped cosine, and on the Runge function.

## [1.735.0] - 2026-09-12

### Documentation
- README SVI section now documents the jump-wing parameterization
  (`raw_to_jumpwing`, `jumpwing_to_raw`) with a worked example, noting the map is
  closed-form and round-trips to machine precision.

## [1.734.0] - 2026-09-12

### Added
- `svi_jumpwing.py`: the Gatheral-Jacquier jump-wing parameterization of an SVI
  slice (`SVIJumpWing`, `raw_to_jumpwing`, `jumpwing_to_raw`) -- trader-friendly ATM
  variance, ATM skew, and left/right wing slopes, with the closed-form map to and
  from the raw `(a, b, rho, m, s)` parameters at a fixed expiry. Cross-checked: the
  jump-wing ATM variance equals the raw ATM total variance, the skew matches a
  finite-difference of the total-variance curve, the wing slopes match the raw
  asymptotics, and raw -> jump-wing -> raw reproduces the smile to machine precision
  across expiries and both skew signs.

## [1.733.0] - 2026-09-12

### Documentation
- README hypothesis-test section now documents the multiple-testing corrections
  (`bonferroni`, `holm`, `benjamini_hochberg`, `benjamini_yekutieli`) with a worked
  example and guidance on family-wise vs false-discovery control.

## [1.732.0] - 2026-09-12

### Added
- `multiple_testing.py`: p-value corrections for multiple hypotheses --
  `bonferroni` and `holm` (family-wise error rate), `benjamini_hochberg` (FDR under
  independence/positive dependence) and `benjamini_yekutieli` (FDR under arbitrary
  dependence). Each returns monotone adjusted p-values aligned with the input.
  Cross-checked: Benjamini-Hochberg matches the known worked example, every method
  gives an adjusted p at least the raw one, the conservativeness ordering holds
  (Bonferroni >= Holm, Benjamini-Yekutieli >= Benjamini-Hochberg), and the input
  order is preserved for unsorted p-values.

## [1.731.0] - 2026-09-12

### Documentation
- README hypothesis-test section now documents the power and sample-size routines
  (`two_sample_t_power`, `two_sample_t_sample_size`, `proportion_sample_size`, and
  the one-sample analogues) with a worked example, noting the Monte-Carlo agreement
  and textbook sample sizes.

## [1.730.0] - 2026-09-12

### Added
- `power.py`: statistical power and sample-size calculations (normal approximation)
  for the two-sample and one-sample mean tests (Cohen's d) and the two-proportion
  test -- `two_sample_t_power`/`two_sample_t_sample_size`,
  `one_sample_z_power`/`one_sample_z_sample_size`,
  `proportion_power`/`proportion_sample_size`. Cross-checked: the computed power
  equals an empirical Monte-Carlo rejection rate, the sample sizes reproduce
  textbook values (~64 per group for d=0.5 at 80% power; ~170 per group for a
  0.50 vs 0.65 proportion), each solved sample size delivers at least the target
  power, and a zero effect gives power equal to the significance level.

## [1.729.0] - 2026-09-12

### Documentation
- README hypothesis-test section now documents the four proportion confidence
  intervals (`wald_interval`, `wilson_interval`, `agresti_coull_interval`,
  `clopper_pearson_interval`) with a worked example and guidance on which to use.

## [1.728.0] - 2026-09-12

### Added
- `proportion_ci.py`: four confidence intervals for a binomial proportion --
  `wald_interval` (normal approximation), `wilson_interval` (score),
  `agresti_coull_interval`, and the exact `clopper_pearson_interval` (by inverting
  the binomial CDF). Cross-checked: Clopper-Pearson matches the textbook `(2, 10)`
  interval `[0.0252, 0.5561]`, every interval contains the point estimate and stays
  in `[0, 1]`, Wald converges to Wilson for large `n`, the exact interval is the
  widest, and its Monte-Carlo coverage is at least the nominal level.

## [1.727.0] - 2026-09-12

### Documentation
- README hypothesis-test example now includes `one_sample_t_test`, `paired_t_test`
  and `mann_whitney_u`, noting the paired-t / one-sample-on-differences equivalence
  and that Mann-Whitney is the distribution-free alternative.

## [1.726.0] - 2026-09-12

### Added
- `hypothesis.py` gains three more tests: `one_sample_t_test` (mean vs a reference),
  `paired_t_test` (dependent samples), and `mann_whitney_u` (distribution-free
  rank-sum with a tie- and continuity-corrected normal approximation).
  Cross-checked: the one-sample t matches a by-hand calculation and is zero at the
  sample mean, the paired t equals the one-sample t on the within-pair differences,
  and the Mann-Whitney U gives complete separation (U = 0) for disjoint groups,
  matches a known worked example, and is non-significant for identical samples.

## [1.725.0] - 2026-09-12

### Documentation
- README rank-dependence section now documents `pearson_r` and
  `pearson_correlation_test` (the linear correlation with a t-test and Fisher-z
  interval) alongside the rank-based Kendall/Spearman measures.

## [1.724.0] - 2026-09-12

### Added
- `correlation_test.py`: `pearson_r` and `pearson_correlation_test` -- the sample
  correlation with a two-sided t-test of ``rho = 0`` and a Fisher-z confidence
  interval. Cross-checked: ``r`` matches a by-hand calculation, the t statistic
  equals ``r sqrt((n-2)/(1-r^2))``, the Fisher interval covers the true correlation
  at its nominal rate in a Monte Carlo (94.8% at the 95% level), and a strong
  signal gives a vanishing p-value while noise does not.

## [1.723.0] - 2026-09-12

### Documentation
- README OLS section now documents the new inferential outputs (`p_values`,
  `conf_int`, `f_pvalue`, `confidence` argument) and notes they come from the
  library's t and F distributions, with the single-regressor `F = t^2` identity.

## [1.722.0] - 2026-09-12

### Changed
- `ols_fit` now also reports inferential statistics built on the distribution CDFs:
  two-sided coefficient `p_values` and `conf_int` (from the t distribution at a new
  `confidence` argument, default 0.95) and the overall `f_pvalue`. Existing keys are
  unchanged. Cross-checked: for a single regressor the overall F equals the slope
  t-squared and its p-value matches the slope's, the intervals are symmetric about
  the coefficient and widen with the confidence level, a strong linear signal gives
  a vanishing slope p-value, and noise gives a non-significant one.

## [1.720.0] - 2026-09-12

### Added
- `hypothesis.py`: classical hypothesis tests built on the distribution CDFs, each
  returning `(statistic, p_value)` -- `chi_square_gof_test`,
  `chi_square_independence_test`, `one_way_anova` (F test), `two_sample_t_test`
  (pooled and Welch), and the exact `binomial_test`. Cross-checked: the chi-square
  goodness-of-fit statistic matches a by-hand calculation, the two-group ANOVA F
  equals the pooled t-squared with an identical p-value, the pooled and Welch t
  agree for equal sizes and variance, and the two-sided binomial p-value matches
  the symmetric-tail sum.

## [1.719.0] - 2026-09-12

### Documentation
- README gains a "Probability distributions" section documenting the gamma,
  chi-square, Poisson, F and binomial routines with a worked example (textbook
  chi-square/F critical values, tail p-value, Poisson/binomial CDFs), and a note
  that the discrete CDFs use the gamma/beta identities for large-parameter
  stability. TOC regenerated.

## [1.718.0] - 2026-09-12

### Added
- `distributions.py`: public gamma, chi-square, Poisson, F and binomial
  distributions built on the incomplete gamma and beta functions -- CDFs, densities
  or masses, and quantiles (`gamma_cdf`/`gamma_pdf`/`gamma_ppf`,
  `chi2_cdf`/`chi2_sf`/`chi2_ppf`, `poisson_pmf`/`poisson_cdf`, `f_cdf`/`f_ppf`,
  `binomial_cdf`/`binomial_pmf`). The Poisson and binomial CDFs use the gamma/beta
  identities (`Q(k+1, lam)`, `I_{1-p}(n-k, k+1)`) for stability at large parameters.
  Cross-checked: the gamma reduces to the exponential CDF, the chi-square and F
  quantiles match textbook tables (`chi2_ppf(0.95, 10) = 18.307`,
  `f_ppf(0.95, 1, 10) = 4.965`), the Poisson and binomial CDFs match direct mass
  sums, masses sum to one, and every quantile round-trips its CDF.

## [1.717.0] - 2026-09-12

### Documentation
- README numerics section now documents the public special functions
  (`gammainc`, `gammaincc`, `betainc`, `digamma`, `erfinv`) with a worked example
  and the distribution relationships they provide (gamma/chi-square/Poisson via the
  incomplete gamma, Student-t/F/binomial via the incomplete beta, normal quantiles
  via `erfinv`).

## [1.716.0] - 2026-09-12

### Added
- `special.py`: public special functions that underpin the distribution routines --
  `gammainc`/`gammaincc` (regularized lower/upper incomplete gamma), `betainc`
  (regularized incomplete beta), `digamma`, and `erfinv`. Cross-checked against
  identities and known values: `P + Q = 1`, `P(1/2, x) = erf(sqrt(x))`,
  `P(1, x) = 1 - e^{-x}`, the incomplete-beta symmetry `I_x(a,b) = 1 - I_{1-x}(b,a)`
  and `I_x(1,1) = x`, `digamma(1) = -gamma` with the recurrence
  `psi(x+1) = psi(x) + 1/x`, and `erf(erfinv(y)) = y` to machine precision.

## [1.715.0] - 2026-09-12

### Documentation
- README numerics section now documents Gauss-Laguerre quadrature
  (`gauss_laguerre_integral`, `gauss_laguerre_nodes_weights`) with a worked
  example, noting the `rate` should match the integrand's exponential decay and
  that algebraic-decay integrands will not converge.

## [1.714.0] - 2026-09-12

### Added
- `gauss_laguerre.py`: `gauss_laguerre_nodes_weights` and
  `gauss_laguerre_integral` for integrals over the half-line against the ``e^{-x}``
  weight, exact for polynomials up to degree ``2n - 1``. Nodes and weights are the
  Golub-Welsch eigen-decomposition of the Laguerre Jacobi matrix (diagonal
  ``2k + 1``, off-diagonal ``k``); the wrapper integrates a general
  ``integral_0^inf g(x) dx`` by factoring out an exponential rate. Cross-checked:
  the moments reproduce the factorials ``m!``, and ``integral_0^inf e^{-2x} = 1/2``,
  ``integral_0^inf x^2 e^{-3x} = 2/27``, ``integral_0^inf e^{-x^2} = sqrt(pi)/2``,
  and ``integral_0^inf e^{-x}/(1+x) dx`` all match their closed forms.

## [1.713.0] - 2026-09-12

### Documentation
- README numerics section now documents Gauss-Hermite quadrature
  (`gauss_hermite_expectation`, `gauss_hermite_nodes_weights`) with a worked
  example, noting it is exact for polynomials up to degree `2n-1` and best for
  smooth integrands.

## [1.712.0] - 2026-09-12

### Added
- `gauss_hermite.py`: `gauss_hermite_nodes_weights` and
  `gauss_hermite_expectation` for Gaussian-weighted expectations
  `E[g(X)], X ~ N(mu, sigma^2)`. Nodes and weights use the probabilists'
  convention (weight the standard normal density, weights summing to one) and are
  generated by Golub-Welsch as the eigen-decomposition of the Hermite Jacobi
  matrix, exact for polynomials up to degree `2n - 1`. Cross-checked: weights sum
  to one, the standard-normal even moments (1, 3, 15) and odd moments (0) are
  reproduced, and the moment-generating function and the lognormal mean match their
  closed forms. Mirror pairs are folded to keep the nodes and weights exactly
  symmetric (removing eigensolver float error that would otherwise leak into odd
  moments).

## [1.711.0] - 2026-09-12

### Documentation
- README numerics section now shows `romberg` next to the other quadrature rules,
  with a note that it reaches machine precision in a few halvings on smooth
  integrands (and tanh-sinh for singular endpoints).

## [1.710.0] - 2026-09-12

### Added
- `quadrature.py`: `romberg` integration -- Richardson extrapolation on the
  trapezoid rule via the Romberg tableau, cancelling successive even powers of the
  step so a smooth integrand reaches machine precision in a few halvings, with
  early stopping when the diagonal converges. Cross-checked to 1e-12 against closed
  forms (`sin` over a half period, `exp`, a degree-7 polynomial, `arctan'` giving
  `pi/4`, a Gaussian against `erf`) and against composite Simpson on a smooth
  damped-cosine integrand.

## [1.709.0] - 2026-09-12

### Documentation
- README numerics section now shows `tanh_sinh` next to the other quadrature
  rules, with a note that it is the one to use for integrable endpoint
  singularities where Simpson and Gauss-Legendre lose accuracy.

## [1.708.0] - 2026-09-12

### Added
- `quadrature.py`: `tanh_sinh` double-exponential quadrature. The
  `x = tanh((pi/2) sinh(t))` change of variables clusters abscissae toward the
  endpoints and decays the weights super-fast, so it integrates functions with
  integrable endpoint singularities (`1/sqrt(x)`, `ln x`) that Simpson and
  Gauss-Legendre handle poorly. Cross-checked against closed forms: `integral_0^1
  x^{-1/2} = 2`, `integral_0^1 -ln x = 1`, `integral_0^1 ln(x)/sqrt(x) = -4`, the
  semicircle `integral_-1^1 sqrt(1-x^2) = pi/2`, and a Gaussian against `erf`.

## [1.707.0] - 2026-09-12

### Documentation
- README sizing section now documents the ruin module (`gamblers_ruin_probability`,
  `risk_of_ruin_units`, `ruin_probability_gbm`) with a worked example, alongside the
  Kelly sizing it complements.

## [1.706.0] - 2026-09-12

### Added
- `ruin.py`: risk-of-ruin and first-passage drawdown probabilities.
  `gamblers_ruin_probability` gives the classic unit-stake ruin probability,
  `risk_of_ruin_units` wraps it as a bankroll-in-units measure, and
  `ruin_probability_gbm` gives the chance a drifting log-equity ever falls by a
  given fraction, `(1 - loss)^{2 mu / sigma^2}`. Cross-checked: the gambler's-ruin
  formula matches a direct solve of the ruin recursion (and the fair-game linear
  case), and the geometric-Brownian-motion formula matches its first-passage
  derivation and a simulated lower bound.

## [1.705.0] - 2026-09-12

### Documentation
- README performance-metrics section now documents the Cornish-Fisher pair
  `cornish_fisher_var` and `cornish_fisher_expected_shortfall` with a worked
  example, noting both reduce to the Gaussian figures for a normal series and the
  expected shortfall never falls below the VaR.

## [1.704.0] - 2026-09-12

### Added
- `perfmetrics.py`: `cornish_fisher_expected_shortfall` gives the skew- and
  kurtosis-adjusted expected shortfall as the tail mean of the Cornish-Fisher
  expanded quantile, completing the pair with the existing `cornish_fisher_var`.
  Cross-checked: for a normal series it reduces to the Gaussian expected shortfall,
  it never falls below the Cornish-Fisher VaR, it matches an empirical tail average
  under mild non-normality to within 2%, and a fat left tail raises it.

## [1.703.0] - 2026-09-12

### Documentation
- README structural-credit section now documents `physical_distance_to_default`
  and `physical_default_probability` with a worked example, noting `mu = r`
  recovers the risk-neutral figures and a risky firm's physical default
  probability sits below the risk-neutral one.

## [1.702.0] - 2026-09-12

### Added
- `structural_credit.py`: `physical_distance_to_default` and
  `physical_default_probability` give the Merton default measures under the firm's
  real-world asset drift `mu` rather than the risk-free rate -- the Moody's-KMV
  distance to default. Cross-checked: setting `mu = r` reproduces the risk-neutral
  figures exactly, the physical probability matches a geometric-Brownian-motion
  Monte Carlo, and a higher drift (`mu > r`) puts the physical default probability
  below the risk-neutral one.

## [1.701.0] - 2026-09-12

### Documentation
- README portfolio-credit section now documents `vasicek_loss_pdf` and
  `vasicek_loss_expected_shortfall` with a worked example, noting the density
  integrates to one with mean `pd` and the expected shortfall sits at or above the
  VaR and rises with confidence and asset correlation.

## [1.700.0] - 2026-09-12

### Added
- `copula.py`: `vasicek_loss_pdf` and `vasicek_loss_expected_shortfall` complete the
  large-homogeneous-pool loss distribution alongside the existing CDF and quantile.
  The density is the closed-form derivative of the CDF; the expected shortfall is
  Tasche's (2002) closed form `Phi_2(Phi^{-1}(pd), -Phi^{-1}(q); sqrt(rho))/(1-q)`.
  Cross-checked: the density integrates to one with mean `pd` and matches a
  finite-difference of the CDF, and the expected shortfall matches a tail average of
  the quantile function and a single-factor Monte Carlo, staying at or above the VaR.

## [1.699.0] - 2026-09-12

### Documentation
- README reserving section now documents `mack_standard_error` alongside the
  chain-ladder, Bornhuetter-Ferguson and Cape-Cod methods: a worked example and
  the note that on Mack's Taylor-Ashe triangle it reproduces the reserve of
  18,680,856 and standard error of 2,447,095 to the dollar.

## [1.698.0] - 2026-09-12

### Added
- `mack.py`: `mack_standard_error` computes Mack's (1993) distribution-free
  standard error of chain-ladder reserves -- the per-accident-year and total
  mean-squared error of prediction, combining process and estimation error with
  the between-year correlation term, plus coefficients of variation. Validated
  against the Taylor-Ashe triangle from Mack's paper: total reserve 18,680,856
  and standard error 2,447,095 to the dollar, with the per-year standard errors
  matching his table.

## [1.697.0] - 2026-09-12

### Documentation
- README Structured-notes section now shows a `phoenix_autocall_mc` example
  (memory coupons, autocall and protection barriers) alongside the note-price note
  on the memory feature and coupon barrier.

## [1.696.0] - 2026-09-12

### Added
- `structured.py`: `phoenix_autocall_mc` Monte Carlo prices a Phoenix autocallable
  note -- conditional coupons paid when the spot is above a coupon barrier (with an
  optional memory/snowball feature), early par redemption at an autocall barrier,
  and down-and-in downside at maturity below a protection barrier. Cross-checked:
  the price is positive, the memory feature raises it, a lower coupon barrier pays
  more often, a higher coupon raises the value, and dropping the protection barrier
  never lowers it.

## [1.695.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `exponential_tail_factor` and
  `chain_ladder_with_tail`.

## [1.694.0] - 2026-09-12

### Added
- `chain_ladder.py`: `exponential_tail_factor` extrapolates a tail development
  factor by fitting exponential decay to the age-to-age factors' excess over one,
  and `chain_ladder_with_tail` applies it to capture development beyond the
  triangle. Cross-checked: decaying factors give a finite tail above 1, flat
  factors give 1, the tail scales the ultimates proportionally, a unit tail
  reproduces plain chain-ladder, and a non-decaying factor pattern is rejected.

## [1.693.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with the triangle utilities
  `incremental_to_cumulative`, `cumulative_to_incremental`, and `paid_to_date`.

## [1.692.0] - 2026-09-12

### Added
- `chain_ladder.py`: claims-triangle utilities -- `incremental_to_cumulative`,
  `cumulative_to_incremental`, and `paid_to_date`. Cross-checked: the two
  conversions round-trip, a cumulative triangle from positive increments is
  monotone, paid-to-date is the latest diagonal, and each row's increments sum to
  its cumulative total.

## [1.691.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `cape_cod` reserving.

## [1.690.0] - 2026-09-12

### Added
- `chain_ladder.py`: `cape_cod` (Stanard-Buhlmann) reserving estimates the expected
  loss ratio from the triangle and earned premium (total losses over used-up
  premium) rather than assuming an a-priori, then reserves like Bornhuetter-
  Ferguson. Cross-checked: the ELR equals losses over premium-weighted-by-
  development, a fully-developed year has zero reserve, reserves equal
  ``premium x ELR x undeveloped``, and Cape Cod matches BF with a
  ``premium x ELR`` a-priori.

## [1.689.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `bornhuetter_ferguson`.

## [1.688.0] - 2026-09-12

### Added
- `chain_ladder.py`: `bornhuetter_ferguson` reserving blends the chain-ladder
  `development_pattern` with an a-priori ultimate: reserve = a-priori x (1 - %
  developed). Cross-checked: the development pattern matches the reciprocal
  cumulative factors, a fully-developed year has zero reserve, reserves equal the
  a-priori times the undeveloped fraction, and using the chain-ladder ultimate as
  the a-priori reproduces the chain-ladder reserves exactly.

## [1.687.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `chain_ladder` loss reserving.

## [1.686.0] - 2026-09-12

### Added
- `chain_ladder.py`: `development_factors` and `chain_ladder` implement
  chain-ladder loss reserving -- volume-weighted age-to-age factors from a
  cumulative-claims triangle, projected to ultimate losses and IBNR reserves.
  Cross-checked on a hand-computed triangle: the factors, ultimates, and reserves
  match exactly, a fully-developed accident year has zero reserve, and every
  ultimate is at least the latest paid.

## [1.685.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `buhlmann_premium`,
  `credibility_factor`, and `buhlmann_straub_premium`.

## [1.684.0] - 2026-09-12

### Added
- `credibility.py`: Buhlmann and Buhlmann-Straub credibility (experience rating) --
  `buhlmann_k` (stiffness EPV/VHM), `credibility_factor` (Z = n/(n+k)),
  `buhlmann_premium`, and `buhlmann_straub_premium` (unequal exposures).
  Cross-checked: Z stays in [0,1], the premium lies between the individual and
  collective means, Z -> 1 as data grows and -> 0 with none, more within-risk noise
  lowers Z while more between-risk spread raises it, and Buhlmann-Straub matches
  Buhlmann on equal exposures.

## [1.683.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `aggregate_var` and
  `aggregate_tvar`.

## [1.682.0] - 2026-09-12

### Added
- `panjer.py`: `aggregate_var` and `aggregate_tvar` compute the Value-at-Risk
  (loss quantile) and Tail-VaR / CTE (expected loss beyond VaR) of a discrete
  aggregate-loss distribution. Cross-checked: the VaR is the exact grid quantile,
  TVaR is at least the VaR, both are non-decreasing in confidence, and TVaR at
  near-zero confidence equals the mean.

## [1.681.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `wang_premium` and
  `proportional_hazard_premium`.

## [1.680.0] - 2026-09-12

### Added
- `distortion.py`: distortion (spectral) risk pricing on a discrete loss
  distribution -- `wang_premium` (Wang transform, ``g(u) = Phi(Phi^{-1}(u) +
  lambda)``) and `proportional_hazard_premium` (``g(u) = u^{1/rho}``), plus
  `expected_loss`. Cross-checked: both reduce to the expected loss at zero
  distortion (``lambda = 0`` / ``rho = 1``), a positive load raises the premium
  monotonically and loads the tail, and a negative Wang ``lambda`` discounts below
  the mean.

## [1.679.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `panjer_negative_binomial`.

## [1.678.0] - 2026-09-12

### Added
- `panjer.py`: `panjer_negative_binomial` extends the Panjer recursion to
  compound negative-binomial claim counts (over-dispersed frequency, capturing
  contagion). Cross-checked: the distribution sums to 1, the mean and variance
  match the compound-NB identities, and at a matched mean it is more dispersed
  (heavier-tailed) than the Poisson aggregate.

## [1.677.0] - 2026-09-12

### Documentation
- README: extended the Actuarial section with `panjer_poisson`,
  `stop_loss_premium`, and `layer_expected_loss`.

## [1.676.0] - 2026-09-12

### Added
- `panjer.py`: `panjer_poisson` computes the compound-Poisson aggregate-loss
  distribution exactly by Panjer's recursion, with `aggregate_mean`,
  `stop_loss_premium`, and `layer_expected_loss` (excess-of-loss reinsurance layer
  cost). Cross-checked: the distribution sums to 1, the mean and variance match the
  compound-Poisson identities ``lam E[X]`` and ``lam E[X^2]``, ``g_0`` and ``g_1``
  match their closed forms, the stop-loss premium at zero retention equals the mean
  and decreases with retention, and a full-width layer equals the mean.

## [1.675.0] - 2026-09-12

### Documentation
- README: extended the Weather derivatives section with `degree_day_digital`.

## [1.674.0] - 2026-09-12

### Added
- `weather.py`: `degree_day_digital` prices a binary degree-day option paying a
  fixed amount on a strike breach -- the discounted breach probability times the
  payout under the normal index model. Cross-checked: the call and put digitals sum
  to the discounted payout, an ATM digital is worth half, a deep-ITM call pays the
  full discounted amount, the value is monotone in the strike, and zero vol gives a
  step function.

## [1.673.0] - 2026-09-12

### Documentation
- README: extended the Commodities section with `crack_spread_option`.

## [1.672.0] - 2026-09-12

### Added
- `commodity.py`: `crack_spread_option` prices a refinery crack-spread option
  (weighted refined-product basket minus crude, e.g. 3:2:1) via the normal-model
  Bachelier spread, aggregating the product slate into one leg. Cross-checked: a
  single unit-weight product reduces exactly to `bachelier_spread_option`, a 3:2:1
  crack is positive, put-call parity holds, and a higher product volatility raises
  the price.

## [1.671.0] - 2026-09-12

### Documentation
- README: extended the Fixed income section with `fra_forward_rate` and
  `fra_value`.

## [1.670.0] - 2026-09-12

### Added
- `fra.py`: `fra_forward_rate` and `fra_value` price a forward rate agreement off a
  discount curve -- the simple forward rate over the accrual period and the
  discounted settlement value. Cross-checked: the value is zero at the fair
  (forward) rate, a payer gains when the forward exceeds the contract rate,
  receiver equals negative payer, the notional scales it linearly, and an upward
  curve gives a higher forward rate further out.

## [1.669.0] - 2026-09-12

### Documentation
- README: extended the Fixed income section with `vanilla_swap_value`,
  `single_curve_par_swap_rate`, and `swap_annuity`.

## [1.668.0] - 2026-09-12

### Added
- `swap.py`: single-curve vanilla interest-rate swap valuation --
  `vanilla_swap_value`, `single_curve_par_swap_rate`, and `swap_annuity` (PV01),
  using the telescoping float-leg PV off one discount curve. Cross-checked: the swap
  is worth zero at the par rate, receiver equals negative payer, the value rises as
  the fixed rate falls, the notional scales it linearly, a forward-starting swap is
  zero at its par rate, and the annuity is positive.

## [1.667.0] - 2026-09-12

### Documentation
- README: extended the Fixed income section with `par_yield` and `par_bond_price`.

## [1.666.0] - 2026-09-12

### Added
- `par_yield.py`: `par_yield` computes the par coupon rate (= par swap rate) for a
  maturity from a discount curve, and `par_bond_price` prices a coupon bond off the
  curve. Cross-checked: on a flat curve the par yield equals the flat rate at every
  maturity, a bond bearing the par coupon prices to exactly par, the one-year par
  equals the spot rate, an upward curve gives a par yield between the short and long
  zeros, and a higher coupon yields a premium bond.

### Reverted
- Dropped a `cev_implied_sigma` attempt: the CEV pricer's noncentral chi-square
  series hangs at the large volatilities a bisection bracket sweeps, so it cannot be
  safely inverted without first bounding the pricer.

## [1.665.0] - 2026-09-12

### Documentation
- README: extended the Bachelier / displaced-diffusion section with
  `displaced_diffusion_implied_vol`.

## [1.664.0] - 2026-09-12

### Added
- `displaced.py`: `displaced_diffusion_implied_vol` inverts the displaced-diffusion
  price for the volatility by bisection (the price is monotone in ``sigma``).
  Cross-checked: it round-trips with the pricer across a sigma x shift grid,
  coincides with the Black-Scholes implied vol at zero shift, round-trips on the
  put side, and rejects out-of-band quotes.

## [1.663.0] - 2026-09-12

### Documentation
- README: extended the Bachelier section with `black_to_normal_vol` and
  `normal_to_black_vol` (lognormal <-> normal vol conversion).

## [1.662.0] - 2026-09-12

### Added
- `vol_convert.py`: `black_to_normal_vol` and `normal_to_black_vol` convert between
  Black-76 (lognormal) and Bachelier (normal) implied volatilities by matching the
  option price under both models. Cross-checked: the ATM conversion recovers the
  leading-order ``sigma_N ~ sigma_B * F``, the round-trip is exact, both models
  reproduce the same price at the converted vol, and OTM/put conversions round-trip.

## [1.661.0] - 2026-09-12

### Documentation
- README: extended the Covariance shrinkage section with the EWMA (RiskMetrics)
  covariance and correlation matrices.

## [1.660.0] - 2026-09-12

### Added
- `ewma_cov.py`: `ewma_covariance_matrix` and `ewma_correlation_matrix` estimate an
  exponentially-weighted (RiskMetrics, default ``lambda = 0.94``) covariance and
  correlation from a multi-asset return panel. Cross-checked: the covariance is
  symmetric and positive definite, the correlation has a unit diagonal and stays in
  ``[-1, 1]``, and a shared-factor pair shows high correlation while an independent
  asset shows ~0.

## [1.659.0] - 2026-09-12

### Documentation
- README: extended the Numerical utilities section with `levenberg_marquardt`
  nonlinear least-squares calibration.

## [1.658.0] - 2026-09-12

### Added
- `levenberg.py`: `levenberg_marquardt` fits parameters by nonlinear least squares,
  interpolating between Gauss-Newton and gradient descent via an adaptive damping
  factor, with the residual Jacobian taken by finite differences (only the model
  is needed). Cross-checked: it recovers the parameters of an exponential model
  exactly and from a poor starting point, fits noisy data close to the truth,
  solves a linear model exactly, and drives the residual to zero.

## [1.657.0] - 2026-09-12

### Documentation
- README: extended the Numerical utilities section with `gradient`, `hessian`, and
  `jacobian`.

## [1.656.0] - 2026-09-12

### Added
- `numdiff.py`: central-difference numerical differentiation -- `gradient` and
  `hessian` of a scalar function and the `jacobian` of a vector function, with
  per-coordinate steps scaled by the argument magnitude. Cross-checked against
  exact derivatives: the gradient and Hessian of a quadratic are recovered
  exactly, the Hessian is symmetric, the gradient vanishes at a minimum, and the
  Jacobian of a known vector map matches.

## [1.655.0] - 2026-09-12

### Documentation
- README: extended the Numerical utilities section with the 1-D minimizers
  `golden_section_min` and `brent_min`.

## [1.654.0] - 2026-09-12

### Added
- `minimize1d.py`: one-dimensional minimizers -- `golden_section_min` (robust
  derivative-free bracket shrinking) and `brent_min` (parabolic interpolation with
  a golden-section safeguard). Cross-checked: both locate the minima of a
  quadratic, cosine, quartic, and Gaussian to high precision and agree with each
  other; an invalid bracket raises.

## [1.653.0] - 2026-09-12

### Documentation
- README: extended the Exotic options section with `curran_asian` (sharper
  arithmetic Asian).

## [1.652.0] - 2026-09-12

### Added
- `curran_asian.py`: `curran_asian` prices an arithmetic-average Asian option by
  Curran's (1994) geometric-conditioning approximation -- conditioning the average
  on the exactly-lognormal geometric mean and integrating analytically, sharper
  than the two-moment Turnbull-Wakeman match. Cross-checked: it reduces to Black
  for a single fixing, matches a Monte Carlo arithmetic Asian to ~0.001 at 12
  fixings (closer than Turnbull-Wakeman), satisfies put-call parity, and orders
  ITM above OTM.

## [1.651.0] - 2026-09-12

### Documentation
- README: extended the Two-asset options section with `rainbow_option_mc` for
  n-asset best-of / worst-of options.

## [1.650.0] - 2026-09-12

### Added
- `rainbow_n.py`: `rainbow_option_mc` prices an ``n``-asset best-of / worst-of
  call or put by correlated-GBM Monte Carlo (Cholesky of the correlation matrix),
  generalizing the two-asset rainbow options. Cross-checked: the two-asset case
  matches the closed-form best-of/worst-of prices, ``worst-of <= single-asset <=
  best-of`` holds, and a three-asset worst-of is below its best-of.

## [1.649.0] - 2026-09-12

### Documentation
- README: extended the Two-asset options section with `levy_basket_option` for
  baskets of more than two assets.

## [1.648.0] - 2026-09-12

### Added
- `levy_basket.py`: `levy_basket_option` prices an ``n``-asset basket option by
  Levy's lognormal moment matching (matching the basket forward's first two
  moments to a single lognormal), generalizing the existing two-asset
  ``basket_option``. Cross-checked: it reduces to Black-Scholes for a single asset,
  matches a Monte Carlo basket price to within ~0.05, satisfies put-call parity
  exactly, and rises with volatility.

## [1.647.0] - 2026-09-12

### Documentation
- README: added a "kth-to-default basket" section documenting
  `basket_default_distribution` and `kth_to_default_probability`, the 120th
  documented section; TOC regenerated to 120 entries.

## [1.646.0] - 2026-09-12

### Added
- `basket_default.py`: `basket_default_distribution` gives the full distribution of
  the number of defaults in a homogeneous basket under the one-factor Gaussian
  copula (conditional binomial integrated over the common factor), and
  `kth_to_default_probability` gives the ``P(>= k defaults)`` trigger of a
  kth-to-default swap. Cross-checked: the distribution sums to 1, the trigger is
  monotone decreasing in ``k``, the ``rho = 0`` first-to-default matches
  ``1 - (1-pd)^n``, correlation clusters defaults (senior triggers rise, first-to-
  default falls with ``rho``), and the mean equals ``n * pd``.

## [1.645.0] - 2026-09-12

### Documentation
- README: extended the Hull-White section with `hw_swaption` (Jamshidian
  swaption).

## [1.644.0] - 2026-09-12

### Added
- `hull_white.py`: `hw_swaption` prices a European payer/receiver swaption under
  Hull-White via the Jamshidian decomposition -- solve for the critical short rate
  at which the underlying coupon bond is at par, then sum options on each
  zero-coupon cashflow struck at that rate. Cross-checked: payer and receiver are
  positive, ``payer - receiver`` equals the forward swap value (put-call parity),
  the ATM payer equals the ATM receiver, and a higher volatility raises the price.

## [1.643.0] - 2026-09-12

### Documentation
- README: extended the Hull-White section with `hw_cap` and `hw_floor`.

## [1.642.0] - 2026-09-12

### Added
- `hull_white.py`: `hw_cap` / `hw_floor` (and the `hw_caplet` / `hw_floorlet`
  building blocks) price interest-rate caps and floors under Hull-White as
  portfolios of zero-coupon bond options off the initial curve. Cross-checked: the
  cap equals the sum of its caplets, ``cap - floor`` equals the fixed-vs-float swap
  value (put-call parity) to 1e-9, the price rises with volatility and falls with
  the strike.

## [1.641.0] - 2026-09-12

### Documentation
- README: extended the Hull-White section with `hw_bond_option` (analytic
  zero-coupon bond option).

## [1.640.0] - 2026-09-12

### Added
- `hull_white.py`: `hw_bond_option` prices a European option on a zero-coupon bond
  under Hull-White in closed form (Jamshidian), using the initial discount curve
  and the analytic forward-bond-price volatility. Cross-checked: put-call parity
  holds exactly (``c - p = P(Tb) - K P(To)``), ATM call equals put, a higher
  volatility raises the price, a deep-in-the-money call approaches its intrinsic
  value, and a zero-vol ATM option is worthless.

## [1.639.0] - 2026-09-12

### Documentation
- README: added a "Hull-White (fitted to a curve)" section documenting
  `hw_zero_from_curve` and `hw_B`; TOC regenerated to 119 entries.

## [1.638.0] - 2026-09-12

### Added
- `hull_white.py`: `hw_zero_from_curve` prices a zero-coupon bond under the
  Hull-White (extended Vasicek) model fitted to an arbitrary initial discount
  curve, so the model reprices the input term structure exactly; `hw_B` gives the
  ``B(t, T)`` factor. Cross-checked: at ``t = 0`` it refits flat, upward, and
  Ho-Lee-limit curves to 1e-4, the ``B`` factor tends to ``T - t`` as mean
  reversion vanishes, and fitted discount factors stay in ``(0, 1)``.

## [1.637.0] - 2026-09-12

### Documentation
- README: added a "Black-Karasinski short rate" section documenting
  `bk_zero_coupon_bond`; TOC regenerated to 118 entries.

## [1.636.0] - 2026-09-12

### Added
- `black_karasinski.py`: `bk_zero_coupon_bond` prices a zero-coupon bond under the
  Black-Karasinski log-normal short-rate model (``d ln r = kappa(theta - ln r)dt +
  sigma dW``, so rates stay positive) on a Hull-White-style trinomial tree in the
  log rate. Cross-checked: the price stays in ``(0, 1]`` and near the flat
  discount factor, falls as the rate rises, approaches 1 at short maturity,
  converges across step counts, and falls as volatility rises (Jensen raises the
  expected rate).

## [1.635.0] - 2026-09-12

### Documentation
- README: extended the Decision stump section with the random forest
  (`fit_random_forest`, `predict_random_forest`).

## [1.634.0] - 2026-09-12

### Added
- `random_forest.py`: `fit_random_forest` trains an ensemble of CART trees on
  bootstrap resamples and `predict_random_forest` classifies by majority vote.
  Cross-checked: a separable set is classified at 100%, runs are reproducible per
  seed, the forest generalizes at least as well as a single tree on held-out noisy
  data (0.69 vs 0.63), and a one-tree forest degenerates to a single tree.

## [1.633.0] - 2026-09-12

### Documentation
- README: extended the Decision stump section with the full CART tree
  (`fit_decision_tree`, `predict_decision_tree`, `tree_depth`).

## [1.632.0] - 2026-09-12

### Added
- `decision_tree.py`: `fit_decision_tree` grows a CART classification tree by
  recursively splitting on the best Gini-reducing feature/threshold up to a max
  depth or min node size, `predict_decision_tree` walks it to a leaf, and
  `tree_depth` reports its depth. Cross-checked: a depth-1 tree matches the
  decision stump, a depth-2+ tree solves XOR (which a stump cannot), a separable
  three-class set is fit perfectly, a deeper tree fits at least as well, and a pure
  node becomes a leaf.

## [1.631.0] - 2026-09-12

### Documentation
- README: added a "Decision stump" section documenting `fit_decision_stump`,
  `predict_decision_stump`, and `gini_impurity`; TOC regenerated to 117 entries.

## [1.630.0] - 2026-09-12

### Added
- `decision_stump.py`: `gini_impurity` and a one-split classifier --
  `fit_decision_stump` finds the feature and threshold minimizing the
  size-weighted Gini of the two children, and `predict_decision_stump` applies it.
  Cross-checked: Gini is 0 for a pure node and 0.5 for a 50/50 split, a
  single-feature-separable set is split perfectly (Gini 0), noisy data still
  reduces Gini below the parent with >90% accuracy, and identical rows fall back
  to the global majority.

## [1.629.0] - 2026-09-12

### Documentation
- README: extended the Feature scaling section with categorical encoding
  (`fit_label_encoder`, `label_encode`, `one_hot_encode`).

## [1.628.0] - 2026-09-12

### Added
- `encoding.py`: categorical feature encoding -- `fit_label_encoder`,
  `label_encode` / `label_decode` (integer codes), and `one_hot_encode` (0/1
  indicator columns). Cross-checked: categories are sorted and encode to
  contiguous indices, codes round-trip back to labels, unseen categories map to
  ``-1`` / ``None`` (label) or an all-zero row (one-hot), and one-hot rows sum to
  1 for known categories.

## [1.627.0] - 2026-09-12

### Documentation
- README: added a "Gaussian naive Bayes" section documenting `fit_gaussian_nb`,
  `predict_gaussian_nb`, and `predict_proba_gaussian_nb`; TOC regenerated to 116
  entries.

## [1.626.0] - 2026-09-12

### Added
- `naive_bayes.py`: Gaussian naive Bayes -- `fit_gaussian_nb` estimates each
  class's prior and per-feature mean/variance, `predict_gaussian_nb` picks the
  max-log-posterior class, and `predict_proba_gaussian_nb` returns softmaxed class
  probabilities. Cross-checked: the fitted class means and priors match the
  generating distribution, separable classes are classified at ~100% accuracy,
  the posterior sums to 1, and it favors the nearby class.

## [1.625.0] - 2026-09-12

### Documentation
- README: added a "k-nearest neighbors" section documenting `knn_classify` and
  `knn_regress`; TOC regenerated to 115 entries.

## [1.624.0] - 2026-09-12

### Added
- `knn.py`: `knn_classify` (majority vote of the ``k`` nearest neighbors, ties
  broken by nearest total distance) and `knn_regress` (mean target of the ``k``
  nearest). Cross-checked: ``k=1`` reproduces the training labels exactly and the
  nearest target for regression, separable classes and a linear target are
  predicted correctly, and out-of-range ``k`` or mismatched inputs raise.

## [1.623.0] - 2026-09-12

### Documentation
- README: added a "Hierarchical clustering" section documenting `linkage` and
  `fcluster`; TOC regenerated to 114 entries.

## [1.622.0] - 2026-09-12

### Added
- `agglomerative.py`: `linkage` builds an agglomerative cluster tree (single /
  complete / average linkage) with monotone merge distances, and `fcluster` cuts it
  into a target number of clusters or below a distance threshold. Cross-checked:
  there are ``n-1`` merges with non-decreasing distances, cutting into 3 recovers
  three separated blobs purely, the distance-threshold cut agrees, the extreme cuts
  give 1 and ``n`` clusters, and the complete-linkage final merge is at least the
  single-linkage one.

## [1.621.0] - 2026-09-12

### Documentation
- README: extended the K-means clustering section with `silhouette_score`
  (label-free cluster-quality scoring for choosing `k`).

## [1.620.0] - 2026-09-12

### Added
- `silhouette.py`: `silhouette_score` and `silhouette_samples` measure clustering
  quality -- for each point, ``(b - a) / max(a, b)`` comparing its own-cluster
  cohesion to the nearest other cluster. Cross-checked: the correct ``k`` scores
  high (~0.93) and beats a wrong ``k``, values stay in ``[-1, 1]``, a singleton
  cluster scores 0, tight well-separated clusters approach 1, and random labels
  score near 0.

## [1.619.0] - 2026-09-12

### Documentation
- README: added a "K-means clustering" section documenting `kmeans`; TOC
  regenerated to 113 entries.

## [1.618.0] - 2026-09-12

### Added
- `kmeans.py`: `kmeans` clusters points by Lloyd's algorithm with k-means++
  seeding, returning labels, centroids, inertia, and the iteration count.
  Cross-checked: three well-separated blobs are recovered as pure clusters with
  centroids at the true centers, ``k=1`` gives the global mean, runs are
  reproducible per seed, and inertia decreases as ``k`` grows.

## [1.617.0] - 2026-09-12

### Documentation
- README: extended the Feature scaling section with `polynomial_features` (monomial
  and interaction expansion feeding the regression/classification models).

## [1.616.0] - 2026-09-12

### Added
- `poly_features.py`: `polynomial_features` expands a feature matrix into all
  monomials up to a given degree (squares, cubes, cross-products), with optional
  bias and an ``interaction_only`` mode. Cross-checked: degree 1 with a bias is the
  identity plus a ones column, degree 2 includes the squares and pairwise products,
  the feature count equals ``C(p+d, d)``, and ``interaction_only`` drops the pure
  powers.

## [1.615.0] - 2026-09-12

### Documentation
- README: added a "Feature scaling" section documenting `fit_standardize`,
  `fit_min_max`, `fit_robust`, `scale_transform`, and `scale_inverse_transform`;
  TOC regenerated to 112 entries.

## [1.614.0] - 2026-09-12

### Added
- `scaling.py`: feature scalers with fit / transform / inverse discipline --
  `fit_standardize` (z-score), `fit_min_max` ([0, 1]), `fit_robust` (median / IQR),
  applied by `scale_transform` and undone by `scale_inverse_transform` (exported
  under those names to avoid clashing with the copula ``inverse_transform``).
  Cross-checked: standardized columns have mean 0 and std 1, min-max maps to exactly
  [0, 1], both round-trip through the inverse, and the robust scaler's median center
  is resistant to an injected outlier.

## [1.613.0] - 2026-09-12

### Documentation
- README: added a "Cross-validation" section documenting `k_fold_indices`,
  `train_test_split`, and `cross_val_score`; TOC regenerated to 111 entries.

## [1.612.0] - 2026-09-12

### Added
- `cross_validation.py`: model-agnostic resampling -- `k_fold_indices` (disjoint
  folds, optionally shuffled), `train_test_split` (single fractional split), and
  `cross_val_score` (run a caller-supplied fit/score over the folds). Cross-checked:
  the test folds tile the data exactly and are disjoint, train/test are
  complementary within each fold, fold sizes differ by at most one, a seeded
  shuffle is reproducible and seed-sensitive, and cross_val_score returns one score
  per fold.

## [1.611.0] - 2026-09-12

### Documentation
- README: added a "Classification metrics" section documenting `roc_auc`,
  `precision_recall_f1`, `confusion_matrix`, `log_loss`, and `brier_score`; TOC
  regenerated to 110 entries.

## [1.610.0] - 2026-09-12

### Added
- `classification_metrics.py`: binary-classification evaluation -- `roc_auc` (the
  Mann-Whitney rank AUC), `confusion_matrix`, `precision_recall_f1`, `log_loss`
  (cross-entropy), and `brier_score`. Cross-checked: a perfect ranking gives AUC 1
  and a reversed one 0, random scores give ~0.5, ties count as half, the confusion
  matrix and precision/recall/F1 are exact on separable data, log-loss is ~0 for
  confident-correct predictions, and the Brier score is 0 for exact probabilities
  and 0.25 for all-0.5 guesses.

## [1.609.0] - 2026-09-12

### Documentation
- README: added a "Logistic regression" section documenting `fit_logistic` and
  `predict_proba`; TOC regenerated to 109 entries.

## [1.608.0] - 2026-09-12

### Added
- `logistic.py`: `fit_logistic` fits a binary logistic regression by iteratively
  reweighted least squares (Newton-Raphson), with a small ridge term for stability
  under separation, and `predict_proba` returns fitted class probabilities -- the
  standard default-probability / classification model. Cross-checked: it recovers
  the generating coefficients ([0.55, 2.07, -1.00] for a [0.5, 2, -1] logit), the
  probabilities lie in (0, 1), classification accuracy beats the baseline, perfectly
  separable data is classified ~perfectly, and a positive coefficient makes the
  probability monotone in that feature.

## [1.607.0] - 2026-09-12

### Documentation
- README: extended the OLS regression section with `ridge_regression` (L2-penalized
  regression for collinear or numerous regressors).

## [1.606.0] - 2026-09-12

### Added
- `ridge.py`: `ridge_regression` fits L2-penalized (Tikhonov) least squares,
  ``beta = (X'X + lambda I)^{-1} X'y``, shrinking the slope coefficients toward
  zero (the intercept is left unpenalized). Cross-checked: ``lambda = 0`` reproduces
  OLS exactly, a larger ``lambda`` monotonically shrinks the slopes (with the
  intercept converging to the mean of ``y``), and ridge stays solvable under
  perfect collinearity where OLS is singular.

## [1.605.0] - 2026-09-12

### Documentation
- README: added an "OLS regression" section documenting `ols_fit`, placed before
  the finance-specific factor-model regression; TOC regenerated to 108 entries.

## [1.604.0] - 2026-09-12

### Added
- `ols.py`: `ols_fit` fits a multivariate OLS regression by the normal equations
  and returns the standard diagnostics -- coefficients, standard errors,
  t-statistics, R-squared, adjusted R-squared, the F-statistic, and residuals (with
  an optional intercept column). Cross-checked: an exact line gives R-squared 1 and
  zero residuals, a noisy multivariate model recovers its coefficients with
  significant t-stats and a large F, adjusted R-squared never exceeds R-squared,
  and the single-regressor slope matches the closed form.

## [1.603.0] - 2026-09-12

### Documentation
- README: added a "Forecast accuracy" section documenting `mae`, `rmse`, `mape`,
  `smape`, and `mase`; TOC regenerated to 107 entries.

## [1.602.0] - 2026-09-12

### Added
- `forecast_metrics.py`: forecast-accuracy metrics -- `mae`, `rmse`, `mape`,
  `smape` (symmetric, bounded), and `mase` (mean absolute scaled error vs a
  seasonal-naive benchmark). Cross-checked: a perfect forecast scores zero on all,
  known inputs give the exact hand values, RMSE >= MAE, MAPE is scale-invariant,
  sMAPE is bounded in [0, 2], and MASE equals 1 for a naive-equivalent forecast and
  drops below 1 for a perfect one.

## [1.601.0] - 2026-09-12

### Documentation
- README: added an "Exponential smoothing (Holt / Holt-Winters)" section
  documenting `holt_linear` and `holt_winters_add`; TOC regenerated to 106 entries.

## [1.600.0] - 2026-09-12

### Added
- `holt_winters.py`: exponential smoothing -- `holt_linear` (double smoothing with
  a level and trend; straight-line forecast) and `holt_winters_add` (triple
  smoothing with an additive seasonal component of a given period). Cross-checked:
  Holt recovers a pure linear trend and forecasts a straight line; Holt-Winters
  reproduces a repeating seasonal pattern exactly and recovers an added linear
  trend under seasonality.

## [1.599.0] - 2026-09-12

### Documentation
- README: extended the Autocorrelation section with `select_ar_order` (automatic
  AR order selection by AIC / BIC).

## [1.598.0] - 2026-09-12

### Added
- `ar_select.py`: `select_ar_order` picks the AR order that minimizes AIC or BIC
  over a range of candidates, and `ar_information_criteria` returns the criteria of
  a single fit (from the concentrated Gaussian log-likelihood). Cross-checked: BIC
  recovers the true order of simulated AR(1) and AR(2) series, the BIC-selected
  order is never larger than the AIC one (heavier penalty), and AIC drops at the
  true order.

## [1.597.0] - 2026-09-12

### Documentation
- README: extended the Autocorrelation section with `fit_ar_yule_walker` and
  `ar_forecast` (fit and project an AR model after identifying its order).

## [1.596.0] - 2026-09-12

### Added
- `ar_model.py`: `fit_ar_yule_walker` fits an AR(p) model by the Yule-Walker /
  Durbin-Levinson method (returning coefficients, intercept, and innovation
  variance), and `ar_forecast` iterates the deterministic AR recursion forward.
  Cross-checked: it recovers the coefficients of simulated AR(1) and AR(2) series
  (and the AR(1) intercept and unit noise variance), the forecast mean-reverts
  toward the process mean, and the one-step forecast matches the recursion exactly.

## [1.595.0] - 2026-09-12

### Documentation
- README: added an "Autocorrelation (ACF / PACF)" section documenting `acf` and
  `pacf`; TOC regenerated to 105 entries.

## [1.594.0] - 2026-09-12

### Added
- `acf.py`: `acf` (autocorrelation function) and `pacf` (partial autocorrelation
  via the Durbin-Levinson recursion) for ARMA model identification. Cross-checked
  against theory: the ACF of an AR(1) decays geometrically as ``phi^k`` while its
  PACF cuts off after lag 1; the ACF of an MA(1) cuts off after lag 1 (matching
  ``theta/(1+theta^2)``); white noise is flat beyond lag 0.

## [1.593.0] - 2026-09-12

### Documentation
- README: extended the Performance metrics section with `kappa_ratio` and
  `upside_potential_ratio` (generalized downside-risk-adjusted return).

## [1.592.0] - 2026-09-12

### Added
- `downside_ratios.py`: generalized downside-risk performance measures --
  `kappa_ratio` (Kaplan-Knowles Kappa of any order about a target; order 1 is the
  Omega-Sharpe ratio, order 2 the per-period Sortino), `upside_potential_ratio`,
  and the underlying `lower_partial_moment`. Cross-checked: Kappa-2 times the
  annualization factor reproduces the existing Sortino ratio exactly, Kappa
  decreases with order (deeper shortfalls penalized more), a higher target lowers
  it, a symmetric mean-zero series gives ~0, and the upside-potential ratio is
  positive for a positive-drift series.

## [1.591.0] - 2026-09-12

### Documentation
- README: extended the Performance metrics section with `drawdown_analytics`
  (depth, time-to-recovery, longest underwater).

## [1.590.0] - 2026-09-12

### Added
- `drawdown.py`: `drawdown_analytics` reports the depth-and-duration summary of a
  return series' drawdowns -- the deepest peak-to-trough drop with its peak /
  trough / recovery indices, the time to recovery (``None`` if never), and the
  longest underwater stretch. Cross-checked: the depth matches
  :func:`max_drawdown`, peak precedes trough precedes recovery, a monotone-up
  series has no drawdown, a terminal crash never recovers, and the longest
  underwater count is exact.

## [1.589.0] - 2026-09-12

### Documentation
- README: extended the Gaussian-copula sampling section with
  `student_t_copula_sample` and its tail-dependence behavior.

## [1.588.0] - 2026-09-12

### Added
- `t_copula_sample.py`: `student_t_copula_sample` draws dependent uniforms from a
  Student-t copula -- same rank correlation as the Gaussian copula but with
  symmetric tail dependence that grows as the degrees of freedom fall, the standard
  fix for underestimated joint-crash risk. Cross-checked: the rank correlation
  tracks the target, each margin is uniform, a lower ``df`` gives stronger lower
  tail dependence, a large ``df`` reproduces the Gaussian copula's tail behavior,
  and an identity target gives ~0 rank correlation.

## [1.587.0] - 2026-09-12

### Documentation
- README: added a "Gaussian-copula sampling" section documenting
  `gaussian_copula_sample` and `inverse_transform`; TOC regenerated to 104 entries.

## [1.586.0] - 2026-09-12

### Added
- `copula_sample.py`: `gaussian_copula_sample` draws dependent uniforms from a
  Gaussian copula with a target correlation matrix (correlate standard normals via
  the Cholesky factor, then map back through the normal CDF), and
  `inverse_transform` pushes those uniforms through any inverse-CDF to build
  correlated draws from arbitrary marginals. Cross-checked: the sample rank
  correlation matches the target (0.69 vs 0.7), each margin is uniform (mean ~0.5,
  in (0,1)), independence gives ~0 rank correlation, a negative target gives a
  negative rank correlation, and a non-positive-definite matrix is rejected.

## [1.585.0] - 2026-09-12

### Documentation
- README: added a "Goodness of fit (Jarque-Bera / KS)" section documenting
  `jarque_bera_test` and `ks_two_sample`; TOC regenerated to 103 entries.

## [1.584.0] - 2026-09-12

### Added
- `gof_tests.py`: goodness-of-fit tests -- `jarque_bera_test` returns the
  Jarque-Bera normality statistic with its chi-square(2) p-value, and
  `ks_two_sample` is the two-sample Kolmogorov-Smirnov test (max gap between the
  empirical CDFs, with the asymptotic Kolmogorov p-value). Cross-checked: a normal
  sample is not rejected while a heavy-tailed one is (JB ~ 1200); KS does not reject
  two samples from the same distribution but rejects both a location shift and a
  scale difference.

## [1.583.0] - 2026-09-12

### Documentation
- README: added a "Serial correlation (Ljung-Box / Durbin-Watson)" section
  documenting `ljung_box`, `box_pierce`, and `durbin_watson`; TOC regenerated to
  102 entries.

## [1.582.0] - 2026-09-12

### Added
- `serial_correlation.py`: portmanteau serial-correlation tests -- `ljung_box`
  (with its small-sample refinement), `box_pierce`, and `durbin_watson`, with a
  self-contained chi-square survival function (regularized incomplete gamma) for
  the p-values. Cross-checked: the chi-square SF matches reference points, white
  noise fails to reject (p > 0.05, DW ~ 2), an AR(1) rejects strongly (Q ~ 790,
  p < 1e-6, DW < 1), Ljung-Box exceeds Box-Pierce, and a negative-AR series gives
  DW > 3.

## [1.581.0] - 2026-09-12

### Documentation
- README: extended the Bootstrap and jackknife section with
  `moving_block_bootstrap_ci` alongside the stationary block bootstrap.

## [1.580.0] - 2026-09-12

### Added
- `resample.py`: `moving_block_bootstrap_ci` adds the Kunsch fixed-length
  moving-block bootstrap confidence interval for serially-correlated data
  (overlapping blocks of length ``block``, wrapping at the end). Cross-checked: the
  interval brackets the sample mean, it is wider than the IID `bootstrap_ci` for a
  positively autocorrelated mean (which the IID bootstrap under-covers), it agrees
  with the IID width when ``block = 1``, and it matches the IID interval on
  genuinely independent data.

## [1.579.0] - 2026-09-12

### Documentation
- README: extended the rank-dependence section with tail dependence
  (`upper_tail_dependence`, `lower_tail_dependence`, `exceedance_correlation`).

## [1.578.0] - 2026-09-12

### Added
- `tail_dependence.py`: empirical tail dependence and exceedance correlation --
  `upper_tail_dependence` / `lower_tail_dependence` estimate ``P(U > q | V > q)``
  (and the lower analogue) on the rank/uniform scale, and `exceedance_correlation`
  is the Pearson correlation on joint-tail observations only. Cross-checked: a
  comonotone pair has tail dependence ~1, independent margins match the ``1 - q``
  null and shrink toward 0 at more extreme thresholds, a common-shock pair shows
  positive tail dependence, and the exceedance correlation is ~1 for a comonotone
  pair.

## [1.577.0] - 2026-09-12

### Documentation
- README: added a "Rank dependence (Kendall / Spearman)" section documenting
  `kendall_tau`, `spearman_rho`, and `pseudo_observations`; TOC regenerated to 101
  entries.

## [1.576.0] - 2026-09-12

### Added
- `copula_stats.py`: rank-based dependence measures -- `kendall_tau` (concordant
  minus discordant pairs), `spearman_rho` (Pearson correlation of ranks, with tie
  averaging), and `pseudo_observations` (scaled-rank empirical-copula transform).
  Cross-checked: a strictly increasing relationship gives tau = rho = 1, a
  decreasing one gives -1, Spearman is invariant to a monotone (exp) transform of a
  margin, independent samples give ~0, and the pseudo-observations stay in the open
  unit interval with correct rank scaling and tie handling.

## [1.575.0] - 2026-09-12

### Documentation
- README: added an "Entropy pooling (views on scenarios)" section documenting
  `entropy_pooling_mean` and `relative_entropy`, the 100th documented section; TOC
  regenerated to 100 entries.

## [1.573.0] - 2026-09-12

### Documentation
- README: added a "Spectral analysis" section documenting `periodogram`,
  `dominant_frequency`, and `spectral_energy`; TOC regenerated to 99 entries.

## [1.572.0] - 2026-09-12

### Added
- `spectral.py`: spectral analysis of a real signal -- `dft` (direct discrete
  Fourier transform), `periodogram` (one-sided power spectrum), `dominant_frequency`
  (largest non-DC peak), and `spectral_energy`. Cross-checked: a pure sine peaks
  exactly at its frequency (period recovered to machine precision), Parseval's
  theorem holds (spectral energy equals the time-domain sum of squares), a constant
  signal has power only at DC, white noise spreads across frequencies with no
  dominant peak, and a two-tone signal picks the stronger tone.

## [1.571.0] - 2026-09-12

### Documentation
- README: added a "Structural breaks (CUSUM / Chow)" section documenting
  `cusum_mean`, `cusum_break_detected`, and `chow_test`; TOC regenerated to 98
  entries.

## [1.570.0] - 2026-09-12

### Added
- `structural_break.py`: structural-break diagnostics -- `cusum_mean` /
  `cusum_break_detected` (the standardized CUSUM of deviations from the mean,
  compared to Kolmogorov Brownian-bridge critical values) and `chow_test` (an F
  test for a mean break at a known point). Cross-checked: a stable series stays in
  the band on the large majority of runs (~5% false-positive rate), an injected
  level shift breaches the band (max |CUSUM| ~ 6.6 vs a 1.36 band), the Chow F is
  huge at the true break (>1000) and small on a stable series.

## [1.569.0] - 2026-09-12

### Documentation
- README: added a "Cointegration (ADF / Engle-Granger)" section documenting
  `adf_test` and `engle_granger`, placed before the OU calibration to reflect the
  pairs-trade pipeline; TOC regenerated to 97 entries.

## [1.568.0] - 2026-09-12

### Added
- `cointegration.py`: `adf_test` runs the augmented Dickey-Fuller unit-root test
  (constant, no trend, optional lagged differences) and `engle_granger` runs the
  two-series cointegration test -- regress ``y`` on ``x`` and ADF-test the residual
  spread, returning the hedge ratio and a cointegration verdict. Cross-checked: a
  stationary AR(1) rejects the unit root (statistic ~ -18), a random walk does not
  (verified on the majority of independent draws), a synthetic cointegrated pair is
  detected with the correct hedge ratio, and two independent walks are not
  cointegrated. Feeds directly into `fit_ornstein_uhlenbeck` for the spread.

## [1.567.0] - 2026-09-12

### Documentation
- README: added an "Ornstein-Uhlenbeck calibration" section documenting
  `fit_ornstein_uhlenbeck`; TOC regenerated to 96 entries.

## [1.566.0] - 2026-09-12

### Added
- `ou_fit.py`: `fit_ornstein_uhlenbeck` calibrates the mean-reverting OU process
  ``dX = kappa(theta - X)dt + sigma dW`` from a sampled path in closed form -- it
  fits the exact discrete AR(1) by least squares and inverts the OU relations to
  recover ``(kappa, theta, sigma)`` plus the half-life ``ln(2)/kappa``.
  Cross-checked: it recovers the parameters of a simulated OU path to within ~10%,
  the half-life relation holds exactly, faster true reversion yields a shorter
  fitted half-life, a random walk fits a near-zero kappa (very long half-life), and
  an anti-persistent (oscillating) series is rejected.

## [1.565.0] - 2026-09-12

### Documentation
- README: added a "Variance-ratio test" section documenting `variance_ratio` and
  `variance_ratio_zstat`; TOC regenerated to 95 entries.

## [1.564.0] - 2026-09-12

### Added
- `variance_ratio.py`: `variance_ratio` computes the Lo-MacKinlay (1988) variance
  ratio ``VR(q)`` from overlapping ``q``-period returns (with the unbiased scaling
  factors), and `variance_ratio_zstat` gives the heteroskedasticity-robust ``z``
  statistic for the random-walk null. Cross-checked: white noise gives ``VR ~ 1``
  with ``|z|`` small, a mean-reverting AR(1) gives ``VR < 1`` with ``z`` rejecting
  below -1.96, and a trending AR(1) gives ``VR > 1`` with ``z`` rejecting above
  1.96.

## [1.563.0] - 2026-09-12

### Documentation
- README: added a "Hurst exponent (long memory)" section documenting
  `hurst_exponent` and `rescaled_range`; TOC regenerated to 94 entries.

## [1.562.0] - 2026-09-12

### Added
- `hurst.py`: `hurst_exponent` estimates the Hurst exponent by rescaled-range
  (R/S) analysis -- the slope of ``log(R/S)`` against ``log(window)`` -- plus the
  `rescaled_range` of a single window. Cross-checked against the canonical regimes:
  white noise gives ``H ~ 0.5``, a cumulated random walk gives ``H ~ 1``, a
  mean-reverting AR(1) with negative coefficient gives ``H < 0.5``, and a
  persistent AR(1) gives ``H > 0.5``.

## [1.561.0] - 2026-09-12

### Documentation
- README: extended the Realized volatility section with
  `two_scale_realized_variance` and `realized_variance_naive` (microstructure-noise
  bias correction), with a snippet verified against the package.

## [1.560.0] - 2026-09-12

### Added
- `two_scale_rv.py`: `two_scale_realized_variance` implements the Zhang-Mykland-
  Ait-Sahalia (2005) two-scale estimator that removes the microstructure-noise
  bias of realized variance by combining a subsampled slow scale with the
  bias-estimable fast scale; also `realized_variance_naive` and
  `noise_variance_estimate`. Cross-checked: on a noise-free path it reproduces the
  integrated variance, under i.i.d. noise it is far less biased than the naive RV
  (which is inflated by ``2 n Var(eps)``, confirmed to 20%), and the noise-variance
  estimate recovers ``Var(eps)`` to ~6% in the noise-dominated regime.

## [1.559.0] - 2026-09-12

### Documentation
- README: extended the Realized volatility section with the intraday jump
  decomposition (`realized_variance_from_returns`, `bipower_variation`,
  `jump_variation`), with a snippet verified against the package.

## [1.558.0] - 2026-09-12

### Added
- `realized.py`: jump-robust realized volatility from intraday returns --
  `realized_variance_from_returns` (total quadratic variation),
  `bipower_variation` (Barndorff-Nielsen-Shephard integrated variance, robust to
  jumps), `jump_variation` (``max(RV - BV, 0)``), and
  `realized_volatility_signature`. Cross-checked: on a continuous path RV and BV
  both recover the integrated variance and the jump component is negligible; a
  single injected jump inflates RV by roughly the jump squared while BV stays
  within 15%, and the jump variation recovers the jump magnitude.

## [1.557.0] - 2026-09-12

### Documentation
- README: added a "Robust scale and location" section documenting
  `median_absolute_deviation`, `interquartile_range`, `winsorize`, and
  `trimmed_mean`, with a snippet verified against the package; TOC regenerated to
  93 entries.

## [1.556.0] - 2026-09-12

### Added
- `robust_stats.py`: robust scale and location estimators --
  `median_absolute_deviation` (50%-breakdown scale, 1.4826-scaled to match the
  standard deviation under normality), `interquartile_range` (with a
  normal-consistent scaling), `winsorize` (clip tails to percentiles), and
  `trimmed_mean`. Cross-checked: MAD and the scaled IQR both recover the true
  Gaussian sigma, they are exactly zero on constant data, MAD stays bounded under
  outliers where the standard deviation explodes, winsorize clips extremes while
  preserving length, and the trimmed mean shrugs off a gross outlier that moves
  the ordinary mean by 100x.

## [1.555.0] - 2026-09-12

### Documentation
- README: added a "Theil-Sen robust regression" section documenting `theil_sen`,
  with a snippet verified against the package; TOC regenerated to 92 entries.

## [1.554.0] - 2026-09-12

### Added
- `theil_sen.py`: `theil_sen` computes the Theil-Sen robust regression slope (the
  median of all pairwise slopes) and intercept, with a ~29% breakdown point.
  Cross-checked: it is exact on noiseless linear data, tracks OLS and the true
  slope on clean Gaussian data, and holds the true slope under 20% gross-outlier
  contamination where OLS is dragged away; pairs sharing an ``x`` value are
  skipped.

## [1.553.0] - 2026-09-12

### Documentation
- README: extended the Kalman-filter section with `kalman_regression_beta`
  (time-varying regression slope / dynamic hedge ratio), with a snippet verified
  against the package.

## [1.552.0] - 2026-09-12

### Added
- `kalman_beta.py`: `kalman_regression_beta` filters a time-varying regression
  slope (a dynamic hedge ratio / factor loading) as a random walk observed through
  ``y_t = beta_t x_t + v_t`` -- a scalar Kalman filter with a time-varying
  observation loading. Cross-checked: with ``Q = 0`` and a diffuse prior it equals
  the OLS slope ``sum(x y) / sum(x^2)``, the posterior variance shrinks with data,
  a positive ``Q`` tracks a regime switch in the true slope (1.0 -> 3.0), and the
  static filter settles at a blend rather than following the late regime.

## [1.551.0] - 2026-09-12

### Documentation
- README: added a "Newey-West HAC variance" section documenting
  `newey_west_variance`, `newey_west_mean_se`, and `autocorrelation`, with a
  snippet verified against the package; TOC regenerated to 91 entries.

## [1.550.0] - 2026-09-12

### Added
- `hac.py`: `newey_west_variance` estimates the Newey-West (1987) HAC long-run
  variance of a series via Bartlett-weighted autocovariances (the weighting that
  guarantees a non-negative estimate), plus `newey_west_mean_se` for the HAC
  standard error of the sample mean and `autocovariance` / `autocorrelation`
  helpers. Cross-checked: lag 0 reduces to the sample variance, the estimate is
  always non-negative, white noise stays near the sample variance, and on an
  AR(1) it climbs toward the analytic long-run variance ``sigma^2 / (1 - phi)^2``
  as the lag grows while the lag-1 autocorrelation recovers ``phi``.

## [1.549.0] - 2026-09-12

### Documentation
- README: added a "Kalman filter (local level)" section documenting
  `kalman_local_level` and `kalman_steady_state_gain`, with a snippet verified
  against the package; TOC regenerated to 90 entries.

## [1.548.0] - 2026-09-12

### Added
- `kalman.py`: `kalman_local_level` runs the exact scalar Kalman recursion for the
  local-level (random-walk-plus-noise) state-space model, returning the filtered
  level, its posterior variance, and the Kalman gain at each step;
  `kalman_steady_state_gain` gives the closed-form Riccati fixed-point gain.
  Cross-checked: the running gain converges to the steady-state formula,
  ``R -> 0`` drives the gain to 1 and reproduces the data, ``Q -> 0`` with a
  diffuse prior recovers the running mean with gain ``1/(t+1)`` (recursive least
  squares), the posterior variance stays positive, and the steady-state gain is
  monotone in the signal-to-noise ratio.

## [1.547.0] - 2026-09-12

### Documentation
- README: added a "Hodrick-Prescott filter" section documenting `hp_filter`, with
  a snippet verified against the package; TOC regenerated to 89 entries.

## [1.546.0] - 2026-09-12

### Added
- `hp_filter.py`: `hp_filter` performs a Hodrick-Prescott trend/cycle
  decomposition, solving ``(I + lambda D'D) tau = y`` for the smooth trend by an
  O(n) banded LDL^T factorization of the pentadiagonal system (no dense inverse).
  Cross-checked: the banded solve matches a dense inverse to ~1e-13, the trend and
  cycle reconstruct the series exactly, ``lambda -> 0`` returns the data itself,
  ``lambda -> inf`` returns the least-squares linear trend, and a larger
  ``lambda`` yields a strictly smoother (lower-curvature) trend.

## [1.545.0] - 2026-09-12

### Documentation
- README: extended the Kelly-sizing block with `kelly_fractions_multivariate`
  (growth-optimal leverage across correlated assets), with a snippet verified
  against the package.

## [1.544.0] - 2026-09-12

### Added
- `sizing.py`: `kelly_fractions_multivariate` computes the growth-optimal Kelly
  leverage vector across correlated assets, ``f* = Sigma^{-1} mu``, and
  `kelly_growth_rate_multivariate` evaluates the expected log-growth
  ``f . mu - 0.5 f . Sigma f`` at any leverage vector. Cross-checked: reduces to
  the scalar ``mu / sigma^2`` for one asset, to the per-asset fractions when the
  covariance is diagonal, the growth rate is maximized at ``f*`` (perturbations
  lower it), fractional Kelly scales the vector linearly, and positive correlation
  cuts the total leverage below the uncorrelated case. Pairs naturally with the
  shrunk covariance from `ledoit_wolf_shrinkage`.

## [1.543.0] - 2026-09-12

### Documentation
- README: added a "Covariance shrinkage (Ledoit-Wolf)" section documenting
  `ledoit_wolf_shrinkage`, with a deterministic snippet verified against the
  package; TOC regenerated to 88 entries.

## [1.542.0] - 2026-09-12

### Added
- `shrinkage.py`: `ledoit_wolf_shrinkage` estimates a covariance matrix by the
  Ledoit-Wolf (2004) constant-correlation shrinkage -- a data-driven convex blend
  ``delta * F + (1 - delta) * S`` of the sample covariance and a
  constant-correlation target, with the closed-form optimal intensity. Also exports
  `sample_covariance` and `constant_correlation_target`. Cross-checked: ``delta``
  stays in ``[0, 1]`` and decreases as the sample grows (1.0 at n=50 down to 0.002
  at n=5000 on heterogeneous data), the result is the exact convex combination,
  the diagonal is preserved, off-diagonals are pulled toward the target, and the
  estimate is symmetric and positive definite.

## [1.541.0] - 2026-09-12

### Documentation
- README: extended the "Double-barrier knock-out" section with the knock-in
  counterpart (`double_knockin_call`) and the in-out parity, with a snippet
  verified against the package.

## [1.540.0] - 2026-09-12

### Added
- `double_barrier.py`: `double_knockin_call` prices a double-barrier knock-in call
  -- alive only if the spot touches either barrier before expiry -- via the in-out
  parity ``knock-in + knock-out = vanilla``. Cross-checked: parity holds to 1e-9,
  the value is bounded by the vanilla call, it approaches zero as the barriers move
  far away (a breach becomes rare), rises toward the vanilla as the corridor
  tightens, and matches a fine-grid Monte Carlo in the continuous-monitoring limit
  (discrete MC undercounts breaches and converges up onto it).

## [1.539.0] - 2026-09-12

### Documentation
- README: added a "Range-accrual note" section documenting `range_accrual_note`,
  with snippets verified against the package; TOC regenerated to 87 entries.

## [1.538.0] - 2026-09-12

### Added
- `range_accrual.py`: `range_accrual_note` prices a range-accrual note's coupon
  leg in closed form -- the coupon accrues proportionally to the fraction of
  observation dates on which the index sits inside a band ``[L, U]``, valued as a
  discounted sum of GBM range probabilities. Cross-checked: it collapses to the
  discounted full coupon as the band widens to the whole axis, a wider band raises
  the value, higher volatility lowers it, the notional scales it linearly, and it
  matches a Monte Carlo count to 1e-5.

## [1.537.0] - 2026-09-12

### Documentation
- README: added a "Double-barrier knock-out" section documenting
  `double_knockout_call`, with a snippet verified against the package; TOC
  regenerated to 86 entries.

## [1.536.0] - 2026-09-12

### Added
- `double_barrier.py`: `double_knockout_call` prices a double-barrier knock-out
  call in closed form via the Kunitomo-Ikeda (1992) image series -- the vanilla
  call payoff survives only if the spot stays strictly inside a corridor
  ``(L, U)`` for the whole life. Cross-checked: it never exceeds the vanilla call,
  approaches the vanilla as the barriers move far away, the series converges in a
  handful of terms, a tighter corridor lowers the value, and it matches a
  fine-grid Monte Carlo in the continuous-monitoring limit (discrete MC converges
  down onto the closed form as the step count rises).

## [1.535.0] - 2026-09-12

### Documentation
- README: added an "Installment options" section documenting `installment_call`,
  with snippets verified against the package; TOC regenerated to 85 entries.

## [1.534.0] - 2026-09-12

### Added
- `installment.py`: `installment_call` prices a European installment call on a CRR
  tree. The option is paid for in a stream of premiums; at each installment date
  the holder may lapse (stop paying, forfeit for zero) so the option is kept alive
  only while its continuation value exceeds the next installment -- a compound
  option priced by backward induction. Returns the fair upfront value.
  Cross-checked: zero installment (or empty schedule) reproduces the plain CRR
  European call to 1e-9, a larger installment lowers the upfront value, more
  payment dates lower it further, a prohibitive installment drives it to zero, and
  the value is always non-negative.

## [1.533.0] - 2026-09-12

### Documentation
- README: added a "Shout and ladder options" section documenting `shout_call` and
  `ladder_call`, with snippets verified against the package; TOC regenerated.

## [1.532.0] - 2026-09-12

### Added
- `shout.py`: `ladder_call` prices a ladder call on a CRR tree -- the payoff is
  floored at the highest preset rung the underlying touches before expiry,
  `max(S_T - K, max_touched L_i - K, 0)`, carried as a rung-state variable through
  backward induction. Cross-checked: with no rungs it reproduces the plain CRR
  European call to 1e-9, rungs at or below the strike are ignored, more/higher
  rungs raise the value, and higher volatility raises the value.

## [1.531.0] - 2026-09-12

### Added
- `shout.py`: `shout_call` prices a shout option on a CRR tree -- the holder may
  shout once to lock in the current intrinsic as a floor while keeping the upside.
  Cross-checked: the shout is at least the vanilla call, the in-the-money lock-in
  premium is positive, value rises with volatility, and the tree converges across
  step counts.

## [1.530.0] - 2026-09-12

### Documentation
- README: added a "Leveraged ETFs" section covering the `leveraged_etf` module
  (TOC auto-updated to 83 entries), with runnable examples verified against the
  installed package.

## [1.529.0] - 2026-09-12

### Added
- `leveraged_etf.py`: daily-rebalanced leveraged/inverse ETF path,
  `volatility_drag` (`0.5 L (L-1) sigma^2`), `expected_leveraged_return`, and
  `flat_market_decay`. Cross-checked: 1x matches the underlying, the drag is zero
  at leverage 0 or 1 and larger for higher/inverse leverage, the expected return
  is below naive leverage, and a 3x ETF loses value over a net-flat volatile path.

## [1.528.0] - 2026-09-12

### Documentation
- README: added an "Equity valuation" section covering the `valuation` module
  (TOC auto-updated to 82 entries), with runnable examples verified against the
  installed package.

## [1.527.0] - 2026-09-12

### Added
- `valuation.py`: equity valuation. `capm_cost_of_equity`, `wacc`,
  `gordon_growth_value` (constant-growth DDM), `terminal_value`, and
  `two_stage_dcf`. Cross-checked: WACC lies between the after-tax debt and equity
  costs (and equals `ke` at all-equity), the Gordon and terminal formulas require
  `r > g`, and the two-stage DCF rises with growth and falls with the discount
  rate.

## [1.526.0] - 2026-09-12

### Documentation
- README: added a "Capital budgeting" section covering the `capital_budgeting`
  module (TOC auto-updated to 81 entries), with runnable examples verified against
  the installed package.

## [1.525.0] - 2026-09-12

### Added
- `capital_budgeting.py`: `npv`, `irr` (bisection), `profitability_index`,
  `payback_period` (interpolated), and `mirr`. Cross-checked: NPV is zero at the
  IRR, the IRR recovers a known discount rate, the profitability index is one at
  the IRR and above one below it, and NPV is monotone decreasing in the rate.

## [1.524.0] - 2026-09-12

### Documentation
- README: added a "Money-market yields" section covering the `money_market`
  module (TOC auto-updated to 80 entries), with runnable examples verified against
  the installed package.

## [1.523.0] - 2026-09-12

### Added
- `money_market.py`: money-market yield conventions. `price_from_discount`,
  `bank_discount_yield`, `money_market_yield` (CD-equivalent, actual/360),
  `bond_equivalent_yield` (actual/365), `discount_to_bond_equivalent`, and
  `holding_period_return`. Cross-checked: the discount/price round-trip and the
  yield ordering (bond-equivalent > money-market > bank discount).

## [1.522.0] - 2026-09-12

### Documentation
- README: expanded the "FX forwards" section with cross rates and triangular
  arbitrage, verified against the installed package.

## [1.521.0] - 2026-09-12

### Added
- `cross_rate`, `triangular_arbitrage`, `is_arbitrage_free` (in `fxforward.py`):
  FX cross rates from a common currency, the triangular round-trip product around
  a currency loop, and an arbitrage-free check. Cross-checked: the cross rate,
  the round-trip product equal to one in an arbitrage-free market, a mispriced
  loop signalling profit (>1) or the reverse (<1).

## [1.520.0] - 2026-09-12

### Documentation
- README: added the HAR-RV forecast to the "Realized volatility" section, verified
  against the installed package.

## [1.519.0] - 2026-09-12

### Added
- `fit_har_rv` and `har_rv_forecast` (in `volatility.py`): the HAR-RV (Corsi 2009)
  long-memory realized-variance model -- regress next-day RV on the daily, weekly,
  and monthly average RV, and forecast one step ahead. Cross-checked: the fitted
  persistence (sum of the collinear lag slopes) recovers the true 0.9, the
  in-sample fit correlation exceeds 0.99, and the forecast matches the manual
  linear combination.

## [1.518.0] - 2026-09-12

### Documentation
- README: expanded the "GARCH volatility" section with the GJR-GARCH and EGARCH
  leverage variants, verified against the installed package.

## [1.517.0] - 2026-09-12

### Added
- `EGarchParams`, `egarch_variance`, `egarch_forecast` (in `volatility.py`): the
  Nelson EGARCH(1,1) model on the *log* conditional variance, so the variance is
  positive for any parameters (no constraints) and `gamma < 0` gives the leverage
  effect. Cross-checked: the variance is always positive (even for absurd
  parameters), negative shocks raise it more, it is symmetric at `gamma = 0`, and
  the forecast reverts to the unconditional log-variance level.

## [1.516.0] - 2026-09-12

### Added
- `GJRGarchParams`, `gjr_garch_variance`, `gjr_garch_forecast` (in
  `volatility.py`): the GJR-GARCH(1,1,1) leverage model, where negative shocks
  raise conditional variance more than positive ones (the leverage effect), with
  persistence `alpha + beta + 0.5 gamma`. Cross-checked: down shocks raise
  variance more than up shocks, the model reduces to symmetric GARCH at
  `gamma = 0`, and the forecast reverts to the long-run vol.

## [1.515.0] - 2026-09-12

### Documentation
- README: expanded the "Numerical utilities" section with the `quadrature`
  integration routines, verified against the installed package.

## [1.514.0] - 2026-09-12

### Added
- `quadrature.py`: numerical integration -- `trapezoid`, `simpson` (composite,
  exact for cubics), `gauss_legendre` (orders 2-5, exact to degree 2n-1), and
  `adaptive_simpson` (error-controlled). Cross-checked: Simpson exact on a cubic,
  Gauss-Legendre exact to degree 2n-1, and the adaptive rule matching sin,
  Lorentzian, and Gaussian integrals to 1e-8.

## [1.513.0] - 2026-09-12

### Documentation
- README: added a "Rebalancing" section covering the `rebalance` module (TOC
  auto-updated to 79 entries), with runnable examples verified against the
  installed package.

## [1.512.0] - 2026-09-12

### Added
- `rebalance.py`: portfolio rebalancing analytics. `drift_weights` (buy-and-hold
  weight drift), `turnover` (one-way `0.5 sum |delta w|`), `transaction_cost`
  (round-trip cost drag in bps), and `no_trade_band_rebalance` (only trade
  positions outside a tolerance). Cross-checked: drift renormalizes with the
  winner gaining weight, turnover is zero on target and one for a full swap, the
  cost formula, and the band suppressing small drifts while trading large ones.

## [1.511.0] - 2026-09-12

### Documentation
- README: expanded the "Equity compensation and convertibles" section with the
  `convertible_bond_lattice` equity-tree pricer, verified against the installed
  package.

## [1.510.0] - 2026-09-12

### Added
- `convertible_lattice.py`: `convertible_bond_lattice` prices a convertible bond
  on a Cox-Ross-Rubinstein equity tree by backward induction, with American-style
  conversion, an issuer call, and a holder put, discounting at a credit spread.
  Cross-checked: the value stays at or above both the conversion parity and the
  bond floor, approaches the conversion value deep in the money, is lower when
  callable, and rises with volatility.

## [1.509.0] - 2026-09-12

### Documentation
- README: added an auto-generated table of contents (`docs/gen_toc.py`) between
  `<!-- TOC -->` markers, covering all 78 sections with GitHub-style anchors, and
  a `test_readme_toc_is_up_to_date` check so the TOC cannot drift from the
  headings.

## [1.508.0] - 2026-09-12

### Documentation
- README: added a "Pairs trading" section covering the `pairs` module, with
  runnable examples verified against the installed package.

## [1.507.0] - 2026-09-12

### Added
- `pairs.py`: pairs-trading analytics. `pairs_hedge_ratio` (OLS slope of one leg
  on the other), `spread_series`, `ou_half_life` (mean-reversion half-life from an
  AR(1) fit to the spread), and `spread_zscore`. Cross-checked: the hedge ratio
  recovers a known beta, the half-life falls in the expected range for a
  cointegrated pair, and a constant (non-reverting) spread raises. Named
  `pairs_hedge_ratio` to avoid colliding with the LDI `hedge_ratio`.

### Documentation
- README: expanded the "Trend and momentum signals" section with Bollinger bands
  and the Donchian channel.

## [1.506.0] - 2026-09-12

### Added
- `bollinger_bands`, `average_true_range`, `donchian_channel` (in `signals.py`):
  Bollinger bands (SMA +/- k std), Wilder's average true range, and the Donchian
  breakout channel. Cross-checked: the bands are ordered with the SMA middle and
  bracket ~99% of prices, the ATR is positive, and the Donchian upper is the
  rolling high and never below the lower.

## [1.505.0] - 2026-09-12

### Documentation
- README: added a "Trend and momentum signals" section covering the `signals`
  module, with runnable examples verified against the installed package.

## [1.504.0] - 2026-09-12

### Added
- `signals.py`: trend and momentum signals. `sma`, `ema`, `macd`, `rsi` (Wilder's
  smoothing), `rolling_zscore`, and `time_series_momentum`. Cross-checked: the SMA
  values, the EMA reacting faster than the SMA after a step, MACD positive on an
  uptrend and equal to the fast-minus-slow EMA, RSI in `[0, 100]` (high in an
  uptrend, low in a downtrend), the z-score flagging a spike, and the TSM sign.

## [1.503.0] - 2026-09-12

### Documentation
- README: added a "Volatility targeting" section covering the `vol_target`
  module, with runnable examples verified against the installed package.

## [1.502.0] - 2026-09-12

### Added
- `vol_target.py`: volatility-targeting overlay. `target_leverage`
  (`target_vol/realized_vol`, capped), `vol_targeted_returns` (rolling overlay on a
  return series), and `realized_annualized_vol`. Cross-checked: the leverage
  formula and cap, higher realized vol giving lower leverage, and the overlaid
  series realizing near the target vol where the raw series is much higher.

## [1.501.0] - 2026-09-11

### Documentation
- README: added a "Portfolio insurance (CPPI)" section covering the `cppi` module,
  with runnable examples verified against the installed package.

## [1.500.0] - 2026-09-11

### Added
- `cppi.py`: Constant Proportion Portfolio Insurance. `discounted_floor`,
  `cushion`, `risky_exposure` (multiplier times the cushion, capped at wealth and
  floored at zero), and `cppi_path` (wealth-path simulation). Cross-checked: the
  exposure is the multiplier times the cushion and respects the no-leverage cap,
  higher multipliers take more risk, a positive-return path grows, and a crash
  path stays at the guaranteed floor.

## [1.499.0] - 2026-09-11

### Documentation
- README: expanded the "Portfolio optimization" section with inverse-volatility
  and Hierarchical Risk Parity allocation, with runnable examples verified against
  the installed package.

## [1.498.0] - 2026-09-11

### Added
- `hrp.py`: `inverse_volatility_weights` and `hierarchical_risk_parity` (López de
  Prado HRP -- cluster by correlation distance, quasi-diagonalize, allocate by
  recursive inverse-variance bisection, no matrix inversion). Cross-checked:
  inverse-vol weights sum to one and rank by volatility, HRP weights are positive
  and sum to one, an independent asset gets more than a correlated pair, and
  equal-uncorrelated assets get equal weights.

## [1.497.0] - 2026-09-11

### Documentation
- README: added a "Carry and roll-down" section covering the `carry_rolldown`
  module, with runnable examples verified against the installed package.

## [1.496.0] - 2026-09-11

### Added
- `carry_rolldown.py`: bond carry and roll-down return decomposition off a zero
  curve. `carry_return` (coupon net of financing), `rolldown_return` (the price
  gain purely from the yield rolling down the curve), and `total_carry_rolldown`.
  Cross-checked: the roll-down is exactly zero on a flat curve, positive on an
  upward curve (more when steeper), negative when inverted, and the total is
  carry plus roll-down. The roll-down originally conflated time-value with the
  yield change; isolating the yield-change effect fixed the flat-curve case.

## [1.495.0] - 2026-09-11

### Documentation
- README: added a "Callable bonds and OAS" section covering the `callable_bond`
  module (callable/puttable pricing on a short-rate tree, option-adjusted spread),
  with runnable examples verified against the installed package.

## [1.494.0] - 2026-09-11

### Added
- `callable_bond_price_with_spread` and `option_adjusted_spread` (in
  `callable_bond.py`): reprice a callable/puttable bond with a constant spread
  added to every short-rate node, and solve the option-adjusted spread (the
  constant spread repricing the bond to a market price) by bisection.
  Cross-checked: zero spread matches the base price, a positive spread lowers the
  price, and the OAS recovers a known spread and is zero at the model price
  (using a non-binding call so the price-spread curve is strictly monotone).

## [1.493.0] - 2026-09-11

### Added
- `callable_bond.py`: callable/puttable bond pricing on a Black-Derman-Toy-style
  binomial short-rate tree by backward induction. `callable_bond_price` (call caps
  / put floors the node value), `straight_bond_tree_price`, and `call_option_value`
  (straight minus callable). Cross-checked: callable below and puttable above the
  straight bond, no-optionality equals the straight bond, the call value is
  non-negative, and higher volatility lowers the callable price.

## [1.492.0] - 2026-09-11

### Added
- `gev_cdf`, `gev_return_level`, `gev_fit_block_maxima` (in `evt.py`): the
  Generalized Extreme Value (block-maxima) distribution, the T-block return level,
  and a Gumbel-limit method-of-moments fit. Cross-checked: the CDF is monotone and
  bounded with the correct Gumbel limit, the return level equals the GEV quantile
  and rises with the return period, and the Gumbel fit recovers the location and
  scale.

## [1.491.0] - 2026-09-11

### Documentation
- README: added an "Extreme value theory" section covering the `evt` module (Hill
  tail index, peaks-over-threshold GPD VaR/ES), with runnable examples verified
  against the installed package.

## [1.490.0] - 2026-09-11

### Added
- `evt.py`: extreme value theory tail risk. `hill_estimator` (tail index from the
  top order statistics), `gpd_fit_pot` (method-of-moments Generalized Pareto fit
  to peaks over a threshold), `gpd_var`, and `gpd_expected_shortfall`.
  Cross-checked: the Hill estimator recovers a Pareto tail index (1/3 for
  alpha = 3), the GPD shape is positive for a heavy tail, ES exceeds VaR, and VaR
  rises with confidence.

## [1.489.0] - 2026-09-11

### Documentation
- README: added a "Student-t fat tails" section covering the `student_t` module
  (distribution, fat-tailed VaR/ES, degrees-of-freedom fitting), with runnable
  examples verified against the installed package.

## [1.488.0] - 2026-09-11

### Added
- `fit_df_from_kurtosis` and `fit_student_t` (in `student_t.py`): recover the
  Student-t degrees of freedom from excess kurtosis (`df = 4 + 6/excess`) and fit
  a location-scale t to a return sample by matching the mean, variance, and
  kurtosis. Cross-checked: the df/kurtosis formula, higher kurtosis giving fewer
  degrees of freedom, and the fit recovering a finite df near the generating value.

## [1.487.0] - 2026-09-11

### Added
- `student_t.py`: Student's t distribution and fat-tailed parametric risk.
  `t_pdf`, `t_cdf` (regularized incomplete beta), `t_ppf`, `student_t_var`, and
  `student_t_expected_shortfall` (closed form). Cross-checked: the CDF is
  symmetric and converges to the normal as df grows, the df=5 critical value is
  2.015, the density integrates to one, and the t VaR is fatter than the normal
  VaR (converging to it), with ES at least the VaR.

## [1.486.0] - 2026-09-11

### Documentation
- README: expanded the "Performance metrics" section with the drawdown-pain
  (Ulcer / Martin) and statistical-significance (probabilistic and deflated
  Sharpe) measures, with runnable examples verified against the installed package.

## [1.485.0] - 2026-09-11

### Added
- `deflated_sharpe_ratio` (in `perfmetrics.py`): the Bailey-López de Prado
  deflated Sharpe ratio -- the probabilistic Sharpe ratio against a benchmark
  equal to the expected maximum of `n_trials` Sharpe estimates, correcting for
  the selection bias of testing many strategy variants. Cross-checked: reduces to
  the plain PSR at one trial, deflates for more trials, is monotone decreasing in
  the trial count, and stays in `[0, 1]`.

## [1.484.0] - 2026-09-11

### Added
- `probabilistic_sharpe_ratio` and `minimum_track_record_length` (in
  `perfmetrics.py`): the Bailey-López de Prado probability that the true Sharpe
  ratio beats a benchmark (correcting the estimator's standard error for skewness,
  kurtosis, and sample length) and the track-record length needed to reach a
  confidence. Cross-checked: PSR exceeds 0.5 for a positive edge and falls against
  a higher benchmark, and the minimum track record length is positive and longer
  for a smaller edge.

## [1.483.0] - 2026-09-11

### Added
- `ulcer_index`, `pain_index`, `ulcer_performance_index`, `pain_ratio` (in
  `perfmetrics.py`): drawdown-based risk and return-per-pain ratios -- the RMS
  (Ulcer) and mean (pain) depth of the underwater curve, and the Martin/pain
  ratios of annualized excess return over each. Cross-checked: zero without
  drawdown, Ulcer at least the pain index, both bounded by the max drawdown, and
  the ratios raising when there is no drawdown.

## [1.482.0] - 2026-09-11

### Documentation
- README: added a "Sample risk measures" section covering the `riskmeasures`
  module (VaR, ES, spectral/entropic risk, Euler component allocation), with
  runnable examples verified against the installed package.

## [1.481.0] - 2026-09-11

### Added
- `is_subadditive` and `component_expected_shortfall` (in `riskmeasures.py`): a
  subadditivity diagnostic for expected shortfall (`rho(A+B) <= rho(A) + rho(B)`,
  which the coherent ES always satisfies) and the Euler component-ES allocation
  (each sub-portfolio's mean loss over the total portfolio's tail scenarios).
  Cross-checked: subadditivity holds including for perfectly correlated series,
  and the component contributions sum exactly to the portfolio ES for two and
  three components.

## [1.480.0] - 2026-09-11

### Added
- `riskmeasures.py`: sample risk measures from a P&L series. `value_at_risk`
  (loss quantile), `sample_expected_shortfall` (coherent CVaR tail mean),
  `spectral_risk_exponential` (exponential risk-aversion spectrum), and
  `entropic_risk` (exponential certainty equivalent). Cross-checked: ES is at
  least VaR, VaR rises with confidence, spectral and entropic rise with risk
  aversion, and entropic approaches the mean loss as risk aversion vanishes.
  Exported as `sample_expected_shortfall` to avoid colliding with the smile-based
  `expected_shortfall`.

## [1.479.0] - 2026-09-11

### Documentation
- README: added a "Nelson-Siegel / Svensson curves" section covering the
  `nelson_siegel` module (parametric curves and least-squares calibration), with
  runnable examples verified against the installed package.

## [1.478.0] - 2026-09-11

### Added
- `fit_nelson_siegel` (in `nelson_siegel.py`): least-squares fit of Nelson-Siegel
  parameters to observed zero rates. For each ``tau`` on a grid the three betas
  are solved by OLS (they enter linearly via the factor loadings) and ``tau`` is
  chosen by minimizing the residual sum of squares. Cross-checked: recovers the
  true parameters exactly on noiseless data, reprices the curve, and leaves a
  small residual under noise.

## [1.477.0] - 2026-09-11

### Added
- `nelson_siegel.py`: Nelson-Siegel and Svensson parametric yield curves.
  `nelson_siegel_zero`, `svensson_zero` (adds a second curvature hump),
  `nelson_siegel_discount`, and `nelson_siegel_forward`. Cross-checked: the short
  rate (`t -> 0`) is `beta0 + beta1`, the long level (`t -> inf`) is `beta0`, the
  discount factor is `exp(-z t)` and decreasing, and Svensson reduces to
  Nelson-Siegel when `beta3 = 0`.

## [1.476.0] - 2026-09-11

### Documentation
- README: added a "Numerical utilities" section covering the `interpolation` and
  `rootfind` modules and `SplineZeroCurve`, with runnable examples verified
  against the installed package.

## [1.475.0] - 2026-09-11

### Added
- `rootfind.py`: general-purpose scalar root finders. `bisection` (guaranteed
  convergence on a sign-changing bracket), `brent` (bisection + secant + inverse
  quadratic interpolation), and `newton` (derivative-based with an optional
  bisection safeguard). Cross-checked against known roots (sqrt(2), the cos
  fixed point, a cubic), bracket safeguards, and the no-bracket / zero-derivative
  error cases.

## [1.474.0] - 2026-09-11

### Added
- `SplineZeroCurve` (in `discount_curve.py`): a discount curve that interpolates
  the pillar zero rates with a natural cubic spline (smooth instantaneous
  forwards) rather than log-linear discount factors, with the same
  `df`/`zero_rate`/`forward_rate` interface as `DiscountCurve`. Cross-checked: it
  reprices the pillars exactly, `df = exp(-z T)`, discount factors stay monotone
  for an upward curve, and it flat-extrapolates outside the pillar range.

## [1.473.0] - 2026-09-11

### Added
- `interpolation.py`: cubic curve interpolation. `natural_cubic_spline` (C2
  natural spline via the tridiagonal solve) and `monotone_cubic` (Fritsch-Carlson
  monotone cubic Hermite, no overshoot -- for discount-factor/survival curves).
  Cross-checked: both pass through the knots and are exact on linear data, the
  monotone interpolant preserves monotonicity without overshooting where the
  natural spline does.

## [1.472.0] - 2026-09-11

### Documentation
- README: added a "Matrix utilities" section covering the `linalg` module
  (Cholesky, correlated draws, Higham nearest-correlation, basket Monte Carlo),
  with runnable examples verified against the installed package.

## [1.471.0] - 2026-09-11

### Added
- `basket_option_mc` (in `linalg.py`): a general n-asset basket option Monte Carlo
  under correlated geometric Brownian motion (correlated draws via Cholesky).
  Independently validates the two-asset analytic `basket_option` (agrees within 3%
  at 200k paths), handles any number of assets, is deterministic per seed, and
  prices richer for higher correlation.

## [1.470.0] - 2026-09-11

### Added
- `linalg.py`: matrix utilities for risk work. `cholesky` (`L L^T` factor),
  `is_positive_definite`, `correlated_normals` (turn IID normals into correlated
  draws via Cholesky), and `nearest_correlation` (Higham alternating-projection
  repair of an indefinite correlation matrix to the nearest PSD unit-diagonal
  one). Cross-checked: Cholesky reconstructs the matrix, the correlated transform
  matches `L z`, and Higham repair yields a symmetric unit-diagonal PSD matrix
  while leaving a valid correlation unchanged.

## [1.469.0] - 2026-09-11

### Documentation
- README: added a "Principal component analysis" section covering the `pca`
  module (Jacobi eigendecomposition, variance explained, reconstruction, component
  scenarios), with runnable examples verified against the installed package.

## [1.468.0] - 2026-09-11

### Added
- `reconstruct_covariance` and `pca_scenario` (in `pca.py`): rebuild a covariance
  matrix from its top ``k`` principal components (the best rank-``k`` approximation;
  exact with all components) and generate an ``n``-sigma stress scenario along a
  chosen component (level/slope/curvature shifts for a yield curve). Cross-checked:
  full reconstruction recovers the covariance, the rank-1 trace equals the top
  variance, and the scenario magnitude equals the component standard deviation and
  scales linearly with sigma.

## [1.467.0] - 2026-09-11

### Added
- `pca.py`: principal component analysis via the Jacobi eigenvalue algorithm.
  `jacobi_eigen` (symmetric eigendecomposition), `pca` (sorted variances,
  orthonormal loadings, variance explained -- the yield-curve level/slope/curvature
  decomposition), and `project` (component scores). Cross-checked: correct
  eigenvalues, orthonormal eigenvectors, the eigen equation `Av = lambda v`,
  variances summing to the trace, and norm-preserving full projection.

## [1.466.0] - 2026-09-11

### Documentation
- README: added a "Structural credit (Merton)" section covering the
  `structural_credit` module (distance-to-default, default probability, credit
  spread, KMV asset-value solve), with runnable examples verified against the
  installed package.

## [1.465.0] - 2026-09-11

### Added
- `equity_volatility` and `solve_asset_value_and_vol` (in `structural_credit.py`):
  the Merton equity volatility from the asset volatility (`sigma_E = (V/E) N(d1)
  sigma_V`) and the KMV two-equation solve recovering the unobservable asset value
  and volatility from the observed equity value and equity volatility.
  Cross-checked: equity vol exceeds asset vol under leverage, and the solve
  recovers the true asset value and volatility (round-tripping the equity value).

## [1.464.0] - 2026-09-11

### Added
- `structural_credit.py`: the Merton (1974) structural credit model.
  `equity_value` (a call on the firm's assets struck at the debt face),
  `risk_neutral_default_probability` (`Phi(-d2)`), `distance_to_default` (`d2`),
  `risky_debt_value` (`V - E`), and `credit_spread`. Cross-checked: equity equals
  the BSM call, PD equals `Phi(-distance)`, the firm-value identity `V = E + debt`
  holds, and PD/spread rise with leverage and volatility (near zero for a safe
  firm).

## [1.463.0] - 2026-09-11

### Added
- `cumulative_default_term_structure` and `marginal_default_probabilities` (in
  `markov.py`): the cumulative and per-period default probabilities by horizon
  from a rating-migration matrix with an absorbing default state (the default
  column of `P^n`). Cross-checked: the cumulative curve is monotone and matches
  the n-step default column, riskier start ratings default more, marginals are
  non-negative and sum to the cumulative, and default converges to one.

## [1.462.0] - 2026-09-11

### Documentation
- README: added a "Markov chains" section covering the `markov` module (n-step
  transitions, stationary distribution, hitting times, absorbing-chain
  analytics), with runnable examples verified against the installed package.

## [1.461.0] - 2026-09-11

### Added
- `fundamental_matrix`, `expected_steps_to_absorption`, `absorption_probabilities`
  (in `markov.py`): absorbing-Markov-chain analytics (the credit-rating-migration
  toolkit). The fundamental matrix `N = (I - Q)^{-1}`, expected steps to
  absorption (its row sums), and the probability of ending in each absorbing state
  (`N R`). Cross-checked: `N` inverts `I - Q`, expected steps match the hitting
  time, and absorption probabilities sum to one across absorbing states.

## [1.460.0] - 2026-09-11

### Added
- `markov.py`: finite-state Markov chains. `n_step_transition` (`P^n` by repeated
  squaring), `stationary_distribution` (left eigenvector via power iteration), and
  `expected_hitting_time` (first-passage times by Gaussian elimination).
  Cross-checked: the stationary law matches the two-state closed form and is a
  fixed point of `P`, the n-step matrix stays row-stochastic and converges to the
  stationary rows, and the hitting time matches its recurrence solution.

## [1.459.0] - 2026-09-11

### Documentation
- README: added a "Bootstrap and jackknife" section covering the `resample`
  module (IID / stationary / BCa bootstrap and jackknife CIs), with runnable
  examples verified against the installed package.

## [1.458.0] - 2026-09-11

### Added
- `bca_bootstrap_ci` (in `resample.py`): the bias-corrected accelerated (BCa)
  bootstrap confidence interval, correcting the percentile method for median bias
  (`z0`) and skewness (`a`, from the jackknife). Cross-checked: it brackets the
  point estimate, is close to the plain percentile interval for symmetric data,
  differs for skewed data, and is deterministic per seed.

## [1.457.0] - 2026-09-11

### Added
- `resample.py`: bootstrap and jackknife confidence intervals. `bootstrap_ci`
  (IID percentile bootstrap), `stationary_bootstrap_ci` (Politis-Romano block
  bootstrap for serially-correlated series), and `jackknife_estimate` (delete-one
  SE). All deterministic per seed. Cross-checked: the CIs bracket the point
  estimate, the jackknife SE matches the analytic standard error of the mean, and
  a custom statistic (max) is supported. Default statistics use a ``None``
  sentinel so signatures render deterministically in the API docs.

## [1.456.0] - 2026-09-11

### Documentation
- README: added a "Performance attribution (Brinson)" section covering the
  `brinson` module (allocation/selection/interaction, Cariño linking), with a
  runnable example verified against the installed package.

## [1.455.0] - 2026-09-11

### Added
- `carino_factor`, `linked_active_return`, `carino_linked_effects` (in
  `brinson.py`): Cariño (1999) multi-period linking so per-period Brinson effects
  compound to the geometrically-linked active return. Cross-checked: the linked
  active return is the geometric compounding, the Cariño-smoothed effects sum to it
  (where naive arithmetic summation leaves a residual), and a single period reduces
  to the arithmetic active return.

## [1.454.0] - 2026-09-11

### Added
- `brinson.py`: Brinson-Hood-Beebower performance attribution. `allocation_effect`,
  `selection_effect`, `interaction_effect`, and `brinson_attribution` decompose a
  portfolio's active return versus a benchmark into per-segment allocation,
  selection, and interaction effects. Cross-checked: the three effects sum to the
  active return, an identical portfolio has zero effects, and same-return
  different-weight cases isolate the allocation effect.

## [1.453.0] - 2026-09-11

### Documentation
- README: added "Factor models" and "Dual-currency deposits" sections covering
  the `factor_model` and `dual_currency` modules, with runnable examples verified
  against the installed package.

## [1.452.0] - 2026-09-11

### Added
- `factor_attribution` and `rolling_factor_beta` (in `factor_model.py`): decompose
  a realized return into per-factor contributions plus a residual (alpha + factor
  contributions + residual sum to the total), and a trailing-window single-factor
  beta. Cross-checked: the attribution components sum to the total return, the
  contributions equal beta times factor return, and the rolling beta recovers a
  constant loading.

## [1.451.0] - 2026-09-11

### Added
- `factor_model.py`: multi-factor OLS return regression. `factor_regression`
  returns the alpha, factor betas, R-squared, and residual volatility (solved via
  the normal equations with pure-Python Gaussian elimination); `factor_expected_
  return` builds the expected return from betas and factor premia. Cross-checked:
  exact coefficient recovery and unit R-squared on a noiseless synthetic fit, the
  single-factor beta matching `realized_beta`, and R-squared below one under noise.

## [1.450.0] - 2026-09-11

### Added
- `dual_currency.py`: dual-currency deposit (DCD) analytics. `dcd_enhanced_yield`,
  `dcd_option_premium_rate` (Garman-Kohlhagen premium of the embedded sold
  option), `dcd_maturity_payoff` (with conversion), and `dcd_breakeven_spot`.
  Cross-checked: the enhanced yield exceeds the base deposit rate, higher vol
  raises the premium, conversion bites only past the strike, and at the breakeven
  spot the converted DCD equals a plain deposit.

## [1.449.0] - 2026-09-11

### Documentation
- README: added a "Bond futures" section covering the `bond_future` module
  (conversion factors, delivery basis, cheapest-to-deliver, implied repo, DV01
  hedging), with runnable examples verified against the installed package.

## [1.448.0] - 2026-09-11

### Added
- `bond_future_dv01`, `futures_dv01`, `futures_hedge_ratio` (in `bond_future.py`):
  a bond DV01, the bond-future DV01 geared by the CTD conversion factor
  (`ctd_dv01 / CF`), and the number of futures to hedge a cash bond. Cross-checked:
  the futures DV01 gears up for a sub-one conversion factor, and the hedge ratio's
  futures exactly offset the bond DV01. Exported as `bond_future_dv01` to avoid
  colliding with `bondmath`'s `bond_dv01`.

## [1.447.0] - 2026-09-11

### Added
- `bond_future.py`: bond-futures delivery analytics. `conversion_factor` (price at
  the notional coupon), `invoice_price`, `gross_basis`, `net_basis`,
  `implied_repo_rate`, and `cheapest_to_deliver`. Cross-checked: the conversion
  factor is one at the notional coupon and above/below for higher/lower coupons,
  the CTD bond minimizes the net basis, and the implied repo rises with the
  futures price.

## [1.446.0] - 2026-09-11

### Documentation
- README: added a "Liability-driven investing" section covering the `ldi` module
  (funding ratios, liability duration/convexity, duration hedging, immunization,
  surplus-at-risk), with runnable examples verified against the installed package.

## [1.445.0] - 2026-09-11

### Added
- `liability_convexity`, `surplus_change_under_shock`, `funded_ratio_return` (in
  `ldi.py`): the liability's interest-rate convexity, the second-order surplus
  change under a parallel rate shock, and the funding-ratio return from asset/
  liability returns. Cross-checked: convexity positive, a duration-and-convexity-
  matched book leaves the surplus stable under a shock while a duration-only match
  moves at second order, and the funding-ratio return is level-independent.

## [1.444.0] - 2026-09-11

### Added
- `ldi.py`: liability-driven investing. `liability_pv`, `funding_ratio`,
  `surplus`, `liability_duration`, `hedge_ratio`, `required_hedge_duration`, and
  `surplus_at_risk`. Cross-checked: the funding ratio and surplus, the required
  hedge duration giving a unit hedge ratio, the hedge ratio scaling with asset
  duration, and surplus-at-risk positive and rising with confidence.

## [1.443.0] - 2026-09-11

### Added
- `log_contract` and `log_contract_fair_variance` (in `exotics.py`): the log
  contract paying `ln(S_T / F)` (present value `e^{-rt} * (-0.5 sigma^2 t)`) and
  the variance-swap replication identity `sigma^2 = -2/t * E[ln(S_T/F)]`. The log
  contract is the theoretical basis of variance-swap replication; the fair
  variance recovers the input `sigma^2` exactly in the Black-Scholes world.
  Cross-checked against the closed form, the variance recovery, and horizon/vol
  monotonicity.

## [1.442.0] - 2026-09-11

### Added
- `contingent_premium_option` and `pay_later_option_value` (in `exotics.py`): the
  fair premium of a pay-later (contingent-premium) option -- premium paid at
  expiry only if in the money -- and the holder's value at a contracted premium.
  Cross-checked: the fair premium exceeds the vanilla price, the holder value is
  zero at the fair premium (positive below, negative above), and the premium/
  vanilla ratio approaches one deep in the money.

## [1.441.0] - 2026-09-11

### Documentation
- README: added a "Retirement decumulation" section covering the `retirement`
  module (sustainable withdrawals, depletion, glide paths, Monte Carlo ruin), with
  runnable examples verified against the installed package.

## [1.440.0] - 2026-09-11

### Added
- `withdrawal_stream_pv` and `ruin_probability_mc` (in `retirement.py`): the
  present value of a growing real withdrawal stream and the Monte Carlo
  probability of portfolio ruin under lognormal returns. Cross-checked: the PV
  matches the annuity sum and is monotone in withdrawal/horizon/growth/discount,
  and the ruin probability rises with the withdrawal rate and volatility, is near
  zero for low withdrawals, and is deterministic per seed.

## [1.439.0] - 2026-09-11

### Added
- `retirement.py`: deterministic decumulation tools. `portfolio_depletion_years`
  (annuity exhaustion, infinite when the withdrawal is below the interest),
  `sustainable_withdrawal` (its inverse), `withdrawal_balance_path` (inflation-
  indexed balance path), and `glide_path_equity_weight` (linear target-date glide).
  Cross-checked: depletion inverts the sustainable withdrawal, higher return
  extends the horizon, the 4%-rule limit, and a monotone clamped glide path.

## [1.438.0] - 2026-09-11

### Documentation
- README: added a "Copulas and portfolio credit" section covering the `copula`
  module (Gaussian/Clayton/Gumbel/Frank, tail dependence, credit-default
  applications, Vasicek/CDO loss), with runnable examples verified against the
  installed package.

## [1.437.0] - 2026-09-11

### Added
- `cdo_tranche_expected_loss_mc` (in `copula.py`): a finite-pool Monte Carlo of
  the single-factor default model that independently validates the large-pool
  `cdo_tranche_expected_loss`. Agrees within 10% across equity, mezzanine, and
  senior tranches at 500 names / 60k paths; deterministic per seed.

## [1.436.0] - 2026-09-11

### Added
- `vasicek_loss_cdf`, `vasicek_loss_quantile`, `cdo_tranche_expected_loss` (in
  `copula.py`): the Vasicek large-homogeneous-portfolio loss distribution (the
  single-factor Gaussian-copula limit behind the Basel IRB formula), its loss
  quantile, and the expected loss of a CDO tranche `[attachment, detachment]`.
  Cross-checked: the loss CDF is monotone and bounded, the quantile inverts it and
  rises with confidence/PD/correlation, the full-structure expected loss equals
  the pool PD, and equity tranches lose more than senior tranches.

## [1.435.0] - 2026-09-11

### Added
- `frank_copula`, `gaussian_copula_joint_default`, `first_to_default_probability`
  (in `copula.py`): the Frank copula (symmetric, signed dependence, no tail
  dependence) and Gaussian-copula credit applications -- the joint default
  probability of two names and the first-to-default probability. Cross-checked:
  Frank boundary/independence/symmetry, joint default above the independent
  product for positive correlation, and the FtD bounded by `[max pd, sum pd]` and
  falling as correlation rises.

## [1.434.0] - 2026-09-11

### Added
- `copula.py`: bivariate copulas. The Gaussian copula (via a Drezner-Wesolowsky
  bivariate-normal CDF), the Clayton and Gumbel Archimedean copulas, their lower/
  upper tail-dependence coefficients, and the Kendall's-tau calibration inverses.
  Cross-checked against the copula boundary conditions, the Frechet-Hoeffding
  bounds, the independence limits, closed-form tail dependence, and tau
  round-trips.

## [1.433.0] - 2026-09-11

### Documentation
- README: added "Equity swaps and dispersion" and "Futures/forward convexity"
  sections covering the `equity_swap` and `futures_convexity` modules, with
  runnable examples verified against the installed package.

## [1.432.0] - 2026-09-11

### Added
- `forward_curve_from_futures_strip` and `stub_discount_factors_from_forwards`
  (in `futures_convexity.py`): convert a strip of futures quotes to convexity-
  adjusted forward rates and bootstrap discount factors from the forward strip.
  Cross-checked: every forward sits below its futures rate, the adjustment grows
  down the curve, and the bootstrapped discount factors are decreasing.

## [1.431.0] - 2026-09-11

### Added
- `futures_convexity.py`: interest-rate futures/forward convexity adjustment.
  `ho_lee_convexity_adjustment` (`0.5 sigma^2 t1 t2`),
  `hull_white_convexity_adjustment` (mean-reverting, reducing to Ho-Lee as
  `a -> 0`), and `forward_from_futures` / `futures_from_forward` conversions.
  Cross-checked: the adjustment is non-negative and grows with maturity, the
  forward sits below the futures rate, the conversions round-trip, and mean
  reversion dampens the adjustment.

## [1.430.0] - 2026-09-11

### Added
- `dispersion_trade_pnl` (in `correlation.py`): the variance P&L of a dispersion
  trade (long the weighted member variances, short index variance). Profits when
  realized correlation comes in below the implied strike, loses when it comes in
  higher, and is zero when realized vols match the strikes. Cross-checked across
  all three correlation regimes.

## [1.429.0] - 2026-09-11

### Added
- `variance_swap_payoff`, `vega_notional_to_variance_notional`,
  `volatility_swap_payoff`, `variance_swap_mtm` (in `equity_swap.py`): payoffs and
  mark-to-market for variance and volatility swaps, plus the vega-to-variance
  notional conversion. Cross-checked: the variance swap is zero at the strike, the
  vol swap is linear, the variance swap's convexity makes it pay more than the vol
  swap on both sides of the strike, and the MTM reduces to the discounted
  expected payoff at inception and the realized payoff at expiry.

## [1.428.0] - 2026-09-11

### Added
- `equity_swap.py`: total-return and dividend swaps. `total_return_leg`,
  `financing_leg`, `total_return_swap_value`, and `trs_fair_spread` for TRS;
  `dividend_swap_fair_strike` (PV of expected dividends) and `dividend_swap_value`
  for dividend swaps. Cross-checked: the TRS value is the equity leg minus
  financing, the fair spread zeroes it, and the dividend swap is zero at the fair
  strike and positive when realized dividends exceed it.

## [1.427.0] - 2026-09-11

### Documentation
- README: added "Structured notes" and "Actuarial (life contingencies and cat
  bonds)" sections covering the `structured` and `actuarial` modules, with
  runnable examples verified against the installed package.

## [1.426.0] - 2026-09-11

### Added
- `cat_layer_loss`, `cat_expected_loss`, `cat_bond_spread`, `cat_bond_price` (in
  `actuarial.py`): catastrophe-bond analytics -- the loss ceded to an
  `[attachment, exhaustion]` layer, the expected loss rate over loss scenarios,
  the fair coupon spread (expected loss times a risk load), and the single-period
  cat-bond price. Cross-checked: the layer loss is bounded by the layer width, the
  expected loss falls as the attachment rises, the spread is at least the expected
  loss rate, and the price falls as expected loss rises.

## [1.425.0] - 2026-09-11

### Added
- `gompertz_makeham_hazard`, `gompertz_makeham_survival`,
  `gompertz_makeham_survival_curve`, `curtate_life_expectancy` (in `actuarial.py`):
  the Gompertz-Makeham force of mortality `a + b c^x`, its closed-form survival
  probability, a one-year survival curve for the life-table functions, and the
  curtate life expectancy `e_x = sum kp_x`. Cross-checked: the closed-form
  survival matches a numerical hazard integral, the `c -> 1` limit is the
  exponential (Makeham-only) form, and life expectancy falls with age.

## [1.424.0] - 2026-09-11

### Added
- `temporary_life_annuity_due` and `net_level_premium` (in `actuarial.py`): an
  n-year temporary life annuity-due and the equivalence-principle net level
  premium `P = A / a-due` for whole-life or endowment insurance. Cross-checked:
  the temporary annuity sits below whole life and rises to it with the term, and
  the premium times the annuity recovers the benefit EPV (equivalence principle).

## [1.423.0] - 2026-09-11

### Added
- `actuarial.py`: life-contingent expected present values. `survival_probabilities`
  (cumulative `kp_x` from one-year `p_x`), `life_annuity_due`, `term_insurance`,
  `whole_life_insurance`, `pure_endowment`, and `endowment_insurance`.
  Cross-checked against the annuity/insurance identity `A_x = 1 - d * a-due` on a
  terminating table, term below whole life, and endowment = term + pure endowment.

## [1.422.0] - 2026-09-11

### Added
- `capped_principal_protected_note`, `reverse_convertible_fair_coupon`,
  `buffered_note` (in `structured.py`): a PPN with a capped upside (call spread,
  at or below the uncapped PPN), the coupon that prices a reverse convertible at
  par, and a buffered note whose short put strikes below spot so it absorbs the
  first losses. Cross-checked: capped below uncapped and above the floor, fair
  coupon prices the note at par, and a larger buffer raises value monotonically.

## [1.421.0] - 2026-09-11

### Added
- `structured.py`: structured-note valuation by component decomposition.
  `note_zero_coupon_bond` (the guaranteed bond leg), `principal_protected_note`
  (ZC bond + participation * call, floored at the discounted principal),
  `reverse_convertible` (ZC bond + coupon - short put), and
  `note_embedded_option_value` (note minus the bond leg). Cross-checked against
  the component decompositions, the PPN principal floor, and coupon/participation
  monotonicity. `note_zero_coupon_bond` is named to avoid colliding with the
  Vasicek `zero_coupon_bond`.

## [1.420.0] - 2026-09-11

### Documentation
- README: added "Equity compensation and convertibles" and "Optimal execution"
  sections covering the `equity_comp` and `execution` modules, with runnable
  examples verified against the installed package.

## [1.419.0] - 2026-09-11

### Added
- `twap_schedule`, `vwap_schedule`, `pov_schedule` (in `execution.py`): standard
  execution benchmarks -- equal-size TWAP slices, VWAP slices proportional to a
  volume profile (reducing to TWAP when flat), and a percentage-of-volume schedule
  that trades a fixed fraction of each interval's volume with an optional total-
  shares cap. Cross-checked against equal TWAP slices, VWAP proportionality, the
  flat-profile equivalence, and POV proportionality with truncation.

## [1.418.0] - 2026-09-11

### Added
- `kyle_lambda`, `kyle_impact`, `square_root_impact`, `implementation_shortfall`
  (in `execution.py`): Kyle's linear price-impact coefficient and impact, the
  empirical square-root impact law (concave, scaling like `sqrt(size)`), and a
  decomposition of expected implementation shortfall into permanent, temporary,
  and timing-risk components. Cross-checked: Kyle impact is linear in size, the
  square-root law is concave per share, and the IS components sum to the expected
  cost.

## [1.417.0] - 2026-09-11

### Added
- `execution.py`: the Almgren-Chriss optimal trade-execution model.
  `execution_trajectory` gives the cost/risk-optimal holdings path (linear TWAP at
  zero risk aversion, `sinh` front-loading when risk-averse), `execution_trades`
  the per-interval sizes, `expected_cost` the permanent + temporary impact cost,
  `cost_variance` the timing risk, and `efficient_frontier_point` one cost/variance
  point per risk aversion. Cross-checked: the trajectory runs from the full
  position to zero, reduces to TWAP at `lambda = 0`, front-loads and lowers
  variance as risk aversion rises, and traces the cost/variance frontier.

## [1.416.0] - 2026-09-11

### Added
- `conversion_premium`, `investment_premium`, `convertible_breakeven_years` (in
  `equity_comp.py`): convertible-bond relative-value metrics -- the premium over
  parity, the premium over the bond floor, and the years for the bond's income
  advantage over the shares to recoup the conversion premium. Cross-checked: both
  premiums are non-negative above parity/floor, the conversion premium is zero at
  parity, and the breakeven is infinite when dividends exceed coupon income.

## [1.415.0] - 2026-09-11

### Added
- `conversion_value`, `straight_bond_floor`, `convertible_bond_value` (in
  `equity_comp.py`): the parity value of a convertible, its credit-spread-adjusted
  bond floor, and the component (bond floor + conversion call) valuation.
  Cross-checked: the convertible sits at or above both the bond floor and the
  conversion value, approaches the parity value deep in the money and the bond
  floor deep out of the money, and rises with volatility while falling with the
  credit spread.

## [1.414.0] - 2026-09-11

### Added
- `pv_dividends`, `discrete_dividend_price`, `forward_with_dividends` (in
  `equity_comp.py`): discrete cash-dividend handling via the escrowed-dividend
  (spot-minus-PV) method. Prices a European option on the dividend-adjusted spot
  and the corresponding forward. Cross-checked against the no-dividend vanilla
  limit, the adjusted-spot BSM price, put-call parity with dividends, and the
  exclusion of dividends after expiry.

## [1.413.0] - 2026-09-11

### Added
- `equity_comp.py`: warrants and employee stock options. `dilution_factor` and
  `warrant_price` scale a vanilla call by `M/(M+N)` for the new shares created on
  exercise; `eso_expected_life` gives the Hull-White expected life under a
  post-vest exit rate; `eso_value` prices an ESO at that expected life times the
  pre-vest survival probability (FASB 123R practical model). Cross-checked: the
  warrant sits below the vanilla call by exactly the dilution factor, the ESO
  below the full-term call, and both reduce to the vanilla with no dilution/exit.

## [1.412.0] - 2026-09-11

### Documentation
- README: added "Mortgage-backed securities and CMOs" and "Weather derivatives"
  sections covering the `mbs` and `weather` modules, with runnable examples
  verified against the installed package.

## [1.411.0] - 2026-09-11

### Added
- `pac_schedule` and `pac_support_split` (in `mbs.py`): the planned-amortization-
  class principal schedule as the lower envelope of the pool principal at the two
  ends of a PSA collar, and the allocation of pool principal between the PAC (with
  arrears carried forward) and its support/companion tranche. Cross-checked: the
  schedule is a lower envelope, PAC + support sums to the pool principal, the PAC
  WAL is stable across speeds inside the band, and the support absorbs more WAL
  variability than the PAC.

## [1.410.0] - 2026-09-11

### Added
- `sequential_cmo` and `tranche_wal` (in `mbs.py`): split MBS principal across
  sequential-pay CMO tranches (principal waterfalls strictly in priority order)
  and compute each tranche's weighted-average life. Cross-checked: tranche
  principals sum to their sizes, earlier tranches retire first (shorter WAL), and
  the size-weighted blend of tranche WALs equals the pool WAL.

## [1.409.0] - 2026-09-11

### Added
- `mbs_price_with_spread`, `mbs_zspread`, `mbs_effective_duration`,
  `mbs_effective_convexity` (in `mbs.py`): discount MBS cashflows off a zero curve
  plus a parallel static spread, solve the Z-spread reproducing a price, and
  measure effective duration/convexity from a parallel yield bump (static
  cashflows). Cross-checked against the flat-curve price, Z-spread round-trip on
  flat and sloped curves, and duration matching a central finite difference.

## [1.408.0] - 2026-09-11

### Added
- `mbs_cashflows_psa`, `mbs_price`, `mbs_yield` (in `mbs.py`): MBS cashflows on the
  PSA prepayment ramp (age-varying SMM), and the price/yield of a projected
  cashflow strip at a monthly-compounded yield. Faster PSA shortens the WAL, the
  price is monotone decreasing in yield and equals par at the coupon rate, and the
  yield inverts the price.

## [1.407.0] - 2026-09-11

### Added
- `mbs.py`: mortgage-backed security cashflows and prepayment conventions.
  `monthly_payment` (level fully-amortizing payment), `amortization_schedule`,
  `cpr_to_smm` / `smm_to_cpr` (annual/monthly prepayment conversions), `psa_cpr`
  (the PSA ramp), `mbs_cashflows` (scheduled plus prepaid principal at a constant
  SMM), and `weighted_average_life`. Cross-checked: the balance amortizes to zero,
  principal sums to the original balance with and without prepayment, CPR/SMM
  round-trips, and the WAL shortens as prepayment speeds up.

## [1.406.0] - 2026-09-11

### Added
- `tolling_value` (in `commodity.py`): values a tolling agreement as a strip of
  daily spark-spread call options (the right to run a plant each period), summing
  `spark_spread_option` across the delivery periods at a common heat rate, strike,
  and vols, with optional per-period carbon and a discount-factor override.
  Increasing in the number of run periods; cross-checked against the manual sum,
  the discount override, and the emissions effect.

## [1.405.0] - 2026-09-11

### Added
- `spark_spread_option` (in `commodity.py`): a spark/dark-spread option on a power
  generator's clean margin `power - heat_rate * fuel - emissions_rate * carbon`,
  built on the normal `bachelier_spread_option` (generation margins routinely go
  negative). The heat rate scales both the fuel forward and its volatility, and an
  optional carbon adder charges emissions. Cross-checked against the underlying
  Bachelier price, put-call parity (including the carbon leg), negative margins,
  and the emissions effect.

## [1.404.0] - 2026-09-11

### Added
- `seasonal_mean_temperature`, `expected_temperature`, `temperature_variance` (in
  `weather.py`): the Alaton-Djehiche-Stillberger mean-reverting temperature model.
  A deterministic seasonal trend (linear plus annual sinusoid) is the reversion
  level of an Ornstein-Uhlenbeck temperature; the expected temperature equals the
  current temperature at horizon zero and relaxes to the seasonal mean, while the
  variance rises to the stationary `sigma^2/(2 kappa)`. Cross-checked against the
  end-point limits, mean reversion, and monotone variance.

## [1.403.0] - 2026-09-11

### Added
- `degree_day_swap_rate`, `degree_day_collar`, `degree_day_option_mc` (in
  `weather.py`): the fair degree-day swap strike (the expected index), a
  cap/floor collar (long call, short put; reduces to the discounted forward at a
  common strike), and a Monte Carlo option over simulated daily temperatures that
  independently validates the Bachelier `degree_day_option` (agrees within 10%).
  Cross-checked against the zero-swap-value fair strike, collar/forward parity,
  and MC agreement.

## [1.402.0] - 2026-09-11

### Added
- `weather.py`: weather derivatives. `heating_degree_days` / `cooling_degree_days`
  accumulate daily degree days versus a base temperature; `degree_day_index`
  dispatches on kind; `degree_day_swap_payoff` is the linear tick-value swap; and
  `degree_day_option` prices an option on the accumulated index with the Bachelier
  (normal) model (appropriate for a sum-of-days total), with an optional payoff
  cap priced as a call/put spread. Cross-checked against HDD/CDD complementarity,
  put-call parity, the zero-vol intrinsic, and the cap bound.

## [1.401.0] - 2026-09-11

### Documentation
- README: added a "Counterparty valuation adjustments (XVA)" section covering the
  `xva` module (CVA/DVA/BCVA, FVA/MVA, EPE/PFE profiles, wrong-way risk, CSA
  collateral), with a runnable example verified against the installed package.

## [1.400.0] - 2026-09-11

### Added
- `mva` and `swap_cva` (in `xva.py`): the margin valuation adjustment (funding
  cost of posted initial margin over the trade's life, proportional to the
  funding spread and margin, optionally survival-weighted) and a one-shot par-swap
  CVA convenience that builds the analytic exposure profile and feeds it to `cva`.
  Cross-checked against MVA proportionality and the manual EPE-then-CVA
  composition.

## [1.399.0] - 2026-09-11

### Added
- `collateralized_exposure` and `collateralized_exposure_profile` (in `xva.py`):
  residual exposure under a CSA collateral agreement,
  `max(min(E, threshold + MTA) - independent_amount, 0)`. An infinite threshold
  recovers the uncollateralized exposure, a zero threshold with no MTA leaves only
  the independent-amount offset, and collateralization reduces the CVA. Applied
  across a profile with the vector helper.

## [1.398.0] - 2026-09-11

### Added
- `swap_potential_future_exposure` and `wrong_way_cva` (in `xva.py`): the PFE
  profile of a par swap at a high quantile (`std(t) * Phi^{-1}(q)`, sitting above
  the EPE and rising with the quantile) and a CVA with a linear wrong-way-risk
  tilt of the default buckets toward later dates (renormalized to preserve total
  default probability). `alpha = 0` reduces to `cva`; `alpha > 0` raises the CVA
  when exposure rises with time (wrong-way), `alpha < 0` lowers it (right-way).

## [1.397.0] - 2026-09-11

### Added
- `swap_expected_exposure` and `fva` (in `xva.py`): the analytic expected positive
  exposure profile of a par swap/forward (value diffuses as Brownian motion and
  amortizes linearly with remaining life; `EPE = std/sqrt(2 pi)`, zero at both
  ends and humped in between) and the funding valuation adjustment on an
  uncollateralized exposure (proportional to the funding spread and exposure,
  optionally survival-weighted). Cross-checked against the end-point zeros, the
  hump, the closed-form EPE, and FVA proportionality.

## [1.396.0] - 2026-09-11

### Added
- `xva.py`: counterparty valuation adjustments. `marginal_default_probs` buckets
  a `SurvivalCurve` into per-step default probabilities; `cva` discounts the
  expected loss `LGD * sum EE(t) DF(t) dQ(t)` from counterparty default; `dva` is
  the mirror on our own default; `bcva = cva - dva` is the bilateral adjustment.
  Cross-checked: CVA non-negative, zero without exposure, monotone in hazard and
  exposure, exact LGD scaling, and BCVA zero for identical curves/exposures.

## [1.395.0] - 2026-09-11

### Added
- `schwartz_smith_futures_volatility` (in `commodity.py`): instantaneous return
  volatility of a maturity-`T` future under the Schwartz-Smith two-factor model,
  `sqrt(e^{-2 kappa T} sigma_chi^2 + sigma_xi^2 + 2 e^{-kappa T} rho sigma_chi
  sigma_xi)`. Falls from the front toward the long-term floor `sigma_xi` as
  `T -> inf` (Samuelson effect with a non-zero long-end asymptote) and reduces to
  the one-factor `schwartz_futures_volatility` when `sigma_xi = 0`.

## [1.394.0] - 2026-09-11

### Added
- `schwartz_smith_log_mean`, `schwartz_smith_log_variance`,
  `schwartz_smith_forward` (in `commodity.py`): the Schwartz-Smith (2000)
  two-factor model, decomposing log-spot into a mean-reverting short-term
  deviation and a drifting long-term equilibrium. Equals the spot at `T = 0`; the
  short-term contribution decays while the long-term factor dominates the variance
  (growing like `sigma_xi^2 T`) at long horizons. Reduces exactly to the
  one-factor `schwartz_forward` when the long-term vol, drift, and equilibrium are
  turned off.

## [1.393.0] - 2026-09-11

### Documentation
- README: added "Inflation-linked bonds and derivatives" and "Commodities
  (cost-of-carry and mean reversion)" sections covering the two modules built
  since v1.372, with runnable examples (all snippets verified against the
  installed package).

## [1.392.0] - 2026-09-11

### Added
- `turnbull_wakeman_asian` (in `commodity.py`): Turnbull-Wakeman (1991) two-moment
  matched arithmetic-average Asian option. Matches the exact first two moments of
  the arithmetic average of lognormals to a lognormal (effective variance
  `ln(M2/M1^2)`) and applies Black. More accurate than the fixed 1/3-variance
  `asian_commodity_option` -- closer to a 5-seed 400k-path Monte Carlo mean --
  reduces to the vanilla at `n = 1`, stays above the geometric-average lower
  bound, and satisfies put-call parity.

## [1.391.0] - 2026-09-11

### Added
- `geometric_asian_option` and `asian_commodity_option_mc` (in `commodity.py`):
  the exact discrete geometric-average Asian option (the geometric average of
  lognormals is lognormal, giving a closed form) and a discrete-path Monte Carlo
  Asian pricer (arithmetic or geometric). The geometric closed form is a lower
  bound for the arithmetic `asian_commodity_option` (AM-GM); MC validates the
  geometric closed form (within 2%) and the arithmetic 1/3-variance approximation
  (within 5%). Put-call parity holds on the geometric-average forward.

## [1.390.0] - 2026-09-11

### Added
- `commodity_swap_rate`, `commodity_swap_value`, `asian_commodity_option` (in
  `commodity.py`): the fair fixed price of a commodity swap (DF-weighted average
  of the reset forwards, zeroing the swap PV), the swap value versus a fixed
  price, and an average-price (Asian) commodity option priced Black on the
  average forward with the variance reduced toward the continuous-averaging `1/3`
  limit. Cross-checked against zero-PV at the fair rate, put-call parity, and the
  Asian being cheaper than the vanilla on the same forward.

## [1.389.0] - 2026-09-11

### Added
- `bachelier_spread_option` and `spread_option_mc` (in `commodity.py`): a normal-
  model spread option that prices negative or through-zero spreads where the
  lognormal Kirk struggles (crack/location spreads), and a bivariate-lognormal
  Monte Carlo pricer that independently validates `kirk_spread_option` (agrees
  within 2% at 400k paths). Normal call/put satisfy `C - P = e^{-rT}(F1-F2-K)`;
  the MC stream is deterministic per seed.

## [1.388.0] - 2026-09-11

### Added
- `margrabe_exchange_option` and `kirk_spread_option` (in `commodity.py`): the
  exact Margrabe (1978) option to exchange one forward for another (zero-strike
  spread) and the Kirk (1995) approximation for a struck spread option on two
  forwards (crack/spark spreads). Kirk collapses exactly to Margrabe at zero
  strike, call/put satisfy `C - P = e^{-rT}(F1 - F2 - K)`, and both fall as
  correlation rises (the spread vol shrinks). Cross-checked on all three.

## [1.387.0] - 2026-09-11

### Added
- `roll_yield`, `carry_roll_yield`, `schwartz_futures_volatility` (in
  `commodity.py`): annualized roll yield `ln(F_near/F_far)/(t_far-t_near)`
  (positive in backwardation), the carry-model roll yield `-net_cost_of_carry`
  (matching the forward-implied roll yield exactly), and the Schwartz futures
  return volatility `sigma e^{-kappa T}` capturing the Samuelson effect (front
  contracts more volatile than deferred). Cross-checked against the carry sign,
  carry/roll consistency, and monotone vol decay.

## [1.386.0] - 2026-09-11

### Added
- `schwartz_option`, `mean_reversion_half_life`, `schwartz_implied_alpha` (in
  `commodity.py`): European spot option under the Schwartz one-factor model
  (Black-form off the model forward and log-spot variance), the mean-reversion
  half-life `ln(2)/kappa`, and calibration of the risk-neutral long-run level
  `alpha*` from a single forward quote. Cross-checked against put-call parity
  `C - P = e^{-rT}(F* - K)`, the zero-vol intrinsic, and calibration round-trip.

## [1.385.0] - 2026-09-11

### Added
- `schwartz_log_mean`, `schwartz_log_variance`, `schwartz_forward` (in
  `commodity.py`): Schwartz (1997) one-factor mean-reverting model for commodity
  forwards. Log-spot follows a risk-neutral Ornstein-Uhlenbeck process; the
  forward is `exp(E[X_T] + 0.5 Var[X_T])`. Equals the spot at `T = 0` and
  converges to the long-run forward `exp(alpha* + sigma^2/(4 kappa))` as
  `T -> inf`, with variance rising monotonically to `sigma^2/(2 kappa)` -- the
  mean-reverting alternative to the constant-carry forward.

## [1.384.0] - 2026-09-11

### Added
- `commodity_calendar_spread`, `convenience_yield_curve`, `seasonal_forward` (in
  `commodity.py`): far-minus-near forward spread (sign tracks the net carry, so
  positive in contango and negative in backwardation), per-tenor convenience
  yields bootstrapped from a forward strip (reprices each quote), and a
  seasonally-scaled carry forward. Cross-checked against the carry sign, exact
  repricing, and the unit-factor limit.

## [1.383.0] - 2026-09-11

### Added
- `commodity.py`: cost-of-carry forward pricing for storable commodities.
  `commodity_forward` prices `S exp((r + u - y) T)` with storage cost `u` and
  convenience yield `y`; `implied_convenience_yield` and `implied_storage_cost`
  invert a market forward for each; `net_cost_of_carry`, `commodity_forward_curve`
  and `is_backwardation` give the net carry rate, the forward curve, and the
  contango/backwardation classification. Cross-checked against the carry parity,
  round-trip inversion, and curve monotonicity.

## [1.382.0] - 2026-09-11

### Added
- `real_zero_curve`, `nominal_zero_curve`, `real_discount_factor` (in
  `inflation.py`): bridge between nominal and real curves via the Fisher relation
  per tenor. `real_zero_curve` strips breakevens off a nominal zero curve,
  `nominal_zero_curve` reconstructs it (composing the two is the identity), and
  `real_discount_factor` grows a nominal DF by the period index ratio so a real
  cashflow at the real DF equals its inflated nominal cashflow at the nominal DF.

## [1.381.0] - 2026-09-11

### Added
- `yoy_caplet_price_normal` and `yoy_caplet_implied_normal_vol` (in
  `inflation.py`): Bachelier (normal-model) year-on-year inflation caplet/floorlet
  and its implied-normal-vol inverse. Arithmetic Brownian dynamics admit zero and
  negative inflation forwards where the lognormal Black-76 breaks down.
  Cross-checked against normal cap/floor parity, the zero-vol intrinsic, the ATM
  closed form `DF*N*sigma*sqrt(T/2pi)`, and implied-vol round-trip (including a
  negative forward).

## [1.380.0] - 2026-09-11

### Added
- `yoy_cap_price` and `yoy_cap_implied_vol` (in `inflation.py`): a year-on-year
  inflation cap/floor as a strip of Black-76 caplets across periods at a common
  strike and flat vol, plus the bisection inverse for the flat implied vol. A
  single-period strip equals the caplet; cap minus floor telescopes to
  `sum_i DF_i * N * (F_i - K)`; implied vol round-trips.

## [1.379.0] - 2026-09-11

### Added
- `yoy_caplet_price` (in `inflation.py`): Black-76 price of a single year-on-year
  inflation caplet or floorlet, modelling the YoY rate as lognormal around its
  forward with volatility to expiry. Cross-checked against cap/floor parity
  `cap - floor = DF*N*(F-K)`, the zero-vol discounted intrinsic, monotonicity in
  vol, and ATM cap/floor symmetry.

## [1.378.0] - 2026-09-11

### Added
- `normalize_seasonal_factors`, `apply_seasonality`, `deseasonalize` (in
  `inflation.py`): CPI seasonal adjustment. Twelve monthly factors are normalized
  to a geometric mean of one so they multiply to one over the year (no trend
  shift); apply/remove multiply/divide the index by the monthly factor and round-
  trip to the identity. Normalization is idempotent and preserves the relative
  month-to-month shape.

## [1.377.0] - 2026-09-11

### Added
- `linker_real_duration`, `linker_real_convexity`, `linker_real_dv01` (in
  `inflation.py`): real-yield risk of an inflation-linked bond. The index ratio
  multiplies the whole price so it cancels in the fractional duration and
  convexity -- these equal `bondmath`'s `modified_duration`/`convexity` on the
  real cashflows, independent of the index level -- while the dollar DV01 scales
  with the ratio. Cross-checked against `bondmath` and central finite differences.

## [1.376.0] - 2026-09-11

### Added
- `reference_cpi` and `index_ratio_interpolated` (in `inflation.py`): daily
  reference index by linear interpolation between two monthly CPI fixings (the
  standard linker daily-indexation rule), and the corresponding intra-month index
  ratio. Equals the month-start fixing on the 1st, the mean of the anchors at
  mid-month, and reduces to the boundary `index_ratio` on day one.

## [1.375.0] - 2026-09-11

### Added
- `inflation_curve_from_zc_swaps` (in `inflation.py`): projected index levels
  `I_0 (1+k_T)^T` implied by a strip of zero-coupon inflation swap rates -- the
  market inflation curve as forward fixings; reinverting reprices the input swaps.
- `forward_inflation_rate`: annualized forward inflation between two curve
  horizons, chaining with the near leg via `(1+spot)^t1 (1+fwd)^(t2-t1) =
  (1+spot_end)^t2`.
- `yoy_swap_value`: value of a year-on-year inflation swap off a projected index
  curve, paying realized annual inflation each period against a fixed rate; zero
  when the fixed rate matches flat realized inflation.

## [1.374.0] - 2026-09-11

### Added
- `deflation_floored_redemption` and `deflation_floor_value` (in `inflation.py`):
  TIPS-style principal redemption floored at par, `face * max(index_ratio, 1)`,
  and the intrinsic value of the embedded deflation floor. Equals the adjusted
  principal in inflation and binds to par under net deflation.
- `yoy_inflation_rate`, `zc_inflation_swap_rate`, `zc_inflation_swap_value`
  (in `inflation.py`): year-on-year inflation between fixings, the par rate of a
  zero-coupon inflation swap `(I_T/I_0)^(1/T) - 1`, and the inflation-leg-receiver
  swap value. Par rate cross-checked against the compounding identity
  `(1+k)^T * I_0 == I_T` and the par swap value being zero.

## [1.373.0] - 2026-09-11

### Added
- `linker_price` and `linker_real_yield` (in `inflation.py`): dirty price of an
  inflation-linked bond from real (constant-money) cashflows discounted at a
  continuously-compounded real yield and inflated by the settlement index ratio,
  and the bisection inverse for the real yield. Because the index ratio multiplies
  every flow the price is degree-one homogeneous in it; stripping the ratio
  recovers a standard real-yield bond, cross-checked against `bondmath`'s
  `bond_price_from_yield` and `yield_to_maturity`.

## [1.372.0] - 2026-09-11

### Added
- `inflation.py`: inflation-linked (TIPS-style) bond and breakeven-inflation
  analytics. `index_ratio` (CPI_settle / CPI_base) and
  `inflation_adjusted_principal` scale a linker's notional by realized inflation.
  `fisher_real_rate` and `fisher_nominal_rate` are the exact Fisher conversions
  `(1+n)/(1+i)-1` and `(1+r)(1+i)-1` (their round-trip is the identity).
  `breakeven_inflation` returns the inflation rate equating a nominal and a real
  yield, `(1+n)/(1+r)-1`, with `real_from_breakeven` its inverse. Cross-checked
  against the Fisher round-trip and the small-rate approximation n - r.

## [1.371.0] - 2026-09-11

### Added
- `dated_accrued_interest` and `dated_clean_price` (in `bondmath.py`): calendar-
  accurate accrued interest and clean quote. Accrued interest accrues the coupon
  by the day-count fraction of the current period elapsed at settlement
  (`year_fraction(prev, settle) / year_fraction(prev, next)`); the clean price
  is the dated dirty price minus that accrued.
- Verified: accrued is half a coupon three months into a semiannual period and
  zero at the coupon date, growing toward the next coupon; the clean price is
  the dated dirty price minus accrued and coincides with it at a coupon date; a
  settlement outside the coupon period raises.

## [1.370.0] - 2026-09-11

### Added
- `dated_bond_price` and `dated_bond_yield` (in `bondmath.py`): price and yield
  of a bond from `dated_bond_cashflows` discounted off an explicit settlement
  date. Each cashflow is discounted by `exp(-y * T)` with `T` the day-count year
  fraction from settle to the pay date; cashflows on or before settle drop out.
  `dated_bond_yield` inverts this by bisection.
- Verified: the dated price matches the uniform `bond_price_from_yield` on a
  30/360 semiannual bond; the yield round-trips; settling after a coupon drops
  it (lower price); higher yields give lower prices; a non-positive price raises.

## [1.369.0] - 2026-09-11

### Added
- `dated_bond_cashflows` (in `bondmath.py`): coupon-bond cashflows on a real
  calendar. Generates the coupon schedule with `generate_schedule`, weights each
  coupon by the period's `year_fraction` under a day-count convention, and pays
  `face * coupon_rate * tau` per period plus the face at maturity -- reflecting
  actual day counts, month-end roll, and business-day adjustment rather than the
  uniform `1/freq` of `bond_cashflows`.
- Verified: a 30/360 semiannual bond gives four periods of tau 0.5 and 2.5
  coupons with the face on the last flow; act/360 accrual exceeds 0.5 over a
  half year; a zero coupon pays only the face; the end-of-month roll lands on
  month ends; a negative coupon raises.

## [1.368.0] - 2026-09-11

### Added
- New module `schedule.py` with `generate_schedule` and `adjust_business_day`
  over `(year, month, day)` dates. Generates periodic period-end dates from a
  start, maturity, and tenor in months (with optional end-of-month roll) and
  adjusts off weekends per `following` / `modified_following` / `preceding` /
  `unadjusted` conventions.
- Verified: a quarterly end-of-month schedule snaps to month ends
  (Apr 30 / Jul 31 / Oct 31 / Jan 31); a 5y semiannual schedule has 10 periods;
  the month-end roll clamps Jan 31 + 1m to Feb 29 in a leap year; the following
  convention rolls a Saturday to Monday, preceding to Friday, and
  modified-following rolls back to stay in-month; unadjusted is identity; bad
  frequency/maturity/convention raise.

## [1.367.0] - 2026-09-11

### Added
- New module `daycount.py` with `year_fraction` and `day_count` over
  `(year, month, day)` dates, covering the standard market conventions:
  `act/360`, `act/365`, `30/360` (US bond basis), `30E/360` (Eurobond), and
  `act/act` (ISDA, split across year boundaries).
- Verified: act/360 over 180 days is 0.5; a full leap year is 366/365 under
  act/365 and exactly 1 under act/act; 30/360 gives a clean year and applies the
  month-end rule (Jan 31 -> Feb 28 = 28 days); 30E/360 caps the day at 30;
  act/act across a year boundary splits the leap/non-leap stubs; reversed dates
  and unknown conventions raise.

## [1.366.0] - 2026-09-11

### Added
- `compounded_rate_with_lookback` and `compounded_rate_with_lockout` (in
  `rates.py`): real-world RFR compounding conventions. Lookback (observation
  shift) uses the fixing observed `k` business days earlier for each accrual;
  lockout freezes the last `k` days at the final observed fixing so the coupon
  is known before period end. Both reduce to `compounded_overnight_rate` at
  `k = 0`.
- Verified: `k = 0` reproduces the base compounded rate; on a rising fixing
  series both lower the coupon; a full lockout uses the first fixing for the
  whole period; negative `k` and `k >= n` raise.

## [1.365.0] - 2026-09-11

### Added
- `compounded_overnight_rate` and `simple_average_rate` (in `rates.py`):
  overnight-rate coupon conventions. The compounded (SOFR/SONIA setting-in-
  arrears) rate multiplies daily growth factors `(1 + r_i tau_i)` and annualizes
  by the total accrual; the simple version (Fed-funds style) is the
  accrual-weighted arithmetic average.
- Verified: the compounded rate exceeds the simple average for positive fixings
  (interest-on-interest); the simple average of a flat series is the rate; a
  single fixing gives both equal to it; the weighted average and the compounded
  product match hand calculations; length mismatches raise.

## [1.364.0] - 2026-09-11

### Added
- `cumulative_return`, `annualized_return`, and `annualized_volatility` (in
  `perfmetrics.py`): the compounded total return `prod(1+r)-1`, the geometric
  annualized return `(prod(1+r))^{ppy/n}-1`, and the sample volatility scaled by
  `sqrt(periods_per_year)`.
- Verified: cumulative return compounds (+10% then -10% gives -1%); the
  annualized return is zero for a flat series and matches the compound formula;
  annualized volatility equals the scaled sample stdev and is zero for a
  constant series; empty/short series raise.

## [1.363.0] - 2026-09-11

### Added
- `historical_var_series` and `historical_cvar` (in `perfmetrics.py`): empirical
  Value-at-Risk and conditional VaR (expected shortfall) from a raw return
  series -- the `1 - confidence` order-statistic quantile (negated to a loss)
  and the mean of the returns at or below it. Distribution-free, distinct from
  the option-book `historical_var` in `risk`; exported as `historical_var_series`
  to avoid the name clash.
- Verified: both are positive losses; CVaR is at least the VaR; both rise with
  confidence; the fraction of returns below the 95% VaR is ~5%; too-short series
  raise.

## [1.362.0] - 2026-09-11

### Added
- `cornish_fisher_var` (in `perfmetrics.py`): skew/kurtosis-adjusted historical
  Value-at-Risk. Expands the lower-tail normal quantile with the sample skewness
  and excess kurtosis (Cornish-Fisher) so the VaR reflects a fat or asymmetric
  return distribution, then scales the mean linearly and the deviation by
  `sqrt(horizon)`. Returned as a positive loss.
- Verified: it reduces to the parametric VaR on a normal sample; is a positive
  loss; a negative-skew (fat left tail) series pushes it above the Gaussian VaR;
  rises with confidence and horizon; zero-variance/too-short series raise.

## [1.361.0] - 2026-09-11

### Added
- `sample_skewness`, `sample_kurtosis`, and `jarque_bera` (in `perfmetrics.py`):
  return-distribution diagnostics. Population-convention third/fourth
  standardized moments (kurtosis excess by default, so a normal reads 0) and the
  Jarque-Bera normality statistic `n/6 (skew^2 + exkurt^2/4)`.
- Verified: a large normal sample gives near-zero skew/kurtosis and a
  Jarque-Bera below the 5% chi-squared(2) critical value; a symmetric series has
  zero skew; a right-skewed series is positive; fat tails give positive excess
  kurtosis; raw minus excess kurtosis is 3; Jarque-Bera rejects a skewed series;
  zero-variance/too-short series raise.

## [1.360.0] - 2026-09-11

### Added
- `rho_discount` (in `bsm.py`): discount-only rho `dPrice/dr` holding the cost
  of carry `b` fixed, which equals `-t * price` for both calls and puts. This is
  the correct rate sensitivity when the carry is independent of the funding rate
  (Black-76 on a future, an FX or commodity forward), where a rate change moves
  only the discount factor -- contrast `rho`, which assumes `b` moves with `r`.
- Verified: it matches a central finite difference of the price in `r` at fixed
  `b = 0` for calls and puts across rates; equals `-t * price`; is negative for
  both; a zero maturity gives zero.

## [1.359.0] - 2026-09-11

### Changed
- Trimmed the slowest Monte Carlo cross-check tests without losing coverage:
  the partial-time barrier MC tests (start/end) now run one representative `t1`
  at half the paths (the limit/parity/monotonicity tests already sweep `t1`),
  and the Vasicek/CIR/Ho-Lee rate-moment MC tests drop from 120k to 60k paths.
  This cuts roughly 190 seconds off the `slow` suite (partial-time start
  ~150s -> ~21s, end ~75s -> ~9s, moments ~44s -> ~13s) while keeping every
  closed-form-vs-MC check.

## [1.358.0] - 2026-09-11

### Documentation
- Added README sections for the modules built out over recent releases: fixed
  income (coupon-bond analytics, key-rate durations, short-rate moments), credit
  (survival curves, CDS, defaultable bonds), FX forwards under covered interest
  parity, portfolio optimization (min-variance / max-Sharpe / risk-parity /
  efficient frontier / Black-Litterman with VaR budgeting), performance metrics,
  and GARCH volatility with term-vol option pricing.
- Verified: every new README example imports and runs against the current API.

## [1.357.0] - 2026-09-11

### Added
- `box_spread_implied_rate` (in `bsm.py`): the continuously-compounded
  financing rate implied by a box-spread price, `r = -ln(box/(K2 - K1))/t`. The
  box has a riskless `K2 - K1` payoff, so its price pins the synthetic-lending
  rate the options market charges, independent of the underlying.
- Verified: it recovers the input rate from a `box_spread` priced at that rate
  (for two rates); a box above its notional implies a negative rate and one at
  notional a zero rate; reversed strikes, zero time, and a non-positive price
  raise.

## [1.356.0] - 2026-09-11

### Added
- `implied_discount_factor` (in `bsm.py`): the discount factor to expiry read
  from two same-expiry call-put pairs, `DF = [(C1 - P1) - (C2 - P2)]/(K2 - K1)`.
  Subtracting the two put-call parity relations cancels the forward, so this is
  model-free -- no rate and no volatility.
- Verified: it recovers `e^{-rt}` from Black-Scholes pairs (with and without
  dividends), is independent of strike order, and together with
  `implied_forward_from_parity` reconstructs each pair's `C - P`; equal strikes
  raise.

## [1.355.0] - 2026-09-11

### Added
- `implied_forward_from_parity` (in `bsm.py`, exported under that name to avoid
  clashing with the chain-based `implied_forward`): the forward implied by a
  same-strike call-put pair, `F = K + e^{rt}(C - P)`, from put-call parity with
  no volatility input.
- Verified: it recovers the cost-of-carry forward from Black-Scholes call/put
  prices with and without dividends; it is strike-independent across three
  strikes; at the forward strike the call and put are equal and the implied
  forward equals that strike.

## [1.354.0] - 2026-09-11

### Added
- `forward_price` and `put_call_parity_residual` (in `bsm.py`): the cost-of-
  carry forward `F = S e^{b t}` and the put-call parity residual
  `(C - P) - e^{-rt}(F - K)`, which is zero for arbitrage-consistent quotes and
  otherwise reports the violation in price terms.
- Verified: the forward matches `S e^{b t}` (spot for a future `b = 0`, and the
  dividend case `b = r - q`); Black-Scholes call/put prices give a zero parity
  residual with and without dividends; a deliberate mispricing surfaces as the
  exact residual.

## [1.353.0] - 2026-09-11

### Added
- `dual_delta` and `dual_gamma` (in `greeks2.py`): strike sensitivities of the
  Black-Scholes price. `dual_delta = d(price)/d(strike)` is `-e^{-rt} N(d2)`
  (call) / `e^{-rt} N(-d2)` (put); `dual_gamma = d^2(price)/d(strike)^2 =
  e^{-rt} phi(d2)/(K sigma sqrt t)`, the Breeden-Litzenberger risk-neutral
  density (same for calls and puts, non-negative).
- Verified: dual delta matches a central finite difference of the price in the
  strike for calls and puts; the negated call dual delta equals the
  cash-or-nothing digital; dual gamma matches the second strike difference and
  is non-negative; call minus put dual delta equals `-e^{-rt}` (parity).

## [1.352.0] - 2026-09-11

### Added
- `fx_forward_from_curves` (in `fxforward.py`): the FX forward from two discount
  curves, `S * DF_base(t) / DF_price(t)` (curve-based covered interest parity).
  Each curve is any `curve.df(t)` or plain callable.
- Verified: it matches the closed-form `fx_forward` on flat curves; a higher
  base-currency rate gives a forward discount; plain callables work; a zero
  tenor returns spot; bad spot/tenor raise.

## [1.351.0] - 2026-09-11

### Added
- New module `fxforward.py` for FX forwards under covered interest parity:
  `fx_forward` (`S exp((r_price - r_base) t)`), `forward_points`,
  `fx_swap_points` (forward-forward), and `implied_base_rate` /
  `implied_price_rate` (invert CIP from a quoted forward).
- Verified: the forward matches the CIP formula and shows a base premium when
  the price-currency rate is higher (discount otherwise); forward points equal
  `F - S`; the implied rates round-trip; swap points are positive for a
  premium; a zero tenor returns spot; bad inputs raise.

## [1.350.0] - 2026-09-11

### Added
- `accrued_interest`, `dirty_price`, and `clean_price` (in `bondmath.py`): bond
  settlement mechanics. Accrued interest is the straight-line share of the
  current coupon (`face * coupon/freq * fraction_elapsed`); the dirty price is
  the full present value of the remaining cashflows; the clean (quoted) price
  is dirty minus accrued.
- Verified: accrued is half a coupon at mid-period and zero at a coupon date;
  the dirty price equals `bond_price_from_yield`; clean equals dirty minus
  accrued and coincides with dirty at a coupon date; bad fraction/freq raise.

## [1.349.0] - 2026-09-11

### Added
- `black_litterman_weights` (in `portopt.py`): optimal portfolio weights from
  the Black-Litterman posterior returns, `w = (lambda C)^{-1} mu`, optionally
  normalized to sum to 1. With no views the normalized (and raw) weights
  reproduce the market weights -- the prior is self-consistent.
- Verified: with no views the weights equal the market weights (raw and
  normalized); a bullish absolute view tilts weight toward that asset and a
  bearish one away; normalized weights sum to 1 for both absolute and relative
  views.

## [1.348.0] - 2026-09-11

### Added
- `implied_equilibrium_returns` and `black_litterman_returns` (in `portopt.py`):
  the reverse-optimized market prior `Pi = lambda C w` and the Black-Litterman
  posterior expected returns blending that prior with linear views `P mu = Q`
  (uncertainty `omega`, default `diag(tau P C P^T)`) via the standard closed
  form.
- Verified: the equilibrium equals `lambda C w`; with no views the posterior is
  the prior; a bullish absolute view raises that asset; a relative view moves
  the posterior spread between the prior and the view; a tighter view
  uncertainty pulls the posterior closer to the view; malformed views/weights
  raise.

## [1.347.0] - 2026-09-11

### Added
- `marginal_var` and `var_budget` (in `portopt.py`): the sensitivity of the
  portfolio VaR to each weight (`z sqrt(horizon) (C w)_i / sigma_p`) and the
  percentage VaR budget (`w_i (C w)_i / (w^T C w)`, summing to 1 and level-
  independent).
- Verified: the weight-dot-marginal Euler sum recovers the total VaR; the
  weight-times-marginal equals the component VaR up to the z-scaling; the
  budget sums to 1 and is equal across assets for a risk-parity portfolio;
  a zero-variance portfolio raises.

## [1.346.0] - 2026-09-11

### Added
- `portfolio_var`, `portfolio_cvar`, and `component_var` (in `portopt.py`):
  parametric (Gaussian) portfolio risk from the covariance matrix.
  `portfolio_var` is `z sigma_p sqrt(horizon) - mu_p horizon` (a positive loss),
  `portfolio_cvar` the expected shortfall `phi(z)/(1-c) sigma_p ...`, and
  `component_var` the per-asset risk contributions `w_i (C w)_i / sigma_p`
  (summing to the portfolio volatility). Includes an Acklam inverse-normal
  quantile.
- Verified: VaR equals `z * sigma_p`; CVaR exceeds VaR; a positive expected
  return lowers VaR; the component VaRs sum to the portfolio standard
  deviation; VaR scales with `sqrt(horizon)` and rises with confidence; a
  zero-variance portfolio raises.

## [1.345.0] - 2026-09-11

### Added
- `max_diversification_weights` and `diversification_ratio` (in `portopt.py`):
  the most-diversified portfolio (weights proportional to `C^{-1} sigma`,
  normalized) and the diversification ratio `(w^T sigma)/sqrt(w^T C w)`.
- Verified: the weights sum to 1 and maximize the diversification ratio versus
  the min-variance, risk-parity, and equal-weight portfolios; the ratio is at
  least 1, near 1 for highly-correlated assets, and exactly 1 for a single
  asset; a singular covariance raises.

## [1.344.0] - 2026-09-11

### Added
- `target_return_weights` and `efficient_frontier` (in `portopt.py`):
  minimum-variance weights that hit an exact expected return via the
  two-constraint Lagrangian (efficient-frontier scalars `A, B, C, D`), and a
  helper that sweeps target returns into `(return, portfolio_std)` frontier
  points.
- Verified: the weights hit the requested return exactly and sum to 1; at the
  global min-variance return they reproduce `min_variance_weights`; the frontier
  standard deviation is lowest at the min-variance return; requested returns are
  echoed with positive std; a target-return portfolio's variance is at least the
  global minimum; a length mismatch raises.

## [1.343.0] - 2026-09-11

### Added
- New module `portopt.py` for mean-variance portfolio optimization from a
  covariance matrix (pure-Python Gauss-Jordan inverse): `min_variance_weights`
  (`C^{-1}1` normalized), `max_sharpe_weights` (tangency `C^{-1}(mu-rf)`),
  `risk_parity_weights` (equal-risk-contribution by sqrt-damped fixed point),
  plus `portfolio_variance` and `portfolio_return` helpers.
- Verified: min-variance on a diagonal covariance gives inverse-variance
  weights and beats equal-weight variance; risk parity equalizes the per-asset
  risk contributions and reduces to inverse-sigma on a diagonal covariance;
  max-Sharpe weights sum to 1; the helpers match hand calculations; a singular
  covariance raises.

## [1.342.0] - 2026-09-11

### Added
- `up_capture`, `down_capture`, and `downside_beta` (in `perfmetrics.py`):
  benchmark-conditional statistics. The capture ratios are the asset's
  geometric per-period return over the benchmark's, taken over the up
  (`benchmark > 0`) or down (`benchmark < 0`) periods; the downside beta is
  `Cov/Var` restricted to negative-market periods.
- Verified: all three equal 1 for an asset that mirrors the market; a 1.5x
  leveraged asset gives ~1.5 captures and downside beta and a 0.5x defensive
  one gives captures below 1; the downside beta uses only down periods;
  mismatched lengths and empty up/down subsets raise.

## [1.341.0] - 2026-09-11

### Added
- `tracking_error` and `information_ratio` (in `perfmetrics.py`): benchmark-
  relative statistics from a portfolio and benchmark return series. Tracking
  error is the annualized standard deviation of the active (excess) returns;
  the information ratio is the annualized active return over that tracking
  error -- equivalently the Sharpe of the active-return series.
- Verified: tracking error matches the annualized active-return stdev; the
  information ratio equals `sharpe_ratio` of the active series and is positive
  for a positive active mean; a zero-active benchmark gives zero tracking error;
  zero active variance, length mismatch, and too-short series raise.

## [1.340.0] - 2026-09-11

### Added
- `omega_ratio` and `tail_ratio` (in `perfmetrics.py`): the Omega ratio
  (upside over downside area about a threshold) and the tail ratio (right-tail
  magnitude over left-tail, default 95th over 5th percentile).
- Verified: Omega matches the manual gain/loss-sum formula, exceeds 1 for a
  net-positive series, is `inf` without downside and 1 for a symmetric series;
  the tail ratio is 1 for a symmetric series and above 1 for a right-skewed
  one; degenerate inputs and out-of-range percentiles raise.

## [1.339.0] - 2026-09-11

### Added
- `drawdown_curve`, `longest_drawdown_duration`, and `rolling_sharpe` (in
  `perfmetrics.py`): the per-period underwater curve (fractional drop from the
  running peak, whose max is `max_drawdown`), the longest run of consecutive
  underwater periods, and the annualized Sharpe over each trailing window
  (`nan` for a zero-variance window so the series stays aligned).
- Verified: the underwater curve is non-negative, length-matched, its max
  equals `max_drawdown`, and a known series gives `[0, 0.5, 0.4]`; the drawdown
  duration counts the underwater run and resets on recovery; the rolling Sharpe
  has the right length and its first window matches `sharpe_ratio`, with a flat
  window returning `nan`; bad windows raise.

## [1.338.0] - 2026-09-11

### Added
- New module `perfmetrics.py` with track-record statistics from a periodic
  return series: `sharpe_ratio` and `sortino_ratio` (annualized),
  `max_drawdown` of the compounded equity curve, `calmar_ratio` (annualized
  return over max drawdown), `hit_rate`, and `profit_factor`.
- Verified: the Sharpe matches the manual mean/stdev formula; the Sortino
  exceeds the Sharpe on a mixed series; a known -50% peak-to-trough gives a 0.5
  max drawdown and an all-rising series gives 0; the hit rate and profit factor
  match hand counts (with an infinite profit factor when there are no losses);
  the Calmar is positive and raises without a drawdown; degenerate inputs raise.

## [1.337.0] - 2026-09-11

### Added
- `ewma_covariance`, `ewma_correlation`, and `realized_beta` (in
  `correlation.py`): return-series risk estimators. The EWMA covariance/
  correlation use the RiskMetrics decay recursion (default `lam = 0.94`),
  weighting recent observations more; `realized_beta` is the ordinary sample
  `Cov(asset, market) / Var(market)` regression slope.
- Verified: the realized beta recovers a known 1.5 slope and is 1 for the
  market against itself; the EWMA correlation stays in `[-1, 1]`, is 1 for a
  series with itself and -1 against its negation; EWMA self-covariance is
  positive; mismatched lengths, a bad decay, and too-short series raise.

## [1.336.0] - 2026-09-11

### Added
- `garch_option_price` (in `volatility.py`): a Black-Scholes price that uses the
  GARCH term (average) volatility for the maturity. It feeds
  `garch_term_variance` over `horizon` steps into the BSM formula, with the year
  fraction defaulting to `horizon / periods_per_year` (overridable via `t`). A
  fitted GARCH model then prices options consistently with its own vol
  mean-reversion rather than a flat spot vol.
- Verified: it matches BSM priced at the GARCH term vol for calls and puts;
  respects put-call parity; an elevated starting variance makes a short-dated
  option richer than one at the long-run vol; the explicit `t` override works;
  a missing horizon raises.

## [1.335.0] - 2026-09-11

### Added
- `garch_term_variance` (in `volatility.py`): the annualized GARCH(1,1) term
  (average) volatility over the next `horizon` steps -- the vol an option of
  that maturity is priced off, not a single-step forecast. Uses the closed-form
  geometric-series average of the mean-reverting per-step variances,
  `avg = LR + (h_1 - LR)/n (1 - p^n)/(1 - p)`.
- Verified: it equals the brute-force average of the per-step
  `garch_forecast` variances; the one-step case matches the point forecast; it
  converges to the long-run vol as the horizon grows and (from an elevated
  start) decreases monotonically toward it; a unit-root process falls back to
  the one-step variance; a zero horizon raises.

## [1.334.0] - 2026-09-11

### Added
- `cds_accrual_on_default` (in `credit.py`): the accrued-premium annuity a
  protection buyer owes on a mid-period default, integrated over each coupon
  interval. `cds_premium_leg` and `cds_par_spread` gain an
  `accrual_on_default` flag that folds this term into the premium annuity (the
  market-standard convention).
- Verified: the accrual factor is positive and a small fraction of the coupon
  annuity; including it lowers the par spread and raises the premium-leg PV; it
  is bounded by the half-period default probability and rises with the hazard.

## [1.333.0] - 2026-09-11

### Added
- `risky_bond_price` and `risky_bond_yield_spread` (in `credit.py`): price a
  defaultable coupon bond under a hazard-rate survival curve. Each cashflow is
  survival-weighted (`CF DF(t) Q(t)`) and a `recovery * face` payment is
  grid-integrated over the default time; `risky_bond_yield_spread` finds the
  flat credit spread over `r` that reproduces that price from the promised
  cashflows.
- Verified: a zero hazard recovers the risk-free bond price (and a zero yield
  spread); a positive hazard prices below risk-free and falls with higher
  hazard; higher recovery raises the price; the fitted spread is positive and
  reprices the bond when the promised cashflows are discounted at `r + spread`.

## [1.332.0] - 2026-09-11

### Added
- `cds_greeks` (in `credit.py`): finite-difference risk sensitivities of a CDS
  mark-to-market -- `credit01` (1bp parallel hazard bump), `ir01` (1bp discount
  rate), `recovery01` (1-point recovery rise), plus the `value` and
  `risky_annuity`. Respects the protection-buyer/seller sign.
- Verified: `value` matches `cds_value`; a protection buyer has positive
  `credit01` (gains on widening) and negative `recovery01`; the seller's
  sensitivities and value are the negatives of the buyer's; the risky annuity is
  positive; `credit01` scales roughly linearly with the bump size.

## [1.331.0] - 2026-09-11

### Added
- `bootstrap_survival_curve` (in `credit.py`): calibrate a piecewise-constant
  hazard `SurvivalCurve` from par CDS quotes. Solves each tenor's forward hazard
  in turn (holding earlier segments fixed) by bisection so the model par spread
  matches the quote, since the par spread is monotone in the current-segment
  hazard.
- Verified: the calibrated curve reprices every input quote to 1e-6; it returns
  a `SurvivalCurve` on the quote pillars with monotone-decreasing survival; an
  upward-sloping spread curve produces rising hazards; a flat spread curve gives
  a nearly constant hazard; mismatched or empty quotes raise.

## [1.330.0] - 2026-09-11

### Added
- New module `credit.py` with reduced-form credit pricing: a
  piecewise-constant-hazard `SurvivalCurve` (survival probability, forward
  hazard, default density) plus CDS analytics -- `risky_annuity`,
  `cds_protection_leg`, `cds_premium_leg`, `cds_par_spread`, and `cds_value`.
  The protection leg integrates the discounted loss over a grid; the par spread
  is protection PV over the risky annuity.
- Verified: survival starts at 1 and decreases, matching `e^{-h t}` for a flat
  hazard; the par spread is near the credit-triangle `h (1 - recovery)`; the CDS
  value is zero at the par spread; a protection buyer profits below par; a
  higher hazard widens the spread; the risky annuity is below the risk-free
  annuity; malformed curves raise.

## [1.329.0] - 2026-09-11

### Added
- `price_from_curve`, `key_rate_durations`, and `effective_duration_from_curve`
  (in `bondmath.py`): curve-based bond risk. `price_from_curve` discounts the
  cashflows off any `curve.df(t)` (or plain callable); `key_rate_durations`
  bumps each pillar zero rate in turn and returns the partial durations aligned
  with the pillars; `effective_duration_from_curve` central-differences a
  parallel shift of the whole zero curve.
- Verified: the sum of the key-rate durations equals the effective (parallel)
  duration to first order; on a flat curve that sum matches the single-yield
  Macaulay duration; a bullet bond's largest key-rate duration is at the final
  pillar and all are non-negative; `price_from_curve` matches a manual discount
  and accepts a plain callable; length-mismatched pillars/rates raise.

## [1.328.0] - 2026-09-11

### Added
- New module `bondmath.py` with single-yield coupon-bond analytics:
  `bond_cashflows` (level-coupon schedule), `bond_price_from_yield`,
  `macaulay_duration`, `modified_duration`, `convexity`, `bond_dv01`, and
  `yield_to_maturity`. All work off an explicit `(time, amount)` schedule and a
  continuously-compounded yield (so modified and Macaulay duration coincide);
  YTM is solved by Newton with a bisection fallback.
- Verified: YTM round-trips the yield from a price; duration and convexity match
  central finite differences of the price/yield curve; duration lies below
  maturity and convexity is positive; DV01 equals `-modified_duration * price *
  1e-4`; a zero-coupon bond's duration equals its maturity; bad price/maturity
  inputs raise.

## [1.327.0] - 2026-09-11

### Added
- `holee_expected_rate` and `holee_rate_variance` (in `holee.py`): analytic
  moments of the Ho-Lee short rate. With constant drift the rate is a drifted
  Brownian motion `r_t = r0 + theta t + sigma W_t`, so the mean is `r0 + theta t`
  and the variance `sigma^2 t` -- both linear in time, with no stationary law
  (no mean reversion). Completes the short-rate moment set alongside Vasicek and
  CIR.
- Verified: the mean and variance match their formulas and an Euler Monte Carlo
  of the SDE; the variance grows without bound; negative time raises.

## [1.326.0] - 2026-09-11

### Added
- `cir_expected_rate`, `cir_rate_variance`, and `cir_stationary_distribution`
  (in `cir.py`): analytic moments of the Cox-Ingersoll-Ross square-root rate.
  The mean `theta + (r0 - theta) e^{-kappa t}` matches Vasicek, but the variance
  `r0 (sigma^2/kappa)(e^{-kappa t} - e^{-2 kappa t}) + theta (sigma^2/(2 kappa))
  (1 - e^{-kappa t})^2` is state-dependent. The stationary law is Gamma with
  shape `2 kappa theta / sigma^2` and scale `sigma^2/(2 kappa)` (Feller
  condition `shape >= 1`).
- Verified: mean and variance match a full-truncation Euler Monte Carlo of the
  CIR SDE; the variance grows to the stationary value and rises with the
  starting rate `r0`; the stationary Gamma's shape*scale and shape*scale^2
  reproduce mean `theta` and variance; `kappa <= 0` or `sigma <= 0` raise.

## [1.325.0] - 2026-09-11

### Added
- `vasicek_expected_rate`, `vasicek_rate_variance`, and
  `vasicek_stationary_distribution` (in `vasicek.py`): analytic moments of the
  Ornstein-Uhlenbeck short rate. Mean `theta + (r0 - theta) e^{-kappa t}`,
  variance `sigma^2/(2 kappa) (1 - e^{-2 kappa t})`, and the long-run Normal law
  `(theta, sigma^2/(2 kappa))`.
- Verified: the mean matches its formula and an Euler Monte Carlo of the OU
  process; the variance grows from 0 to the stationary value; the `kappa -> 0`
  limits give `r0` and `sigma^2 t`; the stationary mean is `theta` and requires
  `kappa > 0`.

## [1.324.0] - 2026-09-11

### Added
- `box_spread`, `synthetic_forward`, and `collar` (in `strategy.py`): three more
  strategy builders. The box (bull call + bear put on the same strikes) is a
  synthetic zero-coupon bond worth `e^{-rt} (K_high - K_low)` with a constant
  payoff and flat Greeks; the synthetic forward (long call, short put at one
  strike) replicates a forward with delta ~1 by put-call parity; the collar
  (long put floor, short call cap) prices its two option legs.
- Verified: the box value is the discounted strike width, its payoff is
  constant across spot, and its delta/gamma/vega are ~zero; the synthetic
  forward matches `S - e^{-rt} K` with delta ~1; the collar matches its legs;
  reversed strikes raise.

## [1.323.0] - 2026-09-11

### Added
- `calendar_spread` and `diagonal_spread` (in `strategy.py`): multi-expiry
  option-strategy builders. A calendar shorts the near expiry and longs the far
  at the same strike; a diagonal does the same with different strikes. Both use
  the per-contract maturity, so `price_book` gives the net debit and net Greeks
  directly. Require `t_near < t_far`.
- Verified: the net value equals the far-leg price minus the near-leg price
  (call and put); a long calendar is a net debit and long vega; the diagonal net
  matches its two legs; a reversed expiry order raises.

## [1.322.0] - 2026-09-11

### Added
- `discrete_fixed_strike_lookback_greeks` (in `lookback.py`): finite-difference
  delta, gamma, vega, theta of the discretely-monitored fixed-strike lookback
  (`discrete_fixed_strike_lookback`), holding the monitoring count `n_fixings`
  fixed.
- Verified: delta matches a central finite difference for calls and puts; the
  `price` field matches the pricer; the call delta is positive, put delta
  negative, and gamma/vega positive.

## [1.321.0] - 2026-09-11

### Added
- `implied_exchange_correlation` (in `multiasset.py`): back out the correlation
  implied by a Margrabe exchange-option price. The price depends on `rho` only
  through the spread vol `sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)`, which
  falls as `rho` rises, so the price is monotone decreasing in `rho` and a
  bisection on `(-1, 1)` recovers it. Complements the spread and
  geometric-basket implied-correlation solvers.
- Verified: round-trips the correlation across `rho in {-0.5, 0, 0.3, 0.7}`
  (with and without dividends); the price is monotone decreasing in `rho`; a
  quote outside the `rho`-range raises.

## [1.320.0] - 2026-09-11

### Added
- `partial_time_start_barrier_call` (in `exotics.py`): partial-time (start)
  single-barrier call (Heynen-Kat 1994), where the down barrier is monitored
  only over `[0, t1]` and inactive afterwards, paying off at `T2 > t1`. Exact
  bivariate-normal closed form coupling the monitoring-end date to expiry
  (`rho = sqrt(t1/T2)`); the reflection terms use `(+e3, +rho)` rather than the
  `(-e3, -rho)` of the end-barrier variant. Down-out priced directly, down-in by
  in-out parity. (Up-barrier partial-time calls raise `ValueError`.)
- Verified: as `t1 -> 0` it approaches the vanilla call and as `t1 -> T2` the
  standard continuously-monitored down-out; a mid window sits strictly between;
  the price is monotone decreasing in `t1`; in-out parity holds exactly; matches
  a path Monte Carlo across `t1 in {0.25, 0.5, 0.75}` to within ~2%.

## [1.319.0] - 2026-09-11

### Added
- `caplet_implied_normal_vol` (in `rates.py`): invert a caplet/floorlet price
  back to its normal (Bachelier) volatility. Divides out the
  `discount * accrual` factor to recover the undiscounted option value, then
  inverts with `bachelier_implied_vol` — the inverse of `caplet_price`.
- Verified: round-trips the input vol to 1e-8 across strikes, vols, and
  cap/floor; works with a negative forward (normal model); rejects an expired
  caplet.

## [1.318.0] - 2026-09-11

### Fixed
- `sabr_vol` (lognormal Hagan) no longer divides by zero when the vol-of-vol
  `nu = 0` (which makes `z = 0` and `x(z) = 0`). The `z / x(z)` ratio now takes
  its limit of 1, matching the `nu -> 0` behaviour — the same fix already
  applied to `sabr_normal_vol` in 1.317.0. Nonzero-`nu` values are unchanged.
- Verified: the `nu = 0` vol matches the `nu -> 0` limit; the smile stays
  skewed by `beta < 1` when `nu = 0`; the ATM `nu = 0` vol sits at
  `alpha / F^{1-beta}` up to the small time correction.

## [1.317.0] - 2026-09-11

### Added
- `sabr_cap_price` and `sabr_floor_price` (in `rates.py`): price an
  interest-rate cap/floor under a single SABR smile (normal model). Each caplet
  is valued at the SABR normal vol read at its own forward and expiry, so one
  calibrated `(alpha, beta, rho, nu)` prices the whole strip consistently across
  the smile rather than a flat per-period `sigma_n`.

### Fixed
- `sabr_normal_vol` no longer divides by zero when the vol-of-vol `nu = 0`. The
  leading factor is rewritten as `alpha (F - K)/((FK)^{beta/2} log(F/K)) *
  z/x(z)` with the `z -> 0` limit `z/x(z) = 1`, so the degenerate case returns
  the finite lognormal-spacing vol (nonzero-`nu` values are unchanged).

### Verified
- SABR cap equals the sum of caplets each priced at its own SABR vol;
  cap-floor parity equals `sum disc * accrual * (F - K)`; a `beta = nu = 0` SABR
  cap is within ~2% of a flat normal cap at `sigma_n = alpha`; prices are
  positive; `nu = 0` no longer crashes.

## [1.316.0] - 2026-09-11

### Added
- `sabr_swaption_price` (in `rates.py`): price a European swaption whose smile
  is a SABR model. Reads the SABR-implied vol at the (forward swap rate, strike,
  expiry) point and feeds it into the matching pricer — `sabr_vol` into
  `black_swaption_price` for `model="black"`, or `sabr_normal_vol` into
  `swaption_price` for `model="normal"`. One calibrated smile prices every
  strike consistently.
- Verified: the resulting price inverts back to exactly the SABR vol it was
  built from (both Black and normal, across strikes); prices are positive; an
  in-the-money payer exceeds the receiver; an unknown model raises.

## [1.315.0] - 2026-09-11

### Added
- `swaption_implied_normal_vol` and `swaption_implied_black_vol` (in
  `rates.py`): invert a swaption price back to its normal (Bachelier) or Black
  (lognormal) volatility. Each divides out the annuity to recover the
  undiscounted option value, then inverts with `bachelier_implied_vol` /
  `implied_volatility` at zero carry — the inverses of `swaption_price` and
  `black_swaption_price`.
- Verified: both round-trip the input volatility to 1e-6 across payer/receiver
  and multiple vols; the Black inversion rejects non-positive rates; both reject
  a zero/expired option.

## [1.314.0] - 2026-09-11

### Added
- `black_swaption_price` and `black_swaption_greeks` (in `rates.py`): the Black
  (lognormal) European swaption, the market-standard counterpart to the normal
  `swaption_price`. The forward swap rate is lognormal with Black vol
  `sigma_b`, and the swaption is `annuity * Black76(swap_rate, strike, expiry,
  sigma_b)`. Payer = call on the rate, receiver = put; requires positive rate
  and strike. Greeks are the Black-76 rate delta/gamma/vega scaled by the
  annuity.
- Verified: payer-minus-receiver equals `annuity * (swap_rate - strike)`
  exactly; ATM Black matches the Bachelier swaption with
  `sigma_n = sigma_b * swap_rate`; rate delta matches a finite difference; gamma
  and vega are positive; payer delta positive and receiver negative; zero expiry
  returns the intrinsic; negative rate/strike raise.

## [1.313.0] - 2026-09-11

### Added
- `DiscountCurve.forward_swap_rate` and `DiscountCurve.forward_annuity` (in
  `discount_curve.py`): the par rate and PV01 of a forward-starting swap that
  begins accruing at a future date. The forward float leg is
  `DF(start) - DF(T_n)`, so `fwd = (DF(start) - DF(T_n)) / sum_i tau_i DF(T_i)`
  over the forward schedule — the underlying rate a swaption is written on.
- Verified: with `start = 0` the forward swap rate reduces to
  `par_swap_rate`; at the forward par rate the fixed leg exactly balances the
  forward float leg; on an upward curve the forward rate exceeds the spot par
  rate; a pay date at or before the start is rejected.

## [1.312.0] - 2026-09-11

### Added
- `sabr_normal_vol` (in `sabr.py`): the Hagan (2002) normal (Bachelier) implied
  volatility for the SABR model — the absolute-vol `sigma_N` quoted in
  interest-rate markets. Uses the Hagan normal expansion (leading factor
  `nu (F - K) / x(z)`, third-order time bracket shared with the lognormal
  `sabr_vol`), with the ATM `F == K` limit handled separately.
- Verified: matches the reference normal vol (SABR-Black vol -> Black-Scholes
  price -> Bachelier implied vol) to within 5e-3 across strikes; the ATM branch
  joins continuously with the just-off-ATM value; `beta = nu = 0` gives exactly
  `alpha`; the normal-vol smile is positive and rises in both wings.

## [1.311.0] - 2026-09-11

### Added
- `displaced_diffusion_smile` (in `displaced.py`): the Black-Scholes implied-vol
  smile a displaced-diffusion (shifted-lognormal) model produces. Prices a call
  at each strike and inverts to BS implied vol, returning `(log_moneyness, vol)`
  on the forward `F = S e^{b t}`. Mirrors `cev_smile` / `merton_smile`.
- Verified: `shift = 0` gives a flat smile at `sigma`; the ATM vol sits near
  `sigma`; a positive shift produces a downward skew, steeper for larger shift;
  output is sorted by strike; each smile vol reprices to the displaced-diffusion
  price it came from.

## [1.310.0] - 2026-09-11

### Added
- `cev_smile` (in `cev.py`): the Black-Scholes implied-vol smile a CEV model
  produces. Prices a call at each strike under CEV and inverts to BS implied
  vol, returning `(log_moneyness, vol)` pairs on the forward
  `F = S e^{(r-q) t}`. Mirrors `merton_smile` / `kou_smile`.
- Verified: the ATM implied vol sits near `sigma` (calibrated to the ATM
  instantaneous vol); the smile has a downward skew for `beta < 1`, steeper for
  smaller `beta`; output is sorted by strike; each smile vol reprices to the CEV
  price it came from.

## [1.309.0] - 2026-09-11

### Added
- `bjerksund_stensland_boundary` (in `american.py`): the Bjerksund-Stensland
  (2002) flat exercise trigger `I` at inception — a call is exercised for
  `S >= I`, a put for `S <= I`. Returns `None` when early exercise is never
  optimal (an American call with `b >= r`). The put trigger follows from the
  same put-call transformation the pricer uses: the transformed call's
  spot-axis trigger `I2t` maps back to `K^2 / I2t`.

### Changed
- Extracted the BS2002 trigger computation into a shared `_bs2002_triggers`
  helper used by both the call pricer and the new boundary function (prices
  unchanged).
- Verified: the call boundary lies above the strike and the put below; the
  BS2002 price equals the exercise intrinsic exactly at the boundary; a spot
  just inside the boundary has continuation value strictly above intrinsic; the
  boundary sits within ~10% of the Barone-Adesi-Whaley one (both flat
  approximations); the refactor leaves the reference prices unchanged.

## [1.308.0] - 2026-09-11

### Added
- `baw_critical_spot` (in `baw.py`): the Barone-Adesi-Whaley early-exercise
  boundary `S*` at inception — a call is exercised for `S >= S*`, a put for
  `S <= S*`. Returns `None` when early exercise is never optimal (an American
  call with `b >= r`).
- `baw_american_greeks` (in `baw.py`): delta, gamma, vega, theta of the BAW
  American price by finite difference.

### Changed
- Refactored the BAW critical-spot Newton solves into reusable `_critical_call`
  and `_critical_put` helpers shared by the pricer and the new boundary
  function (prices unchanged).
- Verified: the call boundary lies above the strike and the put boundary below;
  the American price equals the exercise intrinsic exactly at `S*`; greeks match
  central finite differences for calls and puts across several spots; the put
  delta lies in `(-1, 0)` with positive gamma; the BAW value dominates the
  European price.

## [1.307.0] - 2026-09-11

### Added
- `bachelier_theta` (in `bachelier.py`): analytic calendar theta of the normal
  (Bachelier) model, `theta = r * price - e^{-rt} sigma phi(d) / (2 sqrt(t))`.
  `bachelier_greeks` now uses it instead of a finite difference.
- `bachelier_cash_or_nothing` and `bachelier_asset_or_nothing` (in
  `bachelier.py`): normal-model digitals. `F_T` is Gaussian, so a cash-or-nothing
  call is `cash e^{-rt} N(d)` and an asset-or-nothing call is
  `e^{-rt} (F N(d) + sigma sqrt(t) phi(d))`, with `d = (F - K)/(sigma sqrt t)`.
- Verified: analytic theta matches central finite differences at, above, and
  below the strike for calls and puts; cash-digital call/put parity sums to the
  discount factor; the vanilla Bachelier price decomposes into asset-or-nothing
  minus `K` cash-or-nothing; both digitals match a Gaussian Monte Carlo.

## [1.306.0] - 2026-09-11

### Added
- `partial_time_end_barrier_call` (in `exotics.py`): partial-time (end)
  single-barrier call (Heynen-Kat 1994), where the down barrier is monitored
  only over `[t1, T2]` (inactive before `t1`). Exact bivariate-normal closed
  form coupling the monitoring-start date to expiry (`rho = sqrt(t1/T2)`).
  Down-out is priced directly; down-in follows from in-out parity. (Up-barrier
  partial-time calls have a distinct form and raise `ValueError`.)
- Verified: as `t1 -> 0` it approaches the standard continuously-monitored
  down-out barrier, and as `t1 -> T2` it approaches the vanilla call; a mid
  window sits strictly between the two; in-out parity holds exactly; matches a
  path Monte Carlo across `t1 in {0.25, 0.5, 0.75}` to within ~2%.

## [1.305.0] - 2026-09-11

### Added
- `discrete_barrier_option` (in `exotics.py`): discretely-monitored
  single-barrier option via the Broadie-Glasserman-Kou (1999) continuity
  correction. A barrier checked at `n_fixings` equally-spaced dates is breached
  less often than a continuous one, so the correction shifts the barrier away
  from the spot by `exp(+/- beta sigma sqrt(dt))` (up for up-barriers, down for
  down-barriers, `beta ~ 0.5826`) and prices with the continuous
  `barrier_option`. Applies to all four knock in/out types.
- Verified: matches a path Monte Carlo at `n = 50` across all four barrier
  types (down/up x in/out) to within ~3%; a discrete knock-out sits above and a
  knock-in below the continuous price; in-out parity holds exactly at the
  shifted barrier (KI + KO = vanilla); converges to the continuous barrier on a
  fine grid.

## [1.304.0] - 2026-09-11

### Added
- `discrete_fixed_strike_lookback` (in `lookback.py`): discretely-monitored
  fixed-strike lookback via the Broadie-Glasserman-Kou (1999) continuity
  correction. The realized extreme is sampled at `n_fixings` equally-spaced
  dates; the correction shifts the spot fed to the continuous
  `fixed_strike_lookback` by `exp(-/+ beta sigma sqrt(dt))` (down for a call on
  the max, up for a put on the min, `beta ~ 0.5826`). As `n_fixings` grows the
  shift vanishes and the price converges to the continuous lookback.
- Verified: matches a path Monte Carlo at `n = 50` (call/put) to within ~1.5%;
  the discrete call sits below the continuous lookback; the price converges
  (monotonically, and the correction decays like `1/sqrt(n)`) to the continuous
  value on a fine grid.

## [1.303.0] - 2026-09-11

### Added
- `holder_extendible_put` and `holder_extendible_put_greeks` (in
  `extendible.py`): holder-extendible put (Longstaff 1990), completing the full
  2x2 extendible family (holder/writer x call/put). At the first expiry `t1` the
  holder takes the best of exercising against `K1`, lapsing, or paying a fee `A`
  to extend to `T2` as a put struck at `K2`. The terminal spot splits into
  exercise (`S < I_low`), extend (`I_low <= S <= I_high`), and lapse
  (`S > I_high`); the extended-put strip is built from the call-strip pieces via
  `N(-d) = 1 - N(d)`. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity/fee cases to within ~0.03%; a very large fee collapses the
  strip and recovers the vanilla put to `t1`; a finite fee adds value; a lower
  fee is worth more; delta is negative and vega positive.

## [1.302.0] - 2026-09-11

### Added
- `writer_extendible_call` and `writer_extendible_call_greeks` (in
  `extendible.py`): writer-extendible call (Longstaff 1990), completing the
  extendible family. At the first expiry `t1` the call is exercised if in the
  money (`S_{t1} > K1`); otherwise the writer's obligation extends automatically
  to `T2` as a call struck at `K2` (no fee). Closed form: a vanilla call to `t1`
  plus the extended-call value collected on `S_{t1} <= K1`, via bivariate
  normals coupling `t1` and `T2`. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity cases to within ~0.1%; the automatic extension makes it worth
  more than a plain call to `t1`; delta is positive and vega positive.

## [1.301.0] - 2026-09-11

### Added
- `writer_extendible_put` and `writer_extendible_put_greeks` (in
  `extendible.py`): writer-extendible put (Longstaff 1990). At the first expiry
  `t1` the put is exercised if in the money (`S_{t1} < K1`); otherwise the
  writer's obligation extends automatically to `T2` as a put struck at `K2`
  (no fee). Closed form: a vanilla put to `t1` plus the extended-put value
  collected on `S_{t1} >= K1`, via bivariate normals coupling `t1` and `T2`.
  Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity cases to within ~0.1%; the automatic extension makes it worth
  more than a plain put to `t1`; delta is negative and vega positive.

## [1.300.0] - 2026-09-11

### Added
- `holder_extendible_call` and `holder_extendible_call_greeks` (new module
  `extendible.py`): holder-extendible call (Longstaff 1990). At the first expiry
  `t1` the holder takes the best of exercising against `K1`, lapsing, or paying
  a fee `A` to extend the life to `T2` with strike `K2`. Priced in closed form
  by splitting the terminal spot into three regions (lapse / extend / exercise)
  and valuing the extension strip with bivariate normals coupling `t1` and `T2`;
  the strip boundaries are found by bisection. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity/fee cases to within ~0.15%; a very large fee collapses the
  extension strip and recovers the vanilla call to `t1`; a finite fee adds value
  over that plain call; a lower fee is worth more.

## [1.299.0] - 2026-09-11

### Added
- `complex_chooser_option` and `complex_chooser_option_greeks` (in
  `chooser.py`): complex chooser (Rubinstein 1991), where at the choice date the
  holder keeps whichever is worth more of a call (strike `Kc`, expiry `Tc`) or a
  put (strike `Kp`, expiry `Tp`) — the two legs may differ in both strike and
  maturity. Priced by Rubinstein's bivariate-normal formula, with the critical
  spot (where the two legs are equal at the choice date) found by bisection.
  Greeks by finite difference.
- Verified: equal strikes and maturities reduce exactly to the simple
  `chooser_option`; matches a Monte Carlo at the choice date to within ~0.14%;
  the chooser is worth more than either leg valued outright today; delta matches
  a finite difference.

## [1.298.0] - 2026-09-11

### Changed
- `supershare_greeks` now computes `delta` and `gamma` analytically instead of
  by finite difference. With price `A S (N(d1_lo) - N(d1_hi))`, `A = e^{(b-r)t}/K_low`,
  the spot derivatives are
  `delta = A (dN + dphi/(sigma sqrt t))` and
  `gamma = A/(S sigma sqrt t) (dphi - (d1_lo phi(d1_lo) - d1_hi phi(d1_hi))/(sigma sqrt t))`,
  where `dN`, `dphi` are the differences of the two `N`/`phi` at the corridor
  edges. `vega` and `theta` remain finite differences.
- Verified: analytic delta and gamma match central finite differences across
  four strike/vol/maturity cases.

## [1.297.0] - 2026-09-11

### Changed
- `gap_option_greeks` now computes `delta` and `gamma` analytically instead of
  by finite difference. Using the identity
  `S carry phi(d1) = K_trigger disc phi(d2)`, the spot sensitivities collapse to
  closed forms in `d1, d2` at the trigger, with a cash-driven correction term
  proportional to the strike gap `K_trigger - K_payoff`. `vega` and `theta`
  remain finite differences.
- Verified: analytic delta and gamma match central finite differences across
  three strike/vol/maturity cases for both calls and puts; gamma is identical
  for calls and puts; equal trigger/payoff strikes recover the vanilla
  Black-Scholes delta and gamma.

## [1.296.0] - 2026-09-11

### Added
- `range_binary_greeks` (in `exotics.py`): analytic `delta` and `gamma` for the
  range binary plus finite-difference `vega`/`theta`. Delta is
  `cash e^{-rt} (phi(d2_lo) - phi(d2_hi))/(S sigma sqrt t)`; gamma carries the
  double-sided pin risk near either corridor edge.
- `supershare_greeks` (in `exotics.py`): finite-difference `delta`, `gamma`,
  `vega`, `theta` for the supershare option.
- Verified: analytic range-binary delta and gamma match central finite
  differences across three strike/vol/maturity cases; both `price` fields match
  their pricers; a forward-centered corridor has positive calendar theta (the
  terminal mass concentrates inside as expiry nears).

## [1.295.0] - 2026-09-11

### Added
- `range_binary` (in `exotics.py`): range binary / double digital paying `cash`
  iff `K_low <= S_T <= K_high` at expiry. Exactly the difference of two
  cash-or-nothing calls: `cash e^{-rt} (N(d2(K_low)) - N(d2(K_high)))`.
- `supershare` (in `exotics.py`): Hakansson (1976) supershare paying
  `S_T / K_low` iff `K_low <= S_T <= K_high`. A scaled difference of two
  asset-or-nothing calls.
- Verified: each equals the corresponding digital difference to machine
  precision; both match a terminal Monte Carlo to within ~0.4%; a wide corridor
  range binary approaches the full discounted cash; the zero-vol limit pays iff
  the forward lands inside the corridor.

## [1.294.0] - 2026-09-11

### Added
- `implied_geometric_basket_correlation` (in `multiasset.py`): backs out the
  uniform pairwise correlation implied by a geometric-basket quote, assuming an
  equicorrelation matrix (`corr[i][j] = rho` off-diagonal). The basket
  log-variance rises with `rho`, so both call and put prices increase in `rho`
  (parity `C - P = disc(F - K)` is strike-independent), and a bisection over
  `[-1/(n-1), 1]` — the PSD-preserving range — recovers it.
- Verified: round-trips the correlation for calls and puts across
  `rho in {-0.3, 0, 0.4, 0.8}` and for a two-asset basket; raises when the quote
  lies outside the `rho`-range; requires at least two assets.

## [1.293.0] - 2026-09-11

### Added
- `geometric_basket_option` and `geometric_basket_greeks` (in `multiasset.py`):
  weighted geometric-average basket option on `prod_i S_i^{w_i}` for any number
  of assets. Unlike the arithmetic `basket_option` (a Levy moment-match), the
  weighted geometric average of correlated lognormals is itself lognormal, so
  this is exact: `log B` is Gaussian with variance
  `t sum_ij w_i w_j corr[i][j] sigma_i sigma_j`, priced by a Black-Scholes
  formula on the basket forward. Greeks return per-asset delta/gamma lists plus
  total vega and theta.
- Verified: a single unit-weight asset recovers the vanilla Black-Scholes
  price; put-call parity equals the discounted basket forward minus strike;
  matches a Cholesky Monte Carlo (3 correlated assets) to within ~0.2%;
  per-asset deltas and vega are positive.

## [1.292.0] - 2026-09-11

### Added
- `average_strike_arithmetic_asian` and `average_strike_arithmetic_asian_greeks`
  (in `exotics.py`): average-strike (floating-strike) discrete arithmetic Asian,
  payoff `max(S_T - A, 0)`. Priced as an exchange option between `S_T` (exactly
  lognormal) and the arithmetic average `A` (Levy two-moment lognormal). The
  cross-moment `E[S_T A] = (S^2/n) sum_i exp(b(t+t_i) + sigma^2 t_i)` is exact,
  giving the log-space covariance and the Margrabe spread variance. Greeks by
  finite difference.
- Verified: put-call parity equals the discounted exchange forward
  `E[S_T] - E[A]` exactly; the average-strike arithmetic call is cheaper than
  the geometric one (AM-GM raises the strike); matches an antithetic Monte
  Carlo (n=12) to ~1% (a small stable Levy approximation bias); explicit times
  match the `n_fixings` grid.

### Changed
- `average_strike_geometric_asian_greeks` test now asserts the exact
  homogeneity property (gamma zero, delta = price/S) instead of a fragile
  gamma-sign check, since the average-strike payoff is homogeneous of degree 1
  in spot.

## [1.291.0] - 2026-09-11

### Added
- `average_strike_geometric_asian` and `average_strike_geometric_asian_greeks`
  (in `exotics.py`): average-strike (floating-strike) discrete geometric Asian.
  The strike is the realized geometric average, so a call pays
  `max(S_T - G, 0)`. `S_T` and `G` are jointly lognormal, making this an
  exchange option with an exact Margrabe-style closed form on the spread
  variance `sigma^2 t + v_G - 2 sigma^2 mean(t_i)`. Greeks by finite
  difference.
- Verified: put-call parity equals the discounted exchange forward
  `E[S_T] - E[G]` exactly; matches an antithetic Monte Carlo (n=12, call/put)
  to within ~0.1%; explicit equally-spaced times match the `n_fixings` grid;
  gamma and vega are positive.

## [1.290.0] - 2026-09-11

### Added
- `seasoned_geometric_asian_greeks` and `seasoned_arithmetic_asian_greeks` (in
  `exotics.py`): finite-difference Greeks (`delta`, `gamma`, `vega`, `theta`)
  for the seasoned (in-progress) Asian pricers. The calendar bump scales the
  remaining fixing schedule with `t` so a `theta` difference stays within the
  option's life and preserves the fixing shape.
- Verified: with no observations the Greeks match the fresh
  `discrete_geometric_asian_greeks` / `discrete_arithmetic_asian_greeks`; the
  `price` field matches the pricer; partial-window delta/gamma/vega are
  positive and theta negative for a call; a partially-seasoned option has
  smaller delta and vega than a fresh one (less of the average is still
  stochastic).

## [1.289.0] - 2026-09-11

### Added
- `seasoned_arithmetic_asian` (in `exotics.py`): prices an in-progress discrete
  arithmetic-average Asian where some fixings are already observed. Writing the
  average as `A = (Q + sum remaining)/n` with `Q = sum(observed)` known, a call
  payoff is an arithmetic-average option on the remaining fixings with the
  shifted strike `K' = nK - Q`, scaled by `1/n`. The remaining sum's two exact
  moments are matched to a lognormal (Levy 1992). Deep-in/out branches (`K' <=
  0`) are priced exactly as `disc*(E[A] - K)` / zero.
- Verified: with no observations it reduces to `discrete_arithmetic_asian`;
  all-observed gives the deterministic arithmetic intrinsic; the `K' <= 0`
  branch equals the discounted expected-average payoff; a partial window
  matches a path Monte Carlo (n=12, 6 observed) to within 1%; the seasoned
  arithmetic price sits above the seasoned geometric one (AM-GM).

## [1.288.0] - 2026-09-11

### Added
- `seasoned_geometric_asian` (in `exotics.py`): prices an in-progress discrete
  geometric-average Asian where some fixings have already been observed. The
  observed prices contribute a known constant `A = sum log(S_obs)` and the `k`
  remaining log-prices stay jointly Gaussian, so `log G` is still
  `Normal(m, v)` with `m = (A + k log S + (b - sigma^2/2) sum tau_j)/n` and
  `v = (sigma^2/n^2) sum_ij min(tau_i, tau_j)` — an exact Black-Scholes-style
  closed form on the remaining window.
- Verified: with no observations it reduces exactly to
  `discrete_geometric_asian`; with all fixings observed the payoff is the
  deterministic discounted geometric intrinsic; a partial-window case matches a
  path Monte Carlo (n=12, 6 observed) to within 1%; locking in high observed
  fixings raises the call versus at-the-money observations.

## [1.287.0] - 2026-09-11

### Added
- `discrete_arithmetic_asian` and `discrete_arithmetic_asian_greeks` (in
  `exotics.py`): discretely-monitored arithmetic-average Asian by Levy (1992)
  two-moment lognormal matching. The average's first two moments over the
  fixing dates are exact — `M1 = (S/n) sum_i exp(b t_i)` and
  `M2 = (S^2/n^2) sum_ij exp(b(t_i+t_j) + sigma^2 min(t_i,t_j))` — and are
  matched to a lognormal priced with a Black-Scholes formula on `M1` with
  variance `log(M2/M1^2)`. Accepts `n_fixings` or explicit `fixing_times`.
  Greeks by finite difference.
- Verified: a single fixing recovers the vanilla Black-Scholes price; the
  arithmetic Asian sits above the exact discrete geometric Asian (AM-GM) and
  below the vanilla call; explicit equally-spaced times match the `n_fixings`
  grid; matches a geometric-control-variate Monte Carlo (n=12, call/put) to
  within 1%.

## [1.286.0] - 2026-09-11

### Added
- `discrete_geometric_asian` and `discrete_geometric_asian_greeks` (in
  `exotics.py`): exact closed form for a discretely-monitored geometric-average
  Asian option. `log G` over the fixing dates is Gaussian with mean
  `log S + (b - sigma^2/2) mean(t_i)` and variance
  `(sigma^2/n^2) sum_ij min(t_i, t_j)`, giving a Black-Scholes-style price on
  the lognormal average. Accepts equally-spaced `n_fixings` or an explicit
  `fixing_times` sequence. Greeks by finite difference.
- Verified: a single fixing recovers the vanilla Black-Scholes price; the price
  converges (monotonically, from above) to the continuous Kemna-Vorst
  `geometric_asian` as `n_fixings` grows; explicit equally-spaced fixing times
  match the `n_fixings` grid; a discrete geo Asian call is cheaper than the
  vanilla call; matches a path Monte Carlo (n=12, call/put) to within 1%.

## [1.285.0] - 2026-09-11

### Added
- `powered_option` and `powered_option_greeks` (in `exotics.py`): the powered
  option, payoff `max(S_T - K, 0)**p` (call) or `max(K - S_T, 0)**p` (put) for
  a positive integer power `p`. Distinct from `power_option` (payoff
  `max(S_T**p - K, 0)`): here the option payoff itself is raised to a power.
  Closed form by binomial expansion of the polynomial payoff into truncated
  moments `E[S_T**j] N(d_j)` (Esser 2003; Heynen-Kat 1996). Greeks by finite
  difference.
- Verified: `p = 1` recovers the vanilla Black-Scholes price and Greeks;
  matches an antithetic Monte Carlo for `p = 2` (call/put) and `p = 3` (call)
  to within statistical error; zero-vol limit equals the discounted powered
  intrinsic; delta/gamma/vega positive for the call.

## [1.284.0] - 2026-09-11

### Added
- `ultima` (in `greeks2.py`): the third-order vega Greek,
  `d(vomma)/d(sigma) = d^3(price)/d(sigma)^3`, useful for the convexity of a
  volga hedge. Analytic: `ultima = -(vega/sigma^2) * [d1 d2 (1 - d1 d2) + d1^2 + d2^2]`.
  Same for calls and puts.
- Verified: matches a central finite difference of `vomma` in sigma to
  6-7 digits; matches a central third difference of the BSM price in sigma;
  call equals put; ATM ultima is negative (volga concave in sigma there).

## [1.283.0] - 2026-09-11

### Added
- `double_knock_in_call_greeks` (in `exotics.py`): Greeks of a double-barrier
  knock-in call by in-out parity. Spot/vol/time Greeks are the vanilla
  Black-Scholes Greek minus the double knock-out Greek; barrier sensitivities
  are the negatives of the knock-out's (the vanilla has no barrier dependence).
- Verified: `knock_in + knock_out` Greeks sum to the vanilla Greek; `dV/dL` and
  `dV/dU` negate the knock-out's; a wider corridor lowers the knock-in
  (`dV/dL > 0`, `dV/dU < 0`); the price field matches the pricer.

## [1.282.0] - 2026-09-11

### Added
- `double_knock_out_call_greeks` (in `exotics.py`): finite-difference Greeks of
  the Ikeda-Kunitomo double knock-out call — `delta`, `gamma`, `vega`, `theta`,
  plus the two barrier sensitivities `dV/dL` and `dV/dU`.
- Verified: `delta` matches an independent bump; `vega < 0` (short vol);
  widening the corridor raises value (`dV/dL < 0`, `dV/dU > 0`), both matching
  direct bumps; the price field matches the pricer.

## [1.281.0] - 2026-09-11

### Added
- `surface_arbitrage_report` and `surface_is_arbitrage_free` (in `rnd.py`):
  combine the butterfly (per-slice density-sign) and calendar (total-variance
  monotonicity) static-arbitrage checks into one surface-level verdict. Takes
  smiles in log-moneyness (one per expiry) and returns
  `{"butterfly": {t: [strikes]}, "calendar": [(k, t_lo, t_hi)]}`.
- Verified: an arbitrage-free SSVI surface passes (agreeing with SSVI's own
  check); a butterfly-violating slice is reported under its expiry; a
  calendar-violating term structure is reported; a flat surface is free.

## [1.280.0] - 2026-09-11

### Added
- `calendar_arbitrage_violations` and `surface_is_calendar_arbitrage_free` (in
  `rnd.py`): model-free calendar-arbitrage detection across a smile term
  structure. Flags `(k, t_lo, t_hi)` points where total implied variance
  `w(k, t) = sigma^2 t` decreases with maturity — a horizontal-spread arbitrage.
  Works for any smiles expressed in log-moneyness (SVI, SABR, vanna-volga, raw).
- Verified: a flat-vol term structure is arbitrage-free; a dropping total
  variance is flagged; the verdict agrees with SSVI's own calendar check;
  length-mismatch and non-increasing expiries raise.

## [1.279.0] - 2026-09-11

### Added
- `smile_arbitrage_violations` and `smile_is_arbitrage_free` (in `rnd.py`):
  model-free butterfly-arbitrage detection for any implied-vol smile. Scans a
  log-moneyness grid and flags strikes where the Breeden-Litzenberger density is
  negative (a negative-cost butterfly). Works for SVI, SABR, vanna-volga, or raw
  quotes.
- Performance: the extended-Sobol geometric-Asian direction-number test now runs
  at `n_rand=8` instead of 24 across the 2..12 step sweep (`4*SE` still catches a
  broken direction number), cutting ~6s off the fast gate.
- Verified: flat and convex-SVI smiles are flagged arbitrage-free; a
  butterfly-violating SVI slice is flagged; the density-sign verdict agrees with
  SVI's analytic g-function across slices.

## [1.278.0] - 2026-09-11

### Added
- `double_knock_in_call` (in `exotics.py`): double-barrier knock-in call priced
  by in-out parity, `vanilla - double_knock_out_call`. A knock-in and knock-out
  with the same strike and corridor partition every path, so they sum to the
  vanilla call.
- Verified: `DKI + DKO = vanilla` exactly; wide barriers give a near-zero
  knock-in; a tighter band raises the knock-in value; the price stays in
  `[0, vanilla]`.

## [1.277.0] - 2026-09-11

### Added
- `double_knock_out_call` (in `exotics.py`): the Ikeda-Kunitomo (1992)
  double-barrier knock-out call — payoff `max(S_T - K, 0)` paid only if the
  continuously-monitored spot stays inside a corridor `(L e^{delta1 s},
  U e^{delta2 s})`. Prices via the truncated Ikeda-Kunitomo image series
  (`delta1 = delta2 = 0` gives flat barriers).
- Verified: wide barriers recover the vanilla call; a tighter band and higher
  vol both lower the value; the price stays in `[0, vanilla]`; and it matches a
  fine-step (6000-step) continuously-monitored Monte Carlo (which sits just above
  the closed form by the usual discrete-monitoring bias).

## [1.276.0] - 2026-09-11

### Added
- `double_no_touch_greeks` (in `exotics.py`): finite-difference Greeks of a
  double-no-touch — `delta`, `gamma`, `vega`, `theta`, plus the two barrier
  sensitivities `dV/dL` and `dV/dU`.
- Verified: `delta` matches an independent bump; `vega < 0` and `gamma < 0`
  (short vol, peaked inside the band); `theta > 0` (less time to knock);
  widening the band raises value (`dV/dL < 0`, `dV/dU > 0`), and both barrier
  sensitivities match direct bumps.

## [1.275.0] - 2026-09-11

### Added
- `double_no_touch` and `double_one_touch` (in `exotics.py`): continuously
  monitored double-barrier binaries paying `cash` if the spot stays inside
  `(L, U)` to expiry (DNT), or if it touches either barrier (DOT). DNT uses the
  Fourier eigenfunction expansion of the driftful survival probability in a
  strip.
- Verified: `U -> infinity` collapses DNT to the lower `no_touch`, `L -> 0` to
  the upper; `DNT + DOT = cash e^{-rt}`; a narrower band and higher vol both
  lower survival; and DNT matches a fine-step (4000-step) continuously-monitored
  Monte Carlo.

## [1.274.0] - 2026-09-10

### Added (tests)
- Cross-model surface density parity: an SSVI surface calibrated to SABR smiles
  at two expiries reproduces the SABR Breeden-Litzenberger risk-neutral density
  to ~4% at the money and ~6% in the wings at a fitted expiry; the fitted surface
  is arbitrage-free and both densities integrate to one. Complements the
  single-slice SVI/SABR parity check.

## [1.273.0] - 2026-09-10

### Added
- `volatility_swap_bounds_from_smile` (in `varswap.py`): brackets the fair
  volatility-swap strike between the at-the-money-forward implied vol (lower --
  the Carr-Lee zero-correlation proxy) and `sqrt(K_var)` (upper -- the Jensen
  bound from the variance-swap strike). The bracket width is the convexity /
  vol-of-vol premium the smile implies.
- Verified: a flat smile collapses the bracket; a convex smile keeps
  `lower < upper`; the upper bound equals `sqrt` of the variance-swap strike and
  the lower equals the ATMF vol; the premium widens with SVI curvature.

## [1.272.0] - 2026-09-10

### Changed
- Trimmed the fast gate: the RQMC knock-in/knock-out partition-identity tests
  (Parisian, barrier, barrier-digital) now run at `n_rand=6` instead of 24. The
  identity is exact path-by-path at a shared seed, so it holds for any `n_rand`;
  this drops ~10s off the non-slow suite with no loss of coverage.

### Added (tests)
- Cross-model density parity: an SVI slice calibrated to a SABR smile reproduces
  the SABR Breeden-Litzenberger risk-neutral density to ~2% at the money and ~5%
  in the wings, both densities are non-negative, and both integrate to one -- two
  independent constructions agreeing end to end.

## [1.271.0] - 2026-09-10

### Added
- `VannaVolgaSmile.vol_at_delta`, `.risk_reversal`, and `.butterfly` (in
  `vannavolga.py`): query the FX smile in delta space. `vol_at_delta` solves the
  fixed point `sigma = vol(strike_from_delta(sigma))`; `risk_reversal` and
  `butterfly` recompute the delta-consistent skew and convexity at any delta.
- Verified: at 25-delta the risk reversal and butterfly recover the inputs used
  to build the smile, and the pillar vols are returned exactly; a ~50-delta
  option sits near the ATM vol; a negative RR makes the 25d put richer than the
  25d call; a convex smile has a larger 10-delta butterfly than 25-delta.

## [1.270.0] - 2026-09-10

### Added
- `corrado_su_implied_vol` and `corrado_su_smile` (in `gramcharlier.py`): the
  Black-Scholes implied vol of a Corrado-Su price at a strike, and the smile over
  a strike grid. Turns Gram-Charlier `(sigma, skew, excess_kurt)` into an
  implied-vol curve for plotting or seeding an SVI/SABR fit.
- Verified: zero moments give a flat smile at `sigma`; negative skew lifts the
  low-strike put wing (monotone-decreasing smile); positive excess kurtosis lifts
  both wings above the at-the-money level; the implied vol reprices the
  Corrado-Su value; the smile round-trips through `calibrate_corrado_su`.

## [1.269.0] - 2026-09-10

### Added
- `calibrate_corrado_su` (in `gramcharlier.py`): fits the Corrado-Su
  `(sigma, skew, excess_kurt)` to a set of market call prices by Nelder-Mead
  least squares, with a smooth reparametrization keeping `sigma > 0`. Returns the
  three parameters plus the price RMSE.
- Verified: recovers known `(sigma, skew, kurt)` from synthetic Corrado-Su
  prices to 1e-2; a flat Black-Scholes surface calibrates to zero skew and
  excess kurtosis at the input vol; reprices every strike within RMSE.

## [1.268.0] - 2026-09-10

### Added
- `implied_spread_correlation_bs` (in `multiasset.py`): the correlation implied
  by a spread-option price under the Bjerksund-Stensland (2014) model, recovered
  by bisection on the price (monotone decreasing in `rho`). Companion to the Kirk
  `implied_spread_correlation`.
- Verified: round-trips a known `rho` for both calls and puts; the price is
  monotone decreasing in correlation; a quote above the `rho = -1` maximum
  raises.

## [1.267.0] - 2026-09-10

### Added
- `spread_option_bs_greeks` (in `multiasset.py`): finite-difference Greeks of the
  Bjerksund-Stensland (2014) spread option — two spot deltas, own-gammas,
  cross-gamma `d2V/dS1 dS2`, and `corr_vega` — matching the `spread_greeks`
  layout.
- Verified: delta1 > 0, delta2 < 0, positive gammas, negative cross-gamma and
  correlation sensitivity; delta1 matches an independent bump; call-minus-put
  deltas are `(+1, -1)`; the cross-gamma tracks `-sqrt(gamma1 gamma2)`; gammas
  are identical for call and put.

## [1.266.0] - 2026-09-10

### Added
- `spread_option_bs` (in `multiasset.py`): the Bjerksund-Stensland (2014)
  three-`d` closed-form spread-option approximation for `max(S1 - S2 - K, 0)`,
  generally more accurate than the Kirk approximation at wide strikes, high
  volatility, or dispersed leg vols, and exact (Margrabe) at `K = 0`.
- Verified: matches Margrabe at `K = 0`; satisfies put-call parity on the
  spread; agrees with 200k-400k-path Latin-hypercube Monte Carlo within a few
  standard errors; sits at least as close to MC as Kirk at a wide, dispersed-vol
  strike; the call decreases in correlation.

## [1.265.0] - 2026-09-10

### Added
- `risk_neutral_var_from_smile` and `risk_neutral_cvar_from_smile` (in `rnd.py`):
  risk-neutral Value-at-Risk and Conditional VaR (expected shortfall) of the
  terminal simple return `L = 1 - S_T/S0`, from an implied-vol smile.
  `VaR_alpha = 1 - Q(1-alpha)/S0` via the inverse smile CDF; CVaR integrates the
  truncated first moment `E[S_T ; S_T <= K]` against the Breeden-Litzenberger
  density.
- Verified against closed-form lognormal VaR/expected-shortfall on a flat smile;
  CVaR >= VaR; VaR rises with confidence; an equity-skew smile raises both tail
  metrics.

## [1.264.0] - 2026-09-10

### Added
- `risk_neutral_cdf_from_smile` and `risk_neutral_quantile_from_smile` (in
  `rnd.py`): the Breeden-Litzenberger risk-neutral CDF `F(K) = 1 + e^{rt} dC/dK`
  of the terminal spot from an implied-vol smile, and its inverse (bisection over
  a forward-standard-deviation bracket).
- Verified: a flat smile recovers the Black-Scholes `N(-d2)`; the CDF is monotone
  in `[0, 1]`; the quantile inverts the CDF; the CDF equals the integral of the
  Breeden-Litzenberger density; the median matches the lognormal
  `F e^{-sigma^2 t / 2}`; an equity-skew smile fattens the left tail.

## [1.263.0] - 2026-09-10

### Added
- SSVI surface smile consumers at a fitted expiry (in `ssvi.py`):
  `ssvi_variance_swap_strike`, `ssvi_vix`, `ssvi_svix`, `ssvi_density`
  (Breeden-Litzenberger), and `ssvi_bkm_moments`. Each maps a strike to the
  surface's Black vol `implied_vol(ln(K/F), t)` on the forward and feeds the
  existing smile machinery.
- Verified: on a flat-ATM (`theta_t = sigma^2 t`) surface the var-swap and
  VIX/SVIX recover `sigma`; the density integrates to 1 and has mean equal to
  the forward and is non-negative on an arbitrage-free slice; BKM skewness
  tracks the sign of the surface `rho`. Requesting an unfitted expiry raises.

## [1.262.0] - 2026-09-10

### Added
- `svi_density` (in `svi.py`): Breeden-Litzenberger risk-neutral density `g(K)`
  implied by a raw-SVI slice, `g(K) = e^{rt} d^2C/dK^2` with the call priced at
  the slice's smile vol on the forward.
- Verified: a flat slice matches the closed-form lognormal density; the density
  is non-negative on a butterfly-arbitrage-free slice, integrates to 1, and has
  mean equal to the forward; a wing-ok but butterfly-violating slice produces a
  negative density exactly where `svi_g < 0`.

## [1.261.0] - 2026-09-10

### Added
- `sabr_bkm_moments` (in `sabr.py`): risk-neutral variance, skewness, and excess
  kurtosis implied by a SABR smile via Bakshi-Kapadia-Madan moment replication.
- Verified: a near-flat slice is symmetric; negative `rho` gives negative
  risk-neutral skew (positive `rho` positive); higher vol-of-vol raises the
  excess kurtosis.

## [1.260.0] - 2026-09-10

### Added
- `svi_bkm_moments` (in `svi.py`): risk-neutral variance, skewness, and excess
  kurtosis implied by a raw-SVI slice via the Bakshi-Kapadia-Madan moment
  replication of the slice's smile.
- Verified: a flat slice is near-symmetric (skew ~ 0, excess kurtosis ~ 0);
  negative SVI `rho` gives negative risk-neutral skew and positive `rho`
  positive skew; a convex slice has positive excess kurtosis.

## [1.259.0] - 2026-09-10

### Added
- `sabr_variance_swap_strike` and `sabr_vix` (in `sabr.py`): the fair
  variance-swap strike (annualized variance) and the VIX-style index implied by
  a SABR smile, each pricing the strike chain at `sabr_vol(F, K, ...)` and
  replicating via `variance_swap_from_smile` / `vix_from_smile`.
- Verified: at `beta = 1`, `nu -> 0` the strike is `alpha^2` and the VIX is
  `100 * alpha`; a skewed slice's variance exceeds the ATM variance; higher
  vol-of-vol raises the strike.

## [1.258.0] - 2026-09-10

### Added
- `svi_vix` and `svi_svix` (in `svi.py`): the VIX-style and Martin SVIX indices
  implied by a raw-SVI slice, mapping each strike to `p.implied_vol(ln(K/F), t)`
  and replicating via `vix_from_smile` / `svix_from_smile` (both reported as
  `100 * index`).
- Verified: a flat slice returns `100 * sigma` for both; a skewed slice lifts the
  VIX above `100 *` the ATM vol; higher convexity (`b`) raises the index.

## [1.257.0] - 2026-09-10

### Added
- `svi_variance_swap_strike` (in `svi.py`): fair variance-swap strike
  (annualized variance) replicated model-consistently from a raw-SVI slice --
  maps each strike to `p.implied_vol(ln(K/F), t)` and feeds the smile to
  `variance_swap_from_smile`.
- Verified: a flat slice (`b = 0`) returns the flat variance `sigma^2`; a
  skewed/convex slice returns a variance above the ATM variance (the convexity
  premium); higher `b` (more convexity) raises the strike.

## [1.256.0] - 2026-09-10

### Added
- `DiscountCurve` analytics methods: `instantaneous_forward` (`-d ln DF/dT` by
  central FD), `annuity` (fixed-leg PV01), `swap_value` (unit-notional payer/
  receiver), and `swap_dv01` (value change for a 1bp parallel curve drop via a
  shifted-curve wrapper).
- Verified: on a flat curve the instantaneous forward equals the zero rate; a
  par-struck swap has zero value; the annuity matches a manual sum; a payer's
  DV01 is negative on a rate drop and a receiver's positive, of order
  annuity x 1bp.

## [1.255.0] - 2026-09-10

### Added
- Ho-Lee caps and floors: `holee_caplet`, `holee_floorlet`, `holee_cap`,
  `holee_floor`, each a strip via the bond-put/call identity. Completes cap/floor
  coverage across every short-rate model (Vasicek, CIR, Ho-Lee, Cheyette) plus
  G2++.
- Verified: `cap - floor` equals the underlying swap value (parity) to 1e-8; the
  cap equals the sum of its caplets; positive and monotone in strike.

## [1.254.0] - 2026-09-10

### Added
- Vasicek and CIR caps and floors: `vasicek_caplet`/`vasicek_floorlet`/
  `vasicek_cap`/`vasicek_floor` and `cir_caplet`/`cir_floorlet`/`cir_cap`/
  `cir_floor`, each a strip of caplets/floorlets via the bond-put/call identity
  on the model's exact bond option.
- Verified: for both models `cap - floor` equals the underlying swap value
  (parity) to 1e-8, and the cap equals the sum of its caplets; positive and
  monotone in strike.

## [1.253.0] - 2026-09-10

### Added
- `cheyette_floorlet`, `cheyette_cap`, `cheyette_floor` (in `cheyette.py`):
  Cheyette (Hull-White) floorlet (bond-call identity) and cap/floor as strips of
  caplets/floorlets.
- Verified: `cap - floor` equals the underlying swap value (parity) to 1e-8; cap
  and floor equal the sums of their legs; both positive and monotone in strike.

## [1.252.0] - 2026-09-10

### Added
- `g2pp_floorlet`, `g2pp_cap`, `g2pp_floor` (in `g2pp.py`): G2++ floorlet (via
  the bond-call identity) and cap/floor as strips of caplets/floorlets over a
  schedule of `(t_i, P(0, t_i))` dates.
- Verified: `cap - floor` equals the underlying swap value (put-call parity) to
  1e-8; the cap and floor equal the sums of their caplets/floorlets; both
  positive; a higher strike lowers the cap and raises the floor.

## [1.251.0] - 2026-09-10

### Added
- `sabr_option_greeks` now also returns the total `gamma` (`d2Price/dF2`),
  computed by a central difference of the SABR-repriced surface so it captures
  the smile backbone's curvature -- not just the vol-fixed Black gamma.
- Verified: the total gamma matches a finite difference that recomputes the SABR
  vol at each bumped forward, and is positive.

## [1.250.0] - 2026-09-10

### Added
- `sabr_option_greeks` (in `sabr.py`): Greeks of a Black-76 option priced at the
  Hagan SABR smile vol. The total `delta = black_delta + vega * dsigma/dF`
  includes the smile backbone (the exact AD `dsigma/dF` from
  `sabr_sensitivities`), so it differs from the vol-fixed Black delta. Returns
  `price`, `vol`, `delta` (backbone-adjusted), `black_delta`, and `vega`, with an
  optional `discount` to scale the forward figures to present value.
- Verified: the total delta matches a finite difference that recomputes the SABR
  vol at each bumped forward; it differs measurably from the Black delta; vega
  positive; `vol` equals `sabr_vol`; price and delta scale linearly in `discount`.

## [1.249.0] - 2026-09-10

### Added
- `cir_coupon_bond_option` and `cir_swaption` (in `cir.py`): exact European
  coupon-bond option and swaption under CIR by Jamshidian decomposition. The CIR
  bond is monotone in `r0`, so the critical-rate `r*` split into per-cashflow CIR
  zero-coupon-bond options (`cir_bond_option`) is exact; a payer swaption is a put
  on the fixed-leg coupon bond struck at the notional. Completes the full option
  chain (bond -> coupon-bond option -> swaption) for all three affine short-rate
  models: Vasicek, Ho-Lee, CIR.
- Verified: a single cashflow reduces to the scaled `cir_bond_option`; swaption
  parity `payer - receiver = annuity (swap_rate - strike)` holds to 1e-7; both
  legs positive; a higher strike lowers the payer.

## [1.248.0] - 2026-09-10

### Added
- `cir_bond_option` (in `cir.py`): exact European option on a CIR zero-coupon
  bond (Cox-Ingersoll-Ross 1985), using the noncentral chi-square CDF already in
  the library. The call is `P(0,t_bond) X2(...) - strike P(0,t_option) X2(...)`
  with critical rate `r* = ln(A/strike)/B`; puts follow from put-call parity.
- Verified: put-call parity `C - P = P(0,t_bond) - strike P(0,t_option)` to
  1e-9; the call matches a fine-step CIR Monte Carlo (0.02260, slow test); both
  prices positive; a higher strike lowers the call.

## [1.247.0] - 2026-09-10

### Added
- `holee_bond_option`, `holee_coupon_bond_option`, `holee_swaption` (in
  `holee.py`): exact Ho-Lee zero-coupon-bond option (Black-style with bond vol
  `sigma (t_bond - t_option) sqrt(t_option)`), coupon-bond option by Jamshidian
  decomposition, and European swaption via the coupon-bond-option identity (payer
  = put on the fixed-leg coupon bond struck at notional).
- Verified: a single cashflow reduces to the scaled zero-coupon bond option;
  the bond option satisfies put-call parity; swaption parity
  `payer - receiver = annuity (swap_rate - strike)` holds to 1e-8; both swaption
  legs positive.

## [1.246.0] - 2026-09-10

### Added
- `vasicek_swaption` (in `vasicek.py`): exact European swaption under Vasicek via
  the coupon-bond-option identity. The fixed leg plus notional is a coupon bond;
  a payer swaption is a put on it struck at the notional and a receiver a call,
  both priced by the Jamshidian `coupon_bond_option` -- no normal-model
  approximation.
- Verified: swaption parity `payer - receiver = annuity (swap_rate - strike)`
  holds to 1e-8; at the forward swap rate payer equals receiver; a higher strike
  lowers the payer; both legs positive.

## [1.245.0] - 2026-09-10

### Added
- `vasicek_coupon_bond_option` (in `vasicek.py`): European option on a
  coupon bond under Vasicek by Jamshidian's decomposition. Solves for the
  critical short rate `r*` where the bond value at expiry equals the strike, then
  prices the option as the `c_i`-weighted sum of zero-coupon-bond options struck
  at `K_i = P(t_option, t_i | r*)` -- exact, no simulation.
- Verified: a single cashflow reduces to the scaled zero-coupon `bond_option`;
  put-call parity `C - P = sum c_i P(0,t_i) - K P(0,t_option)` holds to 1e-8;
  prices positive; a higher strike lowers the call.

## [1.244.0] - 2026-09-10

### Added
- `cms_adjustment_greeks` (in `cms.py`): sensitivities of the standard-model CMS
  convexity adjustment -- `d_forward` and `d_sigma` by central finite differences
  on `cms_adjustment_standard`.
- Verified: both match finite differences; `d_sigma > 0` (more vol means more
  convexity, a larger adjustment); a zero-vol adjustment has zero sensitivity.

## [1.243.0] - 2026-09-10

### Added
- `swaption_greeks` (in `rates.py`): analytic Greeks of a normal-model European
  swaption. The value is `annuity * Bachelier(swap_rate, strike, expiry, 0,
  sigma_n)`, so `rate_delta`, `rate_gamma`, and `vega` are the Bachelier Greeks
  scaled by the annuity (also returned).
- Verified: `rate_delta` matches a finite difference of `swaption_price`; a payer
  swaption's rate delta is positive and a receiver's negative; the ATM payer
  delta is exactly `annuity/2` (Bachelier ATM delta 0.5); vega positive.

## [1.242.0] - 2026-09-10

### Added
- `caplet_greeks`, `cap_greeks`, `floor_greeks` (in `rates.py`): analytic Greeks
  of normal-model caplets/floorlets and their strips. Each caplet's value is
  `discount * accrual * Bachelier(F, K, ...)`, so its `rate_delta`, `rate_gamma`,
  and `vega` are the Bachelier Greeks in the forward rate scaled by the same
  factor; `cap_greeks`/`floor_greeks` sum them across periods.
- Verified: the cap `rate_delta` matches a finite difference of `cap_price`; a
  cap's rate delta is positive and a floor's negative; both vegas positive; the
  cap Greeks equal the sum of the caplet Greeks.

## [1.241.0] - 2026-09-10

### Added
- `dual_swap_dv01` (in `dualcurve.py`): risk of a dual-curve (OIS-discounted,
  projection-forward) swap -- the exact fixed-leg `pv01` (annuity), the total
  `dv01` for a 1bp parallel drop of both curves, and the split `ois_dv01` /
  `proj_dv01` from shifting only the discount or only the projection curve.
  Curve shifts use a light `_ShiftedCurve` wrapper (`df(T) e^{-dr T}`).
- Verified: `pv01` equals `-dV/d(fixed_rate)` for a payer to machine precision;
  the total `dv01` is the sum of the per-curve DV01s to first order; a payer's
  `dv01` is negative on a rate drop and a receiver's positive.

## [1.240.0] - 2026-09-10

### Added
- `bermudan_swaption_g2pp_greeks` (in `bermudan_swaption.py`): Greeks of a G2++
  Bermudan swaption by common-random-number bumps on `bermudan_swaption_g2pp` --
  `d_fixed` (dV/d fixed_rate) and `curve_dv01` (value change for a 1bp parallel
  curve drop). Same-seed repricing shares the G2++ state paths for low-variance
  differences.
- Verified: `d_fixed` matches a CRN finite difference; it is negative for a payer
  (paying a higher fixed rate is worth less) and positive for a receiver;
  reproducible under a fixed seed.

## [1.239.0] - 2026-09-10

### Added
- `vasicek_bond_greeks` (in `vasicek.py`) and `holee_bond_greeks` (in
  `holee.py`): exact rate sensitivities of the Vasicek and Ho-Lee zero-coupon
  bonds (`rho_r`, `gamma_r`, `duration`, `convexity`). Vasicek uses the affine
  `B(t)`; Ho-Lee is linear in `r0` inside the exponent, so its duration is
  exactly the maturity `t` and convexity `t^2`.
- Verified: both `rho_r` match finite differences of their bond prices;
  `convexity = duration^2`; the Ho-Lee duration equals the maturity; zero-
  maturity bonds are flat.

## [1.238.0] - 2026-09-10

### Added
- `cir_bond_greeks` (in `cir.py`): exact rate sensitivities of a CIR
  zero-coupon bond. Since `P = A(t) e^{-B(t) r0}`, the short-rate delta is
  `rho_r = -B P`, gamma `gamma_r = B^2 P`, rate duration `B`, and convexity
  `B^2` -- all closed form.
- Verified: `rho_r` and `gamma_r` match finite differences of
  `cir_zero_coupon_bond`; `rho_r < 0`, `gamma_r > 0`, `convexity = duration^2`,
  `duration = -rho_r / price`; a zero-maturity bond is flat.

## [1.237.0] - 2026-09-10

### Added
- `cheyette_bond_option_greeks` (in `cheyette.py`): exact discount-factor deltas
  `delta_T` (`dV/dP0T = N(d1)`, underlying-bond delta) and `delta_S` (`dV/dP0S`),
  plus the vol sensitivity `vega`, of a Cheyette (Hull-White) zero-coupon-bond
  option.
- Verified: both deltas match finite differences of `cheyette_bond_option`; the
  call's `delta_T` is in (0,1) and `delta_S` negative; vega positive; the put's
  `delta_T` negative.

## [1.236.0] - 2026-09-10

### Added
- `vasicek_bond_option_greeks` (in `vasicek.py`): short-rate sensitivities
  `rho_r` (dV/dr0) and `gamma_r` (d2V/dr0^2), plus the vol sensitivity `vega`
  (dV/dsigma) of a Vasicek zero-coupon-bond option by central finite differences
  on `bond_option`.
- Verified: `rho_r` matches a finite difference; a bond call's `rho_r` is
  negative (it loses value as the short rate rises) and a put's is positive;
  vega is positive.

## [1.235.0] - 2026-09-10

### Added
- `g2pp_bond_option_greeks` (in `g2pp.py`): Greeks of a G2++ zero-coupon-bond
  option -- the exact discount-factor deltas `delta_T` (`dV/dP0T = N(d1)`, the
  underlying-bond delta) and `delta_S` (`dV/dP0S`), plus the two factor-vol
  vegas `vega_sigma` and `vega_eta` by finite difference.
- Verified: both deltas match finite differences of `g2pp_bond_option` (delta_T
  equals the analytic `N(d1)`); the call's `delta_T` is in (0,1) and `delta_S`
  negative; both vegas positive; the put's `delta_T` is negative.

## [1.234.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~31s -> ~20s. The two-asset LSM-Greeks and
  rBergomi-Greeks tests each run 5+ re-prices; the `test_reproducible`
  determinism checks (which price twice) dropped to tiny sizes (n_steps 6-10,
  n_paths 800-1500), and the sign/ordering checks (min-put `delta2 < delta1`,
  basket `delta1 > delta2`, both-deltas-negative) to 10 steps / 3000-4000 paths
  -- all robust across seeds at the smaller size. Accuracy cross-checks stay
  under `-m slow`.

## [1.234.0] - 2026-09-10

### Added
- `andreasen_huge_strike_greeks` (in `andreasenhuge.py`): strike-space Greeks of
  the arbitrage-free Andreasen-Huge call surface -- the `dual_delta` (`dC/dK`,
  discounted; equals `-e^{-rT}` times the risk-neutral exceedance probability)
  and the `rnd` (Breeden-Litzenberger risk-neutral density `e^{rT} d2C/dK2`) at
  each interior strike.
- Verified: the density is non-negative for any positive local vols (the
  scheme's convexity), the dual delta is monotone in `[-e^{-rT}, 0]`, and the
  density integrates to ~1 -- rising toward 1 as the strike grid widens
  (0.89 -> 0.994 -> 0.9997), confirming the tail-truncation deficit rather than a
  normalization error.

## [1.233.0] - 2026-09-10

### Added
- `rbergomi_greeks_cv` (in `rbergomi.py`): delta, gamma, and `vega_xi0`
  (forward-variance-level sensitivity) of a rough-Bergomi call by
  common-random-number bumps on the conditional control-variate estimator
  `rbergomi_price_cv`. Same-seed repricing shares the volatility-driving Brownian
  paths, so the finite differences are low-variance.
- Verified: delta matches a common-random-number finite difference of
  `rbergomi_price_cv` to machine precision; call delta in (0,1) with positive
  gamma and `vega_xi0`; reproducible under a fixed seed.

## [1.232.0] - 2026-09-10

### Added
- `rough_heston_greeks` (in `rough_heston.py`): spot Greeks (`delta`, `gamma`)
  and the initial-variance sensitivity `vega_v0` of a rough-Heston option by
  central finite differences on `rough_heston_price`. Each re-price runs the
  O(n_grid^2) fractional-Riccati solve, so it is comparatively slow; the default
  `n_grid` matches the pricer's (small `H` needs a fine grid to stay stable).
- Verified: at `H = 0.5` the delta matches a finite difference of the classical
  Heston price (`xi = kappa * nu`); at `H = 0.3` (rough) the delta matches a
  finite difference of `rough_heston_price`, with call delta in (0,1), positive
  gamma and `vega_v0`.

## [1.231.0] - 2026-09-10

### Added
- `double_heston_greeks` (in `double_heston.py`): spot Greeks (`delta`,
  `gamma`) plus a per-factor initial-variance sensitivity `vega_v01` and
  `vega_v02` of a double-Heston option by central finite differences on
  `double_heston_price` (one-sided at each variance floor).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma; both variance vegas are positive.

## [1.230.0] - 2026-09-10

### Added
- `meixner_greeks` (in `meixner.py`): spot Greeks (`delta`, `gamma`, `theta`)
  plus the asymmetry/skew sensitivity `d_b` (dV/db) of a Meixner option by
  central finite differences on `meixner_price` (the `d_b` bump is clipped to
  keep `b` in `(-pi, pi)`). Completes the Levy-model Greek family (Merton, Kou,
  VG, NIG, CGMY, Meixner).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the skew sensitivity is finite.

## [1.229.0] - 2026-09-10

### Added
- `cgmy_greeks` (in `cgmy.py`): spot Greeks (`delta`, `gamma`, `theta`) plus the
  tail-activity sensitivity `d_Y` (dV/dY) of a CGMY option by central finite
  differences on `cgmy_price` (`d_Y` falls back to a one-sided bump near `Y = 2`).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the tail sensitivity is finite.

## [1.228.0] - 2026-09-10

### Added
- `bates_greeks` (in `bates.py`): spot Greeks (`delta`, `gamma`), the
  initial-variance sensitivity `vega_v0`, and the jump-intensity sensitivity
  `d_lambda` of a Bates (Heston + Merton jumps) option by central finite
  differences on `bates_price` (one-sided at the `v0`/`lam` floors).
- Verified: at `lam = 0` the delta matches a finite difference of the exact
  Heston price; with jumps the delta matches a finite difference of `bates_price`;
  call delta in (0,1) with positive gamma and `vega_v0`.

## [1.227.0] - 2026-09-10

### Added
- `nig_greeks` (in `nig.py`): spot Greeks (`delta`, `gamma`, `theta`) plus
  process-parameter sensitivities (`d_alpha` tail steepness, `d_beta` skew) of a
  NIG option by central finite differences on `nig_price`. The `d_beta` bump is
  clipped to keep `|beta| < alpha` on both sides.
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the parameter sensitivities are finite.

## [1.226.0] - 2026-09-10

### Added
- `variance_gamma_greeks` (in `variancegamma.py`): delta, gamma, vega
  (Brownian-vol sensitivity), and `theta_greek` (calendar decay) of a
  Variance-Gamma option by central finite differences on `variance_gamma_price`.
  The calendar Greek is named `theta_greek` to avoid clashing with the VG skew
  *parameter* `theta`.
- Verified: a small `nu` gives a delta close to the Black-Scholes delta (the
  `nu -> 0` limit); the delta matches a finite difference; call delta in (0,1)
  with positive gamma/vega, put delta negative.

## [1.225.0] - 2026-09-10

### Added
- `kou_greeks` (in `kou.py`): delta, gamma, vega (diffusion-vol sensitivity),
  and theta of a Kou double-exponential jump-diffusion option by central finite
  differences on `kou_price`.
- Verified: at `lam = 0` (no jumps) the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; with jumps the delta matches a finite difference, the
  call delta is in (0,1) with positive gamma/vega, and the put delta is negative.

## [1.224.0] - 2026-09-10

### Added
- `bachelier_greeks` (in `bachelier.py`): bundles the Bachelier (normal-model)
  Greeks -- `delta`, `gamma`, `vega` from the existing exact closed forms plus
  `theta` (calendar decay) by central finite difference.
- Verified: delta/gamma/vega equal their standalone closed forms; theta matches
  a finite difference; the ATM call delta is exactly 0.5 at zero rate; the put
  delta is negative.

## [1.223.0] - 2026-09-10

### Added
- `merton_jump_greeks` (in `merton.py`): delta, gamma, vega (diffusion-vol
  sensitivity), and theta of a Merton jump-diffusion option by central finite
  differences on `merton_jump_price`.
- Verified: at `lam = 0` (no jumps) the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; with jumps the delta matches a finite difference, the
  call delta is in (0,1) with positive gamma/vega, and the put delta is negative.

## [1.222.0] - 2026-09-10

### Added
- `cev_greeks` (in `cev.py`): delta, gamma, vega, and theta of a CEV option by
  central finite differences on `cev_price`.
- Verified: delta matches a finite difference across `beta` in {0.3, 0.5, 0.7};
  the call delta is in (0,1) with positive gamma/vega and the put delta is
  negative; a lower `beta` (steeper local-vol skew) lifts the ATM call delta.
  (The `beta -> 1` Black-Scholes limit is *not* asserted -- the Schroder
  noncentral-chi-square parameters are numerically unstable as `beta -> 1`, so
  the checks stay in the well-conditioned interior.)

## [1.221.0] - 2026-09-10

### Added
- `cliquet_greeks` (in `forwardstart.py`): delta, gamma, vega, and theta of a
  cliquet (ratchet) by central finite differences on `cliquet_price` (theta
  shifts every reset date together).
- Verified: every strike scales with the spot (`alpha * S_reset`), so the
  cliquet is homogeneous of degree 1 in `S` -- `delta` equals `price / S`
  exactly and `gamma` is zero. Delta and vega also match finite differences of
  `cliquet_price`.

## [1.220.0] - 2026-09-10

### Added
- `forward_start_greeks` (in `forwardstart.py`): Greeks of a forward-start
  option. Since the Rubinstein price `S e^{(b-r) t_start} u` is exactly linear in
  the spot (the unit price `u` does not depend on `S`), `delta = e^{(b-r) t_start} u`
  is constant in the spot and `gamma = 0` -- a forward-start has no spot gamma
  until its strike is fixed. `vega` and `theta` are central finite differences of
  the closed form.
- Verified: delta matches a finite difference and is constant across spot levels;
  gamma is exactly zero; vega matches a finite difference; theta is zero at
  `b = r` (only the maturity gap matters) and nonzero under a carry.

## [1.219.0] - 2026-09-10

### Added
- `displaced_diffusion_greeks` (in `displaced.py`): delta, gamma, vega, and
  theta of a displaced-diffusion option by central finite differences on
  `displaced_diffusion_price`.
- Verified: at `shift = 0` the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; the delta matches a finite difference at a positive
  shift; call delta in (0,1) with positive gamma/vega, put delta negative.

## [1.218.0] - 2026-09-10

### Added
- `gap_option_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a gap
  option (Reiner-Rubinstein) by central finite differences on `gap_option`.
- Verified: setting `K_trigger = K_payoff` reduces the Greeks exactly to the
  vanilla Black-Scholes Greeks; the gap delta matches a finite difference for
  both call (trigger above payoff) and put.

## [1.217.0] - 2026-09-10

### Added
- `compo_option_greeks` (in `quanto.py`): Greeks of a composite (compo) FX
  option. The price is a Black-Scholes price on the domestic-currency asset with
  the combined vol `sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)`
  and carry `r_domestic - q`, so `delta`/`gamma` are exact BSM Greeks; `vega`,
  `fx_vega`, and `corr_vega` are finite differences.
- Verified: delta matches a finite difference; at `sigma_fx = 0` the delta
  equals the plain BSM delta at carry `r_domestic - q`; a compo is long FX
  volatility (`fx_vega > 0`), unlike a quanto; call gamma and vega positive.

### Fixed
- The `fx_vega` finite difference in `quanto_option_greeks` /
  `compo_option_greeks` probed a negative `sigma_fx` when `sigma_fx` sat below
  the bump size (e.g. `sigma_fx = 0`), raising `ValueError`. It now falls back to
  a one-sided difference at the volatility floor.

## [1.216.0] - 2026-09-10

### Added
- `quanto_option_greeks` (in `quanto.py`): Greeks of a quanto option. The price
  is a Black-Scholes price on the foreign asset with the quanto-adjusted carry
  `b_q = r_foreign - q_asset - rho sigma_asset sigma_fx`, so the spot enters only
  through that BSM price -- `delta` and `gamma` are the exact BSM Greeks at `b_q`
  (no finite difference). `vega` (asset vol), `fx_vega` (FX vol), and `corr_vega`
  (correlation) are central finite differences of the closed form.
- Verified: delta matches a finite difference; at `rho = 0` the delta equals the
  plain BSM delta at carry `r_foreign - q`; `corr_vega` is negative for a
  positive-rho call (higher correlation lowers the quanto carry); call gamma and
  asset vega are positive.

## [1.215.0] - 2026-09-10

### Added
- `power_option_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a
  power option (`max(S_T^power - K, 0)`) by central finite differences on the
  closed-form `power_option`.
- Verified: at `power = 1` the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; at `power = 2` delta matches a finite difference and the
  call gamma/vega are positive; the put delta is negative.

## [1.214.0] - 2026-09-10

### Added
- `chooser_option_greeks` (in `chooser.py`): delta, gamma, and vega of a simple
  chooser option, exact by decomposition. The chooser is exactly a call to `T`
  plus a put struck at the discounted-forward level expiring at the choice date,
  and both legs are Black-Scholes prices in `S` and `sigma`, so the Greeks are
  the exact sums of the two legs' BSM Greeks -- no finite difference.
- Verified: delta, gamma, and vega match finite differences of `chooser_option`
  to step precision; gamma and vega are positive (long both a call and a put).

## [1.213.0] - 2026-09-10

### Added
- `compound_option_greeks` (in `compound.py`): delta, gamma, vega, and theta of
  a Geske compound option by central finite differences on `compound_option`
  (all four kinds: call/put on call/put). Theta shifts both expiries together.
- Verified: delta matches a finite difference of the closed form; a
  call-on-call has positive delta, gamma, and vega; a put-on-call has negative
  delta. Noted (via cross-check) that vega is *not* universally positive -- a put
  on an option is short the compound optionality, so `put-on-call` and
  `put-on-put` can have negative vega; only the call compounds are asserted
  positive.

## [1.212.0] - 2026-09-10

### Added
- `geometric_asian_greeks` (in `exotics.py`): Greeks of the continuously-
  monitored geometric-average Asian. Since the Kemna-Vorst price is exactly a
  Black-Scholes price at the adjusted vol `sigma_A = sigma/sqrt(3)` and carry
  `b_A = (b - sigma^2/6)/2`, the spot enters only through that BSM price, so
  `delta` and `gamma` are the exact BSM Greeks there (no finite difference);
  `vega`, `theta`, and `rho` are central differences of the exact closed form.
- Verified: delta, gamma, and vega match finite differences of `geometric_asian`
  to machine / step precision; call delta in (0,1), gamma and vega positive; put
  delta negative; the `price` field equals `geometric_asian`.

## [1.211.0] - 2026-09-10

### Added
- `bermudan_basket_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American basket option by common-random-number bumps on
  `bermudan_basket_lsm`. Same-seed repricing shares the Brownian shocks so the
  finite differences are low-variance; the regression is re-fit at each bump.
  Completes the two-asset LSM Greek family (max-call, spread, min-put, basket).
- Verified: without dividends the two spot deltas match the European
  `basket_greeks` deltas within noise; the call deltas are positive and the
  larger-weight asset carries the larger delta; the put deltas are negative;
  reproducible.

## [1.210.0] - 2026-09-10

### Added
- `sobol_arithmetic_asian_rqmc` (in `sobol.py`): randomized-QMC arithmetic Asian
  with a geometric control variate -- the two strongest variance-reduction
  techniques together. With `control_variate=True` each path's estimator is
  `arith - geo + E[geo]`, where `E[geo]` is the exact discrete-geometric closed
  form over the same `n_steps` dates and the two averages (nearly perfectly
  correlated) share the path. Bridge construction + per-dimension
  Cranley-Patterson rotation.
- Verified: matches the control-variate `arithmetic_asian_mc` for call and put;
  the control cuts the standard error to ~0.11x the RQMC-only estimator (on top
  of the QMC gain), so ~9x below RQMC alone.

## [1.209.0] - 2026-09-10

### Changed
- Extended the Sobol generator from 6 to 12 dimensions, adding the canonical
  Joe-Kuo primitive-polynomial coefficients and direction-number seeds for dims
  7-12. Validated end-to-end: the RQMC geometric-Asian price matches the exact
  discrete-geometric closed form at every `n_steps` from 2 to 12 (a wrong
  polynomial/seed would break the low-discrepancy property and bias it),
  regression-tested by `test_extended_sobol_dims_match_closed_form`. This lets
  every `sobol_*_rqmc` routine use up to 12 monitoring dates.

### Added
- `sobol_parisian_rqmc` (in `sobol.py`): randomized-QMC Parisian barrier option
  with an honest standard error -- knock-out/knock-in triggered only after the
  spot spends a *consecutive* `window` on the barrier's far side. Bridge
  construction + per-dimension Cranley-Patterson rotation; the discrete analogue
  of `parisian_barrier_mc`.
- Verified: matches `parisian_barrier_mc` at equal `n_steps` for down-out and
  down-in; a knock-in plus its knock-out sum to the vanilla at the same seed; a
  longer required window raises the knock-out value.

## [1.208.0] - 2026-09-10

### Added
- `sobol_geometric_asian_rqmc` (in `sobol.py`): randomized-QMC geometric-average
  Asian option (call `max(G - K, 0)`, put `max(K - G, 0)`) with an honest
  standard error. Bridge construction + per-dimension Cranley-Patterson
  rotation. Because the discrete geometric average is exactly lognormal it has a
  closed form (`_discrete_geometric_asian`), so this is cross-checked against a
  *deterministic* reference -- the tightest possible.
- Verified: matches the exact closed form for call and put (also under a
  dividend carry) within the tiny RQMC standard error (< 0.02).

## [1.207.0] - 2026-09-10

### Added
- `bermudan_basket_lsm` (in `lsm.py`): American basket option on
  `w1 S1 + w2 S2` (call or put) by Longstaff-Schwartz. Two correlated GBMs; the
  continuation value is regressed on a quadratic basis in both spots plus the
  basket `{1, S1, S2, S1^2, S2^2, S1 S2, w1 S1 + w2 S2}` over the in-the-money
  paths.
- Verified: without dividends the call matches the European moment-matched
  `basket_option`; a 6% dividend on both assets gives a positive early-exercise
  premium on the call (~0.38); the put carries an early-exercise premium even
  without dividends (~0.26); reproducible.

## [1.206.0] - 2026-09-10

### Added
- `sobol_average_strike_rqmc` (in `sobol.py`): randomized-QMC average-strike
  Asian option with an honest standard error. The strike is the realized
  arithmetic average of the path (call `max(S_T - A, 0)`, put `max(A - S_T, 0)`).
  Bridge construction + per-dimension Cranley-Patterson rotation; the discrete
  analogue of `average_strike_asian_mc`.
- Verified: matches `average_strike_asian_mc` for call and put within MC error,
  at roughly a seventh of the plain-MC standard error at equal points; prices
  positive; SE < 0.02.

## [1.205.1] - 2026-09-10

### Tests
- Trimmed the two `bermudan_min_put_lsm_greeks` sign-check fast tests (each
  re-prices five times) from 20-25 steps / 15000-30000 paths to 12-15 / 5000-8000;
  the delta signs and the `delta2 < delta1` ordering are robust across seeds at
  the smaller size. Fast gate ~60s -> ~30s.

## [1.205.0] - 2026-09-10

### Added
- `bermudan_min_put_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American min-put (worst-of protective put) by
  common-random-number bumps on `bermudan_min_put_lsm`. Same-seed repricing
  shares the Brownian shocks so the finite differences are low-variance; the LSM
  regression is re-fit at each bump. Completes the two-asset LSM Greek family
  (max-call, spread, min-put).
- Verified: both spot deltas are negative (a higher spot lifts the min and
  shrinks the put); the lower-starting asset carries the larger-magnitude delta
  (it is more often the min); the price exceeds the European
  `worst_of_put_closed` by the early-exercise premium; reproducible.

## [1.204.0] - 2026-09-10

### Added
- `sobol_cliquet_rqmc` (in `sobol.py`): randomized-QMC capped cliquet (ratchet)
  note with an honest standard error -- the same product as `capped_cliquet_mc`
  (per-period returns clipped to `[local_floor, local_cap]`, running sum clipped
  to `[global_floor, global_cap]`). Brownian motion at the reset times comes from
  one Sobol point via a bridge on the reset grid; the per-period standardized
  shock is the bridge increment over `sqrt(dt_i)`, randomized by a per-dimension
  Cranley-Patterson rotation.
- Verified: matches `capped_cliquet_mc` with full caps and with uncapped-local /
  no-global settings within MC error; a tighter local cap lowers the value; the
  across-randomization SE is tight (< 0.01).

## [1.203.0] - 2026-09-10

### Added
- `sobol_double_knockout_rqmc` (in `sobol.py`): randomized-QMC double-knockout
  (corridor) option with an honest standard error. Pays the vanilla payoff only
  if the spot stays strictly inside `(lower, upper)` at every monitoring date,
  else the cash `rebate`. Bridge construction + per-dimension Cranley-Patterson
  rotation; the discrete analogue of `double_knockout_mc`.
- Verified: matches `double_knockout_mc` within MC error; a very wide corridor
  approaches the vanilla call; a narrower corridor is worth strictly less.

## [1.202.0] - 2026-09-10

### Added
- `sobol_autocallable_rqmc` (in `sobol.py`): randomized-QMC autocallable
  structured note with an honest standard error -- the same product as
  `autocallable_mc` (early redemption with accrued coupons at each observation,
  down-and-in protection at maturity). Each path's Brownian values at the
  observation dates come from one Sobol point via a Brownian bridge built
  directly on the (possibly non-uniform) observation grid (`_bridge_on_times`),
  randomized by a per-dimension Cranley-Patterson rotation.
- Verified: matches `autocallable_mc` with and without a protection barrier
  within MC error; the price rises with the coupon; the across-randomization SE
  is tight (< 0.01).

## [1.201.0] - 2026-09-10

### Added
- `sobol_barrier_digital_rqmc` (in `sobol.py`): randomized-QMC barrier-contingent
  cash-or-nothing digital with an honest standard error. Pays `cash` iff the
  option finishes in the money AND the barrier condition holds over the `n_steps`
  monitoring dates (`up-in`/`down-in`/`up-out`/`down-out`). Bridge construction +
  per-dimension Cranley-Patterson rotation; the discrete analogue of
  `barrier_digital_mc`.
- Verified: matches `barrier_digital_mc` across all four barrier types within MC
  error; a knock-in plus its knock-out sum to the plain `cash_or_nothing` digital
  at the same seed; the value scales linearly in `cash`.

## [1.200.1] - 2026-09-10

### Tests
- Trimmed the two `bermudan_spread_lsm_greeks` sign-check fast tests (each
  re-prices five times) from 20 steps / 15000 paths to 12 / 5000; they assert
  only delta signs, which are robust at the smaller size. Fast gate ~35s -> ~20s.

## [1.200.0] - 2026-09-10

### Added
- `bermudan_spread_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American spread option by common-random-number bumps on
  `bermudan_spread_lsm`. Same-seed repricing shares the Brownian shocks so the
  finite differences are low-variance; the LSM regression is re-fit at each bump.
- Verified: without dividends the two spot deltas match the European Kirk
  `spread_greeks` deltas within noise (slow test); the spread-call long-leg
  delta is positive and the short-leg delta negative, with the signs flipping
  for the put; reproducible under a fixed seed.

## [1.199.0] - 2026-09-10

### Added
- `sobol_barrier_rqmc` (in `sobol.py`): randomized-QMC discretely-monitored
  single-barrier vanilla (down/up, in/out, with rebate) and an honest standard
  error. Normals come from an `n_steps`-dim Sobol point through the Brownian
  bridge, randomized by a per-dimension Cranley-Patterson rotation; the discrete
  analogue of `barrier_mc(brownian_bridge=False)`.
- Verified: matches the discretely-monitored `barrier_mc` across all four
  barrier types within MC error; a knock-in plus its knock-out sum to the
  vanilla call at the same seed (they partition every path); the knock-out is
  below the vanilla.

## [1.198.0] - 2026-09-10

### Added
- `two_asset_digital_greeks` (in `multiasset.py`): Greeks of a two-asset
  correlated digital by finite differences on the exact `two_asset_digital`
  closed form (no Monte Carlo noise) -- the two spot deltas, own-gammas,
  cross-gamma, and the correlation sensitivity `corr_vega`.
- Verified: a both-above digital has positive spot deltas and a positive
  `corr_vega` (the two in-the-money events move together), a mixed above/below
  digital has negative `corr_vega`, and summed over the four exhaustive quadrants
  the correlation sensitivity is zero (total probability is `rho`-independent);
  Greeks scale linearly in `cash`; a `below` condition flips that spot delta's
  sign.

## [1.197.0] - 2026-09-10

### Added
- `bermudan_min_put_lsm` (in `lsm.py`): American put on the minimum of two assets
  `max(K - min(S1_T, S2_T), 0)` by Longstaff-Schwartz -- the worst-of protective
  put. Two correlated GBMs; the continuation value is regressed on a quadratic
  basis in both spots plus the running min
  `{1, S1, S2, S1^2, S2^2, S1 S2, min(S1,S2)}` over the in-the-money paths.
- Verified: exceeds the European Stulz `worst_of_put_closed` by a positive
  early-exercise premium (puts carry early-exercise value even without
  dividends, ~0.42 here); is worth at least its intrinsic on the min when deep
  in the money; reproducible under a fixed seed.

### Audit
- Swept every test for Monte Carlo functions called without a `seed` whose
  result feeds an assertion (the class of the `test_pde2d_american` flake fixed
  in 1.196.1). All remaining seedless calls are `pytest.raises(ValueError)`
  argument-validation checks (no sampling) or pass the seed via `**kw`; no
  further nondeterministic references found.

## [1.196.1] - 2026-09-10

### Fixed
- `test_pde2d_american.py` compared the deterministic two-asset ADI price to the
  *Monte Carlo* `best_of_call` / `worst_of_call` (seed unset), so the reference
  itself was random and its noise occasionally breached the tolerance -- the
  source of the rare "1 failed" seen in full-suite runs. Switched both to the
  exact Stulz `best_of_call_closed` / `worst_of_call_closed`; the checks are now
  deterministic (and the worst-of one no longer needs `-m slow`).

## [1.196.0] - 2026-09-10

### Changed
- The test suite now runs across all cores by default via `pytest-xdist`
  (`addopts = "-n auto"` in `pyproject.toml`; `pytest-xdist>=3` added to the
  `test` extra). The ~1700 tests are independent with no shared state, so
  distribution is safe; the fast gate drops from ~55s to ~12s on a 24-core box
  and CI gets a proportional speedup. Run serially with `-n0`.

### Tests
- Reduced the six RQMC lookback structural checks (floating and fixed strike)
  from 4096 paths / 24 randomizations to 2048 / 12; the sign and ordering
  assertions and the SE bound (actual SE ~0.006 vs a 0.05 threshold) hold
  comfortably at the smaller size.

## [1.195.0] - 2026-09-10

### Added
- `sobol_fixed_lookback_rqmc` (in `sobol.py`): randomized-QMC discretely-
  monitored fixed-strike lookback with an honest standard error. Call payoff
  `max(max_i S_{t_i} - K, 0)`, put payoff `max(K - min_i S_{t_i}, 0)` over `S_0`
  and the `n_steps` monitoring dates; normals from an `n_steps`-dim Sobol point
  through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
  rotation.
- Verified: the discrete price is below the continuously-monitored Conze-
  Viswanathan `fixed_strike_lookback` for call and put, and rises toward it as
  monitoring frequency grows (n=2 -> 6: 12.2 -> 14.6 vs 19.2 continuous); the
  call dominates the plain vanilla.

## [1.194.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~66s -> ~50s by shrinking the heavy two-asset LSM/RQMC
  structural checks (each stores full two-asset paths and re-prices several
  times): the max-call cross-gamma check 15 steps / 8000 paths -> 12 / 5000 (its
  sign is robust across seeds there), the delta-interval check 12 / 3000 ->
  10 / 2500, and the American-spread lower-bound check 25 / 15000 -> 15 / 6000.
  The accuracy cross-checks stay under `-m slow`.

## [1.194.0] - 2026-09-10

### Added
- `sobol_lookback_rqmc` (in `sobol.py`): randomized-QMC discretely-monitored
  floating-strike lookback with an honest standard error. Call payoff
  `S_T - min_i S_{t_i}`, put payoff `max_i S_{t_i} - S_T` over the `n_steps`
  monitoring dates; normals come from an `n_steps`-dim Sobol point through the
  Brownian bridge, randomized by a per-dimension Cranley-Patterson rotation so
  `n_rand` shifts give a genuine SE.
- Verified: the discrete price is below the continuously-monitored
  Goldman-Sosin-Gatto `floating_strike_lookback` (fewer sampling dates see less
  extreme highs/lows), for both call and put, and rises toward it as the
  monitoring frequency increases (n=2 -> 6: 11.9 -> 13.8 vs 17.2 continuous);
  the across-randomization SE is tight (< 0.05).

## [1.193.0] - 2026-09-10

### Added
- `bermudan_spread_lsm` (in `lsm.py`): American spread option
  `max(S1_T - S2_T - K, 0)` (call) or `max(K - (S1_T - S2_T), 0)` (put) by
  Longstaff-Schwartz. Two correlated GBMs; the continuation value is regressed
  on a quadratic basis in both spots plus the spread
  `{1, S1, S2, S1^2, S2^2, S1 S2, S1 - S2}` over the in-the-money paths.
- Verified: without dividends it matches the European Kirk `spread_option`
  (11.83); an 8% dividend on the long leg produces a positive early-exercise
  premium (~0.75); it is never cheaper than the European; call and put prices
  are positive and reproducible.

## [1.192.0] - 2026-09-10

### Added
- `bermudan_max_call_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American max-call by common-random-number bumps on
  `bermudan_max_call_lsm`. Repricing at bumped spots with the same seed shares
  the Brownian shocks, so the finite differences are low-variance; the LSM
  regression is re-fit at each bump. Deltas are reliable; the gammas (second
  differences over a re-fit regression) are indicative and need many paths.
- Verified: without dividends the two spot deltas match the exact European
  `rainbow_greeks` deltas within noise; both deltas stay in `(0,1)`; the
  cross-gamma is robustly negative across seeds (the two spots are substitutes
  in a max payoff); reproducible under a fixed seed.

## [1.191.1] - 2026-09-10

### Tests
- Trimmed the `bermudan_max_call_lsm` single-asset-lower-bound fast test from 30
  steps / 20000 paths to 20 / 8000; the bound has a ~7-point margin, so the
  smaller sample is comfortably safe. Fast gate back to ~51s.

## [1.191.0] - 2026-09-10

### Added
- `bermudan_max_call_lsm` (in `lsm.py`): American call on the maximum of two
  assets by Longstaff-Schwartz -- the classic two-asset early-exercise
  benchmark. Simulates two correlated GBMs and regresses the continuation value
  on a quadratic basis in both spots plus the running max
  `{1, S1, S2, S1^2, S2^2, S1 S2, max(S1,S2)}` over the in-the-money paths at
  each date.
- Verified: without dividends it matches the European Stulz
  `best_of_call_closed` (little early-exercise value); with a 6% dividend on
  both assets it exceeds the European by a positive early-exercise premium; it
  dominates a single-asset vanilla and is reproducible under a fixed seed.

## [1.190.1] - 2026-09-10

### Tests
- Trimmed the fast suite: the three `heston_mc_greeks` structural checks (each
  runs 9 QE re-prices per Greek set, ~14s combined) drop to smaller sizes -- the
  reproducibility check to 12 steps / 1500 paths and the two positivity/interval
  checks to 15 / 3000. They only assert determinism, a sign, or a `(0,1)` range,
  not accuracy; the Fourier-FD accuracy cross-checks stay under `-m slow`.

## [1.190.0] - 2026-09-10

### Added
- `rainbow_greeks` (in `multiasset.py`): Greeks of a best-of/worst-of rainbow
  option by finite differences on the exact Stulz closed forms (no Monte Carlo
  noise). Returns the two spot deltas, the two own-gammas, the cross-gamma
  `d2V/dS1 dS2`, and the correlation sensitivity `dV/drho`, for calls and puts.
- Verified: the best-of and worst-of call deltas in each asset sum to the
  single-asset Black-Scholes delta (differentiating the Stulz identity
  `C_max + C_min = c1 + c2`), to 1e-4; the max-call `corr_vega` is negative and
  the min-call's positive with the two summing to zero; own-gammas positive;
  worst-of put deltas negative.

## [1.189.0] - 2026-09-10

### Added
- `heston_mc_greeks` (in `heston_mc.py`): Heston Greeks by common-random-number
  finite differences on the QE Monte Carlo. Repricing at bumped inputs with the
  same seed shares the random draws, so the bumped difference is dominated by the
  true sensitivity rather than MC noise. Returns `price`, `delta`, `gamma`,
  `vega_v0` (initial-variance), `vega_theta` (long-variance), `volvol` (dV/dxi),
  and `rho_sens` (dV/drho) -- the sensitivities that matter for a stochastic-vol
  book.
- Verified: at 120000-180000 paths every Greek matches a finite difference of the
  exact Fourier `heston_price` (delta 0.703, gamma 0.0191, vega_v0 47.6,
  vega_theta 48.6, volvol -2.08 vs -2.17, rho 0.20); delta stays in (0,1), gamma
  and both variance vegas are positive.

## [1.188.0] - 2026-09-10

### Added
- `two_asset_gap_option` (in `multiasset.py`): a gap option on asset 1 gated by
  asset 2, separating the trigger strike from the payoff strike. The asset-1 gap
  payoff `(S1_T - K_payoff) 1[S1_T > K_trigger]` (call) fires only if asset 2
  clears its barrier. Decomposes into the two two-asset digitals with the
  trigger strike setting the asset-1 condition and the payoff strike scaling the
  cash leg, so it is an exact bivariate-normal closed form.
- Verified: matches a correlated-GBM Monte Carlo for call and put; reduces to
  `correlation_option` when `K_payoff = K_trigger`; and a larger payoff strike
  lowers the call.

## [1.187.0] - 2026-09-10

### Added
- `correlation_option` (in `multiasset.py`): a two-asset correlation option -- a
  vanilla on asset 1 that pays only if asset 2 satisfies a barrier condition
  (`above`/`below` K2). Decomposes exactly into the two two-asset digitals:
  `call = AoN(S1>K1, cond2) - K1 CoN(S1>K1, cond2)` and
  `put = K1 CoN(S1<K1, cond2) - AoN(S1<K1, cond2)`, where `AoN` is
  `two_asset_asset_or_nothing` and `CoN` is `two_asset_digital`, so it is a
  closed form.
- Verified: matches a correlated-GBM Monte Carlo for call, put, and a below
  barrier; reduces to the plain Black-Scholes vanilla on S1 when the barrier is
  always met; and the `above` + `below` gated calls sum to that vanilla (1e-9).

## [1.186.0] - 2026-09-10

### Added
- `two_asset_asset_or_nothing` (in `multiasset.py`): asset-or-nothing digital
  paying the first asset's terminal value `S1_T` iff both conditions hold (each
  `above`/`below` its strike). Priced under the asset-1 (share) measure, where
  asset 1's drift gains `sigma1^2` and asset 2's shock inherits an extra
  `rho sigma1 sqrt(t)`: `S1 e^{-q1 t} M(s1 a1, s2 a2; s1 s2 rho)`.
- Verified: all four quadrants match a correlated-GBM Monte Carlo; they sum to
  the discounted forward `S1 e^{-q1 t}` (asset 1 is always delivered on some
  quadrant) to 1e-8, including under a dividend yield; and at `rho = 0` the
  price factorizes into the single-asset `asset_or_nothing` on S1 times the
  risk-neutral probability that S2 clears K2.

## [1.185.0] - 2026-09-10

### Added
- `two_asset_digital` (in `multiasset.py`): exact closed form for a cash-or-
  nothing digital on two correlated assets. Pays `cash` iff both single-asset
  conditions hold (each `above` or `below` its strike). Price is
  `cash e^{-rt} M(s1 d1, s2 d2; s1 s2 rho)` with `di` the usual `d2`,
  `si = +/-1` for above/below, and `M` the standardized bivariate-normal CDF;
  flipping a condition flips that `d`'s sign and the correlation.
- Verified: all four quadrant prices match a correlated-GBM Monte Carlo, sum to
  `e^{-rt}` (exhaustive) to 1e-9, factorize into the product of two single-asset
  `cash_or_nothing` digitals at `rho = 0`, and scale linearly in `cash`.

## [1.184.0] - 2026-09-10

### Added
- `best_of_put_closed` / `worst_of_put_closed` (in `multiasset.py`): exact
  closed forms for rainbow puts on the maximum / minimum of two assets. By
  rainbow put-call parity `P = C - disc E[chosen] + K e^{-rt}`, built on the
  call closed forms and `disc E[min] = S2 e^{-q2 t} - exchange_option(S2, S1)`
  (`min(a,b) = b - max(b-a,0)`), with `disc E[max]` the two forwards less the
  discounted expected min.
- Verified: the put identity `P_min + P_max = p(S1) + p(S2)` holds to 1e-9;
  both match their Monte Carlo counterparts within MC error; put-on-min exceeds
  put-on-max, and put-on-max is below each single-asset vanilla put.

## [1.183.0] - 2026-09-10

### Added
- `best_of_call_closed` / `worst_of_call_closed` (in `multiasset.py`): exact
  Stulz (1982) closed forms for rainbow calls on the maximum / minimum of two
  assets, replacing Monte Carlo (`best_of_call` / `worst_of_call`) with an
  analytic price. Built on the spread vol
  `sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)` and the standardized
  bivariate-normal CDF; the call-on-min is priced directly and the call-on-max
  from the Stulz identity `C_max + C_min = c(S1) + c(S2)`.
- Verified: the identity holds to 1e-9; both match their Monte Carlo
  counterparts within MC error (best 17.150 vs 17.16, worst 4.574 vs 4.58), also
  under negative correlation with dividends; the max-call dominates each vanilla
  and the min-call.

## [1.182.0] - 2026-09-10

### Added
- `basket_option_lhs_mc` (in `montecarlo.py`): two-asset basket option
  `max(w1 S1 + w2 S2 - K, 0)` by Latin hypercube sampling -- the analogue of
  `spread_option_lhs_mc` for a weighted-sum payoff. Each driving normal is
  stratified into `n_paths` equiprobable bins with independently permuted bin
  orders, mapped through the inverse normal CDF; asset 2's shock is correlated
  by `corr z1 + sqrt(1 - corr^2) z2`.
- Verified: agrees with the Levy moment-matched (approximate) `basket_option`
  for call and put (to a few cents), reduces to a Black-Scholes call when all
  weight is on one asset, and -- wrapped in `replicated_mc` for an honest
  across-seed SE -- places the moment-match price within a few SEs. As with any
  LHS estimator the per-run `std_error` overstates the true error; documented.

## [1.181.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~50s -> ~47s. The Heston-CV standard-error-ordering
  check drops from 30000 paths / 100 steps to 8000 / 30 (a robust inequality,
  not an accuracy claim), and the dividend-yield Fourier match moved under
  `-m slow` with the other `heston_cv_mc` accuracy cross-checks. Together these
  removed ~19s of QE simulation from the fast gate.

## [1.181.0] - 2026-09-10

### Added
- `heston_cv_mc` (in `heston_mc.py`): Heston QE Monte Carlo with the underlying
  as a control variate. Andersen's QE step is martingale-corrected, so the
  discounted terminal spot `Y = e^{-rt} S_T` has the known mean `S0 e^{-qt}`;
  the controlled estimator `X - beta (Y - E[Y])` with the regression-optimal
  `beta = Cov(X,Y)/Var(Y)` cuts the standard error at no bias. Same QE variance
  step, `K0..K4` asset constants, and antithetic draws as `heston_qe_mc`.
- Verified: matches the Fourier `heston_price` for call and put (also under a
  dividend yield) within MC error; at equal path count the standard error is
  ~0.56x the plain `heston_qe_mc` (the spot/payoff correlation is weaker under
  stochastic vol than in Black-Scholes, so the reduction is milder than the
  lognormal `european_cv_mc`).

## [1.180.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~50s -> ~46s. The two RQMC-Asian cross-checks use a
  120000-path (was 400000) control-variate reference and n_rand=16 (was 24); the
  vol-surface example builds its SVI calibration once via a module-scoped fixture
  instead of three times; and the local-vol-MC positivity/smoke checks drop from
  100 steps / 20000 paths to 40 / 8000. Accuracy cross-checks stay under `-m slow`.

## [1.180.0] - 2026-09-10

### Fixed
- `arithmetic_asian_mc` control variate was biased. It used the *continuous*
  geometric-Asian closed form (`geometric_asian`, `sigma/sqrt(3)`) as the
  control, but the Monte Carlo averages over `n_steps` *discrete* dates, so the
  control's known expectation did not match the simulated geometric payoff. The
  CV price was biased low -- e.g. 5.76 vs a true ~6.55 at `n_steps=6` -- while
  reporting a tiny SE. Existing tests used `n_steps=50`/`250` (near-continuous)
  and a loose `3*plain_SE` band, so they missed it. Added
  `_discrete_geometric_asian` (exact closed form of the discretely-monitored
  geometric average, `E[ln G] = ln S0 + (b-sig^2/2) dt (n+1)/2`,
  `Var[ln G] = sig^2 dt (n+1)(2n+1)/(6n)`) and use it as the control. CV and
  no-CV prices now agree at every `n_steps`; a tight `n_steps=6` regression test
  guards it.

### Added
- `sobol_asian_rqmc` (in `sobol.py`): randomized-QMC arithmetic Asian with an
  honest standard error -- the multi-dimensional analogue of
  `sobol_european_rqmc`. Each path's normals come from an `n_steps`-dimensional
  Sobol point through the Brownian bridge; a per-dimension Cranley-Patterson
  rotation randomizes the point set, so `n_rand` shifts give i.i.d. QMC estimates
  whose spread is a genuine SE. Cross-checks the (now-fixed) control-variate
  `arithmetic_asian_mc` for call and put; ~0.06x the plain-MC SE at equal points.

## [1.179.0] - 2026-09-10

### Added
- `sobol_european_rqmc` (in `sobol.py`): randomized-QMC European price with an
  honest standard error. Plain Sobol QMC returns a single deterministic number
  with no error bar; this applies a Cranley-Patterson rotation -- shifting the
  whole Sobol point set by a random `U ~ Uniform[0,1)` mod 1 -- which preserves
  the low discrepancy but makes each estimate an unbiased draw. `n_rand`
  independent shifts give i.i.d. QMC estimates whose spread is a genuine SE.
  Returns the mean price, the across-randomization SE, and `n_paths` = total
  points (`n_rand * n_paths`).
- Verified: matches the closed-form Black-Scholes call and put (also under a
  dividend carry); at equal total points the RQMC standard error is ~0.05x the
  plain `european_mc` for this smooth 1-D integral.

## [1.178.0] - 2026-09-10

### Added
- `replicated_mc` (in `montecarlo.py`): batched-replication honest standard
  error for variance-reduced estimators. Stratified, Latin hypercube, and QMC
  estimators draw dependent samples, so the plain i.i.d. `std_error` they report
  is not the truth (see the note on `spread_option_lhs_mc`). Running the whole
  estimator `n_batches` times with distinct seeds gives independent batch means
  whose spread is an unbiased standard error. Takes any callable
  `seed -> MCResult | float`; returns the mean price, the across-batch SE
  `s/sqrt(n_batches)`, and `n_paths = n_batches`.
- Verified: for a plain `european_mc` the across-batch SE matches a single-run
  i.i.d. SE / sqrt(n_batches) (0.032 vs 0.031, calibrated); the mean is unbiased
  vs Black-Scholes for a stratified estimator; and it exposes the LHS spread
  reduction (LHS honest SE < 0.6x a plain two-asset MC's honest SE) that the
  naive per-run SE hides.

## [1.177.0] - 2026-09-10

### Added
- `spread_option_lhs_mc` (in `montecarlo.py`): two-asset spread option
  `max(S1 - S2 - K, 0)` by Latin hypercube sampling. Each of the two driving
  normals is stratified into `n_paths` equiprobable bins with one draw per bin,
  the two dimensions' bin orders are independently permuted, and the stratified
  uniforms are mapped through the inverse normal CDF; the second asset's shock
  is correlated by `rho z1 + sqrt(1 - rho^2) z2`.
- Verified: matches the Kirk `spread_option` and, at `K = 0`, the exact Margrabe
  `exchange_option`. The LHS samples are not independent, so the reported
  `std_error` (plain i.i.d. formula) overstates the true error; the genuine gain
  is in the across-seed spread -- at 4000 paths the RMSE vs Kirk is ~0.10 versus
  ~0.27 for a plain two-asset draw (~2.6x), documented on the function.

## [1.176.0] - 2026-09-10

### Added
- `european_stratified_mc` (in `montecarlo.py`): European price by stratified
  sampling of the terminal normal. The driving normal is split into `n_strata`
  equiprobable strata in probability space; `n_per` uniforms are drawn within
  each and mapped through the inverse normal CDF, spreading draws evenly and
  removing the clustering that inflates plain MC variance. With equiprobable
  strata the estimate is the average of the per-stratum means and its variance
  is `(1/n_strata^2) sum_i s_i^2 / n_per`.
- Verified: matches the closed-form Black-Scholes call and put (also deep OTM
  and under a dividend carry); at 10000 paths (400 strata x 25) the standard
  error is ~0.05x the plain `european_mc` for the smooth call payoff.

## [1.175.0] - 2026-09-10

### Added
- `european_is_adaptive_mc` (in `montecarlo.py`): importance sampling with a
  pilot-tuned optimal shift. The strike-centring shift of `european_is_mc` is
  near-optimal for a digital but not for a vanilla, whose payoff keeps growing
  past the strike and pulls the best shift further OTM. A change of measure
  gives the shifted second moment from plain draws,
  `M(mu) = E_0[payoff(z)^2 exp(-mu z + mu^2/2)]`, so one `N(0,1)` pilot sample
  scores every candidate `mu` on a grid at negligible cost; the main run samples
  at the variance-minimising `mu` and stays unbiased via the likelihood ratio.
- Verified: matches the closed-form Black-Scholes call and put (deep OTM, ATM,
  under carry); for a K=160 deep-OTM call the tuned shift gives ~0.87x the
  standard error of the fixed strike-centring shift at equal main-run paths.

## [1.174.0] - 2026-09-10

### Added
- `digital_is_mc` (in `montecarlo.py`): cash-or-nothing digital price by
  importance sampling. A deep-OTM digital is even harder to simulate plainly
  than a vanilla -- the 0/`cash` indicator has relative SE that blows up like
  `sqrt((1-p)/p)` for a small hit probability `p`. Shifting the terminal normal
  to `N(mu, 1)` and reweighting by `L(z) = exp(-mu z + mu^2/2)` moves paths into
  the money unbiased; the default shift puts the mean draw on the strike
  boundary so about half the paths pay.
- Verified: matches the closed-form `cash_or_nothing` (deep OTM call and put,
  ATM, cash scaling, dividend carry); for a K=160 deep-OTM digital the standard
  error is ~0.17x a plain indicator estimator at equal paths.

## [1.173.0] - 2026-09-10

### Added
- `european_is_mc` (in `montecarlo.py`): European price by importance sampling,
  for deep out-of-the-money options where a plain simulation wastes almost every
  path. Draws the terminal normal from a shifted mean `N(mu, 1)` to push mass
  into the money and corrects with the likelihood ratio
  `L(z) = exp(-mu z + mu^2/2)`. The default shift centres the terminal log-spot
  on the strike, `mu* = (ln(K/S0) - (b - sig^2/2) t) / (sig sqrt(t))`.
- Verified: matches the closed-form Black-Scholes call and put (deep OTM and
  ATM, under a dividend carry, and with an explicit shift); for a K=160 deep-OTM
  call the standard error is ~0.086x the plain `european_mc` at equal paths.

## [1.172.0] - 2026-09-10

### Added
- `lr_digital_greeks` (in `mc_greeks.py`): delta, vega, and gamma of a
  cash-or-nothing digital by the likelihood-ratio method. The digital payoff is
  discontinuous, so the pathwise method is undefined for every Greek; the LR
  method differentiates the log-normal density, so the same one-step
  Black-Scholes score weights as `lr_greeks` (delta `Z/(S0 sig sqrt t)`, vega
  `(Z^2-1)/sig - Z sqrt t`, gamma `(Z^2 - Z sig sqrt t - 1)/(S0^2 sig^2 t)`)
  apply unchanged to the digital.
- Verified: delta and gamma match the analytic `digital_greeks`, vega matches a
  sigma-bump of `cash_or_nothing`, within MC error for both call and put (e.g.
  call delta 0.0274 vs 0.0274, vega -0.478 vs -0.479, gamma -0.00048 vs
  -0.00048); call delta positive, put delta negative.

## [1.171.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~49s -> ~35s. The LSM-Greeks sign/price-field checks
  drop from 25 steps / 8000 paths to 20 / 4000, and the Asian pathwise-vega
  positivity check from 40 steps / 40000 paths to 20 / 15000. These assert
  only a sign or an exact price-field match, not accuracy; the accuracy
  cross-checks stay under `-m slow`, unchanged.

## [1.171.0] - 2026-09-10

### Added
- `barrier_lr_delta` (in `mc_greeks.py`): delta of a discretely-monitored
  single-barrier option by the likelihood-ratio method. The knock-out/knock-in
  payoff is discontinuous in the spot, so the pathwise method fails; the LR
  method uses the score of the path density. The initial spot enters only
  through the mean of the first log-increment, so the score is
  `Z_1 / (S0 sig sqrt(dt))` and `delta = E[disc_payoff * Z_1 / (S0 sig sqrt(dt))]`.
- Verified: matches a common-random-number finite difference of the
  discretely-monitored `barrier_mc` (brownian_bridge=False) within MC error for
  a down-out call (0.784 vs 0.788) and an up-in put (0.089 vs 0.092); down-out
  call delta is positive.

## [1.170.0] - 2026-09-10

### Added
- `european_cv_mc` (in `montecarlo.py`): European Monte Carlo combining both
  variance-reduction techniques -- antithetic sampling and a control variate.
  The discounted terminal spot `Y = e^{-rt} S_T` has the known mean
  `S0 e^{(b-r)t}` and is correlated with the payoff, so the estimator
  `X - beta (Y - E[Y])` with the regression-optimal `beta = Cov(X,Y)/Var(Y)`
  cuts the variance. Antithetic pairs are averaged into one sample before the
  control is applied.
- Verified: matches the closed-form Black-Scholes call and put within MC error
  (also under a dividend carry `b = r - q`); at equal path count the standard
  error is ~0.24x the plain `european_mc` for the ITM call and ~0.42x for the
  put, and under half at the money.

## [1.169.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~49s -> ~35s: the two Heston pathwise-delta
  sign/range checks run at n_steps=40 / n_paths=6000 instead of 50/20000 (1.8s
  each -> under 0.7s). The Fourier-FD accuracy cross-checks stay under `-m slow`.

## [1.169.0] - 2026-09-10

### Added
- `heston_pathwise_delta` (in `heston_mc.py`): Heston delta by the pathwise
  method. In the Andersen-QE simulation the initial spot enters only as the
  additive `ln S0` in the terminal log-price, so `S_T = S0 e^Y` with `Y`
  independent of `S0` and the pathwise delta is exact:
  `delta = e^{-rt} E[1_{S_T>K} S_T/S0]` (put: `-1_{S_T<K}`).
- Verified: the call and put deltas match a finite-difference of the exact
  Fourier `heston_price` within Monte Carlo error; the call delta is in (0, 1)
  and the put delta is negative.

## [1.168.0] - 2026-09-10

### Added
- `smoothed_digital_delta` (in `mc_greeks.py`): pathwise delta of a
  cash-or-nothing digital via a call-spread smoothing of the discontinuous
  indicator (ramp of relative width `eps_rel`), making the payoff Lipschitz so
  the pathwise derivative is well-defined. A single-pass, model-agnostic
  alternative to the likelihood-ratio digital delta.
- Verified: it matches the analytic digital delta; a narrower spread lowers the
  smoothing bias (converging as `eps_rel -> 0`, at the cost of higher variance
  from the `1/eps` ramp -- the classic bias/variance tradeoff, documented); call
  delta positive, put delta negative.

## [1.167.1] - 2026-09-10

### Tests
- Trimmed the fast suite: the two fast LSM-Greeks checks (put-delta sign, price
  field consistency) run at n_steps=25 / n_paths=8000 instead of 40/40000 (9s
  each -> under 2s), cutting ~14s off the fast run. The delta/gamma accuracy
  cross-checks against the binomial tree still run under `-m slow`.

## [1.167.0] - 2026-09-10

### Added
- `bermudan_lsm_greeks` (in `lsm.py`): delta and gamma of a Bermudan/American
  LSM price by common-random-number bumps -- the option is repriced at
  `S(1 +/- h)` on the same seeded path stream, so the finite differences are
  low-variance. Returns price, delta and gamma.
- Verified: delta matches a 2000-step binomial tree (~1e-2); gamma is positive
  and in a loose band around the tree gamma (a second difference over a re-fit
  regression is only indicative and needs many paths -- documented). Delta
  converges to the tree value as paths grow (-0.40 -> -0.412).

## [1.166.0] - 2026-09-10

### Added
- `mixed_gamma` (in `mc_greeks.py`): European gamma by the mixed
  pathwise-likelihood-ratio estimator. It differentiates the pathwise delta
  payoff `D = e^{-rt} 1_{S_T>K} S_T/S0` w.r.t. `S0` -- through both the density
  (LR weight `Z/(S0 sigma sqrt t)`) and the explicit `1/S0` factor -- giving
  `gamma = E[D (Z/(S0 sigma sqrt t) - 1/S0)]`. Well-defined despite the
  indicator (pure pathwise gamma is not) and lower variance than a double
  likelihood-ratio.
- Verified: matches the Black-Scholes gamma; its standard error is ~4x smaller
  than the LR gamma at equal paths; put gamma equals call gamma; positive.

## [1.165.0] - 2026-09-10

### Added
- Multi-level Monte Carlo (new `mlmc.py`): `mlmc_asian` prices a fixed-strike
  arithmetic Asian by the Giles (2008) telescoping estimator
  `E[P_L] = E[P_0] + sum_l E[P_l - P_{l-1}]`, with each correction from coupled
  fine/coarse paths sharing Brownian increments (each coarse step sums `M` fine
  ones). Deeper levels use fewer paths since their corrections have lower
  variance, cutting the cost to a target accuracy.
- Verified: the level corrections decay geometrically; the estimate converges
  to the Turnbull-Wakeman value as levels grow (err ~0.01 at 6 levels) and more
  levels reduce the discretisation bias; call and put prices are positive.

## [1.164.0] - 2026-09-10

### Added
- `asian_pathwise_vega` (in `mc_greeks.py`): pathwise vega of a fixed-strike
  arithmetic-average Asian option -- differentiating the payoff along each path
  w.r.t. `sigma`, with `dS_i/dsigma = S_i (W_i - sigma t_i)` and
  `dA/dsigma = mean_i dS_i/dsigma`. A lower-variance alternative to bumping for
  this path-dependent, Lipschitz payoff.
- Verified: it matches a central-difference bump of the arithmetic-Asian Monte
  Carlo within MC error, and agrees with a finite-difference of the
  Turnbull-Wakeman continuous-average closed form; call and put vega are
  positive.

## [1.163.0] - 2026-09-10

### Added
- `density_var_es` (new `density_var.py`): risk-neutral Value-at-Risk and
  Expected Shortfall of an option position whose horizon P&L is `pnl(S_T)`,
  computed by integrating the smile-implied density into a P&L distribution and
  taking the loss quantile and tail mean. VaR/ES returned as positive losses.
- Verified: a long-forward VaR matches the lognormal quantile; ES >= VaR; a long
  call's VaR is capped at the premium (max loss); a short call has a positive
  VaR with a larger ES; higher confidence gives a larger VaR.
- Bug fixed during validation: the ES accumulation broke on the first
  zero-density deep-tail atom (returning ~0); it now fills the whole tail mass
  regardless of individual atom probabilities.

## [1.162.0] - 2026-09-10

### Added
- `wasserstein_smiles` (in `density_metrics.py`): the 1-Wasserstein
  (earth-mover) distance between two smile-implied densities, computed as the L1
  gap between their CDFs `integral |F_p - F_q| dK` on a shared strike grid.
  Unlike the KL divergence it is a true metric (symmetric, triangle inequality)
  measured in price units.
- Verified: identical smiles give ~0; it is symmetric; positive for differing
  vols and growing with the vol gap; a skewed smile is a positive distance from
  a flat one; the triangle inequality holds.

## [1.161.0] - 2026-09-10

### Added
- `kl_divergence_smiles` (in `density_metrics.py`): the Kullback-Leibler
  divergence `KL(g_p || g_q) = integral g_p ln(g_p/g_q) dK` between two
  smile-implied risk-neutral densities (the second interpolated onto the first's
  grid). Measures how far one implied distribution sits from another -- two
  dates, two models, or implied vs a reference.
- Verified: identical smiles give ~0; differing vols give a positive,
  asymmetric, non-negative divergence that grows with the vol gap; a skewed
  smile diverges positively from a flat one.

## [1.160.0] - 2026-09-10

### Added
- Risk-neutral density tail/shape metrics (new `density_metrics.py`):
  `tail_probability` (`Q(S_T < L)` / `> U`, the undiscounted digital price),
  `density_entropy` (differential entropy of the terminal-spot density) and
  `expected_shortfall` (`E^Q[S_T | S_T` in a tail]), all from the
  Breeden-Litzenberger smile density.
- Verified: the upper-tail probability matches the lognormal `N(d2)`; the two
  tails sum to 1; a downward skew fattens the left tail; entropy increases with
  vol; the tail expected-shortfalls sit beyond their thresholds; an empty tail
  returns NaN.

### Fixed
- `density_grid_from_smile` could evaluate a call at a non-positive strike at the
  low end of a wide/high-vol grid (the finite-difference step exceeded the first
  strike). The FD step is now capped at `0.5 * K` so `K - h > 0` everywhere.

## [1.159.0] - 2026-09-10

### Added
- Breeden-Litzenberger risk-neutral density from a smile (new `rnd.py`):
  `risk_neutral_density_from_smile` gives `g(K) = e^{r t} d^2 C/dK^2`,
  `density_grid_from_smile` returns the density on a strike grid, and
  `price_payoff_from_density` prices any European payoff model-free by
  integrating it against the density.
- Verified: the density integrates to 1; it reprices calls, puts and digitals
  to Black-Scholes; the payoff `K` recovers the spot (discounted forward); the
  ATM density is positive; and a concave (arbitraging) smile produces a negative
  density in the wings.

## [1.158.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime from ~45s to ~30s: the variance-term-structure
  tests build the strip at n_strikes=101, the LSM local-vol American>=European
  check runs at n_steps=25/n_paths=8000, the rBergomi auto/forced identity uses
  150 paths, and the Levy VG parameter recovery is marked slow. All accuracy
  claims still run under `-m slow`.

## [1.158.0] - 2026-09-10

### Added
- `variance_term_structure` (in `varswap.py`): the term structure of
  variance-swap strikes across expiries plus the forward (instantaneous)
  variance curve between them. Returns `(spot_var, forward_var)` where
  `forward_var[i] = (K_i t_i - K_{i-1} t_{i-1})/(t_i - t_{i-1})` by
  total-variance additivity.
- Verified: a flat term structure gives a constant curve (to the ~1e-6 strip
  truncation level); a rising vol term structure gives positive, increasing
  forward variances; the forward curve matches
  `forward_variance_swap_from_smile`; and `forward_var[0]` equals the first spot
  strike.

## [1.157.0] - 2026-09-10

### Added
- `moment_risk_premia` (new `moment_premium.py`): skewness- and
  kurtosis-risk premia -- the risk-neutral (BKM) skewness/excess-kurtosis the
  option smile implies versus what realized in the price history. Returns
  realized, implied and premium (implied - realized) for both moments, wiring
  `bkm_moments_from_smile` to the realized-moment estimators.
- Verified: a downward-skewed smile against a near-symmetric GBM history gives a
  negative skew premium (crash-protection demand) and a positive kurtosis
  premium; a flat smile gives a near-zero skew premium; the components are
  internally consistent.

## [1.156.0] - 2026-09-10

### Added
- `equity_premium_lower_bound` (in `vix.py`): Martin's (2013) model-free lower
  bound on the expected equity excess return, `Rf * SVIX^2`, computed from the
  option smile via the simple-variance index. Under the negative-correlation
  condition the annualized expected market excess return is bounded below by
  this quantity.
- Verified: a flat 20% vol gives a bound near `Rf * sigma^2`; a higher vol
  raises it; and it equals `Rf` times the SVIX variance exactly.

## [1.155.0] - 2026-09-10

### Added
- `bermudan_lsm_local_vol` (in `lsm.py`): Longstaff-Schwartz American/Bermudan
  pricing under a local-volatility surface. Same backward-induction regression
  as `bermudan_lsm`, but each Euler step uses `local_vol_fn(S, tau)` -- so it
  prices early-exercise options directly on a calibrated Dupire/SVI local-vol
  surface. Carry `b = r - q`.
- Verified: a flat `local_vol_fn` reproduces the constant-vol LSM price and a
  binomial tree; a genuinely skewed local vol matches an independent
  Crank-Nicolson American PDE (~0.03).

## [1.154.0] - 2026-09-10

### Added
- Monte Carlo Greeks without bumping (new `mc_greeks.py`): `lr_greeks` returns
  European delta, gamma and vega by the likelihood-ratio (Malliavin-flavoured)
  method -- `E[payoff * weight]` with density-derivative weights -- which works
  even for discontinuous payoffs; `pathwise_delta` is the lower-variance
  pathwise estimator for smooth payoffs; and `lr_digital_delta` gives the delta
  of a cash-or-nothing digital, where the pathwise method fails.
- Verified: LR delta/gamma/vega match Black-Scholes within Monte Carlo error;
  the pathwise delta matches and has a smaller standard error than the LR delta;
  and the LR digital delta matches the analytic (finite-difference) digital
  delta.

## [1.153.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime from ~44s to ~32s: the Kim-Greeks structural
  checks (price equals the direct price, sign checks) run at n_steps=50 instead
  of 120; the Bermudan-swaption sign/monotonicity checks at n_paths=6000 instead
  of 20000; and the Heston-calibration parameter-validity check is marked slow.
  All accuracy claims still run under `-m slow`.

## [1.153.0] - 2026-09-10

### Added
- `svix_from_smile` (in `vix.py`): Martin's (2013) simple-variance index (SVIX).
  Unlike the VIX log-contract (`1/K^2` weights), the simple variance swap weights
  the OTM strip by the constant `1/F^2`, corresponding to the payoff
  `(S_T - F)^2 / F^2` with no log approximation -- jump-robust and a genuine
  lower bound on the equity premium. Reported as `100 * SVIX`.
- Verified: a flat smile returns ~`100 * sigma`; `SVIX = 100 * sqrt(variance)`;
  it is positive; and it differs from the VIX under a skew (the two coincide only
  to leading order).

## [1.152.0] - 2026-09-10

### Added
- Variance risk premium (new `vrp.py`): `realized_variance` computes annualized
  realized variance from a close series, and `variance_risk_premium` combines it
  with a supplied implied (variance-swap) variance into the additive premium
  (realized - implied), the ratio, and the vol-point premium
  (implied vol - realized vol).
- Verified: realized variance recovers `sigma^2` from a GBM path; implied above
  realized gives a negative VRP with a positive vol premium and ratio < 1;
  realized above implied flips the signs; the components are internally
  consistent.

## [1.151.0] - 2026-09-10

### Added
- Bakshi-Kapadia-Madan risk-neutral moments (new `bkm.py`):
  `bkm_moments_from_smile` extracts the model-free risk-neutral variance,
  skewness and excess kurtosis of the log-return from a smile `vol_fn(K)` via the
  quadratic/cubic/quartic option-spanning contracts (strike-weighted OTM
  strips). `skew_swap_from_smile` returns the risk-neutral skewness (a skew
  swap's fair value).
- Verified: a flat smile gives variance `sigma^2 t`, zero skewness and zero
  excess kurtosis; a downward (equity) skew gives negative skewness and positive
  excess kurtosis; an upward skew gives positive skewness; a steeper skew is more
  negative.

## [1.150.0] - 2026-09-10

### Added
- CBOE VIX-style fair volatility index (new `vix.py`): `vix_from_chain`
  implements the exact discrete CBOE formula -- the `dK/K^2`-weighted OTM strip
  minus the `(F/K0 - 1)^2` forward-correction term, reported as
  `100 * sqrt(variance)` -- and `vix_from_smile` builds the chain from a smile
  `vol_fn(K)`.
- Verified: a flat smile returns `VIX = 100 * sigma` exactly; a downward skew
  lifts the index above the ATM level (the fear premium); a manually-built
  flat-vol chain matches; and `VIX = 100 * sqrt(variance)`.

## [1.149.0] - 2026-09-10

### Added
- `forward_variance_swap_from_smile` (in `varswap.py`): fair forward-start
  variance-swap strike over `[t1, t2]` from the two expiries' smiles. Total
  variance is additive in time, so the accrued variance is
  `(K_var(t2) t2 - K_var(t1) t1)/(t2 - t1)` with each spot-starting leg replicated
  from its smile.
- Verified: a flat term structure returns `sigma^2`; `t1 = 0` recovers the spot
  variance swap exactly; an upward-sloping vol term structure gives a forward
  variance above the front; and the additivity identity holds to 1e-6.

## [1.148.0] - 2026-09-10

### Added
- `gamma_swap_from_smile` (in `varswap.py`): fair gamma-swap (price-weighted
  variance) strike from a smile `vol_fn(K)`. A gamma swap accrues `(S_t/S0)`-
  weighted realized variance, so by Carr-Lewis it is replicated by an option
  strip weighted `1/K` (vs the variance swap's `1/K^2`), scaled by
  `2 e^{rt}/(S0 t)`, with no forward remainder (the price-weighted log contract
  has none).
- Verified: a flat smile returns `sigma^2`; a denser strip reduces the error; a
  downward skew makes the gamma swap worth less than the same smile's variance
  swap (spot-weighting down-weights the high-vol low-strike puts).

### Fixed
- An initial gamma-swap formula carried a spurious forward "drift" term (copied
  from the variance-swap log contract); the price-weighted log contract has no
  such remainder. Removing it makes the flat case equal `sigma^2` and restores
  the correct gamma < variance ordering under a downward skew.

## [1.147.0] - 2026-09-10

### Added
- `corridor_variance_swap_from_smile` (in `varswap.py`): fair corridor
  variance-swap strike from a smile `vol_fn(K)`. A corridor variance swap
  accrues realized variance only while the spot is in `[lower, upper]`; by
  Carr-Lewis static replication this restricts the `1/K^2`-weighted option strip
  to strikes inside the corridor. Each strip option is priced at its smile vol
  with Black-Scholes.
- Verified: a narrower corridor accrues less variance; nested corridors are
  monotone at fixed strike density; a wide corridor recovers the full
  variance-swap strike (~5e-3); the strike is positive.

## [1.146.0] - 2026-09-10

### Added
- `variance_swap_from_smile` (in `varswap.py`): fair variance-swap strike
  replicated directly from a volatility smile `vol_fn(K)`. Builds the OTM
  option strip (puts below the forward split, calls above) by pricing each
  strike at its smile vol with Black-Scholes and delegates to
  `variance_swap_strike`. Convenient for marking a swap off a fitted smile
  (SVI/SABR/vanna-volga).
- Verified: a flat smile returns its variance `sigma^2` (to ~3e-3, converging as
  the strip widens/densifies -- the residual is finite-strip truncation); a
  downward skew adds convexity variance over the same-grid flat level; and it
  agrees with a manually-built flat-vol chain.

## [1.145.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime back from ~45s to ~37s: the smooth `K -> 0`
  Asian-PDE sanity check runs on a modest grid (it needs no strike-kink
  resolution), the ATM/OTM Asian-PDE put check is marked slow, and the
  Heston-calibration parameter-validity check uses a short optimiser run (it
  tests the reparametrization, not fit accuracy). Full accuracy coverage is
  unchanged under `-m slow`.

## [1.145.0] - 2026-09-10

### Added
- `asian_pde_price` (new `pde_asian.py`): continuously-averaged arithmetic Asian
  option by an augmented-state 2D PDE in spot and running integral. The integral
  state is pure transport (`dI = S dt`, no diffusion), so each backward step is
  operator-split: a semi-Lagrangian transport in `I` (interpolate along the
  characteristic) then a Crank-Nicolson diffusion in `S` on each `I`-line. The
  `I`-grid is bounded by the average's reachable range so it resolves the strike
  kink.
- Verified: a `K -> 0` average-price call equals the discounted expected average
  (~0.05); prices agree with Turnbull-Wakeman to ~0.25 and converge toward a
  continuous-monitoring Monte Carlo as the `I`-grid refines.

## [1.144.0] - 2026-09-10

### Added
- `adi_two_asset_american` (in `pde2d.py`): two-asset options with optional
  early exercise on the Peaceman-Rachford ADI grid. When `american=True` the
  value grid is floored at the immediate-exercise payoff after each time step
  (the 2D analogue of the vanilla PSOR floor), for American best-of / worst-of /
  spread payoffs.
- Verified: the European limit (`american=False`) matches the best-of and
  worst-of closed forms and converges in the grid; the American value is at
  least the European; a dividend yield produces a positive early-exercise
  premium; and the American value is floored at intrinsic deep in the money.

## [1.143.0] - 2026-09-10

### Added
- `adi_two_asset_cs` (in `pde2d.py`): the Craig-Sneyd ADI scheme for two-asset
  options. Peaceman-Rachford is only first-order in time when a mixed
  (correlation) derivative is present; Craig-Sneyd restores second order with a
  Douglas predictor followed by a corrector that re-applies the explicit cross
  term at the predicted value. The directional operators are solved implicitly
  (Thomas), the cross term explicitly.
- Verified: matches the Margrabe exchange closed form; at few time steps it is
  markedly more accurate than Peaceman-Rachford (0.0004 vs 0.043 at n_time=10),
  its purpose; it converges in space; and the zero-correlation case matches
  Margrabe.

## [1.142.0] - 2026-09-10

### Added
- Two-asset ADI PDE solver (new `pde2d.py`): `adi_two_asset` prices a European
  two-asset option for any terminal payoff by the Peaceman-Rachford
  alternating-direction-implicit scheme in log-prices -- each time step is two
  half-steps (implicit in x1, then x2) with the mixed correlation derivative
  explicit, each half a set of tridiagonal Thomas solves. `adi_spread_option`
  wraps it for `max(S1 - S2 - K, 0)`.
- Verified: the zero-strike spread matches the Margrabe exchange closed form
  (~2e-3, converging in the grid), a non-zero-strike spread matches the Kirk
  approximation (~3e-3), higher asset correlation lowers the spread price
  monotonically, and a generic max-of-two payoff prices sanely.

## [1.141.0] - 2026-09-10

### Added
- `VannaVolgaSmile.vol_cm` (in `vannavolga.py`): the Castagna-Mercurio analytic
  vanna-volga implied-vol approximation from the three delta pillars. `order=1`
  is the vega/vanna/volga-weighted first-order average
  `sigma_atm + sum_i x_i (sigma_i - sigma_atm)`; `order=2` adds the standard
  second-order convexity correction
  `(-sigma_atm + sqrt(sigma_atm^2 + d1 d2 (2 sigma_atm D1 + D2)))/(d1 d2)`.
- Verified: both orders are exact at the three pillars; the second order matches
  the price-corrected vol off-pillar to ~1e-3; the two orders differ off-pillar;
  and a negative risk reversal gives a downward skew.

## [1.140.0] - 2026-09-10

### Added
- Digital and no-touch binaries via the Crank-Nicolson PDE (`pde.py`):
  `crank_nicolson_digital` prices a cash-or-nothing digital (digital-specific
  boundary conditions, Rannacher damping on by default for the payoff jump), and
  `crank_nicolson_no_touch` prices a no-touch binary as a knock-out of a
  constant cash payoff with an absorbing barrier (a node placed exactly on `H`).
  A pay-at-expiry one-touch is `cash * e^{-rt} - no_touch`. Constant `sigma` or
  `local_vol_fn`, carry `b = r - q`.
- Verified: the digital matches the closed-form `cash_or_nothing` (~4e-3),
  digital call + put equals the discount factor exactly, the no-touch matches
  its closed form (~5e-3), and the one-touch complement matches the analytic
  pay-at-expiry one-touch.

## [1.139.0] - 2026-09-10

### Added
- FX delta-space quoting conventions (new `fxdelta.py`): `strike_from_delta` /
  `delta_from_strike` convert between strike and delta in the forward- or
  spot-delta convention, with an optional premium adjustment; `atm_dns_strike`
  gives the delta-neutral-straddle ATM strike `F exp(0.5 sigma^2 t)`; and
  `rr_bf_to_pillars` turns (ATM, 25d risk reversal, 25d butterfly) quotes into
  the put/ATM/call pillar strikes and vols.
- Verified: strike<->delta round-trips to ~1e-16 (forward, spot and
  premium-adjusted); the ATM DNS strike makes a straddle delta-neutral; the
  pillars reproduce the RR and BF quotes exactly and are strike-ordered; a
  negative risk reversal makes the put vol exceed the call vol.

## [1.138.0] - 2026-09-10

### Added
- `calibrate_ssvi` gains a `vega_weighted` option: each quote's total-variance
  error is weighted by an approximate Black vega
  (`sqrt(w) exp(-d1^2/2)`), so near-the-money quotes -- largest vega, deepest
  liquidity -- dominate the fit rather than the deep wings.
- Verified: a clean surface is still recovered exactly; under wing noise the
  vega-weighted fit has a smaller at-the-money vol error than the unweighted
  fit; the fitted surface stays valid and arbitrage-free; and the default
  (unweighted) behaviour is unchanged.

## [1.137.0] - 2026-09-10

### Added
- `calibrate_svi_from_prices` (in `svi.py`): calibrate a raw SVI slice directly
  from market call *prices* rather than pre-inverted vols. Inverts each call to
  its Black-Scholes implied vol, converts to total variance, and fits raw SVI,
  vega-weighting the quotes by default (so a price error in the deep wings,
  where vega is tiny, does not dominate the fit). Returns
  `(params, iv_rmse, price_rmse)`.
- Verified: recovers a synthetic SVI slice from clean prices exactly (iv/price
  RMSE ~0); fits noisy prices; vega-weighted and unweighted both fit clean
  prices; and the fitted slice reprices each strike to within the reported price
  RMSE.

## [1.136.0] - 2026-09-10

### Added
- Rannacher time-stepping in the Crank-Nicolson solver (`pde.py`): the first
  `rannacher` steps out of expiry are taken fully implicit (backward Euler)
  before switching to Crank-Nicolson, damping the spurious oscillation the
  non-smooth payoff kink induces in CN (which otherwise corrupts gamma/theta
  near the strike). Exposed as a `rannacher` argument on `crank_nicolson_price`
  and `crank_nicolson_greeks` (default 2; `0` recovers pure CN).
- The theta-scheme was generalised to a per-step weight so the same code runs
  backward Euler and Crank-Nicolson. Verified: prices still match Black-Scholes;
  `rannacher=0` recovers pure CN; and on a coarse time grid the Rannacher gamma
  near the strike is no worse than (here better than) pure CN.

## [1.135.0] - 2026-09-10

### Added
- `crank_nicolson_barrier` (in `pde.py`): continuously-monitored single-barrier
  options by a Crank-Nicolson PDE with an absorbing boundary -- the value is
  pinned to the rebate on the dead side of the barrier at every time step. The
  grid is anchored so a node lands exactly on `H` (otherwise the barrier is
  applied at an offset node). `down-out`/`up-out` are solved directly;
  `down-in`/`up-in` come from in + out = vanilla. Supports a constant `sigma` or
  a `local_vol_fn` and a carry `b = r - q`.
- Verified: matches the Reiner-Rubinstein closed form on all four barrier kinds
  and a put (within the barrier's O(ds) discretisation bias on a fine grid);
  knock-in + knock-out reproduces the vanilla; a knock-out is below the vanilla;
  spot at the barrier returns the rebate.

## [1.134.0] - 2026-09-10

### Added
- `crank_nicolson_greeks` (in `pde.py`): delta, gamma and theta read straight
  off the Crank-Nicolson grid -- no extra solves. Delta and gamma are node-level
  central differences interpolated to the spot (so the grid offset does not bias
  them), and theta is the difference between the `t=0` grid and the grid one
  time step earlier. The solver was refactored into a shared `_cn_solve` used by
  both the price and Greeks entry points.
- Verified against Black-Scholes: delta to ~1e-4, gamma to ~5e-6, theta to
  ~3e-3; the put delta matches; American-put Greeks have the right signs.

## [1.133.0] - 2026-09-10

### Added
- Crank-Nicolson PDE solver (new `pde.py`): `crank_nicolson_price` solves the
  Black-Scholes-Merton PDE on a spot grid with the (unconditionally stable,
  second-order) Crank-Nicolson scheme. Supports a constant `sigma` or a
  `local_vol_fn(S, t)`, a cost-of-carry `b = r - q`, and American exercise via
  projected SOR (PSOR) enforcing `V >= payoff`. The interior tridiagonal system
  is solved by the Thomas algorithm.
- Verified: European calls and puts match Black-Scholes; the error converges
  `O(1/n^2)` in the grid; the dividend-carry and flat-local-vol cases match BS;
  the American put matches a 3000-step binomial tree to ~1e-2 and exceeds the
  European value.

## [1.132.0] - 2026-09-10

### Added
- Arbitrage-free SSVI calibration (in `ssvi.py`): `calibrate_ssvi` gains an
  `arb_weight` that penalises butterfly and calendar no-arbitrage violations,
  and `calibrate_ssvi_arbitrage_free` ramps that weight (warm-starting each fit
  from the last) until the surface passes `ssvi_is_arbitrage_free`. Trades a
  little fit RMSE for a guaranteed no-arbitrage surface.
- Verified: an intentionally arbitraging SSVI surface is flagged and the plain
  fit reproduces it, while the penalised fit returns an arbitrage-free surface
  (pulling `eta` 8.0 -> 1.8); on a clean surface the arb-free fit matches the
  plain fit's RMSE and removing arbitrage never improves the fit to an
  arbitraging market.

## [1.131.0] - 2026-09-10

### Added
- LSV leverage-function calibration (new `lsv.py`): `calibrate_lsv_leverage`
  calibrates a local-stochastic-volatility leverage surface `L(S, t)` by the
  Guyon-Henry-Labordere particle method -- Gyongy's condition
  `L(K,t)^2 = sigma_Dupire(K,t)^2 / E[V_t | S_t = K]`, with the conditional
  expectation estimated by simulating the Heston variance and binning `V_t` by
  the spot level, marched forward with Euler sub-steps between the calibration
  expiries. Returns the leverage grid and an interpolating `lev_fn(spot, t)`.
- Verified: the front-expiry ATM leverage equals `sigma_local / sqrt(v0)`
  exactly (Gyongy at `t -> 0`, `E[V|S] = v0`); leverage is positive everywhere;
  a higher target vol scales it up; and the flat-target front-expiry leverage is
  uniform at `sigma_local / sqrt(v0)`.

## [1.130.0] - 2026-09-10

### Added
- `calibrate_heston` (new `heston_calib.py`): fit the five Heston parameters
  `(v0, kappa, theta, xi, rho)` to an implied-vol surface by least squares on
  Black vol over `(expiry, strike, vol)` quotes, pricing each candidate with the
  Fourier `heston_price` and optimising with Nelder-Mead under a smooth
  reparametrization. An optional `feller_weight` penalises Feller-condition
  violations (`2 kappa theta < xi^2`) to prefer a strictly-positive variance
  process.
- Verified: recovers a synthetic Heston surface's parameters (rmse ~0); the
  Feller penalty pushes `2 kappa theta` up to the `xi^2` boundary; all fitted
  parameters stay in their valid regions.

## [1.129.0] - 2026-09-10

### Added
- `calibrate_double_heston` (new `double_heston_calib.py`): fit all ten
  double-Heston parameters (two variance factors) to an implied-vol surface by
  least squares on Black vol over `(expiry, strike, vol)` quotes, pricing each
  candidate with the Fourier `double_heston_price` and optimising with
  Nelder-Mead under a smooth reparametrization (variances/vol-of-vols positive,
  correlations in `(-1, 1)`). A two-scale (fast + slow factor) seed is used by
  default.
- Verified: calibrating to a synthetic double-Heston surface (4 expiries x 5
  strikes) fits to ~2e-5 RMSE with all fitted parameters in their valid regions.

## [1.128.1] - 2026-09-10

### Fixed
- `rough_heston_price` now raises a clear error when the fractional Riccati
  diverges (a too-coarse `n_grid` for small `H`) instead of silently returning
  NaN -- found while trimming the test suite, where an H=0.2 parity check on a
  60-step grid produced NaN.

### Tests
- Trimmed the fast-suite runtime from ~66s to ~28s: marked the heaviest
  Fourier/calibration cross-checks `slow` (rough-Heston Hurst-reduction and
  short-skew, Levy NIG/Meixner/cross-model recovery, the Levy-surface calendar
  scan, the 256-step rBergomi FFT identity) and cut redundant path counts on the
  FFT-identity checks (numerical-identity, not statistical). Full coverage is
  unchanged; the heavy checks still run under `-m slow`.

## [1.128.0] - 2026-09-10

### Added
- SABR butterfly-arbitrage detection and repair (in `sabr.py`): `sabr_density`
  is the Breeden-Litzenberger implied density of a SABR smile (second strike
  derivative of the Black call at the SABR vol); `sabr_butterfly_arbitrage`
  returns the strikes where that density goes negative (the Hagan expansion is
  not guaranteed arb-free in the wings); `sabr_is_arbitrage_free` is the
  boolean; and `sabr_repair_butterfly` shrinks the vol-of-vol `nu` until the
  density is non-negative, preserving alpha/beta/rho.
- Verified: a benign smile's density integrates to 1 and is arbitrage-free; an
  extreme vol-of-vol flags negative-density strikes; repair restores arbitrage
  freedom (`nu` 2.5 -> 1.66) while leaving a clean smile untouched.



### Added
- CMS convexity adjustment (new `cms.py`): `cms_adjustment_standard` gives the
  linear-TSR / Hagan standard-model adjustment `G * Var_A(S_T)` (using the exact
  lognormal variance), `cms_rate_convexity_replication` computes it by model-free
  static replication over a swaption strip with a caller-supplied volatility
  smile (Carr-Madan second-moment replication times the level factor `G`), and
  `cms_rate` returns the convexity-adjusted expected CMS rate.
- Verified: zero vol gives no adjustment; the adjustment is positive; the
  flat-smile replication matches the standard model to ~1e-5; both match a
  linear-TSR Monte Carlo to ~2e-4; a convex smile raises the adjustment; and a
  payment lag increases it.



### Added
- Swaption volatility cube (new `volcube.py`): `VolCube` stores a calibrated
  SABR smile at each `(expiry, tenor)` node and interpolates across all three
  market axes -- SABR handles strike/smile analytically, while expiry and tenor
  use bilinear interpolation of the total variance `sigma^2 * expiry` (the
  no-calendar-arbitrage-friendly choice). `VolCube.fit` builds it from market
  smiles; `vol(expiry, tenor, strike)` queries any point.
- Verified: each node reprices its own SABR smile, an interpolated `(3, 7)` node
  lands between its neighbours, expiry interpolation stays in range, the SABR
  smile shape is preserved, and beyond the grid the nearest node is used.

## [1.125.0] - 2026-09-10

### Added
- OIS/LIBOR dual-curve discounting (new `dualcurve.py`): value swaps under the
  post-2008 convention where forward rates are projected off a forward (LIBOR)
  curve while every cashflow is discounted on the collateral (OIS) curve.
  `dual_forward_rate`, `dual_float_leg_value`, `dual_par_swap_rate` and
  `dual_swap_value` take separate OIS and projection `DiscountCurve`s plus an
  additive basis; `dual_calibrate_basis` solves the constant basis spread that
  reprices a set of par swaps.
- Verified: with a single shared curve the dual-curve par rate collapses to the
  single-curve one exactly; the flat-curve forward equals `e^r - 1`; a known
  basis is recovered to ~1e-6; a payer swap is worth zero at par and
  payer + receiver = 0; adding basis raises the float leg; and a projection
  curve above OIS raises the par rate.

## [1.124.0] - 2026-09-10

### Added
- Discount curve (new `discount_curve.py`): `DiscountCurve` builds a log-linear
  (piecewise-constant instantaneous-forward) discount curve from pillar
  `(T, DF)` points or continuously-compounded zero rates, with zero-rate,
  forward-rate and par-swap-rate accessors. `bootstrap_from_swaps` bootstraps
  pillar discount factors from par swap rates (each pillar solves a linear
  equation given the shorter ones). This is exactly the initial-curve object the
  Gaussian short-rate models (`cheyette`, `g2pp`) consume.
- Verified: a flat zero-rate curve reproduces `e^{-zT}`; interpolation is exact
  at pillars; a consecutive-annual-tenor bootstrap reproduces its par rates to
  1e-10; and G2++ reprices the bootstrapped curve's own bonds exactly.

## [1.123.0] - 2026-09-10

### Added
- Rough-Heston model (new `rough_heston.py`): `rough_heston_price` and
  `rough_heston_smile`. The variance is driven by a fractional kernel with Hurst
  `H in (0, 0.5]`; its characteristic function comes from the fractional Riccati
  equation, solved here by the fractional Adams predictor-corrector
  (Diethelm-Ford-Freed) and inverted by the same Gil-Pelaez two-probability
  Gauss-Legendre integral as Heston.
- Verified: `H = 0.5` reduces to classical Heston (with `xi = kappa * nu`, the
  El Euch-Rosenbaum convention) to <5e-3 and converges in the grid size; the
  characteristic function is a martingale; put-call parity holds; and `H = 0.1`
  gives a much steeper short-dated skew (-0.88) than `H = 0.5` (-0.42) -- the
  rough-vol signature.

## [1.122.0] - 2026-09-10

### Added
- Double-Heston two-factor stochastic-volatility model (new `double_heston.py`):
  `double_heston_price` and `double_heston_smile`. Two independent Heston
  variance factors (a fast- and a slow-reverting one) let the short- and
  long-dated skew move more independently than single-factor Heston allows.
  Since the factors are independent, the log-spot characteristic function is the
  product of the two single-factor pieces; priced by the same Gil-Pelaez
  two-probability Gauss-Legendre integral as Heston. The `xi -> 0` factor limit
  is handled analytically (deterministic-variance contribution).
- Verified: zeroing the second factor recovers single-factor Heston exactly,
  put-call parity holds, two factors are worth more than one, negative
  correlations give a downward skew, and the Fourier price matches an
  independent two-factor QE Monte Carlo to ~0.6 SE.

## [1.121.0] - 2026-09-10

### Added
- `bermudan_swaption_g2pp` (new `bermudan_swaption.py`): Bermudan swaption
  pricing by Longstaff-Schwartz on the G2++ state. Simulates the correlated
  ``(x, y)`` factors on the exercise schedule (exact OU transitions), carries a
  discretely-compounded money-market numeraire from the one-period G2++ bonds,
  values the co-terminal swap in closed form from the state, and regresses the
  discounted continuation on a quadratic basis in ``(x, y)`` for the
  exercise decision.
- Verified: price is positive; a multi-date Bermudan is worth at least the
  single-exercise value; payer and receiver are both positive; a payer's value
  falls as the fixed rate rises; payer/receiver swap values are antisymmetric;
  the simulated factors have ~zero mean.

## [1.120.0] - 2026-09-10

### Added
- Two-factor G2++ Gaussian short-rate model (new `g2pp.py`): two correlated
  mean-reverting factors, `r = x + y + phi(t)`. `g2pp_zero_bond` is the
  exponential-affine bond off the initial curve (with the closed-form variance
  `g2pp_V`), `g2pp_bond_option` the exact Gaussian bond-option price, and
  `g2pp_caplet` via the bond-put identity.
- Verified: the bond option matches a forward-measure Monte Carlo to ~1e-5,
  put-call parity holds, turning the second factor off (`eta -> 0`, `rho = 0`)
  reproduces the single-factor Cheyette price exactly, factor correlation raises
  the bond vol and option price monotonically, and `sigma = eta = 0` gives the
  discounted intrinsic.

## [1.119.0] - 2026-09-10

### Added
- Single-factor Cheyette / quasi-Gaussian short-rate model (new `cheyette.py`):
  the Markovian HJM representation with state `(x, y)` and short rate
  `r = f(0,t) + x`. `cheyette_zero_bond` reconstitutes `P(t,T)` off the initial
  curve via the `G(t,T)` function, and `cheyette_bond_option` / `cheyette_caplet`
  give the exact prices for a constant short-rate vol (where the factor
  coincides with Hull-White), the caplet using the bond-put identity.
- Verified: the bond option matches a forward-measure Monte Carlo to ~1e-5,
  put-call parity holds, the caplet equals its `(1 + K tau)` bond-put
  composition exactly, `kappa -> 0` stays finite, and `sigma -> 0` collapses to
  the discounted intrinsic.

## [1.118.0] - 2026-09-10

### Added
- Andreasen-Huge single-step arbitrage-free local-vol smile (new
  `andreasenhuge.py`): `andreasen_huge_prices` solves one implicit Dupire
  finite-difference step (a tridiagonal M-matrix system, Thomas algorithm) to
  produce call prices that are monotone-decreasing and convex in strike -- hence
  arbitrage-free -- for any positive local-vol grid. `andreasen_huge_smile`
  inverts them to implied vols, and `andreasen_huge_calibrate` bootstraps the
  per-strike local vols that reprice a market smile.
- Verified: the prices are monotone and convex for both flat and steeply skewed
  local vols; calibration reproduces a market skew to <5e-3 RMSE (interior
  strikes to ~1e-4); the calibrated smile is downward-sloping.

## [1.117.0] - 2026-09-10

### Added
- `kim_put_greeks` (in `kim.py`): delta, gamma and theta of a Kim American put.
  The early-exercise boundary is spot-independent, so it is solved once (the
  expensive step) and the spot bumps only re-run the cheap European-plus-premium
  evaluation on that fixed boundary; theta uses a maturity bump.
- Verified against binomial bump Greeks: delta matches to ~5e-3, gamma to ~2e-3
  (stable in the grid size -- the small residual is Kim's boundary
  approximation, not quadrature noise), put delta is negative, gamma positive,
  and theta negative.

## [1.116.0] - 2026-09-10

### Added
- `baw_american` (new `baw.py`): the Barone-Adesi-Whaley (1987) quadratic
  approximation for American options. Splits the price into the European value
  plus an early-exercise premium (a power of spot anchored at the critical price
  ``S*``, found by a 1-D Newton solve of the value-matching condition). Fast
  closed-form, dividends via ``b = r - q``; a no-dividend call returns the
  European value.
- Verified: within ~0.1 of a 4000-step binomial tree for puts and dividend
  calls (BAW's known approximation accuracy), within ~0.05 of the accurate Kim
  integral-equation price, and correctly floored at intrinsic.

## [1.115.0] - 2026-09-10

### Added
- American pricing via Kim's (1990) integral equation (new `kim.py`):
  `kim_american_put`, `kim_american_call`, `kim_exercise_boundary`. The
  early-exercise boundary solves a Volterra integral equation by backward
  marching from expiry, bisecting the value-matching condition at each step; the
  price is the European value plus the early-exercise premium integrated at
  spot. The call uses the McDonald-Schroder put-call symmetry
  `C(S,K,r,q) = P(K,S,q,r)`.
- Verified: puts and dividend calls match a 4000-step binomial tree to ~1e-3 at
  200 grid steps; the boundary rises to K at expiry; a no-dividend call equals
  the European value.

### Fixed
- Two bugs caught against the binomial reference: the premium integral indexed
  the boundary backward (`B[n-m]` instead of `B[m]`, the boundary `s` years
  ahead), and the `s -> 0` integrand endpoint was forced to zero when on the
  boundary it is `rK/2 - qS/2` (from `N(0)=1/2`). With both fixed the price
  converges to the tree instead of ~1.8% low.

## [1.114.0] - 2026-09-10

### Added
- Sobol sequence + Brownian-bridge QMC (new `sobol.py`): a `Sobol` generator
  built from primitive-polynomial direction numbers (Joe-Kuo initial values) via
  the Gray-code recurrence, `brownian_bridge_path` which loads a path's dominant
  variance onto the leading (most uniform) Sobol coordinates, and two pricers --
  `sobol_european` (1-D) and `sobol_asian` (Brownian-bridge, multi-step).
- Verified: Sobol points cover the unit cube evenly (each quarter of the 1-D
  coordinate within 2%), `sobol_european` matches Black-Scholes and its error is
  materially smaller than pseudo-random MC at matched N (~0.007 vs ~0.09 at
  N=8192), the Brownian bridge has the right terminal moments (mean ~0, variance
  ~t), and `sobol_asian` matches a discrete-monitoring pseudo-random MC to ~0.01.

## [1.113.1] - 2026-09-10

### Changed
- The rBergomi conditional (turbocharged) estimator's `I1`/`QV` path loop now
  shares the same cached-kernel FFT Volterra convolution as `rbergomi_paths`
  (factored into module-level `_build_far_kernel` / `_far_sums_fft`), with the
  same `fast="auto"` switch at `n_steps >= 200`. Verified the FFT and direct
  `I1`/`QV` agree to ~1e-14 and the CV price/tests are unchanged.

## [1.113.0] - 2026-09-10

### Added
- `LevySurface` (new `levysurface.py`): a unified implied-vol surface generated
  by a single exponential-Levy parameter set (`vg`, `nig`, `meixner`, `cgmy`).
  `implied_vol(k, t)` / `total_variance(k, t)` evaluate any point, `grid(...)`
  builds an expiry x log-moneyness surface with one Carr-Madan FFT strip per
  maturity, and `calendar_violations` / `is_calendar_arbitrage_free` check that
  total variance is non-decreasing in maturity (the same test as the SVI
  `VolSurface`).
- `levy_psi(model, params)` factory (in `levycalib.py`): returns the
  characteristic exponent for a named model + raw parameters, for pricing or
  surface-building outside calibration.
- Verified: the surface reprices the NIG and VG per-model smiles to ~1e-3, a
  genuine Levy law is calendar-arbitrage-free across five expiries, total
  variance grows with maturity, and `beta < 0` gives a downward skew at each
  expiry.

## [1.112.0] - 2026-09-10

### Changed
- rBergomi path generation (`rbergomi_paths`) now evaluates the hybrid-scheme
  Volterra convolution with a cached-kernel FFT instead of the `O(n^2)` inner
  double loop. The far-cell contribution `F[i] = sum_j dW[j] g[i-j]` is a
  discrete convolution of the fixed weight kernel with the Brownian increments;
  the kernel's FFT is precomputed once and reused across every path (one forward
  + one inverse transform of the increments per path), giving `O(n log n)`. A
  new `fast` argument (`"auto"` default) uses the FFT when `n_steps >= 200` and
  the direct loop below that (where pure-Python FFT overhead dominates); the two
  agree to ~1e-12. Measured ~2.3x faster at 512 steps, growing with `n_steps`.

## [1.111.0] - 2026-09-10

### Added
- `calibrate_levy_smile` (new `levycalib.py`): fit an exponential-Levy model
  (`vg`, `nig`, `meixner`, `cgmy`) to a one-expiry market implied-vol smile by
  least squares on vol. Each candidate smile is priced with a single Carr-Madan
  FFT strip (interpolated to the market strikes) rather than an integral per
  strike, and a per-model smooth parameter transform keeps the optimiser
  unconstrained while the raw parameters stay in their valid region.
- Verified: recovers synthetic VG, NIG and Meixner parameters to <1e-3 RMSE;
  CGMY (weakly identified from few strikes) is checked on fit quality (RMSE
  < 5e-3); a cross-model NIG fit to a VG smile is close (RMSE ~1e-3).

## [1.110.0] - 2026-09-10

### Added
- `variance_gamma_smile` (in `variancegamma.py`): the Black-Scholes implied-vol
  smile the VG model produces, `theta < 0` giving a downward skew.

### Changed
- `variance_gamma_price` refactored onto the shared `carrmadan.levy_price`
  engine via the closed-form VG exponent
  `psi(u) = -(1/nu) log(1 - i theta nu u + 0.5 sigma^2 nu u^2)`, replacing its
  own two-probability Gauss-Legendre loop. Same prices, and it now inherits the
  COS-method cross-check like the other Levy models. Verified: Carr-Madan and
  COS agree to ~4.5e-6 on three parameter sets, the existing VG tests still
  pass, parity holds, and small `nu` recovers Black-Scholes.

## [1.109.0] - 2026-09-10

### Added
- `cos_greeks` (in `carrmadan.py`): analytic delta and gamma of the COS-method
  price for any Levy model. Spot enters the cosine sum only through the
  characteristic function's `e^{i u (x + mu)}` factor with `x = ln(S/K)`, so
  differentiating term-by-term gives delta (`i u / S` per term) and gamma
  (`i u (i u - 1)/S^2`) with no re-pricing and no finite differences; the put
  follows by parity.
- Verified: matches Black-Scholes delta/gamma to 1e-6 for the GBM exponent, and
  matches finite-difference COS Greeks on CGMY and NIG to ~1e-4.

## [1.108.0] - 2026-09-10

### Added
- Meixner Levy model (new `meixner.py`): `meixner_price` and `meixner_smile`.
  Schoutens' Meixner process has the analytic exponent
  `psi(u) = 2 d (ln cos(b/2) - ln cosh((a u - i b)/2))` (scale `a`, asymmetry
  `b in (-pi, pi)`, activity `d`) and is priced through the shared Carr-Madan
  engine, so it also gets the COS cross-check for free.
- Verified: Carr-Madan and COS prices agree to ~1.7e-6 on three parameter sets,
  put-call parity holds, `b < 0` gives a downward skew and `b > 0` an upward
  skew, and more activity `d` raises the price.

## [1.107.0] - 2026-09-10

### Added
- COS method (in `carrmadan.py`): `cos_price` implements the Fang-Oosterlee
  (2008) Fourier-cosine expansion of the risk-neutral log-return density on a
  cumulant-based truncation range, with closed-form call payoff coefficients
  (chi/psi). Exponentially convergent and independent of the Carr-Madan
  transform, so it is a genuine cross-check for the whole Levy family.
- Verified: matches Black-Scholes to 1e-8 (GBM exponent), and matches the
  Carr-Madan pricer to ~1e-6/1e-4 on CGMY and three NIG parameter sets; parity
  holds and the error shrinks with the term count.

### Fixed
- The COS characteristic function omitted the `x = ln(S/K)` shift, so only ATM
  strikes priced correctly; and the fourth-difference c4 cumulant estimate blew
  up on non-smooth exponents (CGMY's Gamma(-Y) power law), ballooning the
  truncation range until the density undersampled. c4 now uses a larger step and
  is clamped to a sane multiple of c2^2 (it only fine-tunes the range).

## [1.106.0] - 2026-09-10

### Added
- Carr-Madan FFT strip (in `carrmadan.py`): `carr_madan_strip` prices a whole
  log-strike grid of European calls in a single transform, and
  `carr_madan_smile_strip` returns the implied-vol smile over a log-moneyness
  window from that one pass. Backed by a pure-Python radix-2 Cooley-Tukey `_fft`
  (bit-reversal + butterflies, no NumPy) and Simpson-weighted frequency
  sampling; works for any Levy model via its characteristic exponent (CGMY, NIG,
  ...). The log-strike spacing is `2 pi / (n_fft * eta)`.
- Verified: `_fft` matches a naive DFT and round-trips to ~1e-12; on the native
  FFT grid the CGMY and NIG strip prices match the per-strike Gauss-Legendre
  pricer to ~2e-6, and the smile strip is sorted, in-window, and shows the
  expected skew.

## [1.105.0] - 2026-09-10

### Added
- Shared Carr-Madan engine (new `carrmadan.py`): `levy_price` / `carr_madan_call`
  price any exponential-Levy model from its characteristic exponent `psi(u)`,
  applying the martingale correction `omega = -psi(-i)` and the alpha-damped
  Fourier inversion over the shared Gauss-Legendre nodes.
- Normal Inverse Gaussian model (new `nig.py`): `nig_price` and `nig_smile`.
  Barndorff-Nielsen's NIG has the analytic exponent
  `psi(u) = delta (sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta + i u)^2))`
  (tail `alpha`, asymmetry `beta`, scale `delta`) and is priced through the
  shared engine.
- Verified NIG against an independent Gil-Pelaez inversion to ~1e-6 on three
  parameter sets, with put-call parity, `beta < 0` -> downward skew and
  `beta > 0` -> upward skew.

### Changed
- `cgmy_price` now routes through the shared `carrmadan.levy_price` instead of
  its own inlined Fourier loop (identical prices; the `_cgmy_char_logspot`
  helper is kept for the Gil-Pelaez cross-check).

## [1.104.0] - 2026-09-10

### Added
- CGMY tempered-stable Levy model (new `cgmy.py`): `cgmy_price` and `cgmy_smile`.
  A pure-jump process with a stable Levy density tempered independently on each
  tail (`C` activity, `G`/`M` down/up tempering, fine-structure `Y < 2`), whose
  closed-form characteristic exponent uses `Gamma(-Y)` and complex powers.
  Priced by Carr-Madan Fourier inversion of the damped call transform (damping
  `alpha`, requiring `alpha + 1 < M`) over the shared Gauss-Legendre nodes, with
  a martingale drift correction `omega = -psi(-i)`.
- Verified against an independent Gil-Pelaez inversion of the same
  characteristic function to ~1e-6 across three regimes (`Y = 0.5, 0.8, 1.2`);
  put-call parity holds, `G < M` produces a downward skew and `G = M` a
  symmetric smile, and more activity (`C`) raises the price.

## [1.103.0] - 2026-09-10

### Added
- `leisen_reimer_american_accel` (in `leisen_reimer.py`): the Broadie-Detemple
  (1996) two-point Richardson extrapolation of the American Leisen-Reimer price,
  `V_ext = 2*V(2n) - V(n)`, which cancels the leading `1/n` error that the
  American tree carries (its exercise boundary breaks the smooth-payoff
  assumption behind Peizer-Pratt). Measured against a 6000-step CRR tree it is
  1.6-20x more accurate than a single 101-step LR tree for the same work,
  reaching a few mils; `b = r - q` supports dividends.

## [1.102.0] - 2026-09-10

### Added
- Leisen-Reimer binomial tree (new `leisen_reimer.py`): `leisen_reimer_price`
  and `leisen_reimer_greeks`. Centres the tree on the strike and sets the
  up-move and probability from a Peizer-Pratt inversion of the Black-Scholes
  d1/d2, giving smooth `O(1/n^2)` convergence instead of CRR's slow oscillation.
  Handles European and American exercise and a continuous dividend yield via
  `b = r - q`; the step count is forced odd so the strike sits at the tree
  centre.
- Verified: European prices match Black-Scholes to <2e-3 at ~51 steps and beat
  CRR at the same step count; American prices with dividends reach a 3000-step
  CRR tree to ~1 cent (the American convergence is slower than the European
  `O(1/n^2)`, since the smooth-payoff assumption breaks at the exercise
  boundary, so a larger step count is used).

## [1.101.0] - 2026-09-10

### Added
- Kou (2002) double-exponential jump-diffusion (new `kou.py`): `kou_price` and
  `kou_smile`. The jump size is an asymmetric double exponential (up-jump tail
  rate `eta1 > 1`, down-jump tail rate `eta2 > 0`, up-probability `p`), giving
  fatter, asymmetric tails than Merton's Gaussian jumps with a still-analytic
  Levy characteristic function. Priced with the same two-probability
  Gauss-Legendre Fourier integral as the Heston/Bates pricers, with the jump
  drift compensator keeping the discounted spot a martingale.
- Verified: `lambda = 0` recovers Black-Scholes to 1e-14, put-call parity holds,
  jumps raise the price, asymmetric jumps produce a downward skew, and the
  Fourier price matches an independent double-exponential-jump Monte Carlo to
  <0.9 standard errors on three parameter sets.

## [1.100.0] - 2026-09-10

### Added
- `svi_local_variance` and `svi_surface_local_vol` (in `svi.py`): analytic
  Dupire local volatility for raw SVI, mirroring the SSVI local-vol pair. The
  strike derivatives w_k, w_kk come in closed form from the SVI parametrization
  (reusing `_svi_derivs`); `svi_local_variance` takes a caller-supplied dw/dt
  (one slice carries no maturity), while `svi_surface_local_vol` interpolates
  total variance linearly across a term structure of slices to supply dw/dt.
  Both raise on a non-positive Dupire denominator (butterfly-arbitrage flag).
- Verified: the surface local vol matches a finite-difference Dupire on the same
  slice term structure to ~1e-7 across strikes and maturities.

## [1.99.0] - 2026-09-10

### Added
- Bates (1996) model (new `bates.py`): Heston stochastic volatility plus Merton
  lognormal jumps. `bates_price` reuses the Heston two-probability Fourier
  integral and Gauss-Legendre machinery, multiplying the Heston characteristic
  function by the independent compound-Poisson jump factor with a martingale
  drift compensator. `bates_smile` returns the implied-vol smile. `lambda = 0`
  recovers the Heston price exactly.
- Verified: `lambda = 0` matches Heston to 1e-8, put-call parity holds, jumps
  raise the price vs Heston, a negative mean jump steepens the downward skew,
  and -- the decisive check -- the Fourier price matches an independent
  Andersen-QE variance simulation with a compound-Poisson jump overlay to <0.35
  standard errors on three parameter sets.

### Fixed
- The jump characteristic factor initially used a `+-1/2` measure shift and an
  `i*phi` compensator; the Monte Carlo cross-check exposed a 12-86 sigma error.
  The correct Heston-decomposition shift is `s = 1` (share measure, argument
  `phi - i`) for P1 and `s = 0` for P2, with the compensator `-a*lambda*t*k`;
  after the fix the MC agreement is <0.35 sigma. Parity and the `lambda = 0`
  limit alone did not catch this -- only the independent jump-aware MC did.

## [1.98.1] - 2026-09-10

### Tests
- Cross-check the Turnbull-Wakeman arithmetic-Asian Greeks against an
  independent method: bump Greeks of the arithmetic-Asian Monte Carlo under
  common random numbers. Delta, gamma and vega agree to ~1-2% (TW is a
  two-moment approximation), and the arithmetic and geometric Greeks converge at
  low vol. The prior Asian-Greeks tests only finite-differenced the same closed
  form; these confirm the analytic values against a separate pricer.

## [1.98.0] - 2026-09-10

### Added
- `ssvi_local_vol_fn` and `ssvi_reprice_mc` (in `ssvi.py`): close the calibrate
  -> local-vol -> reprice loop. `ssvi_local_vol_fn` turns a fitted SSVI surface
  into a `(spot, tau) -> local vol` callable (mapping spot to SSVI log-moneyness
  on the forward), and `ssvi_reprice_mc` simulates that surface through
  `local_vol_mc`. A correct local-vol construction reprices the SSVI *implied*
  smile, and it does: the Monte Carlo price matches the closed-form price at the
  SSVI implied vol across the smile to within ~1.2 standard errors, provided the
  fitted surface includes short maturities.

### Fixed
- `ssvi_local_vol_from_params` / `ssvi_local_vol_fn` froze the ATM variance below
  the shortest fitted expiry, which biased the local-vol Monte Carlo more as the
  step count grew (the path integrates `tau` from 0 but saw a frozen short-end
  vol). The short end now linearly extrapolates `theta -> 0` toward the origin,
  matching SSVI's small-time behaviour; the bias-vs-steps drift is gone.

## [1.97.0] - 2026-09-10

### Added
- `ssvi_local_variance` and `ssvi_local_vol_from_params` (in `ssvi.py`): the
  Dupire local variance of an SSVI surface in fully analytic form. The total
  variance's strike derivatives (w_k, w_kk) and its theta-derivative come from
  the SSVI parametrization in closed form -- no finite differences -- and feed
  Gatheral's total-variance Dupire formula. `ssvi_local_vol_from_params` builds
  theta(t) and theta'(t) by piecewise-linear interpolation of the fitted ATM
  variances and returns the local vol at any (k, t) in range; it also raises on
  a non-positive Dupire denominator (a butterfly-arbitrage flag).
- Verified: the analytic local variance matches a finite-difference Dupire on
  the same surface to ~1e-9 across strikes and maturities.

## [1.96.0] - 2026-09-10

### Added
- Surface SVI (new `ssvi.py`): the Gatheral-Jacquier (2014) arbitrage-free
  whole-surface parametrization, tying every expiry's smile together through a
  shared skew function and the ATM total-variance term structure `theta_t`.
  `calibrate_ssvi` fits the global `(rho, eta, gamma)` power-law skew plus one
  `theta` per expiry to a `(t, k, iv)` market surface; `ssvi_butterfly_free`,
  `ssvi_calendar_free` and `ssvi_is_arbitrage_free` implement the sufficient
  no-static-arbitrage conditions (density positivity within a slice, total
  variance non-decreasing in maturity across slices).
- Verified: recovers a synthetic surface's parameters to ~1e-8 RMSE, the
  butterfly condition flags an excessive-skew slice, the calendar condition
  flags a decreasing-`theta` term structure, and `phi` decays with maturity.

## [1.95.0] - 2026-09-10

### Added
- `rbergomi_price_cv` and `rbergomi_smile_cv` (in `rbergomi.py`): the conditional
  ("turbocharged") rough Bergomi estimator of McCrickerd & Pakkanen (2018).
  Conditioning on the volatility-driving Brownian motion makes the terminal
  log-spot Gaussian, so each path contributes a smooth Black-Scholes conditional
  price (effective spot `S exp(rho I1 - rho^2 QV/2)`, effective variance
  `(1 - rho^2) QV`) instead of a noisy indicator payoff. Same price as the plain
  estimator within Monte Carlo error, with the standard error cut severalfold --
  ~2x at rho=-0.7 up to ~6-7x as `|rho|` shrinks (the integrated-out noise
  fraction is `1 - rho^2`).

## [1.94.0] - 2026-09-10

### Added
- `rbergomi_price` and `rbergomi_smile` (new `rbergomi.py`): the rough Bergomi
  stochastic-volatility model (Bayer-Friz-Gatheral 2016), whose variance is
  driven by a rough fractional process with Hurst exponent `H < 1/2`. Simulated
  with the Bennedsen-Lunde-Pakkanen (2017) hybrid scheme (kappa=1): the singular
  cell nearest each step is sampled exactly from the joint law of the Brownian
  increment and its kernel integral, the far cells use the optimally-placed
  discretised kernel. `rbergomi_smile` reprices every strike on common terminal
  spots and inverts to a Black-Scholes vol.
- Verified: the discounted spot is a martingale (`E[S_T] e^{-rt} = S` to <5e-3),
  `eta = 0` collapses to a flat smile at `sqrt(xi0)`, `rho < 0` gives a downward
  skew, and -- the point of rough vol -- `H = 0.1` yields a steeper short-dated
  ATM skew (-0.58) than the `H = 0.5` diffusive case (-0.41).

### Fixed
- Hybrid-scheme weights divided by `alpha = H - 1/2`, blowing up at `H = 1/2`;
  that case now returns the constant-kernel (standard Brownian) weights.

## [1.93.0] - 2026-09-10

### Added
- `VannaVolgaSmile.price` and `.vol_price_corrected` (in `vannavolga.py`): the
  exact second-order Castagna-Mercurio vanna-volga construction. Starts from the
  flat-ATM Black-Scholes price and adds the three pillar options' market-minus-
  ATM price gaps, weighted (via a 3x3 vega/vanna/volga solve) so the hedging
  portfolio matches the target's vega, vanna and volga. Unlike the existing
  quadratic vol interpolation (`.vol`), this reprices the three market
  instruments *exactly in price*, the standard FX smile pricer.
- Verified: reprices all three market pillars to <1e-10, price-corrected vol
  equals the pillar vols, put-call parity holds to 1e-16, a negative risk
  reversal puts the put wing above the call wing, and a flat market gives a flat
  smile.

## [1.92.0] - 2026-09-10

### Added
- `calibrate_sabr_lm` (in `sabr.py`): a Levenberg-Marquardt SABR calibrator that
  reuses the exact `sabr_jacobian` from v1.91.0 instead of the derivative-free
  Nelder-Mead. Solves the damped 3x3 normal equations (new `_solve3` Gaussian
  elimination), adapts the damping to guarantee a downhill step, and box-clamps
  (alpha, rho, nu) to the valid region. Returns `(params, rmse, n_iter)`.
- Recovers synthetic (alpha, rho, nu) to ~1e-12 RMSE in ~4 iterations (vs
  Nelder-Mead's ~2e-9) and lands on the same optimum as `calibrate_sabr` on a
  noisy market; also verified for the beta=1 lognormal smile.

## [1.91.0] - 2026-09-10

### Added
- `sabr_sensitivities` and `sabr_jacobian` (in `sabr.py`): exact partial
  derivatives of the Hagan SABR implied vol via a small forward-mode dual-number
  type, with no finite-difference truncation error. `sabr_sensitivities` returns
  the vol plus its partials w.r.t. F, K, alpha, rho and nu; `sabr_jacobian`
  stacks the (alpha, rho, nu) columns into the calibration Jacobian a
  Gauss-Newton / Levenberg-Marquardt step (and the parameter covariance) needs.
- Verified every partial against central finite differences across strikes to a
  max error of ~1e-10.

### Fixed
- The dual-number vol reused a scalar `log(F/K)`, which zeroed the F/K partials;
  `logFK` is now carried as a dual so the smile backbone/skew slopes propagate.

## [1.90.0] - 2026-09-10

### Added
- `heston_qe_mc` (new `heston_mc.py`): Heston Monte Carlo via Andersen's (2008)
  Quadratic-Exponential scheme. The CIR variance is advanced by moment-matching
  to a shifted squared-Gaussian (low vol-of-vol) or an exponential-with-atom
  (high vol-of-vol), so variances stay non-negative by construction; the
  log-asset step uses Andersen's K0..K4 constants with a martingale correction.
  Includes an Acklam inverse-normal `_norm_ppf` for the QE branch draw.
- Verified against the Fourier `heston_price` on four parameter sets (including
  a long-dated, Feller-violating xi=1 case and deep OTM) to within ~1.3 SE, the
  discounted spot is a martingale (forward recovered to 8e-5), and puts match.

## [1.89.0] - 2026-09-10

### Added
- `barrier_mc` (in `montecarlo.py`): Monte Carlo for single-barrier vanilla
  options with a Brownian-bridge crossing correction, the natural cross-check
  for the `barrier_option` closed form including a continuous dividend yield
  (`b = r - q`). Each step contributes the exact conditional probability the
  bridge between its endpoints touched the barrier, removing the discrete-
  monitoring bias that otherwise over-prices knock-outs. Agrees with the
  Reiner-Rubinstein closed form on all four barrier kinds (with q) to a few
  basis points, and knock-in + knock-out reproduces the vanilla.

## [1.88.0] - 2026-09-10

### Added
- `displaced_implied_shift` (in `displaced.py`): calibrate the displaced-diffusion
  shift that reproduces an observed Black-Scholes vol skew. For each trial shift
  the local `sigma` is re-solved to match the ATM quote exactly, so the shift is
  driven purely by the off-ATM skew and is not biased by the local-vs-implied
  vol convention; golden-section search minimises the RMS vol error. Recovers a
  known shift exactly on synthetic quotes (40 -> 40, 100 -> 100, 0 -> 0).

## [1.87.0] - 2026-09-10

### Added
- `merton_smile` (in `merton.py`): the Black-Scholes implied-vol smile a Merton
  jump-diffusion produces — prices calls across strikes and inverts each to a
  BSM vol. Flat at `sigma` with no jumps; symmetric jumps lift the wings into a
  smile and a negative mean jump tilts it into a downward skew.

## [1.86.1] - 2026-09-10

### Added
- Test hardening: a barrier in-out parity property test (`KI + KO = vanilla`
  across all four barrier kinds, calls/puts, and many strike/barrier
  combinations, with and without dividends) and a Bachelier-to-Black-Scholes
  low-vol ATM convergence test.

## [1.86.0] - 2026-09-10

### Added
- `heston_smile` (in `heston.py`): the Black-Scholes implied-vol smile a Heston
  model produces — prices calls across strikes and inverts each to a BSM vol,
  returning `(log_moneyness, vol)` pairs. Flat at `sqrt(v0)` when the vol-of-vol
  is zero; a negative `rho` gives a downward equity skew, positive an upward
  one.

## [1.85.1] - 2026-09-10

### Added
- Benchmark suite now reports implied-vol solver iterations — Newton with the
  Corrado-Miller seed (~4-5 per quote) vs pure bisection (~29), quantifying the
  seed's benefit. Covered by a smoke test.

## [1.85.0] - 2026-09-10

### Added
- `digital_greeks` (in `exotics.py`): delta and gamma of a cash-or-nothing
  digital by finite differences, demonstrating the pin-risk delta spike as
  expiry nears the strike (the motivation for the call-spread over-hedge).
- Regression guard test: the Kirk spread at strike 0 matches the exact Margrabe
  exchange price.

## [1.84.0] - 2026-09-10

### Added
- `implied_spread_correlation` (in `multiasset.py`): back out the correlation
  implied by a spread-option market price, bisecting the Kirk price (monotone
  decreasing in rho) over `(-1, 1)`. Round-trips exactly and raises for quotes
  outside the rho-spanned price range.

## [1.83.0] - 2026-09-10

### Added
- `svi_repair_butterfly` (in `svi.py`): repair a single SVI slice's butterfly
  arbitrage by geometrically shrinking the wing angle `b` (flattening the smile,
  lifting the g-function) until `svi_is_butterfly_free` passes, preserving the
  level/skew/shift/curvature. A clean slice is returned unchanged.

## [1.82.0] - 2026-09-10

### Added
- `strategy_report` (in `strategy.py`): summarize any strategy `Book`'s expiry
  P&L — net premium, max profit / max loss over a terminal-spot grid, whether
  the profit/loss tail is unbounded, and the break-even spots. Works for every
  builder (verticals, straddles, condors, ratio/backspreads, ...).

## [1.81.0] - 2026-09-10

### Added
- `holee.py`: the Ho-Lee (1986) short-rate model — `holee_zero_coupon_bond` and
  `holee_zero_coupon_yield` via the affine closed form (constant-drift case),
  the simplest no-mean-reversion model. Yield is linear in drift with a
  `sigma^2 t^2 / 6` convexity pull-down; MC-verified.

## [1.80.0] - 2026-09-10

### Added
- `basket_greeks` (in `multiasset.py`): Greeks of a two-asset basket option by
  finite differences — both spot deltas, own-gammas, cross-gamma, and
  correlation sensitivity. Both deltas positive, correlation vega positive
  (higher correlation raises the basket vol), and weight rebalancing shifts the
  delta split.

## [1.79.0] - 2026-09-10

### Added
- `bjerksund_stensland_1993` (in `american.py`): the single-flat-boundary
  Bjerksund-Stensland (1993) American approximation — simpler and slightly less
  accurate than the 2002 two-region version, within a few cents of the binomial
  tree. Calls direct, puts via the exact transformation.

## [1.78.0] - 2026-09-10

### Added
- `lee_wing_slopes` and `lee_bounds_ok` (in `svi.py`): the asymptotic total-
  variance wing slopes of an SVI slice, `b(1-rho)` (left) and `b(1+rho)`
  (right), and a check that both satisfy Lee's moment bound (slope <= 2). Agrees
  with `SVIParams.is_arbitrage_free_wings` and matches the empirical
  large-|k| slope.

## [1.77.0] - 2026-09-10

### Added
- `spread_greeks` (in `multiasset.py`): Greeks of a Kirk spread option by finite
  differences — the two spot deltas, own-gammas, cross-gamma, and correlation
  sensitivity. delta1 > 0 / delta2 < 0, cross-gamma near minus the own-gamma,
  correlation vega negative.

## [1.76.0] - 2026-09-10

### Added
- `cir.py`: the Cox-Ingersoll-Ross (1985) short-rate model — `cir_zero_coupon_bond`
  and `cir_zero_coupon_yield` via the affine closed form with the square-root
  diffusion that keeps the rate non-negative. Short yield equals the short rate,
  the long yield approaches the CIR limit `2 kappa theta / (gamma + kappa)`, and
  bond prices match a floored-Euler Monte Carlo.

## [1.75.0] - 2026-09-10

### Added
- `book_bump_greeks` (in `bookgreeks.py`): net book delta/gamma/vega/theta by
  bumping the shared spot/vol/time across every leg and repricing via
  `price_book` — model-free, so it works for any instrument in the book, not
  only ones with analytic Greeks. Matches the analytic net Greeks on a plain
  option book.

## [1.74.0] - 2026-09-10

### Added
- `income.py`: covered-call and cash-secured-put income analytics —
  `covered_call` and `cash_secured_put` return an `IncomeMetrics` with the
  premium (model or supplied), static and annualized yield, the if-assigned
  return, and the breakeven (`S - premium` / `K - premium`).

## [1.73.0] - 2026-09-10

### Added
- `vasicek.py`: the Vasicek (1977) short-rate model — `zero_coupon_bond` and
  `zero_coupon_yield` (affine closed form) and `bond_option` (Jamshidian
  closed-form European option on a zero-coupon bond). The short yield equals the
  short rate, the long yield approaches the mean-reversion level less the
  convexity term, and bond-option parity holds; MC-verified.

## [1.72.0] - 2026-09-10

### Added
- `VolSurface.vol_grid` and `VolSurface.strike_vol_grid`: export a
  `(expiries x strikes)` grid of implied vols from the surface for charting —
  indexed by log-moneyness or by strike (using the carry-implied forward per
  expiry). Returns the grid plus its axes.

## [1.71.0] - 2026-09-10

### Added
- `double_knockout_mc` (in `montecarlo.py`): Monte Carlo a double-knockout
  (corridor) barrier option that pays the vanilla payoff only if the spot stays
  inside `(lower, upper)` for the whole path, else the cash rebate. Wide
  barriers approach the vanilla; a tighter corridor is cheaper.

## [1.70.0] - 2026-09-10

### Added
- `autocallable_mc` (in `montecarlo.py`): Monte Carlo an autocallable note —
  early redemption with accrued coupon when the spot is at or above the autocall
  barrier at an observation date, otherwise notional at maturity with a
  down-and-in downside below the protection barrier. Higher coupons raise the
  value; a protection barrier lowers it.

## [1.69.0] - 2026-09-10

### Added
- `exchange_greeks` (in `multiasset.py`): Greeks of a Margrabe exchange option
  by finite differences — the two spot deltas, own-gammas, the cross-gamma
  (`d2V/dS1 dS2`), and the correlation sensitivity. Cross-gamma equals minus the
  own-gamma (degree-1 homogeneity) and correlation vega is negative.

## [1.68.0] - 2026-09-10

### Added
- `carry_roll_pnl` (in `attribution.py`): the roll-down / carry-roll P&L of an
  option over a horizon at constant vol — rolls the spot to its forward and
  reprices at the shorter maturity, returning a `CarryRoll` with the value now,
  the rolled value, and the roll P&L (theta bleed net of carry drift). Negative
  for a long option, positive for a short.

## [1.67.0] - 2026-09-10

### Added
- `implied_vol_smile` (in `implied.py`): invert a whole option chain to an
  implied-vol smile in one call, returning `(log_moneyness, vol)` pairs sorted
  by strike and dropping any quote outside the no-arbitrage band.

### Changed
- Marked several heavy Monte Carlo tests `slow`, cutting the default
  (`pytest -m "not slow"`) run back to ~12s.

## [1.66.0] - 2026-09-10

### Added
- `VolSurface.fit_arbitrage_free`: fits each expiry with SVI, then walks from
  the short end up and lifts each slice's level just enough to keep total
  variance non-decreasing in maturity — repairing calendar arbitrage. The
  result passes `is_calendar_arbitrage_free`; a consistent surface is left as
  fitted.

## [1.65.0] - 2026-09-09

### Added
- `barrier_rebate` (in `exotics.py`): the standalone rebate cashflow on a
  barrier — a knock-out rebate (pays cash on breach, at-hit or at-expiry) is a
  one-touch, a knock-in rebate (pays at expiry if never breached) is a
  no-touch. Knock-out (expiry) + knock-in sums to the discounted cash.

## [1.64.0] - 2026-09-09

### Added
- `lookback_greeks` (in `lookback.py`): delta, gamma, vega, and theta of a
  floating- or fixed-strike lookback by central finite differences on the
  closed form (`kind="floating"`/`"fixed"`, with the running extreme).

## [1.63.0] - 2026-09-09

### Added
- `ratio_spread` and `backspread` (in `strategy.py`): a ratio spread (long 1,
  short `ratio` at a further strike; net short options) and its mirror the
  backspread (short 1, long `ratio`; net long options). Both return a leg
  `Book`, so net price/Greeks and the payoff diagram come from the engine.

## [1.62.0] - 2026-09-09

### Added
- GARCH(1,1) volatility forecasting (in `volatility.py`): `fit_garch` estimates
  the variance model by Gaussian quasi-MLE (stationary reparametrization,
  Nelder-Mead), and `garch_forecast` gives the annualized vol `horizon` steps
  ahead, mean-reverting to the long-run level. Recovers persistence on synthetic
  GARCH data.

## [1.61.0] - 2026-09-09

### Added
- `compo_option` (in `quanto.py`): composite (compo) FX option on a foreign
  asset converted at the *floating* exchange rate, so the effective vol combines
  the asset and FX vols with their correlation
  (`sqrt(sa^2 + sfx^2 + 2 rho sa sfx)`). Complements the fixed-FX quanto; higher
  correlation raises the price, and it matches a combined-lognormal Monte Carlo.

## [1.60.0] - 2026-09-09

### Added
- `variancegamma.py`: `variance_gamma_price` prices European options under the
  Variance-Gamma (Madan-Carr-Chang) pure-jump model via its characteristic
  function, integrated with the shared Gauss-Legendre quadrature (no SciPy).
  `nu` controls kurtosis and `theta` the skew; `nu -> 0` recovers Black-Scholes.
  Monte-Carlo verified.

## [1.59.0] - 2026-09-09

### Added
- `average_strike_asian_mc` (in `montecarlo.py`): Monte Carlo an average-strike
  Asian option, where the strike is the realized arithmetic average of the path
  (call pays `max(S_T - A, 0)`). Complements the fixed-strike Asians; a single
  monitoring date collapses the payoff to zero.

## [1.58.0] - 2026-09-09

### Added
- `attribution.py`: `attribute_pnl` explains an option position's realized P&L
  over a market move via a second-order Greek expansion (delta / gamma / vega /
  theta / rho P&L), returning a `PnLAttribution` with the true revaluation, the
  explained sum, and the unexplained residual. The residual is tiny for modest
  moves and grows for large ones (higher-order Greeks).

## [1.57.0] - 2026-09-09

### Changed
- `implied_volatility` now seeds the Newton solve with the Corrado-Miller (1996)
  rational approximation instead of the ATM-only Brenner-Subrahmanyam guess.
  It is accurate away from the money too, cutting the Newton iteration count
  several-fold (~5x fewer across a strike/vol grid), and falls back to the old
  seed in the deep wings. The solver's answers are unchanged.

## [1.56.0] - 2026-09-09

### Added
- `local_vol_mc` (in `montecarlo.py`): Monte Carlo a European option under a
  Dupire local-volatility surface `sigma_loc(S, t)`, evolving the spot in
  log-space with a spot/time-dependent vol at each step. A flat local vol
  reproduces the Black-Scholes price; pairs with `dupire_local_vol` /
  `sabr_local_vol` for surface-consistent pricing.

## [1.55.0] - 2026-09-09

### Added
- `parisian_barrier_mc` (in `montecarlo.py`): Monte Carlo for a Parisian barrier
  option, which activates only after the spot stays past the barrier for a
  *consecutive* window (robust to brief spikes) rather than on a single touch.
  Supports all four in/out, up/down kinds; knock-in + knock-out equals the
  vanilla, and a Parisian knock-out is worth more than the instantaneous one.

## [1.54.0] - 2026-09-09

### Added
- `vol_cone` (in `volatility.py`): the realized-volatility cone — for each
  rolling window length it reports the min / 25th / median / 75th / max and the
  current realized vol (annualized), returning `VolConePoint` per window. The
  cone narrows as the window grows and the median tracks the true vol.

## [1.53.0] - 2026-09-09

### Added
- SVI butterfly-arbitrage check (in `svi.py`): `svi_g` evaluates the
  Gatheral-Jacquier g-function of a slice (>= 0 everywhere iff no butterfly /
  density arbitrage), with `svi_butterfly_arbitrage` listing violating strikes
  and `svi_is_butterfly_free` the boolean. Cross-checked against the sign of the
  Breeden-Litzenberger density.

## [1.52.0] - 2026-09-09

### Added
- `theta_carry_report` (in `bookgreeks.py`): decomposes a book's net theta into
  the gamma-rent term (`-0.5 * Gamma * sigma^2 * S^2`) and a residual
  drift/financing carry, returning a `ThetaCarry`. With zero carry the theta is
  pure gamma rent; a non-zero rate produces the financing residual.

## [1.51.0] - 2026-09-09

### Added
- `vannavolga.py`: `VannaVolgaSmile` builds an FX smile from the three market
  quotes (ATM vol, 25-delta risk reversal, 25-delta butterfly), recovering the
  pillar vols/strikes and interpolating the vol at any strike. `pillar_vols`
  exposes the 25P/ATM/25C vols. Exact at the three pillars.

## [1.50.0] - 2026-09-09

### Added
- Kelly-criterion sizing (in `sizing.py`): `kelly_fraction_binary` (optimal
  stake for a binary bet from win probability and odds),
  `kelly_fraction_continuous` (growth-optimal leverage `mu / sigma^2`, with a
  fractional-Kelly multiplier), and `kelly_growth_rate` (expected log-growth at
  a given leverage, maximized at full Kelly).

## [1.49.0] - 2026-09-09

### Added
- `early_exercise_premium` (in `american.py`): decomposes the American price
  (Bjerksund-Stensland) into the European (BSM) value plus the early-exercise
  premium. The premium is zero for a no-dividend American call and positive for
  ITM puts and dividend-paying calls.

## [1.48.0] - 2026-09-09

### Added
- `displaced.py`: `displaced_diffusion_price` prices under Rubinstein's
  displaced-diffusion (shifted-lognormal) model — a Black-Scholes price on
  `S + shift` / `K + shift` with a rescaled vol. `shift = 0` recovers
  Black-Scholes; a positive shift allows negative strikes/spot and moves the
  skew toward normal-model behavior. Monte-Carlo verified.

## [1.47.0] - 2026-09-09

### Added
- `dividend_curve` (in `forward.py`): bootstraps an implied dividend-yield term
  structure from a multi-expiry option chain, applying `implied_forward` per
  expiry and returning sorted `(t, ForwardResult)` points. Recovers a known
  dividend term structure exactly.

## [1.46.0] - 2026-09-09

### Added
- `quanto.py`: `quanto_option` prices a foreign-asset option settled in domestic
  currency at a fixed exchange rate. The quanto adjustment shifts the carry by
  `-rho * sigma_asset * sigma_fx` and discounts at the domestic rate; `rho=0`
  removes it. Monte-Carlo verified.

## [1.45.0] - 2026-09-09

### Added
- `perpetual.py`: `perpetual_american` prices a no-expiry American option in
  exact closed form (Merton 1973), and `perpetual_exercise_boundary` returns the
  flat optimal-exercise spot. Matches the binomial American price at long
  maturity; a call with carry >= rate is never exercised early.

## [1.44.0] - 2026-09-09

### Added
- `gap_option` and `power_option` (in `exotics.py`): closed forms for a gap
  option (separate trigger and payoff strikes; Reiner-Rubinstein) and a power
  option (payoff on `S^power`; adjusted-drift/vol Black-Scholes). Equal gap
  strikes and `power=1` recover the vanilla option; the power form is
  Monte-Carlo verified.

## [1.43.0] - 2026-09-09

### Added
- `correlation_term_structure` (in `correlation.py`): implied correlation at each
  expiry from index and member vol term structures, returning
  `(expiry, rho)` pairs. Recovers a constant or maturity-varying correlation
  exactly.

## [1.42.0] - 2026-09-09

### Added
- `asian_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a Asian
  option by central finite differences on either closed form
  (`average="geometric"` Kemna-Vorst or `"arithmetic"` Turnbull-Wakeman). The
  arithmetic-average delta exceeds the geometric one, matching the price order.

## [1.41.0] - 2026-09-09

### Added
- `compound.py`: `compound_option` prices Geske (1979) compound options — an
  option on an option — for all four kinds (call/put-on-call/put), reusing the
  bivariate-normal CDF and solving the critical-spot exercise boundary by
  bisection. Cross-checked against Monte Carlo; call-on-call collapses to the
  vanilla as the first strike goes to zero.

## [1.40.0] - 2026-09-09

### Added
- `best_of_call` / `worst_of_call` (in `multiasset.py`): rainbow options on the
  max / min of two correlated assets, priced by Monte Carlo on correlated GBM
  (antithetic, seeded). Best-of + worst-of equals the sum of the two single-name
  calls (Stulz identity), which the tests verify.

## [1.39.0] - 2026-09-09

### Added
- `chooser.py`: `chooser_option` prices a simple chooser (Rubinstein 1991) in
  closed form — the holder picks call or put at a future date. Decomposes into a
  call to expiry plus a put on the discounted-forward strike expiring at the
  choice date; equals a straddle when the choice is at expiry.

## [1.38.0] - 2026-09-09

### Added
- `barrier_digital_mc` (in `montecarlo.py`): Monte Carlo for a cash-or-nothing
  digital contingent on a barrier condition (up/down, knock-in/knock-out) — pays
  the cash only if the option finishes in the money AND the barrier condition
  holds over the path. Knock-in + knock-out sums to the plain digital.

## [1.37.0] - 2026-09-09

### Added
- `sabr_local_vol` (in `localvol.py`): Dupire local volatility of a single SABR
  smile — builds the Hagan implied-vol smile on the forward and feeds it through
  the Dupire formula. A flat (nu->0, beta=1) SABR gives a constant local vol =
  alpha; a skewed smile gives the steeper-than-implied local skew.

## [1.36.0] - 2026-09-09

### Added
- `vegabucket.py`: `vega_buckets` groups a multi-expiry book's position-scaled
  vega into maturity buckets defined by upper-edge tenors, returning a
  `VegaBuckets` whose buckets sum to the net book vega — so a desk can see where
  its vol risk sits along the curve.

## [1.35.0] - 2026-09-09

### Added
- `lsm.py`: `bermudan_lsm` prices Bermudan/American options by Longstaff-Schwartz
  least-squares Monte Carlo — backward induction over exercise dates, regressing
  the discounted continuation value on a polynomial basis of spot (normal
  equations solved in pure Python). Converges to the binomial American value as
  the number of exercise dates grows.

## [1.34.0] - 2026-09-09

### Added
- `gramcharlier.py`: Corrado-Su skew/kurtosis-adjusted pricing. `corrado_su_call`
  / `corrado_su_price` add the first skewness and excess-kurtosis corrections to
  Black-Scholes via a Gram-Charlier expansion (skew=kurt=0 recovers BSM;
  kurtosis fattens the tails). `realized_skewness` and
  `realized_excess_kurtosis` estimate those moments from a return series.

## [1.33.0] - 2026-09-09

### Added
- `capped_cliquet_mc` (in `montecarlo.py`): Monte Carlo pricer for a locally-
  and globally-capped cliquet (ratchet) note — sums clipped periodic returns
  and clips the running total, with antithetic variates and a seed. Tightening
  either cap lowers the price; the global cap bounds the payoff.

## [1.32.0] - 2026-09-09

### Added
- `sizing.py`: hedge-quantity helpers. `delta_hedge_shares` zeros a book's net
  delta with the underlying; `neutralize` solves the units of a hedge option to
  move delta/gamma/vega to a target; `vega_neutral_quantity` and
  `gamma_neutral_quantity` are the zero-target shortcuts. All work off the net
  Greeks from `price_book`.

## [1.31.0] - 2026-09-09

### Added
- `VolSurface.forward_variance` / `forward_vol`: the forward (instantaneous-
  average) variance and volatility between two maturities, from the additive
  total-variance surface — the vol of a forward-starting option. Raises on a
  negative forward variance (calendar arbitrage).

## [1.30.0] - 2026-09-09

### Added
- `barrier_greeks` (in `exotics.py`): delta, gamma, vega, and theta (calendar)
  of a single-barrier option by central finite differences on the
  Reiner-Rubinstein price. A far knock-out matches the vanilla Greeks, and
  knock-in + knock-out delta equals the vanilla delta (in-out parity).

## [1.29.0] - 2026-09-09

### Added
- `cev.py`: Constant-Elasticity-of-Variance pricing (Schroder/Hull) for
  `0 <= beta < 1`, with a from-scratch noncentral chi-square CDF
  (`noncentral_chisq_cdf`) and regularized incomplete gamma. The chi-square
  summation starts at the Poisson mode so it is stable at large noncentrality.
  `cev_price` matches Black-Scholes when scaled, satisfies parity, and produces
  the leverage skew (lower beta -> richer downside puts). Cross-checked vs
  Monte Carlo.

## [1.28.0] - 2026-09-09

### Added
- `correlation.py`: index implied correlation. `implied_correlation` inverts
  the index-variance decomposition for the single common correlation consistent
  with a quoted index vol; `index_vol_from_correlation` is the forward map; and
  `dispersion_basket_vol` is the zero-correlation reference. The standard
  dispersion-trading measure.

## [1.27.0] - 2026-09-09

### Added
- `qmc.py`: quasi-Monte Carlo. `halton` generates low-discrepancy points (van
  der Corput radical inverse per prime base) and `european_qmc` prices a
  European option by deterministic Halton integration of the payoff — it
  converges to the Black-Scholes value several times faster than pseudo-random
  Monte Carlo at the same point count.

## [1.26.0] - 2026-09-09

### Added
- `book_second_order` (in `bookgreeks.py`): position-scaled net second-order
  Greeks across a book — vanna, vomma/volga, charm, veta, speed, zomma, color —
  returned as a `BookSecondOrder`. Complements the first-order net Greeks from
  `price_book`.

## [1.25.0] - 2026-09-09

### Added
- Documentation site. `docs/gen_api.py` generates `docs/api.md` by introspecting
  `quantforge.__all__` (signatures + docstrings, grouped by module) with zero
  dependencies; `docs/index.md` and `mkdocs.yml` wire up an mkdocs site.
- CI now runs `docs/gen_api.py --check`, and a test asserts the reference is in
  sync and every public callable/class is documented, so the docs can't drift.

## [1.24.0] - 2026-09-09

### Added
- `epsilon` (dividend rho) in `bsm.py`: analytic sensitivity of the option
  price to the continuous dividend yield, `dPrice/dq`. Negative for calls,
  positive for puts, and verified against finite differences.

## [1.23.0] - 2026-09-09

### Added
- `dv01.py`: `key_rate_dv01` computes bucketed (key-rate) DV01 for any book
  expressed as `price(zero_curve)`. Bumps each tenor independently (central or
  one-sided), reports per-bucket sensitivities normalized to 1bp plus the
  parallel DV01; the buckets sum to the parallel shift.

## [1.22.0] - 2026-09-09

### Added
- `overhedge.py`: super-replicate a cash-or-nothing digital with a tight
  vanilla spread. `digital_call_overhedge` / `digital_put_overhedge` return an
  `Overhedge` (spread cost as a conservative price, the fair digital value, and
  the cushion between them); `overhedge_payoff` gives the spread's terminal
  payoff, which dominates the digital everywhere and converges to it as the
  spread width shrinks.

## [1.21.0] - 2026-09-09

### Added
- `swaption_price` (in `rates.py`): European payer/receiver swaptions on the
  Bachelier model — the swap's PV annuity times a normal-model option on the
  forward swap rate (handles negative rates). Plus `annuity` (PV01) and
  `swaption_parity` (payer - receiver = annuity * (swap_rate - strike)).

## [1.20.0] - 2026-09-09

### Added
- `rates.py`: interest-rate caps, floors, and collars priced as Bachelier
  (normal-vol) caplet/floorlet strips, so they handle negative rates.
  `CapletPeriod` describes each accrual period; `cap_price`/`floor_price` sum
  the strip, `collar_price` is long-cap/short-floor, and `caplet_floorlet_parity`
  gives the check identity. Cap - floor at one strike equals the swap PV.

## [1.19.0] - 2026-09-09

### Added
- `spline.py`: a pure-stdlib natural cubic spline (`CubicSpline`, Thomas-solved
  tridiagonal moments, C2-continuous, clamped outside its range) and
  `SmileSpline`, a strike->implied-vol interpolator with flat extrapolation — a
  model-free alternative to SVI/SABR for a single smile.

## [1.18.0] - 2026-09-09

### Added
- `one_touch` / `no_touch` (in `exotics.py`): continuously-monitored touch
  binaries. `one_touch` pays cash if the barrier is ever reached (immediately
  on hit, the FX convention, or deferred to expiry); `no_touch` pays if it
  never is. Direction (up/down barrier) is inferred from `H` vs `S`.
  Cross-checked against a barrier-crossing Monte Carlo; touch + no-touch (paid
  at expiry) sum to the discounted cash.

## [1.17.0] - 2026-09-09

### Added
- `bjerksund_stensland_greeks`: delta, gamma, vega, theta (calendar), and rho of
  the Bjerksund-Stensland American price by central finite differences (the
  2002 closed form has no simple Greek expressions). For a no-dividend American
  call the Greeks equal the European BSM Greeks, as they must.

## [1.16.1] - 2026-09-09

### Changed
- `calibrate_svi` now uses a deterministic multi-start (several fixed seeds,
  keep the best) so it no longer stalls in the degenerate huge-`b` valley of
  the raw-SVI objective. On the example chain the worst per-expiry fit improved
  from rmse ~9e-4 to ~5e-6 and the assembled surface is calendar-arbitrage free.

### Added
- `examples/vol_surface.py`: end-to-end surface workflow (invert quotes -> fit
  SVI per expiry -> assemble `VolSurface` -> calendar check -> interpolate vol
  and extract Dupire local vol), with a smoke test.

## [1.16.0] - 2026-09-09

### Added
- `localvol.py`: Dupire local volatility. `dupire_local_vol` evaluates the
  Dupire formula from a call-price surface `C(K, T)` by finite differences;
  `local_vol_from_implied` wraps an implied-vol surface via Black-Scholes.
  Recovers a flat implied vol as a constant local vol and matches the analytic
  term-structure local variance `dw/dT`.

## [1.15.0] - 2026-09-09

### Added
- `strategy.py`: multi-leg option-strategy builders returning a `Book` (so net
  price/Greeks come from the existing engine): `vertical_spread`, `straddle`,
  `strangle`, `risk_reversal`, `butterfly`, `iron_condor`. Plus
  `payoff_at_expiry` / `payoff_profile` for the P&L diagram and `break_evens`
  (grid scan + bisection) for the zero-P&L spots.

## [1.14.0] - 2026-09-09

### Added
- `arithmetic_asian` (in `exotics.py`): closed-form arithmetic-average Asian via
  Turnbull-Wakeman moment matching (match the average's first two moments to a
  lognormal, then Black-Scholes). Complements the Monte Carlo pricer and agrees
  with it to a few cents; arithmetic value dominates the geometric Asian.

## [1.13.0] - 2026-09-09

### Added
- `multiasset.py`: two-asset options. `exchange_option` (Margrabe, exact),
  `spread_option` (Kirk approximation, puts via parity), and `basket_option`
  (Levy lognormal moment-match on the two-asset weighted sum). Exchange and
  basket match correlated-GBM Monte Carlo; basket single-asset reduces to
  Black-Scholes exactly.

## [1.12.0] - 2026-09-09

### Added
- `varswap.py`: model-free variance- and volatility-swap fair strikes by static
  option replication (Demeterfi-Derman-Kamani-Zou log-strip). `variance_swap_strike`
  integrates an OTM put/call strip weighted by `1/K^2` around the forward;
  `volatility_swap_strike` returns the `sqrt` proxy. Recovers `sigma^2` from a
  flat-vol Black-Scholes chain to strip-truncation error.

## [1.11.0] - 2026-09-09

### Added
- `density.py`: Breeden-Litzenberger risk-neutral density extraction from a
  call-price curve. `risk_neutral_density` (second strike-derivative, non-
  uniform grid), `risk_neutral_cdf` (first derivative), `price_from_density`
  (integrate any payoff against the recovered density), and
  `density_total_mass` (sanity check ~1). Recovers the lognormal pdf from a
  BSM curve and reprices vanilla and digital payoffs.

## [1.10.0] - 2026-09-09

### Added
- `merton.py`: `merton_jump_price` prices European options under Merton (1976)
  jump-diffusion as a Poisson-weighted sum of Black-Scholes prices with
  jump-adjusted volatility and carry. Collapses to Black-Scholes at zero jump
  intensity; the compensated drift keeps the forward a martingale so put-call
  parity holds exactly. Cross-checked against Monte Carlo.

## [1.9.0] - 2026-09-09

### Added
- `hedgesim.py`: `simulate_delta_hedge` Monte Carlos a discretely delta-hedged
  short option, returning the hedging-error distribution (mean, std, min, max).
  Supports hedging at a different vol than the realized path (`hedge_vol` vs
  `real_vol`) to study vol-mismatch P&L.
- Test suite gains a `slow` marker; run `pytest -m "not slow"` for a ~4s fast
  pass (deep Monte Carlo / tree cross-checks are marked slow).

## [1.8.0] - 2026-09-09

### Added
- `bachelier.py`: the Bachelier (normal) model — `bachelier_price`, analytic
  delta/gamma/vega, and `bachelier_implied_vol` (Newton + bisection). Prices
  options on a forward following arithmetic Brownian motion, so it handles
  negative forwards/strikes (rates and spread options) where the lognormal
  model breaks down. `sigma` here is the normal (absolute) volatility.

## [1.7.0] - 2026-09-09

### Added
- `heston.py`: `heston_price` prices European options under the Heston (1993)
  stochastic-volatility model via its characteristic function (Albrecher
  "little trap" form), integrated with a self-contained 64-point
  Gauss-Legendre rule (no SciPy). Puts follow from parity. Collapses to
  Black-Scholes as the vol-of-vol goes to zero.

## [1.6.0] - 2026-09-09

### Added
- `lookback.py`: continuously-monitored lookback options —
  `floating_strike_lookback` (Goldman-Sosin-Gatto; payoff against the realized
  extreme) and `fixed_strike_lookback` (Conze-Viswanathan; ordinary strike on
  the realized extreme). Both take the running extreme and cost of carry `b`,
  with the `b -> 0` singularity handled by a nudge. Cross-checked against Monte
  Carlo path max/min.

## [1.5.0] - 2026-09-09

### Added
- `hedging.py`: smile-aware delta. `smile_delta` returns the effective delta
  under the sticky-strike rule (equals BS delta) or the sticky-delta /
  sticky-moneyness rule (adds a `-vega * (dsigma/dk) / S` skew term).
  `skew_slope` finite-differences a supplied smile, and
  `smile_delta_from_smile` wires the two together.

## [1.4.0] - 2026-09-09

### Added
- `forwardstart.py`: `forward_start_price` prices forward-start options
  (strike fixed at a future date as a multiple of the then-spot) via
  Rubinstein's closed form, and `cliquet_price` values a cliquet/ratchet as a
  strip of consecutive forward-starts.

## [1.3.0] - 2026-09-09

### Added
- `forward.py`: `implied_forward` extracts the implied forward price and
  discount factor from a call/put chain via a put-call-parity least-squares
  fit (no volatility assumption), and backs out the implied rate and dividend
  yield.

## [1.2.0] - 2026-09-09

### Added
- `trinomial.py`: `trinomial_price` prices American/European options on a Boyle
  trinomial lattice (smoother convergence than the binomial tree), and
  `richardson_american` combines `n`/`2n` solves to cancel the leading O(1/n)
  error for a more accurate American price.

## [1.1.0] - 2026-09-09

### Added
- `surface.py`: `VolSurface` stitches per-expiry SVI smiles into a term
  structure, interpolates total variance linearly in maturity, and reports
  calendar arbitrage (total variance must be non-decreasing in `t` at each
  strike). `VolSurface.fit` calibrates one SVI slice per expiry.

## [1.0.0] - 2026-09-09

First stable release. The public API is now considered stable under SemVer.

### Added
- Single-sourced version: `pyproject.toml` reads `quantforge.__version__`.

### Summary of the 1.0 feature set
- **Pricing** — generalized Black-Scholes-Merton (stock, dividend, Black-76,
  FX via cost-of-carry `b`); European and American exercise.
- **Greeks** — analytic delta, gamma, vega, theta, rho, plus the second-order
  vanna, vomma/volga, charm, veta, speed, zomma, color.
- **Implied volatility** — robust Newton-with-bisection solver with
  arbitrage-band rejection.
- **American options** — Cox-Ross-Rubinstein binomial tree and the
  Bjerksund-Stensland (2002) closed form (~500x faster than the tree).
- **Exotics** — cash/asset-or-nothing digitals, single-barrier options
  (Reiner-Rubinstein, all four kinds, with rebate), geometric-Asian.
- **Volatility surfaces** — Gatheral raw SVI and SABR (Hagan expansion), each
  with calibration.
- **Monte Carlo** — GBM engine with antithetic and control-variate variance
  reduction; standard errors on every estimate.
- **Portfolio** — batch pricing, net Greeks, VaR / Expected Shortfall
  (parametric delta-gamma, historical, full-reprice MC), spot×vol stress grid.
- **Realized volatility** — close-to-close, EWMA, Parkinson, Garman-Klass,
  Rogers-Satchell, Yang-Zhang estimators.
- **Performance** — zero-dependency core (~1.4M prices/sec) with an optional
  NumPy vectorized fast path (~6x on large batches).
- **Tooling** — 182 tests, CI on Python 3.8/3.10/3.12, tag-triggered PyPI
  release via Trusted Publishing, reproducible benchmark suite.

## [0.1.0 - 0.13.0]

Pre-1.0 development. Each minor version added one major capability area:
BSM core (0.1), portfolio/CLI/CI (0.2), SVI surface (0.3), exotics (0.4),
Monte Carlo (0.5), VaR/ES (0.6), realized vol (0.7), second-order Greeks and
the release workflow (0.8), benchmarks (0.9), the NumPy fast path (0.10), the
stress grid (0.11), Bjerksund-Stensland American (0.12), and SABR (0.13).
