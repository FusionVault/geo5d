# Changelog

All notable changes to `geo5d`. Plain semver; tags `vX.Y.Z` at the published commit.

## 0.5.2 — 2026-10-25

- **New:** `relativity::perihelion_advance`, `6πGM/(c²a(1 − e²))` per orbit.
- **Tests:** published values for the existing formulas, each read from its paper: JILA's and
  Wisconsin's clock gradients, the solar redshift seen from the Earth, Cassini's Shapiro delay,
  the lunar de Sitter precession, Mercury's perihelion and S2's Schwarzschild precession.
- `cargo fmt` over the crate, which CI checks.

## 0.5.1

- Gravity Probe B check uses its real orbit (a = 7027.4 km) and the published predictions
  (6606.1 and 39.2 mas/yr, the latter projected on the guide star).

## 0.5.0 — 2026-10-24

- **New:** frame dragging in `relativity`: `lense_thirring_node_rate`, `lense_thirring_perigee_rate`, `geodetic_precession`, `gyroscope_frame_dragging` and `MAS_PER_YEAR`; constants `G_NEWTON` and `EARTH_ANGULAR_MOMENTUM`. Checked against LAGEOS and Gravity Probe B.

## 0.4.1 — 2026-10-24

- **Licence:** now MIT (was Apache-2.0).
- **Docs:** a shorter README with a figure, linked sources and a tested example. It is now also the crate's front page on docs.rs.

## 0.4.0 — 2026-10-23

The spherical-harmonic synthesis now comes from the [`geopotential`](https://crates.io/crates/geopotential) crate, which is re-exported. That is the same Holmes–Featherstone algorithm, with pole-safe first and second derivatives from neighbouring orders throughout. The two engines agreed to ~1e-15 in `V` on the real degree-300 GOCO06s model, and every existing test passes unchanged.

- **New:** `SphericalHarmonicModel::from_coefficients` and `coefficients()`, to build a model from any `geopotential::Coefficients`, e.g. a time-variable ICGEM model read and evaluated at an epoch with the `gravity-field` crate.
- **Breaking:** the serialized form of `SphericalHarmonicModel` (with the `serde` feature) is now a `geopotential::Coefficients` plus `omega`. Deserializing checks that the coefficient arrays match the degree.
- **Breaking:** `from_gfc` refuses `max_degree` above 2700 (`GfcError::Unsupported`), the synthesis limit, instead of failing later.

## 0.3.0 — 2026-10-23

Fixes near the poles and in the `.gfc` reader, from a cross-check against an independent synthesis (the `geopotential` crate) on the real degree-300 GOCO06s model. Away from the poles the two agree to ~1e-15 in `V` and ~1e-6 µGal in `|g|`. With these fixes they agree near the poles too.

- **Fix:** the gravity vector was **NaN within ~1 cm of the rotation axis**, and so was everything built on it (`gravity`, `gravity_magnitude`, `gravity_disturbance`, the tensor). Its horizontal components also lost precision within ~1 m of the axis. Colatitude came from `acos(z/r)`, which rounds to 0 there. It is now computed with `atan2`, and a point exactly on the axis is evaluated 1 nm off it.
- **Fix:** the zonal (m = 0) θ-derivatives are now taken from the order-1 column, `dP̄n0/dθ = −√(n(n+1)/2)·P̄n1`. The old form divided a vanishing difference by `sinθ` and lost about 2 digits per decade of θ approaching the poles; the new one is exact at the poles.
- **Breaking:** `from_gfc` refuses time-variable models (`gfct`/`trnd`/`dot`/`acos`/`asin` records, e.g. GOCO06s) and `norm unnormalized` files with the new `GfcError::Unsupported`. It used to read `gfct` as static and drop the trend and periodic terms without a warning. Evaluate a time-variable model at your epoch first and build it with `SphericalHarmonicModel::new`.

## 0.2.0 — 2026-09-14

Precision fix to the geoid, from a cross-check against GeographicLib and pyshtools on the real EGM96 model. The gravity synthesis was already exact — geo5d matches pyshtools' `|∇V|` to **0.05 µGal** (5×10⁻¹¹) across a global grid — but the geoid formula was missing its zero-degree term.

- **Breaking:** `SphericalHarmonicModel::geoid_undulation` now takes the adopted geoid potential `w0` and applies the zero-degree term `−(W₀ − U₀)/γ` (leaving it out biased the result by ~0.5 m). With the term, geo5d agrees with GeographicLib on the EGM96 geoid to a uniform ~9 cm (a W₀/tide-system convention constant) away from high terrain. It is documented as the **height anomaly / quasigeoid** (data-free); the classical geoid differs by the Bouguer term over mountains.
- **Breaking:** `NormalGravity` gains a `u0` field (the ellipsoid's normal potential), with `WGS84_U0` added to `constants`.

## 0.1.0 — 2026-09-14

Initial release — the geopotential & relativistic dimension of geodesy, on top of `geo4d`/`geo3d`. **Work in progress**: the core is implemented and internally validated; an external cross-check against GeographicLib on a real EGM model and the ns-level TDB−TT series are still to come.

- `gravity` — spherical-harmonic synthesis (`SphericalHarmonicModel` behind a `GravityModel` trait): the Holmes–Featherstone stable recursion with first and second derivatives, the geopotential `W`, the gravity vector, geoid undulation, gravity disturbance and anomaly, deflections of the vertical, the gravity-gradient (Marussi) tensor, GRS80/WGS84 normal fields, and an ICGEM `.gfc` parser. Validated against Somigliana normal gravity, the Laplace-trace identity, and the analytic `(n−1)` anomaly.
- `relativity` — chronometric leveling (clock frequency ratio ↔ geopotential ↔ height, with the correct sign convention and uncertainty propagation) and the relativistic clock-rate layer (`dτ/dTT`, `dτ/dTCG`). Validated against the Tokyo Skytree redshift experiment.
- `transfer` — GNSS relativistic corrections (constant offset, eccentricity, Sagnac, Shapiro), one-/two-way time transfer, and proper-time integration over a `geo4d` `Trajectory` with a `GravityModel`.
- `normal` — Somigliana normal gravity (GRS80/WGS84) with the second-order height reduction.
- `height` — `GeopotentialNumber` and dynamic / normal / orthometric heights.
- `constants` — defining constants including both conventional reference geopotentials `W₀`.
- Re-exports the `geo4d` (and thus `geo3d`) API; optional `serde` feature.
