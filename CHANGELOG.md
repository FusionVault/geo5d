# Changelog

All notable changes to `geo5d`. Plain semver; tags `vX.Y.Z` at the published commit.

## 0.1.0 — 2026-09-14

Initial release — the geopotential & relativistic dimension of geodesy, on top of `geo4d`/`geo3d`. **Work in progress**: the core is implemented and internally validated; an external cross-check against GeographicLib on a real EGM model and the ns-level TDB−TT series are still to come.

- `gravity` — spherical-harmonic synthesis (`SphericalHarmonicModel` behind a `GravityModel` trait): the Holmes–Featherstone stable recursion with first and second derivatives, the geopotential `W`, the gravity vector, geoid undulation, gravity disturbance and anomaly, deflections of the vertical, the gravity-gradient (Marussi) tensor, GRS80/WGS84 normal fields, and an ICGEM `.gfc` parser. Validated against Somigliana normal gravity, the Laplace-trace identity, and the analytic `(n−1)` anomaly.
- `relativity` — chronometric leveling (clock frequency ratio ↔ geopotential ↔ height, with the correct sign convention and uncertainty propagation) and the relativistic clock-rate layer (`dτ/dTT`, `dτ/dTCG`). Validated against the Tokyo Skytree redshift experiment.
- `transfer` — GNSS relativistic corrections (constant offset, eccentricity, Sagnac, Shapiro), one-/two-way time transfer, and proper-time integration over a `geo4d` `Trajectory` with a `GravityModel`.
- `normal` — Somigliana normal gravity (GRS80/WGS84) with the second-order height reduction.
- `height` — `GeopotentialNumber` and dynamic / normal / orthometric heights.
- `constants` — defining constants including both conventional reference geopotentials `W₀`.
- Re-exports the `geo4d` (and thus `geo3d`) API; optional `serde` feature.
