# geo5d

The **geopotential & relativistic dimension** of geodesy, in pure `f64`.

> ⚠️ **Work in progress.** The core is implemented and internally validated — normal gravity, height
> systems, spherical-harmonic gravity synthesis (geoid, anomalies, deflections, gradient tensor),
> chronometric leveling, the relativistic clock-rate layer, and relativistic time/frequency transfer.
> Not yet published; the remaining steps are an external cross-check against GeographicLib on a real
> EGM model and the ns-level TDB−TT series. See the roadmap below.

Where [`geo3d`](https://crates.io/crates/geo3d) answers *where* (frames on an ellipsoid) and
[`geo4d`](https://crates.io/crates/geo4d) answers *when* (a nanosecond atomic clock), `geo5d` answers
**how gravity and relativity warp both** — the field the Earth's mass defines, and the relativity
that couples it to clocks. It builds on `geo4d` (which re-exports `geo3d`), so one dependency brings
frames, time, and this. Everything is exact, data-free formula and `#![forbid(unsafe_code)]`.

## What's in the slice

- **`gravity`** — the gravitational field as a spherical-harmonic expansion (`SphericalHarmonicModel`)
  behind a `GravityModel` trait: the Holmes–Featherstone stable recursion (with first and second
  derivatives), the geopotential `W`, the gravity vector, geoid undulation, gravity disturbance and
  anomaly, deflections of the vertical, the gravity-gradient (Marussi) tensor, GRS80/WGS84 normal
  fields, and an ICGEM `.gfc` parser. Bring your own coefficients — the algorithm is small, the data
  is not.
- **`NormalGravity`** — Somigliana normal gravity on and above a level ellipsoid (GRS80/WGS84), with
  the second-order height reduction.
- **`GeopotentialNumber`** — the gravity-field-true vertical coordinate `C = W₀ − W`, and the
  dynamic / normal / orthometric heights derived from it.
- **`relativity`** — **chronometric leveling** (a clock frequency ratio *is* a geopotential
  difference, hence a height difference — the 1 cm ↔ 10⁻¹⁸ regime of optical-clock geodesy), the
  synchronization-drift helpers, and the relativistic clock-rate layer (`dτ/dTT`, `dτ/dTCG`).
- **`transfer`** — relativistic time & frequency transfer: GNSS relativistic corrections (constant
  offset, eccentricity, Sagnac, Shapiro), one-/two-way time transfer, and proper-time integration of
  a `Trajectory` with a `GravityModel` supplying the potential — the full stack tying together.
- **`constants`** — the defining constants, including **both** conventional reference geopotentials
  `W₀` (the IAU value frozen into `L_G`, and the IHRS height-datum value) with distinct roles.

The spherical-harmonic engine is validated internally (synthesizing the GRS80 normal field from its
even zonal coefficients reproduces closed-form Somigliana normal gravity across latitudes, and the
gravity-gradient tensor satisfies the Laplace identity) **and against independent implementations on
the real EGM96 model**: geo5d's gravity `|∇V|` matches [pyshtools](https://shtools.github.io/) to
**0.05 µGal (5×10⁻¹¹)** across a global grid, and its geoid matches [GeographicLib](https://geographiclib.sourceforge.io/)
to a uniform **~9 cm** (a W₀/tide-system convention constant) away from high terrain.

Feed the engine a real ICGEM `.gfc` model and it gives an absolute `W`, which flows straight into the
chronometric-leveling and clock-rate layers.

> Note: `geoid_undulation` returns the **height anomaly / quasigeoid** (data-free), with the essential
> `−(W₀−U₀)/γ` zero-degree term. The classical geoid differs from it by the topographic (Bouguer)
> term over high terrain — up to ~1 m in the Himalaya — which needs a crustal-density hypothesis.

```rust
use geo5d::{relativity, GeopotentialNumber, NormalGravity};

// Two optical clocks compared over a fibre link differ by 3.2×10⁻¹⁷ in frequency.
// How far apart are they in height?
let g = NormalGravity::GRS80.at_latitude(47.0);
let dh = relativity::height_difference(3.2e-17, g);   // ≈ 0.29 m, higher clock = faster

// A surveyed geopotential number as a dynamic height.
let h = GeopotentialNumber::new(9_806.2).dynamic_height(NormalGravity::GRS80.gamma_45());
assert!((h - 1_000.0).abs() < 0.1);
```

The sign conventions are the physically-correct ones and are unit-tested against the **Tokyo Skytree**
redshift experiment (Takamoto et al. 2020: Δh = 452.6 m ⇒ Δν ≈ 21.2 Hz on the ⁸⁷Sr clock).

## Design & scope

- **Pure formula, bring-your-own field.** An *absolute* geopotential `W` comes from a gravity model
  you supply — the same "pass in the value" pattern `geo3d` uses for the ellipsoid. The crate carries
  the relativity and the normal field, not a global gravity grid.
- **Layered and self-contained**, like the rest of the family (`geo3d ← geo4d ← geo5d`), no new
  external dependencies for the core.

## Roadmap (researched, not yet built)

- A first-class classical-geoid option (height anomaly → geoid via the Bouguer term) for users who
  supply topographic heights, and higher-degree synthesis benchmarks.
- The ns-level periodic `TDB − TT` series and the `TCG/TCB/TDB` rate transforms, layered on `geo4d`'s
  TAI epoch (the linear rates are here in `constants`; the periodic series is the remaining piece).

Deliberately **out of scope**: X-ray pulsar navigation (novel but inherently data-coupled — better as
its own crate) and Conformal Geometric Algebra (elegant, but it linearizes the spherical/Euclidean
layer, not the ellipsoidal one that is geodesy's actual difficulty).

## License

Apache-2.0, see [LICENSE](LICENSE).
