# geo5d

**Folded into [physica](https://github.com/FusionVault/physica).** Version 0.6.0 is the final
release: every item is deprecated with a note naming its new home, and nothing else changed. The
same formulas, tested against the same published values and more, now live in

| geo5d | physica |
|---|---|
| `NormalGravity`, GRS80 and WGS84 constants | [`geopotential::normal`](https://docs.rs/geopotential) |
| `GeopotentialNumber` (heights) | `geopotential::normal::NormalGravity::{dynamic_height, normal_height}` |
| `SphericalHarmonicModel` (gravity, geoid, deflections, gradients) | `geopotential::Coefficients`, `Partials`, `normal::bruns` |
| `relativity` (clocks) and `transfer` | [`gravimetry::clocks`](https://docs.rs/gravimetry) |
| `relativity` (periapsis, de Sitter) | [`gravitational::energy`](https://docs.rs/gravitational) |
| `relativity` (frame dragging) | [`gravitics::gem`](https://docs.rs/gravitics) |
| `constants` | [`metrology`](https://docs.rs/metrology) |

```bash
cargo add geopotential gravimetry
```

The original description follows.

---

Gravity and relativity for geodesy. It covers normal gravity, heights from the geopotential,
gravity from spherical-harmonic models, and the link between clocks and heights: a clock runs
faster the higher it is. Pure Rust, `f64`, no unsafe code.

Part of a family:
- [`geo3d`](https://crates.io/crates/geo3d) says *where* (positions on the ellipsoid);
- [`geo4d`](https://crates.io/crates/geo4d) says *when* (precise time);
- `geo5d` adds gravity and relativity, and re-exports both.

![The geoid from GOCO06s](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/geoid_goco06s.svg)

*The geoid from the GOCO06s model, computed with
[`geopotential`](https://crates.io/crates/geopotential), the engine geo5d uses for gravity.*

## Install

```bash
cargo add geo5d
```

## Example: how much higher is a clock?

```rust
use geo5d::{relativity, GeopotentialNumber, NormalGravity};

// Two optical clocks, compared over a fibre, differ in frequency by 3.2 parts in 10¹⁷.
let g = NormalGravity::GRS80.at_latitude(47.0);
let height = relativity::height_difference(3.2e-17, g);
assert!((0.29..0.30).contains(&height)); // 29 cm: the faster clock is higher

// A surveyed geopotential number, as a height in metres.
let h = GeopotentialNumber::new(9_806.2).dynamic_height(NormalGravity::GRS80.gamma_45());
assert!((h - 1000.0).abs() < 0.1);
```

## What is inside

| Item | What it does |
|---|---|
| `NormalGravity` | Gravity of the reference ellipsoid (GRS80, WGS84) on and above it ([Moritz 2000](https://doi.org/10.1007/s001900050278)) |
| `GeopotentialNumber` | Heights from the geopotential: dynamic, normal and orthometric |
| `SphericalHarmonicModel` | Gravity, geoid, gravity anomalies, deflections and the gradient tensor from a model such as an ICGEM `.gfc` file, through [`geopotential`](https://crates.io/crates/geopotential) |
| `relativity` | Clocks as height meters (chronometric levelling), clock rates in TT and TCG ([IERS Conventions 2010](https://www.iers.org/TN36)), frame dragging: the Lense–Thirring drag of satellite orbits and the gyroscope precessions Gravity Probe B measured ([Everitt et al. 2011](https://arxiv.org/abs/1105.3456)), and the relativistic advance of a periapsis ([Will 2014](https://arxiv.org/abs/1403.7377)) |
| `transfer` | Relativistic corrections for GNSS and time transfer (Sagnac, Shapiro), and proper time along a trajectory |
| `constants` | The defining constants, including both conventional values of W₀ |

`geoid_undulation` returns the height anomaly (the quasigeoid). Over high mountains the classical
geoid differs from it by up to about 1 m, because that needs an assumed crustal density.

## Example: frame dragging

A spinning Earth drags space around with it. The orbit of the LAGEOS satellite turns by about 31
milliarcseconds a year because of it, which LAGEOS and LARES-2 measure to about 0.1%.

```rust
use geo5d::constants::EARTH_ANGULAR_MOMENTUM;
use geo5d::relativity::{MAS_PER_YEAR, lense_thirring_node_rate};

let drag = lense_thirring_node_rate(EARTH_ANGULAR_MOMENTUM, 12_270e3, 0.0045);
println!("{:.2} mas/yr", drag * MAS_PER_YEAR); // 30.67
```

## How it is checked

- **Relativity:** the clock signs and sizes match the Tokyo Skytree experiment: 452.6 m of height
  gives 21.2 Hz on a strontium clock
  ([Takamoto et al. 2020](https://doi.org/10.1038/s41566-020-0619-8)).
- **Gravity:** on the EGM96 model, it matches [pyshtools](https://shtools.github.io/SHTOOLS/) to
  0.05 µGal.
- **Geoid:** it matches [GeographicLib](https://geographiclib.sourceforge.io/) to a constant
  offset of about 9 cm, which comes from convention.
- **GOCO06s:** on the degree-300 model, its synthesis matches `geopotential` to about 1e-15.
- **Frame dragging:** LAGEOS's node drag (30.67 mas/yr) and Gravity Probe B's predicted drifts,
  6,606 mas/yr geodetic and 39.2 mas/yr frame dragging projected on its guide star.
- **Clocks across a sample:** the redshift gradients JILA and Wisconsin measured, 1.09×10⁻¹⁹ per
  mm ([Bothwell et al. 2022](https://arxiv.org/abs/2109.12238)) and 10.9×10⁻¹⁹ per cm
  ([Zheng et al. 2023](https://arxiv.org/abs/2207.07145)); the Sun's redshift seen from the
  Earth, 633.10 m/s ([González Hernández et al. 2020](https://arxiv.org/abs/2009.10558)).
- **Orbits and light:** Mercury's 42.98″ per century and the Moon's de Sitter precession of about
  2″ per century ([Will 2014](https://arxiv.org/abs/1403.7377)); S2's 12.1′ per orbit round
  Sgr A* ([GRAVITY 2020](https://arxiv.org/abs/2004.07187)); Cassini's 262.4 µs round-trip
  Shapiro delay.

**Not included yet:**
- a classical-geoid option from topographic heights;
- the periodic TDB − TT series.
