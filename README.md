# geo5d

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
| `relativity` | Clocks as height meters (chronometric levelling), clock rates in TT and TCG ([IERS Conventions 2010](https://www.iers.org/TN36)), and frame dragging: the Lense–Thirring drag of satellite orbits and the gyroscope precessions Gravity Probe B measured ([Everitt et al. 2011](https://arxiv.org/abs/1105.3456)) |
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
- **Frame dragging:** LAGEOS's node drag (30.67 mas/yr) and Gravity Probe B's spherical-Earth
  drifts (6,638 and 41.1 mas/yr).

**Not included yet:**
- a classical-geoid option from topographic heights;
- the periodic TDB − TT series.
