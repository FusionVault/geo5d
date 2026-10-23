//! The **geopotential & relativistic dimension** of geodesy, in pure `f64`. *(Work in progress — a
//! thin vertical slice.)*
//!
//! Where [`geo3d`](geo4d::geo3d) answers *where* (frames on an ellipsoid) and [`geo4d`] answers
//! *when* (a high-precision atomic clock), `geo5d` answers **how gravity and relativity warp both** —
//! the field the Earth's mass defines and the relativity that couples it to clocks:
//!
//! - [`NormalGravity`] — Somigliana normal gravity on and above a level ellipsoid (GRS80/WGS84).
//! - [`GeopotentialNumber`] — the gravity-field-true vertical coordinate `C = W₀ − W`, and the
//!   dynamic / normal / orthometric heights derived from it.
//! - [`relativity`] — **chronometric leveling** (a clock frequency ratio *is* a geopotential
//!   difference, hence a height difference — the 1 cm ↔ 10⁻¹⁸ regime of optical-clock geodesy) and
//!   the relativistic clock-rate layer (`dτ/dTT`, `dτ/dTCG`).
//! - [`transfer`] — relativistic time & frequency transfer: GNSS relativistic corrections, one-/two-
//!   way time transfer (Shapiro, Sagnac, light time), and proper-time integration of a
//!   [`Trajectory`] with a [`GravityModel`] supplying the potential.
//! - [`constants`] — the defining constants, including **both** conventional reference geopotentials
//!   `W₀` (the IAU value frozen into `L_G`, and the IHRS height-datum value) with distinct roles.
//!
//! Everything here is exact, data-free formula. The one thing it does not carry is a global gravity
//! field: an *absolute* geopotential `W` comes from a model you supply — the same "bring-your-own,
//! behind a value/trait" pattern `geo3d` uses for the ellipsoid. It builds on `geo4d` (and thus
//! `geo3d`), which it re-exports, and is `#![forbid(unsafe_code)]`.
//!
//! ```
//! use geo5d::{relativity, GeopotentialNumber, NormalGravity};
//!
//! // Two optical clocks compared over a fibre link differ by 3.2×10⁻¹⁷ in frequency.
//! // How far apart are they in height? (local gravity ≈ 9.81 m/s²)
//! let g = NormalGravity::GRS80.at_latitude(47.0);
//! let dh = relativity::height_difference(3.2e-17, g);
//! assert!(dh > 0.29 && dh < 0.30); // ≈ 29 cm, the higher clock being the faster one
//!
//! // A surveyed geopotential number, expressed as a dynamic height.
//! let dyn_h = GeopotentialNumber::new(9_806.2).dynamic_height(NormalGravity::GRS80.gamma_45());
//! assert!((dyn_h - 1_000.0).abs() < 0.1);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod constants;
pub mod gravity;
mod height;
mod normal;
pub mod relativity;
pub mod transfer;

pub use constants::{L_G, W0_IAU, W0_IHRS};
pub use gravity::{GravityModel, SphericalHarmonicModel};
pub use height::GeopotentialNumber;
pub use normal::NormalGravity;

// Re-export the layers geo5d builds on, so one dependency brings frames, time, and this.
pub use geo4d;
pub use geo4d::{Duration, Ecef, Ellipsoid, Epoch, Geodetic, StateVector, Trajectory, Vec3, GRS80, WGS84};
/// The spherical-harmonic engine behind [`gravity`], re-exported so its types match.
pub use geopotential;
