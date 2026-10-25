#![doc = include_str!("../README.md")]
#![allow(deprecated)]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

#[deprecated(
    since = "0.6.0",
    note = "moved to the metrology crate (GRS80_F, GRS80_GM, W0_IAU, W0_IHRS, L_G, L_B, EARTH_ANGULAR_MOMENTUM)"
)]
pub mod constants;
#[deprecated(
    since = "0.6.0",
    note = "use geopotential::Coefficients with normal::bruns, Partials::gravity_neu and Partials::deflections"
)]
pub mod gravity;
mod height;
mod normal;
#[deprecated(
    since = "0.6.0",
    note = "clocks moved to gravimetry::clocks; perihelion_advance and geodetic_precession to gravitational::energy; frame dragging to gravitics::gem"
)]
pub mod relativity;
#[deprecated(
    since = "0.6.0",
    note = "moved to gravimetry::clocks (eccentricity_correction, sagnac, shapiro_delay, one_way_transfer, proper_time_offset)"
)]
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
