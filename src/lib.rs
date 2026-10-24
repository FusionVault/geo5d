#![doc = include_str!("../README.md")]

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
