//! Normal gravity on and above a level (equipotential) reference ellipsoid — the closed-form
//! Somigliana–Pizzetti field. This is the reference against which the real field's disturbances are
//! measured, and the gravity value that turns a geopotential difference into a height.

use crate::constants::*;
use geo4d::Geodetic;

/// The normal-gravity field of a level ellipsoid: the Somigliana constants plus what the second-order
/// height reduction needs. Use the [`GRS80`](Self::GRS80) or [`WGS84`](Self::WGS84) constants.
///
/// ```
/// use geo5d::NormalGravity;
/// // Normal gravity increases from equator to pole.
/// let g_eq = NormalGravity::GRS80.at_latitude(0.0);
/// let g_pole = NormalGravity::GRS80.at_latitude(90.0);
/// assert!((g_eq - 9.780_326_7715).abs() < 1e-9);
/// assert!((g_pole - 9.832_186_3685).abs() < 1e-9);
/// // At 45° it equals the dynamic-height reference γ₄₅.
/// assert!((NormalGravity::GRS80.at_latitude(45.0) - 9.806_199_203).abs() < 1e-6);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NormalGravity {
    /// Normal gravity at the equator, m/s².
    pub gamma_e: f64,
    /// Somigliana constant `k`.
    pub k: f64,
    /// First eccentricity squared.
    pub e2: f64,
    /// Semi-major axis, m.
    pub a: f64,
    /// Flattening.
    pub f: f64,
    /// Geocentric gravitational constant, m³/s².
    pub gm: f64,
    /// Angular velocity, rad/s.
    pub omega: f64,
    /// Normal potential on the ellipsoid surface `U₀`, m²/s² — the reference potential the geoid is
    /// measured against.
    pub u0: f64,
}

impl NormalGravity {
    /// GRS80 (the ellipsoid under most modern national datums).
    pub const GRS80: NormalGravity = NormalGravity {
        gamma_e: GRS80_GAMMA_E,
        k: GRS80_K,
        e2: GRS80_E2,
        a: GRS80_A,
        f: GRS80_F,
        gm: GRS80_GM,
        omega: GRS80_OMEGA,
        u0: GRS80_U0,
    };

    /// WGS84 (GPS).
    pub const WGS84: NormalGravity = NormalGravity {
        gamma_e: WGS84_GAMMA_E,
        k: WGS84_K,
        e2: WGS84_E2,
        a: 6_378_137.0,
        f: 1.0 / 298.257_223_563,
        gm: 3.986_004_418e14,
        omega: 7.292_115e-5,
        u0: WGS84_U0,
    };

    /// Normal gravity on the ellipsoid surface at a geodetic latitude (degrees) — the Somigliana
    /// closed form `γ = γ_e·(1 + k·sin²φ)/√(1 − e²·sin²φ)`.
    pub fn at_latitude(&self, lat_deg: f64) -> f64 {
        let s2 = lat_deg.to_radians().sin().powi(2);
        self.gamma_e * (1.0 + self.k * s2) / (1.0 - self.e2 * s2).sqrt()
    }

    /// Normal gravity at a geodetic latitude (degrees) and ellipsoidal height (metres), with the
    /// second-order free-air reduction — accurate to well below a µGal for aircraft/near-surface
    /// heights; at satellite heights use a full model.
    pub fn at(&self, lat_deg: f64, height_m: f64) -> f64 {
        let s2 = lat_deg.to_radians().sin().powi(2);
        let b = self.a * (1.0 - self.f);
        let m = self.omega * self.omega * self.a * self.a * b / self.gm;
        let g0 = self.at_latitude(lat_deg);
        let h = height_m;
        g0 * (1.0 - (2.0 / self.a) * (1.0 + self.f + m - 2.0 * self.f * s2) * h + (3.0 / (self.a * self.a)) * h * h)
    }

    /// Normal gravity at a geodetic position (latitude and ellipsoidal height), the `geo3d`/`geo4d`
    /// convenience.
    #[inline]
    pub fn at_geodetic(&self, position: Geodetic) -> f64 {
        self.at(position.lat_deg, position.height_m)
    }

    /// Normal gravity at 45° on the ellipsoid — the fixed divisor used for dynamic heights.
    #[inline]
    pub fn gamma_45(&self) -> f64 {
        self.at_latitude(45.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn somigliana_hits_the_reference_values() {
        assert!(close(NormalGravity::GRS80.at_latitude(0.0), GRS80_GAMMA_E, 1e-10));
        assert!(close(NormalGravity::GRS80.at_latitude(90.0), GRS80_GAMMA_P, 1e-9));
        assert!(close(NormalGravity::GRS80.gamma_45(), GRS80_GAMMA_45, 1e-6));
        // WGS84 differs from GRS80 only in the last digits — the two fields agree to ~0.1 mGal.
        assert!(close(NormalGravity::WGS84.at_latitude(0.0), WGS84_GAMMA_E, 1e-10));
        assert!(close(
            NormalGravity::WGS84.at_latitude(45.0),
            NormalGravity::GRS80.at_latitude(45.0),
            1e-5
        ));
    }

    #[test]
    fn height_reduction_decreases_gravity_at_the_free_air_rate() {
        let g0 = NormalGravity::GRS80.at_latitude(45.0);
        let g1000 = NormalGravity::GRS80.at(45.0, 1000.0);
        // Free-air gradient ≈ −0.3086 mGal/m ⇒ ≈ −0.3086 m/s² per 1000 m (×1e-5). ~0.3086e-2 m/s² drop.
        let drop = g0 - g1000;
        assert!(close(drop, 0.003_086, 1e-5), "{drop}");
        // Same via at_geodetic.
        assert!(close(
            NormalGravity::GRS80.at_geodetic(Geodetic::new(45.0, 0.0, 1000.0)),
            g1000,
            1e-12
        ));
    }
}
