//! Physical height systems built on the geopotential number — the gravity-field-true vertical
//! coordinate. A geopotential number is not a length; each metric height is it divided by a gravity
//! value, and which gravity you divide by is what distinguishes the systems.

use crate::normal::NormalGravity;

/// A **geopotential number** `C = W₀ − W`, m²/s² — the work per unit mass against gravity from the
/// geoid up to a point. Non-negative above the geoid, path-independent, and the primitive from which
/// every metric height is derived.
///
/// ```
/// use geo5d::{GeopotentialNumber, NormalGravity, W0_IHRS};
/// // A point whose geopotential is 9806.2 m²/s² below the datum → ~1000 m of dynamic height.
/// let c = GeopotentialNumber::new(9_806.199_203);
/// let dyn_h = c.dynamic_height(NormalGravity::GRS80.gamma_45());
/// assert!((dyn_h - 1_000.0).abs() < 1e-6);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GeopotentialNumber {
    /// The geopotential number `C = W₀ − W`, m²/s².
    pub value: f64,
}

impl GeopotentialNumber {
    /// From a geopotential number in m²/s².
    #[inline]
    pub const fn new(value: f64) -> Self {
        GeopotentialNumber { value }
    }

    /// From the geopotential `W` at a point and the datum `W₀` (both m²/s²): `C = W₀ − W`.
    #[inline]
    pub fn from_potential(w: f64, w0: f64) -> Self {
        GeopotentialNumber { value: w0 - w }
    }

    /// The geopotential `W = W₀ − C` at the point, given the datum `W₀`.
    #[inline]
    pub fn potential(self, w0: f64) -> f64 {
        w0 - self.value
    }

    /// **Dynamic height** `H_dyn = C / γ₀`, dividing by a single fixed normal gravity `γ₀` (normal
    /// gravity at 45°). Equal geopotential numbers give equal dynamic heights — it is a true
    /// equipotential label, though not a geometric length.
    #[inline]
    pub fn dynamic_height(self, gamma_45: f64) -> f64 {
        self.value / gamma_45
    }

    /// **Orthometric height** `H = C / ḡ`, dividing by the mean *actual* gravity along the plumb
    /// line (which you supply — e.g. from a Helmert reduction). Referenced to the geoid.
    #[inline]
    pub fn orthometric_height(self, mean_gravity: f64) -> f64 {
        self.value / mean_gravity
    }

    /// **Normal height** `H* = C / γ̄`, dividing by the mean *normal* gravity along the ellipsoidal
    /// normal — no assumption about crustal density. Referenced to the quasigeoid. Since γ̄ depends
    /// on `H*` itself, this iterates to convergence (a few passes).
    pub fn normal_height(self, normal: &NormalGravity, lat_deg: f64) -> f64 {
        let s2 = lat_deg.to_radians().sin().powi(2);
        let g0 = normal.at_latitude(lat_deg);
        let b = normal.a * (1.0 - normal.f);
        let m = normal.omega * normal.omega * normal.a * normal.a * b / normal.gm;
        let coeff = 1.0 + normal.f + m - 2.0 * normal.f * s2;
        let mut h = self.value / g0; // first guess
        for _ in 0..5 {
            // Mean normal gravity from ellipsoid to height h (Heiskanen–Moritz).
            let ratio = h / normal.a;
            let gamma_bar = g0 * (1.0 - coeff * ratio + ratio * ratio);
            let next = self.value / gamma_bar;
            if (next - h).abs() < 1e-9 {
                h = next;
                break;
            }
            h = next;
        }
        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GRS80_GAMMA_45;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn dynamic_height_is_exact_division() {
        let c = GeopotentialNumber::new(GRS80_GAMMA_45 * 1000.0);
        assert!(close(c.dynamic_height(GRS80_GAMMA_45), 1000.0, 1e-9));
        // Potential round trip.
        let c2 = GeopotentialNumber::from_potential(62_636_853.4 - 500.0, 62_636_853.4);
        assert!(close(c2.value, 500.0, 1e-9));
        assert!(close(c2.potential(62_636_853.4), 62_636_853.4 - 500.0, 1e-6));
    }

    #[test]
    fn the_three_heights_agree_to_order_and_scale() {
        // A geopotential number of ~9806 m²/s² is ~1 km of height; the three systems differ by cm-dm.
        let c = GeopotentialNumber::new(9_806.0);
        let dyn_h = c.dynamic_height(NormalGravity::GRS80.gamma_45());
        let norm_h = c.normal_height(&NormalGravity::GRS80, 45.0);
        let orth_h = c.orthometric_height(9.806); // supplied mean gravity
        for h in [dyn_h, norm_h, orth_h] {
            assert!(close(h, 1000.0, 2.0), "{h}"); // all ~1 km, within a couple of metres
        }
        // Normal height uses a smaller mean gravity aloft → slightly larger than dynamic at 45°.
        assert!(norm_h > dyn_h - 0.5);
    }
}
