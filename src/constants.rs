//! Defining physical and geodetic constants, and the two conventional reference geopotentials.

/// Speed of light in vacuum, m/s (exact, SI).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// `c²`, m²/s².
pub const C2: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

/// Standard gravity, m/s² (the conventional `g₀`).
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// `L_G = W₀/c²`, the defining rate difference `1 − dTT/dTCG` (IAU 2000 Resolution B1.9).
pub const L_G: f64 = 6.969_290_134e-10;

/// `L_B`, the defining rate `1 − dTDB/dTCB` (IAU 2006 Resolution B3).
pub const L_B: f64 = 1.550_519_768e-8;

/// Earth's geocentric gravitational constant (WGS84 / IERS), m³/s².
pub const GM_EARTH: f64 = 3.986_004_418e14;

/// Earth's rotation rate (WGS84 / IS-GPS-200), rad/s.
pub const EARTH_ROTATION_RATE: f64 = 7.292_115_146_7e-5;

/// The Sun's gravitational constant (IAU 2009, TDB-compatible), m³/s².
pub const GM_SUN: f64 = 1.327_124_400_18e20;

/// Newtonian constant of gravitation (CODATA 2018), m³ kg⁻¹ s⁻².
pub const G_NEWTON: f64 = 6.674_30e-11;

/// Earth's spin angular momentum `J = Iω`, kg m²/s (the value used by the LAGEOS and LARES
/// frame-dragging analyses).
pub const EARTH_ANGULAR_MOMENTUM: f64 = 5.86e33;

// ── The reference geopotential W₀ — two conventional values, on purpose ──────────────────────────────
//
// There are two, and a careful crate must keep both with distinct roles:
//
//  - `W0_IAU`  = 62 636 856.0 m²/s²: the value frozen into the defining constant L_G (L_G = W0_IAU/c²).
//     It stays the IERS/IAU convention for the TT ↔ TCG relation and relativistic time scales.
//  - `W0_IHRS` = 62 636 853.4 m²/s²: the value the IAG adopted in 2015 (Resolution 1) as the reference
//     potential of the International Height Reference System — the height-datum zero point.
//
// They differ by 2.6 m²/s² (≈ 26 cm of height). Use `W0_IAU` for time-scale / clock-rate work and
// `W0_IHRS` for height-datum work; never silently mix them.

/// The reference geopotential frozen into `L_G` (IERS/IAU), m²/s². Use for time scales and clock rates.
pub const W0_IAU: f64 = 62_636_856.0;

/// The IHRS reference geopotential (IAG 2015 Resolution 1), m²/s². Use for the height datum.
pub const W0_IHRS: f64 = 62_636_853.4;

// ── GRS80 normal-gravity constants (Moritz 2000) ─────────────────────────────────────────────────────

/// GRS80 semi-major axis, m.
pub const GRS80_A: f64 = 6_378_137.0;
/// GRS80 flattening.
pub const GRS80_F: f64 = 1.0 / 298.257_222_101;
/// GRS80 geocentric gravitational constant, m³/s².
pub const GRS80_GM: f64 = 3.986_005e14;
/// GRS80 angular velocity, rad/s.
pub const GRS80_OMEGA: f64 = 7.292_115e-5;
/// GRS80 normal gravity at the equator, m/s².
pub const GRS80_GAMMA_E: f64 = 9.780_326_771_5;
/// GRS80 normal gravity at the pole, m/s².
pub const GRS80_GAMMA_P: f64 = 9.832_186_368_5;
/// GRS80 Somigliana constant `k = (b·γ_p − a·γ_e)/(a·γ_e)`.
pub const GRS80_K: f64 = 0.001_931_851_353;
/// GRS80 first eccentricity squared.
pub const GRS80_E2: f64 = 0.006_694_380_022_90;
/// GRS80 normal potential on the ellipsoid `U₀`, m²/s².
pub const GRS80_U0: f64 = 62_636_860.850;
/// GRS80 normal gravity at 45° latitude, m/s² — the fixed divisor for dynamic heights.
pub const GRS80_GAMMA_45: f64 = 9.806_199_203;

// ── WGS84 normal-gravity constants ───────────────────────────────────────────────────────────────────

/// WGS84 normal gravity at the equator, m/s².
pub const WGS84_GAMMA_E: f64 = 9.780_325_335_9;
/// WGS84 normal gravity at the pole, m/s².
pub const WGS84_GAMMA_P: f64 = 9.832_184_937_8;
/// WGS84 Somigliana constant.
pub const WGS84_K: f64 = 0.001_931_852_652;
/// WGS84 first eccentricity squared.
pub const WGS84_E2: f64 = 0.006_694_379_990_14;
/// WGS84 normal potential on the ellipsoid `U₀`, m²/s².
pub const WGS84_U0: f64 = 62_636_851.714_6;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn l_g_is_w0_iau_over_c2() {
        // The defining constant L_G is exactly W0_IAU/c² — that is why the old W₀ lives on in it.
        assert!((L_G - W0_IAU / C2).abs() < 1e-18, "{}", W0_IAU / C2);
    }

    #[test]
    fn the_two_w0_differ_by_about_26_cm() {
        let dh = (W0_IAU - W0_IHRS) / STANDARD_GRAVITY;
        assert!((dh - 0.265).abs() < 0.01, "{dh}"); // ≈ 0.26 m
    }
}
