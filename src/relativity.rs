//! Chronometric leveling and the relativistic clock-rate layer.
//!
//! General relativity ties the geopotential to the tick rate of a clock: a clock higher in the
//! gravity well (smaller geopotential `W`, in the geodesy convention where `W` is largest at the
//! geoid) runs **faster**. Comparing two clocks' frequencies therefore *measures* their geopotential
//! difference — "chronometric leveling" — and an optical clock at `10⁻¹⁸` resolves about **1 cm** of
//! height. The sign convention here is the physically-correct one and is unit-tested against the
//! Tokyo Skytree redshift experiment.
//!
//! All of this is exact, data-free formula. The one thing you must supply from outside is the actual
//! geopotential `W` (from a gravity model) or the local gravity `g` (from [`NormalGravity`](crate::NormalGravity) or a
//! measurement) — the crate does not carry a global gravity field.

use crate::constants::C2;
use geo4d::{Duration, StateVector};

// ── Chronometric leveling: clock frequency ratio ↔ geopotential ↔ height ─────────────────────────────

/// The geopotential difference `W_B − W_A` (m²/s²) between two clocks from their fractional frequency
/// difference `(f_B − f_A)/f`. Geodesy sign convention: the **higher** clock (smaller `W`) is the
/// **faster** one, so a positive frequency difference gives a negative geopotential difference.
///
/// ```
/// use geo5d::relativity::geopotential_difference;
/// // Clock B is 1×10⁻¹⁶ higher in frequency than clock A ⇒ B sits at smaller geopotential.
/// let dw = geopotential_difference(1e-16);
/// assert!((dw + 8.9876e16 * 1e-16).abs() < 1e3); // ΔW = −c²·(Δf/f) ≈ −8.99 m²/s²
/// ```
#[inline]
pub fn geopotential_difference(frac_freq_diff: f64) -> f64 {
    -C2 * frac_freq_diff
}

/// The height difference `h_B − h_A` (metres) between two clocks from their fractional frequency
/// difference `(f_B − f_A)/f` and the local gravity `g` — the chronometric-leveling result. The
/// higher clock is the faster one.
///
/// ```
/// use geo5d::relativity::height_difference;
/// // The field's headline: 1×10⁻¹⁸ of fractional frequency ≈ 1 cm of height (at g ≈ 9.81).
/// let dh = height_difference(1e-18, 9.806_65);
/// assert!((dh - 0.009_17).abs() < 1e-4); // ≈ 0.92 cm
/// ```
#[inline]
pub fn height_difference(frac_freq_diff: f64, gravity: f64) -> f64 {
    C2 * frac_freq_diff / gravity
}

/// The fractional frequency shift `Δf/f` a clock shows relative to a reference clock, given the
/// geopotential difference `ΔW = W_clock − W_ref` (m²/s²). A clock above the reference (smaller `W`)
/// is blue-shifted (positive).
#[inline]
pub fn fractional_frequency_shift(geopotential_difference: f64) -> f64 {
    -geopotential_difference / C2
}

/// The height uncertainty (metres) a chronometric-leveling measurement carries, propagating the
/// clock frequency-ratio uncertainty `σ(Δf/f)` through `σ_h = (c²/g)·σ(Δf/f)`.
#[inline]
pub fn height_uncertainty(frac_freq_sigma: f64, gravity: f64) -> f64 {
    C2 * frac_freq_sigma / gravity
}

// ── Synchronization / "clocks at different heights drift" (HFT, data centres, fibre) ─────────────────

/// The fractional rate difference `g·Δh/c²` between two clocks a height `Δh` apart near the surface
/// (flat-gravity approximation) — the higher clock runs faster. Over baselines beyond ~10 km use the
/// geopotential form [`fractional_frequency_shift`] with a real `ΔW` instead of flat `g·Δh`.
///
/// ```
/// use geo5d::relativity::sync_fractional_offset;
/// // Two data-centre clocks 100 m apart in height diverge at ~1.1×10⁻¹⁴.
/// let rate = sync_fractional_offset(100.0, 9.806_65);
/// assert!((rate - 1.09e-14).abs() < 1e-16);
/// ```
#[inline]
pub fn sync_fractional_offset(height_diff_m: f64, gravity: f64) -> f64 {
    gravity * height_diff_m / C2
}

/// The clock offset (seconds) accumulated over an elapsed span at a constant fractional rate
/// difference — `rate × elapsed`. Returned in `f64` seconds because for terrestrial rates the offset
/// is sub-nanosecond (below a nanosecond-integer [`Duration`]'s resolution).
#[inline]
pub fn accumulated_offset_seconds(fractional_rate: f64, elapsed: Duration) -> f64 {
    fractional_rate * elapsed.as_seconds()
}

// ── The relativistic clock-rate layer ────────────────────────────────────────────────────────────────

/// The fractional rate offset `dτ/dTT − 1 = −(W − W₀)/c²` of an Earth-fixed clock relative to
/// Terrestrial Time, given the clock's geopotential `W` and the datum `W₀` (both m²/s²). **Zero on
/// the geoid** (`W = W₀`); a clock below the geoid (larger `W`) runs slow, one above runs fast.
///
/// Returned as the small offset itself (not `1 + offset`), because the offset (~10⁻¹⁶) is at the
/// edge of `f64` resolution near 1.
///
/// ```
/// use geo5d::{relativity::fractional_rate_tt, W0_IAU};
/// assert_eq!(fractional_rate_tt(W0_IAU, W0_IAU), 0.0);          // a geoid clock ticks at the TT rate
/// let below = fractional_rate_tt(W0_IAU + 9.8, W0_IAU);         // ~1 m below the geoid
/// assert!(below < 0.0 && (below + 1.09e-16).abs() < 1e-17);     // runs slow by ~1.1×10⁻¹⁶
/// ```
#[inline]
pub fn fractional_rate_tt(geopotential_w: f64, w0: f64) -> f64 {
    -(geopotential_w - w0) / C2
}

/// The fractional rate offset `dτ/dTCG − 1 = −(U + v²/2)/c²` of a clock relative to Geocentric
/// Coordinate Time, from its gravitational potential `U` (m²/s², positive) and its coordinate speed
/// `v` (m/s). This is the general form (gravitational redshift + second-order Doppler) for a moving
/// clock — a satellite, an aircraft, or a ground clock (`v = Ω×r`).
#[inline]
pub fn fractional_rate_tcg(gravitational_potential: f64, speed_mps: f64) -> f64 {
    -(gravitational_potential + 0.5 * speed_mps * speed_mps) / C2
}

/// [`fractional_rate_tcg`] for a `geo4d` [`StateVector`], taking the clock's speed from the state and
/// the gravitational potential at its position from you (a gravity model). The bridge from a 4D state
/// to its relativistic clock rate.
#[inline]
pub fn fractional_rate_tcg_of(state: StateVector, gravitational_potential: f64) -> f64 {
    fractional_rate_tcg(gravitational_potential, state.speed_mps())
}

// ── Frame dragging: the gravitomagnetic field of a spinning Earth ───────────────────────────────────

/// Milliarcseconds per year in one rad/s, to state precession rates the way the literature does.
pub const MAS_PER_YEAR: f64 = 365.25 * 86_400.0 * 180.0 / core::f64::consts::PI * 3.6e6;

/// Lense–Thirring drag of an orbit's node, `2GJ/(c² a³ (1 − e²)^(3/2))`, rad/s, for a body of spin
/// angular momentum `j` (kg m²/s), semi-major axis `a` (m) and eccentricity `e`. Measured for
/// Earth to about one part in a thousand by LARES-2 and LAGEOS (Ciufolini et al., Nature 2026).
///
/// ```
/// use geo5d::{constants::EARTH_ANGULAR_MOMENTUM, relativity::{MAS_PER_YEAR, lense_thirring_node_rate}};
/// // LAGEOS: a = 12 270 km, e = 0.0045.
/// let rate = lense_thirring_node_rate(EARTH_ANGULAR_MOMENTUM, 12_270e3, 0.0045) * MAS_PER_YEAR;
/// assert!((rate - 30.67).abs() < 0.01); // milliarcseconds per year
/// ```
pub fn lense_thirring_node_rate(j: f64, a: f64, e: f64) -> f64 {
    2.0 * crate::constants::G_NEWTON * j / (C2 * a.powi(3) * (1.0 - e * e).powf(1.5))
}

/// Lense–Thirring drag of the perigee, `−6GJ cos i/(c² a³ (1 − e²)^(3/2))`, rad/s, at inclination
/// `inclination` (rad).
pub fn lense_thirring_perigee_rate(j: f64, a: f64, e: f64, inclination: f64) -> f64 {
    -3.0 * inclination.cos() * lense_thirring_node_rate(j, a, e)
}

/// Geodetic (de Sitter) precession of a gyroscope on a circular orbit of radius `a` (m) around a
/// body of gravitational parameter `gm`, orbit-averaged: `(3/2)(GM)^(3/2)/(c² a^(5/2))`, rad/s.
pub fn geodetic_precession(gm: f64, a: f64) -> f64 {
    1.5 * gm.powf(1.5) / (C2 * a.powf(2.5))
}

/// Frame-dragging precession of a gyroscope on a circular polar orbit of radius `a` (m),
/// orbit-averaged: `GJ/(2c²a³)`, rad/s (the effect Gravity Probe B measured).
pub fn gyroscope_frame_dragging(j: f64, a: f64) -> f64 {
    crate::constants::G_NEWTON * j / (2.0 * C2 * a.powi(3))
}

/// The relativistic advance of an orbit's periapsis, `6πGM/(c²a(1 − e²))` radians per orbit, for
/// semi-major axis `a` (m) and eccentricity `e` about a body of gravitational parameter `gm`: the
/// first term of [Will (2014)](https://arxiv.org/abs/1403.7377), Eq. 64, with `γ = β = 1`.
///
/// ```
/// use geo5d::{constants::GM_SUN, relativity::perihelion_advance};
/// // Mercury (JPL's mean elements at J2000): 42.98″ per century (Will 2014, Eq. 65).
/// let (au, arcsec) = (149_597_870_700.0, core::f64::consts::PI / 648_000.0);
/// let (a, e) = (0.387_099_27 * au, 0.205_635_93);
/// let period = 2.0 * core::f64::consts::PI * (a * a * a / GM_SUN).sqrt();
/// let per_century = perihelion_advance(GM_SUN, a, e) * 36_525.0 * 86_400.0 / period / arcsec;
/// assert!((per_century - 42.98).abs() < 0.005);
/// ```
pub fn perihelion_advance(gm: f64, a: f64, e: f64) -> f64 {
    6.0 * core::f64::consts::PI * gm / (C2 * a * (1.0 - e * e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::W0_IAU;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn chronometric_leveling_headline_numbers() {
        // 1 cm of height ↔ 1.091×10⁻¹⁸ fractional frequency (g = 9.80665).
        let frac_per_cm = sync_fractional_offset(0.01, 9.806_65);
        assert!(close(frac_per_cm, 1.0912e-18, 1e-21), "{frac_per_cm}");
        // ...and the inverse: 1×10⁻¹⁸ ↔ ~0.92 cm.
        assert!(close(height_difference(1e-18, 9.806_65), 0.009_17, 1e-5));
        // Geopotential ↔ frequency round trip.
        let dw = geopotential_difference(3.2e-17);
        assert!(close(fractional_frequency_shift(dw), 3.2e-17, 1e-30));
        // Uncertainty propagates linearly: a 1×10⁻¹⁸ clock gives ~1 cm height uncertainty.
        assert!(close(height_uncertainty(1e-18, 9.806_65), 0.009_17, 1e-5));
    }

    #[test]
    fn tokyo_skytree_redshift() {
        // Takamoto et al. 2020: two Sr optical clocks, Δh = 452.6 m on the Tokyo Skytree, measured a
        // redshift of Δν ≈ 21.18 Hz on the ⁸⁷Sr clock transition f₀ = 429.228 004 229 87 THz.
        let f0 = 429_228_004_229_873.0; // Hz
        let frac = sync_fractional_offset(452.6, 9.806_65);
        assert!(close(frac, 4.94e-14, 1e-15), "{frac}");
        let delta_nu = f0 * frac;
        assert!(close(delta_nu, 21.2, 0.3), "{delta_nu} Hz"); // matches the measured ~21.18 Hz
    }

    #[test]
    fn clock_rate_layer() {
        // A clock on the geoid ticks at the TT rate exactly.
        assert_eq!(fractional_rate_tt(W0_IAU, W0_IAU), 0.0);
        // A GPS satellite: gravitational blueshift dominates over the 2nd-order Doppler → net fast.
        // U_sat ≈ GM/r ≈ 3.986e14/2.656e7 ≈ 1.500e7 m²/s²; v ≈ 3874 m/s.
        let sat = fractional_rate_tcg(1.5006e7, 3874.0);
        let geoid = fractional_rate_tcg(6.2636e7, 465.0); // U+centrifugal≈W₀ at the geoid, v=Ω·R_eq
                                                          // The satellite clock runs faster than a geoid clock (the ~+4.5×10⁻¹⁰ GPS offset).
        assert!(sat - geoid > 4.0e-10 && sat - geoid < 5.0e-10, "{}", sat - geoid);
    }

    #[test]
    fn sync_drift_accumulates() {
        // 100 m height difference → ~0.94 ns/day of clock divergence.
        let rate = sync_fractional_offset(100.0, 9.806_65);
        let per_day = accumulated_offset_seconds(rate, Duration::DAY);
        assert!(close(per_day, 9.4e-10, 1e-11), "{per_day}");
    }

    #[test]
    fn clock_gradients_across_a_millimetre_and_a_centimetre() {
        // A clock Δh higher sits at ΔW = −gΔh and runs faster by gΔh/c². Bothwell et al. (2022,
        // arXiv:2109.12238) at g = 9.796 m/s² (USGS): 1.09×10⁻¹⁹ per mm, measured 0.98(23)×10⁻¹⁹.
        let jila = fractional_frequency_shift(-9.796 * 1e-3);
        assert!(
            close(jila, 1.09e-19, 5e-22) && (jila - 0.98e-19).abs() < 2.3e-20,
            "{jila}"
        );
        // Zheng et al. (2023, arXiv:2207.07145) at g = 9.803 m/s²: 10.9×10⁻¹⁹ per cm, measured
        // 12.4 ± 0.7 ± 2.5, 0.57σ away.
        let wisconsin = fractional_frequency_shift(-9.803 * 1e-2);
        assert!(close(wisconsin, 10.9e-19, 5e-21), "{wisconsin}");
        assert!(((12.4e-19 - wisconsin) / 0.7f64.hypot(2.5) / 1e-19 - 0.575).abs() < 0.001);
    }

    #[test]
    fn the_suns_redshift_seen_from_the_earth() {
        // González Hernández et al. (2020, arXiv:2009.10558): GM☉/(cR☉) = 636.31 m/s at the
        // photosphere, 633.35 at 1 AU, 633.10 at the Earth on 2010 November 25 (147.75 million km,
        // less the Earth's own 0.21; the text calls the distance correction +4 cm/s, but it lowers
        // the value); measured 638 ± 6 m/s.
        use crate::constants::{GM_EARTH, GM_SUN, SPEED_OF_LIGHT as C};
        let (r_sun, r_earth, au) = (6.957e8, 6.3781e6, 149_597_870_700.0); // IAU 2015 nominal
        let speed = |w_clock: f64, w_ref: f64| -C * fractional_frequency_shift(w_clock - w_ref);
        assert!(close(speed(GM_SUN / r_sun, 0.0), 636.31, 0.005));
        assert!(close(speed(GM_SUN / r_sun, GM_SUN / au), 633.35, 0.005));
        let seen = speed(GM_SUN / r_sun, GM_SUN / 147.75e9 + GM_EARTH / r_earth);
        assert!(close(seen, 633.10, 0.01) && (638.0 - seen).abs() < 6.0, "{seen}");
    }

    #[test]
    fn precession_of_orbits_and_gyroscopes() {
        use crate::constants::GM_SUN;
        let (au, arcsec) = (149_597_870_700.0, core::f64::consts::PI / 648_000.0);
        let century = 36_525.0 * 86_400.0;
        // The Earth–Moon system as a gyroscope round the Sun: "about 2 arcseconds per century"
        // (Will 2014, after Eq. 70).
        let de_sitter = geodetic_precession(GM_SUN, au) * century / arcsec;
        assert!(
            close(de_sitter, 1.919, 0.001) && close(de_sitter, 2.0, 0.1),
            "{de_sitter}"
        );
        // S2 round Sgr A* (GRAVITY 2020, arXiv:2004.07187): M = 4.261×10⁶ M☉, R₀ = 8246.7 pc,
        // a = 125.058 mas, e = 0.884649 (Table E.1) give 12.16′ per orbit, against their 12.1′
        // (Eq. 1); they measured f_SP = 1.10 ± 0.19 times it.
        let a = 125.058e-3 * 8246.7 * au; // an angle in arcseconds at a distance in parsecs is AU
        let s2 = perihelion_advance(4.261e6 * GM_SUN, a, 0.884_649) / arcsec / 60.0;
        assert!(close(s2, 12.156, 0.001) && close(s2, 12.1, 0.06), "{s2}");
        // A circle advances 6πGM/(c²a), and the advance diverges as e → 1 at fixed a.
        assert!(close(
            perihelion_advance(GM_SUN, au, 0.0),
            6.0 * core::f64::consts::PI * GM_SUN / (C2 * au),
            1e-25
        ));
        assert!(perihelion_advance(GM_SUN, au, 0.99) > 50.0 * perihelion_advance(GM_SUN, au, 0.0));
    }

    #[test]
    fn frame_dragging_lageos_and_gravity_probe_b() {
        use crate::constants::{EARTH_ANGULAR_MOMENTUM as J, GM_EARTH};
        // LAGEOS node: 30.67 mas/yr. Perigee drag is −3 cos i times the node rate.
        assert!(close(
            lense_thirring_node_rate(J, 12_270e3, 0.0045) * MAS_PER_YEAR,
            30.67,
            0.01
        ));
        let i = 109.8f64.to_radians();
        assert!(close(
            lense_thirring_perigee_rate(J, 12_270e3, 0.0045, i) / lense_thirring_node_rate(J, 12_270e3, 0.0045),
            -3.0 * i.cos(),
            1e-12
        ));
        // Gravity Probe B, polar orbit at a = 7027.4 km: 6603.9 mas/yr geodetic and 40.81 mas/yr
        // frame dragging. The published predictions are 6606.1 and 39.2: the latter is the drift
        // projected on the guide star IM Pegasi, 40.81 × cos 16.84° = 39.06 (Everitt et al. 2011;
        // measured −6601.8 ± 18.3 and −37.2 ± 7.2).
        assert!(close(
            geodetic_precession(GM_EARTH, 7_027.4e3) * MAS_PER_YEAR,
            6603.9,
            0.5
        ));
        assert!(close(
            gyroscope_frame_dragging(J, 7_027.4e3) * MAS_PER_YEAR,
            40.81,
            0.02
        ));
        assert!(close(
            geodetic_precession(GM_EARTH, 7_027.4e3) * MAS_PER_YEAR,
            6606.1,
            3.0
        ));
        assert!(close(
            gyroscope_frame_dragging(J, 7_027.4e3) * MAS_PER_YEAR * 16.84f64.to_radians().cos(),
            39.2,
            0.2
        ));
    }
}
