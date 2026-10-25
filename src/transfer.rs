//! Relativistic time & frequency transfer, and the relativistic corrections navigation and clock
//! comparison need — the layer that sits on top of `geo4d`'s TAI clock.
//!
//! Three things live here:
//!
//! - **GNSS relativistic corrections** — the constant factory offset, the periodic eccentricity
//!   term, the Sagnac term of signal propagation in the rotating Earth frame, and the Shapiro delay.
//! - **Time transfer** — the one-way coordinate-time transfer between two points (geometric light
//!   time + Shapiro + Sagnac) and the residual Sagnac term of a two-way link.
//! - **Proper-time integration** — the proper time a clock accumulates along a [`Trajectory`],
//!   integrating `dτ/dt − 1 = −(U + v²/2)/c²` with a [`GravityModel`] supplying `U`.
//!
//! All exact formula; the constants come from IS-GPS-200 / IERS.

use crate::constants::{C2, EARTH_ROTATION_RATE, GM_EARTH, SPEED_OF_LIGHT};
use crate::gravity::GravityModel;
use crate::relativity::fractional_rate_tcg;
use geo4d::{Duration, Ecef, Trajectory};

/// The constant fractional frequency offset built into GPS satellite clocks to compensate the net
/// (gravitational blueshift − second-order-Doppler redshift) rate relative to a geoid clock:
/// **−4.4647×10⁻¹⁰**. The clock is manufactured low by this fraction so it reads 10.23 MHz on orbit.
pub const GPS_CLOCK_RATE_OFFSET: f64 = -4.464_7e-10;

/// `F = −2√μ_E/c²`, the coefficient of the GNSS eccentricity clock correction, s·m^(−1/2)
/// (IS-GPS-200).
pub const GNSS_ECCENTRICITY_F: f64 = -4.442_807_633e-10;

/// The periodic relativistic clock correction `Δt_r = F·e·√a·sin E` (seconds) for a satellite whose
/// orbit has semi-major axis `a` (m) and eccentricity `e`, at eccentric anomaly `E` (rad). Reaches
/// ~±23 ns for GPS. Add it to the satellite clock reading (IS-GPS-200 user algorithm).
#[inline]
pub fn gnss_eccentricity_correction(semi_major_axis_m: f64, eccentricity: f64, eccentric_anomaly_rad: f64) -> f64 {
    GNSS_ECCENTRICITY_F * eccentricity * semi_major_axis_m.sqrt() * eccentric_anomaly_rad.sin()
}

/// The Sagnac correction `(ω/c²)(x_e·y_r − y_e·x_r)` (seconds) for a signal from `emitter` to
/// `receiver` (both ECEF) in the rotating Earth frame — up to ~133 ns for a ground receiver, and
/// sign-sensitive to east/west geometry. Never negligible in a full GNSS solution.
#[inline]
pub fn sagnac_correction(emitter: Ecef, receiver: Ecef) -> f64 {
    EARTH_ROTATION_RATE / C2 * (emitter.x * receiver.y - emitter.y * receiver.x)
}

/// The Shapiro (gravitational) signal-propagation delay past a body of gravitational constant `gm`
/// (m³/s²), for endpoints at distances `r1`, `r2` from the body and separated by `r12` (all metres):
/// `Δt = (2·gm/c³)·ln[(r1+r2+r12)/(r1+r2−r12)]`, seconds.
pub fn shapiro_delay(gm: f64, r1: f64, r2: f64, r12: f64) -> f64 {
    2.0 * gm / (C2 * SPEED_OF_LIGHT) * ((r1 + r2 + r12) / (r1 + r2 - r12)).ln()
}

/// The Shapiro delay near Earth for a signal between two ECEF points, seconds (uses `GM_EARTH`) —
/// a few tens of picoseconds for GNSS, included in precise processing.
#[inline]
pub fn earth_shapiro_delay(a: Ecef, b: Ecef) -> f64 {
    shapiro_delay(GM_EARTH, a.vec().norm(), b.vec().norm(), (a.vec() - b.vec()).norm())
}

/// The one-way coordinate-time transfer from `emitter` to `receiver` (both ECEF), seconds:
/// geometric light time `|x_r − x_e|/c`, plus the Earth [`Shapiro`](earth_shapiro_delay) and
/// [`Sagnac`](sagnac_correction) terms. This is what a clock-comparison or ranging link corrects for.
pub fn one_way_transfer(emitter: Ecef, receiver: Ecef) -> f64 {
    let geometric = (receiver.vec() - emitter.vec()).norm() / SPEED_OF_LIGHT;
    geometric + earth_shapiro_delay(emitter, receiver) + sagnac_correction(emitter, receiver)
}

/// The residual Sagnac term of a **two-way** link between `a` and `b`, seconds. In two-way transfer
/// the geometric and Shapiro terms cancel by reciprocity, but the Sagnac term reverses sign between
/// the two directions and so does **not** cancel — it is the dominant relativistic residual.
#[inline]
pub fn two_way_sagnac_residual(a: Ecef, b: Ecef) -> f64 {
    sagnac_correction(a, b)
}

/// The proper time minus coordinate time (TCG) a clock accumulates along a [`Trajectory`], seconds —
/// integrating `dτ/dt − 1 = −(U + v²/2)/c²` (gravitational redshift + second-order Doppler) with
/// `model` supplying the gravitational potential `U` at each point. Sampled at `step` and integrated
/// trapezoidally; `None` for an empty trajectory.
///
/// The trajectory's velocity is taken as the **GCRS (inertial) coordinate velocity** — as an orbital
/// ephemeris naturally is. For an Earth-fixed clock, add the `Ω×r` rotation velocity so the `v²/2`
/// term becomes the centrifugal potential.
pub fn proper_time_offset(trajectory: &Trajectory, model: &impl GravityModel, step: Duration) -> Option<f64> {
    let (start, end) = (trajectory.start()?, trajectory.end()?);
    let mut total = 0.0;
    let mut t = start;
    let mut prev: Option<(f64, geo4d::Epoch)> = None;
    loop {
        let sample = if t > end { end } else { t };
        let state = trajectory.state_at(sample)?;
        let u = model.potential(state.position);
        let rate = fractional_rate_tcg(u, state.speed_mps());
        if let Some((prev_rate, prev_t)) = prev {
            total += 0.5 * (prev_rate + rate) * (sample - prev_t).as_seconds();
        }
        prev = Some((rate, sample));
        if sample == end {
            break;
        }
        t = t + step;
    }
    Some(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SphericalHarmonicModel;
    use geo4d::{Epoch, StateVector, Vec3};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn gnss_corrections_have_the_right_scale() {
        // Eccentricity term: e = 0.01, a = 2.656e7, sin E = 1 → ~ −22.9 ns.
        let dtr = gnss_eccentricity_correction(2.656_175e7, 0.01, core::f64::consts::FRAC_PI_2);
        assert!(close(dtr, -22.9e-9, 0.5e-9), "{dtr}");
        // Sagnac: a satellite ~20 000 km east of a ground receiver → tens to ~133 ns, sign follows
        // the east/west geometry.
        let sat = Ecef::new(6_378_137.0, 20_000_000.0, 0.0);
        let rcv = Ecef::new(6_378_137.0, 0.0, 0.0);
        let sag = sagnac_correction(sat, rcv);
        assert!(sag.abs() > 10e-9 && sag.abs() < 200e-9, "{sag}");
        // Earth Shapiro between a GNSS satellite and the ground: tens of ps.
        let shap = earth_shapiro_delay(sat, rcv);
        assert!(shap > 1e-12 && shap < 1e-9, "{shap}");
    }

    #[test]
    fn cassini_round_trip_shapiro_delay() {
        // Cassini's 2002 conjunction (Will 2014, arXiv:1403.7377): the spacecraft at 8.43 AU, the
        // ray 1.6 R☉ from the Sun. Twice the one-way delay is 262.4 µs; Will's Eq. 63 shorthand,
        // ½(1 + γ)[240 − 20 ln(d²/r)] µs, gives 263.8.
        let (au, b): (f64, f64) = (149_597_870_700.0, 1.6 * 6.957e8);
        let (r1, r2) = (au, 8.43 * au);
        let r12 = (r1 * r1 - b * b).sqrt() + (r2 * r2 - b * b).sqrt();
        let round_trip = 2.0 * shapiro_delay(crate::constants::GM_SUN, r1, r2, r12) * 1e6;
        let will = 240.0 - 20.0 * (1.6f64.powi(2) / 8.43).ln();
        assert!(
            close(round_trip, 262.42, 0.01) && close(round_trip / will, 1.0, 0.006),
            "{round_trip} µs"
        );
    }

    #[test]
    fn one_way_transfer_is_dominated_by_light_time() {
        let sat = Ecef::new(26_000_000.0, 5_000_000.0, 3_000_000.0);
        let rcv = Ecef::new(6_378_137.0, 0.0, 0.0);
        let total = one_way_transfer(sat, rcv);
        let geometric = (sat.vec() - rcv.vec()).norm() / SPEED_OF_LIGHT;
        // ~0.07 s of light time, with the relativistic terms only ns-scale corrections.
        assert!(close(total, geometric, 1e-6) && geometric > 0.06 && geometric < 0.08);
        // Two-way keeps only the (small, non-cancelling) Sagnac residual.
        assert_eq!(two_way_sagnac_residual(sat, rcv), sagnac_correction(sat, rcv));
    }

    #[test]
    fn proper_time_offset_matches_the_gps_satellite_rate() {
        // A clock at GPS altitude moving at orbital speed: dτ/dTCG − 1 ≈ −2.5×10⁻¹⁰, so over 60 s it
        // gains ≈ −15 ns relative to coordinate time. Use a monopole field for U = GM/r.
        let monopole = SphericalHarmonicModel::new(GM_EARTH, 6_378_137.0, 0.0, 0, vec![1.0], vec![0.0]);
        let t0 = Epoch::from_unix_seconds(1_700_000_000.0);
        let s0 = StateVector::new(t0, Ecef::new(2.656_175e7, 0.0, 0.0), Vec3::new(0.0, 3874.0, 0.0));
        let traj = Trajectory::from_states([s0, s0.at(t0 + Duration::from_seconds(60.0))]).unwrap();
        let offset = proper_time_offset(&traj, &monopole, Duration::from_seconds(5.0)).unwrap();
        // Expected rate: −(GM/r + v²/2)/c² ≈ −2.505e-10; ×60 s ≈ −1.503e-8 s. Within ~2 % (r drifts).
        let expected = -(GM_EARTH / 2.656_175e7 + 0.5 * 3874.0_f64.powi(2)) / C2 * 60.0;
        assert!(close(offset, expected, expected.abs() * 0.02), "{offset} vs {expected}");
        assert!(offset < 0.0);
    }
}
