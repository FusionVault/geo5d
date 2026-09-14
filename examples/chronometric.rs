//! A 5D geodesy walkthrough: normal gravity and a geopotential, chronometric leveling (turning a
//! clock frequency ratio into a height difference), the relativistic clock rate, a proper-time
//! integration over a satellite pass, and the GNSS relativistic corrections.
//!
//! `cargo run --example chronometric`

use geo5d::{
    relativity, transfer, Duration, Ecef, Epoch, Geodetic, GeopotentialNumber, NormalGravity, SphericalHarmonicModel,
    StateVector, Trajectory, Vec3, W0_IHRS,
};

fn main() {
    // ── Normal gravity and the geopotential from a spherical-harmonic model ──────────────────────────
    let field = SphericalHarmonicModel::grs80_normal_field(); // swap for SphericalHarmonicModel::from_gfc(egm2008)
    let site = Geodetic::new(46.9, 7.45, 550.0); // Bern, ~550 m
    let g = NormalGravity::GRS80.at_geodetic(site);
    let w = field.geopotential(geo5d::WGS84.to_ecef(site));
    println!("Bern: normal gravity {:.5} m/s², geopotential W = {:.1} m²/s²", g, w);

    let c = GeopotentialNumber::from_potential(w, W0_IHRS);
    println!(
        "  geopotential number C = {:.1} m²/s²  →  dynamic height {:.2} m",
        c.value,
        c.dynamic_height(NormalGravity::GRS80.gamma_45())
    );

    // ── Chronometric leveling: two optical clocks compared over a fibre link ─────────────────────────
    // Clock B reads 3.1×10⁻¹⁷ higher in frequency than clock A. How much higher is it?
    let frac = 3.1e-17;
    let dh = relativity::height_difference(frac, g);
    let sigma_h = relativity::height_uncertainty(2e-18, g); // a 2×10⁻¹⁸ clock
    println!(
        "\nchronometric leveling: Δf/f = {:.1e}  →  Δh = {:.3} m  (±{:.1} cm at 2×10⁻¹⁸)",
        frac,
        dh,
        sigma_h * 100.0
    );

    // ── The relativistic clock rate relative to Terrestrial Time ─────────────────────────────────────
    let rate = relativity::fractional_rate_tt(w, geo5d::W0_IAU);
    println!(
        "Bern clock vs TT: fractional rate {:+.3e}  →  {:+.2} ns/day",
        rate,
        relativity::accumulated_offset_seconds(rate, Duration::DAY) * 1e9
    );

    // ── Proper time accumulated by a GPS satellite over a 60 s arc ───────────────────────────────────
    let t0 = Epoch::from_gregorian_utc(2026, 3, 14, 0, 0, 0.0);
    let sat = StateVector::new(t0, Ecef::new(2.656_175e7, 0.0, 0.0), Vec3::new(0.0, 3874.0, 0.0));
    let arc = Trajectory::from_states([sat, sat.at(t0 + Duration::from_seconds(60.0))]).unwrap();
    let offset = transfer::proper_time_offset(&arc, &field, Duration::from_seconds(5.0)).unwrap();
    println!(
        "\nGPS satellite proper−coordinate time over 60 s: {:+.2} ns",
        offset * 1e9
    );

    // ── GNSS relativistic corrections for a fix ──────────────────────────────────────────────────────
    let ecc = transfer::gnss_eccentricity_correction(2.656_175e7, 0.012, 1.2);
    let sagnac = transfer::sagnac_correction(Ecef::new(2.0e7, 1.2e7, 0.9e7), geo5d::WGS84.to_ecef(site));
    println!(
        "GNSS corrections: eccentricity {:+.1} ns, Sagnac {:+.1} ns",
        ecc * 1e9,
        sagnac * 1e9
    );

    // ── Geoid & gravity anomaly from a (here synthetic) disturbing field ─────────────────────────────
    // With a real EGM model these are the true geoid undulation and free-air anomaly; here we add one
    // small tesseral term to the normal field so the numbers are illustrative but the machinery — the
    // ICGEM `.gfc` parse, the disturbing potential, Bruns' formula — is exactly what you'd run.
    let perturbed = SphericalHarmonicModel::from_gfc(
        "earth_gravity_constant 3.986005e14\nradius 6378137.0\nmax_degree 3\n\
         gfc 0 0 1.0 0.0\ngfc 2 0 -0.484166774985e-3 0.0\ngfc 3 1 2.0e-6 -1.0e-6\n",
    )
    .unwrap();
    let n = perturbed.geoid_undulation(site, &field, &NormalGravity::GRS80, W0_IHRS);
    let dg = perturbed.gravity_anomaly(site, &field);
    println!(
        "\n(synthetic) geoid undulation {:.2} m, gravity anomaly {:+.1} mGal",
        n,
        dg * 1e5
    );
}
