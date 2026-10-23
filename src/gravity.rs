//! The gravitational field as a spherical-harmonic expansion — the source of an *absolute*
//! geopotential `W`, and of the geoid, gravity anomalies and deflections that physical geodesy is
//! built on.
//!
//! The algorithm is small; the coefficients are enormous (EGM2008 is ~4.7 million of them). So this
//! module ships the **math** — the Holmes–Featherstone stable recursion and the synthesis — and
//! takes the coefficients as a value you provide, exactly as `geo3d` takes an `Ellipsoid` rather
//! than baking one in. Load them from an ICGEM `.gfc` file with [`SphericalHarmonicModel::from_gfc`],
//! or use the built-in GRS80/WGS84 normal fields.
//!
//! The whole engine is validated end-to-end with **no external data**: synthesizing the GRS80 normal
//! field from its even zonal coefficients reproduces Somigliana [`NormalGravity`] to a µGal.

use crate::normal::NormalGravity;
use geo4d::{Ecef, Geodetic, Vec3};

/// A source of the Earth's gravitational potential — the boundary between the pure relativity/height
/// math and a bring-your-own gravity field.
pub trait GravityModel {
    /// Geocentric gravitational constant `GM`, m³/s².
    fn gm(&self) -> f64;
    /// Reference radius `a` of the expansion, m.
    fn reference_radius(&self) -> f64;
    /// The **gravitational** potential `V` at an ECEF position, m²/s² (excludes the centrifugal term).
    fn potential(&self, position: Ecef) -> f64;
    /// The **gravitational** acceleration `∇V` at an ECEF position, m/s² (excludes centrifugal).
    fn acceleration(&self, position: Ecef) -> Vec3;
}

/// The large scaling factor for the Holmes–Featherstone recursion: sectorials near the pole would
/// otherwise underflow `f64` at high degree, so the whole recursion runs pre-multiplied by this and
/// is un-scaled at the end.
const HF_SCALE: f64 = 1e280;
const HF_SCALE_INV: f64 = 1e-280;

/// The potential `V` and its first and second spherical partial derivatives at a point.
#[derive(Clone, Copy, Debug)]
struct Partials {
    v: f64,
    vr: f64,
    vrr: f64,
    vt: f64,
    vtt: f64,
    vrt: f64,
    vl: f64,
    vll: f64,
    vrl: f64,
    vtl: f64,
}

/// Geocentric radius, colatitude and longitude of an ECEF position.
fn spherical(position: Ecef) -> (f64, f64, f64) {
    let rho = position.x.hypot(position.y);
    let r = rho.hypot(position.z);
    // atan2, not acos(z/r), which rounds to θ = 0 within ~1 cm of the axis. A point exactly on the
    // axis is taken 1 nm off it: the field changes by far less than f64 resolution there.
    let theta = rho.max(1e-9).atan2(position.z);
    let lambda = position.y.atan2(position.x);
    (r, theta, lambda)
}

/// A gravitational field as a set of fully-normalized (4π, no Condon–Shortley phase) Stokes
/// coefficients `C̄_nm, S̄_nm` to a maximum degree, with its own `GM`, reference radius `a` and
/// rotation rate `ω`.
///
/// ```
/// use geo5d::gravity::SphericalHarmonicModel;
/// use geo5d::NormalGravity;
/// use geo4d::GRS80;
///
/// // The GRS80 normal field, synthesized, reproduces Somigliana normal gravity.
/// let model = SphericalHarmonicModel::grs80_normal_field();
/// let on_ellipsoid = GRS80.to_ecef(geo4d::Geodetic::new(45.0, 0.0, 0.0));
/// let g = model.gravity_magnitude(on_ellipsoid);
/// assert!((g - NormalGravity::GRS80.at_latitude(45.0)).abs() < 1e-5);
/// ```
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SphericalHarmonicModel {
    gm: f64,
    a: f64,
    omega: f64,
    nmax: usize,
    /// C̄_nm at triangular index `n(n+1)/2 + m`.
    cnm: Vec<f64>,
    /// S̄_nm at the same index.
    snm: Vec<f64>,
}

/// Why a `.gfc` model could not be parsed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GfcError {
    /// A required header key (`earth_gravity_constant`, `radius`, `max_degree`) was missing.
    MissingHeader(&'static str),
    /// A `gfc` coefficient line was malformed.
    BadCoefficientLine,
    /// The file is valid but holds something this static model cannot represent: time-variable
    /// records (`gfct`/`trnd`/`dot`/`acos`/`asin`) or unnormalized coefficients.
    Unsupported(&'static str),
}

impl core::fmt::Display for GfcError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            GfcError::MissingHeader(k) => write!(f, "missing .gfc header key: {k}"),
            GfcError::BadCoefficientLine => f.write_str("malformed .gfc coefficient line"),
            GfcError::Unsupported(what) => write!(f, "unsupported .gfc content: {what}"),
        }
    }
}
impl std::error::Error for GfcError {}

#[inline]
fn idx(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

impl SphericalHarmonicModel {
    /// From explicit coefficients. `cnm`/`snm` are indexed by `n(n+1)/2 + m` and must hold at least
    /// `(nmax+1)(nmax+2)/2` entries.
    pub fn new(gm: f64, a: f64, omega: f64, nmax: usize, cnm: Vec<f64>, snm: Vec<f64>) -> Self {
        let need = (nmax + 1) * (nmax + 2) / 2;
        assert!(
            cnm.len() >= need && snm.len() >= need,
            "coefficient arrays too short for nmax"
        );
        SphericalHarmonicModel {
            gm,
            a,
            omega,
            nmax,
            cnm,
            snm,
        }
    }

    /// The GRS80 normal field (the even zonal coefficients of the level ellipsoid). Synthesizing it
    /// reproduces Somigliana [`NormalGravity::GRS80`]; used as the reference field for the disturbing
    /// potential, geoid and gravity anomaly.
    pub fn grs80_normal_field() -> Self {
        Self::normal_field(
            crate::constants::GRS80_GM,
            crate::constants::GRS80_A,
            crate::constants::GRS80_OMEGA,
            [
                -0.484_166_774_985e-3,
                0.790_303_733_511e-6,
                -0.168_724_961_151e-8,
                0.346_052_468_394e-10,
            ],
        )
    }

    /// The WGS84 normal field.
    pub fn wgs84_normal_field() -> Self {
        Self::normal_field(
            3.986_004_418e14,
            6_378_137.0,
            7.292_115e-5,
            [
                -0.484_166_774_985e-3,
                0.790_303_733_534e-6,
                -0.168_724_961_15e-8,
                0.346_052_468e-10,
            ],
        )
    }

    /// Build a zonal-only field from `C̄_20, C̄_40, C̄_60, C̄_80` (a level ellipsoid's normal potential).
    fn normal_field(gm: f64, a: f64, omega: f64, even_zonals: [f64; 4]) -> Self {
        let nmax = 8;
        let n = (nmax + 1) * (nmax + 2) / 2;
        let mut cnm = vec![0.0; n];
        let snm = vec![0.0; n];
        cnm[idx(0, 0)] = 1.0;
        cnm[idx(2, 0)] = even_zonals[0];
        cnm[idx(4, 0)] = even_zonals[1];
        cnm[idx(6, 0)] = even_zonals[2];
        cnm[idx(8, 0)] = even_zonals[3];
        SphericalHarmonicModel {
            gm,
            a,
            omega,
            nmax,
            cnm,
            snm,
        }
    }

    /// Parse an ICGEM `.gfc` model file: the `earth_gravity_constant`, `radius` and `max_degree`
    /// header keys, then `gfc  n  m  C̄  S̄  [σC σS]` coefficient lines. Rotation rate defaults to the
    /// GRS80 value (override with [`with_omega`](Self::with_omega)).
    ///
    /// Static models only: a time-variable model (GOCO06s, EIGEN-6S4, …) or one with unnormalized
    /// coefficients is refused with [`GfcError::Unsupported`] rather than read wrongly. Evaluate a
    /// time-variable model at your epoch first and build it with [`new`](Self::new).
    pub fn from_gfc(text: &str) -> Result<Self, GfcError> {
        let (mut gm, mut a, mut nmax) = (None, None, None);
        for line in text.lines() {
            let mut it = line.split_whitespace();
            match it.next() {
                Some("earth_gravity_constant") => gm = it.next().and_then(parse_f64),
                Some("radius") => a = it.next().and_then(parse_f64),
                Some("max_degree") => nmax = it.next().and_then(|s| s.parse::<usize>().ok()),
                Some("norm") if it.next() == Some("unnormalized") => {
                    return Err(GfcError::Unsupported("unnormalized coefficients"));
                }
                Some("gfct" | "trnd" | "dot" | "acos" | "asin") => {
                    return Err(GfcError::Unsupported("time-variable model (gfct/trnd/acos/asin)"));
                }
                _ => {}
            }
        }
        let gm = gm.ok_or(GfcError::MissingHeader("earth_gravity_constant"))?;
        let a = a.ok_or(GfcError::MissingHeader("radius"))?;
        let nmax = nmax.ok_or(GfcError::MissingHeader("max_degree"))?;
        let size = (nmax + 1) * (nmax + 2) / 2;
        let mut cnm = vec![0.0; size];
        let mut snm = vec![0.0; size];
        for line in text.lines() {
            let mut it = line.split_whitespace();
            if it.next() != Some("gfc") {
                continue;
            }
            let n: usize = it
                .next()
                .and_then(|s| s.parse().ok())
                .ok_or(GfcError::BadCoefficientLine)?;
            let m: usize = it
                .next()
                .and_then(|s| s.parse().ok())
                .ok_or(GfcError::BadCoefficientLine)?;
            let c = it.next().and_then(parse_f64).ok_or(GfcError::BadCoefficientLine)?;
            let s = it.next().and_then(parse_f64).ok_or(GfcError::BadCoefficientLine)?;
            if n <= nmax && m <= n {
                cnm[idx(n, m)] = c;
                snm[idx(n, m)] = s;
            }
        }
        Ok(SphericalHarmonicModel {
            gm,
            a,
            omega: crate::constants::GRS80_OMEGA,
            nmax,
            cnm,
            snm,
        })
    }

    /// The same model with a different rotation rate (rad/s).
    pub fn with_omega(mut self, omega: f64) -> Self {
        self.omega = omega;
        self
    }

    /// Maximum degree of the expansion.
    #[inline]
    pub fn max_degree(&self) -> usize {
        self.nmax
    }

    /// The geopotential `W = V + ½ω²(x² + y²)` (gravitational plus centrifugal) at an ECEF position,
    /// m²/s². This is the potential a clock senses and the relativity layer consumes.
    pub fn geopotential(&self, position: Ecef) -> f64 {
        self.potential(position) + 0.5 * self.omega * self.omega * (position.x * position.x + position.y * position.y)
    }

    /// The gravity vector `g = ∇W` (gravitational plus centrifugal) in ECEF, m/s².
    pub fn gravity(&self, position: Ecef) -> Vec3 {
        let grav = self.acceleration(position);
        let cent = Vec3::new(
            self.omega * self.omega * position.x,
            self.omega * self.omega * position.y,
            0.0,
        );
        grav + cent
    }

    /// The magnitude of gravity `|∇W|` at an ECEF position, m/s².
    #[inline]
    pub fn gravity_magnitude(&self, position: Ecef) -> f64 {
        self.gravity(position).norm()
    }

    /// The disturbing potential `T = W − U = V − U_grav` at a position, m²/s² — the model's
    /// gravitational potential minus a reference normal field's (the centrifugal parts cancel).
    pub fn disturbing_potential(&self, position: Ecef, normal: &SphericalHarmonicModel) -> f64 {
        self.potential(position) - normal.potential(position)
    }

    /// The height anomaly / quasigeoid `ζ = T/γ − (W₀ − U₀)/γ` at a geodetic position, metres — the
    /// separation between the geopotential surface `W = W₀` and the reference ellipsoid. `T` is the
    /// disturbing potential against the `normal` field, `gravity` the normal gravity (carrying its
    /// `u0`), and `w0` the adopted geoid potential (e.g. [`W0_IHRS`](crate::W0_IHRS)).
    ///
    /// The `−(W₀−U₀)/γ` **zero-degree term** is essential and easy to forget: leaving it out biases
    /// the result by ~0.5 m. This is the data-free Molodensky quantity; the classical geoid `N`
    /// differs from the height anomaly by the topographic (Bouguer) term over high terrain — up to
    /// ~1 m in the Himalaya — which needs a crustal-density hypothesis and is out of scope here.
    pub fn geoid_undulation(
        &self,
        position: Geodetic,
        normal: &SphericalHarmonicModel,
        gravity: &NormalGravity,
        w0: f64,
    ) -> f64 {
        let ecef = geo4d::WGS84.to_ecef(position);
        let gamma = gravity.at_geodetic(position);
        self.disturbing_potential(ecef, normal) / gamma - (w0 - gravity.u0) / gamma
    }

    /// The gravity disturbance `δg = |g_model| − γ` at a geodetic position, m/s² (positive where the
    /// real field is stronger than normal). Convert to mGal by ×10⁵.
    pub fn gravity_disturbance(&self, position: Geodetic, gravity: &NormalGravity) -> f64 {
        let ecef = geo4d::WGS84.to_ecef(position);
        self.gravity_magnitude(ecef) - gravity.at_geodetic(position)
    }

    /// The gravity anomaly `Δg = −∂T/∂r − (2/r)·T` (the spherical/Molodensky approximation) at a
    /// geodetic position, m/s² — the classic reduced quantity of physical geodesy, from the
    /// disturbing potential `T = V − U` against the reference `normal` field. Convert to mGal by ×10⁵.
    pub fn gravity_anomaly(&self, position: Geodetic, normal: &SphericalHarmonicModel) -> f64 {
        let (r, theta, lambda) = spherical(geo4d::WGS84.to_ecef(position));
        let (pm, pn) = (self.synthesize(r, theta, lambda), normal.synthesize(r, theta, lambda));
        let (tval, tr) = (pm.v - pn.v, pm.vr - pn.vr);
        -tr - 2.0 / r * tval
    }

    /// The deflections of the vertical `(ξ, η)` — north and east components, **radians** — at a
    /// geodetic position: the tilt between the true plumb line and the ellipsoid normal, from the
    /// horizontal slope of the disturbing potential. Multiply by 206 265 for arc-seconds. A
    /// longitude-independent (zonal) field gives `η = 0`.
    pub fn deflections(
        &self,
        position: Geodetic,
        normal: &SphericalHarmonicModel,
        gravity: &NormalGravity,
    ) -> (f64, f64) {
        let (r, theta, lambda) = spherical(geo4d::WGS84.to_ecef(position));
        let (pm, pn) = (self.synthesize(r, theta, lambda), normal.synthesize(r, theta, lambda));
        let (tt, tl) = (pm.vt - pn.vt, pm.vl - pn.vl);
        let g = gravity.at_geodetic(position);
        (tt / (g * r), -tl / (g * r * theta.sin()))
    }

    /// The gravity-gradient (Marussi) tensor `∂²V/∂xᵢ∂xⱼ` of the **gravitational** potential at an
    /// ECEF position, in the local East-North-Up frame, s⁻² — the quantity a gradiometer (GOCE,
    /// cold-atom) measures. Symmetric, and trace-free outside the masses (Laplace's equation), which
    /// is a good self-check. Rows/columns are ordered East, North, Up. Exactly on the rotation axis,
    /// where East and North are undefined, it is evaluated 1 nm off the axis toward λ = 0 and holds
    /// to ~1e-5 relative.
    pub fn gradient_tensor_enu(&self, position: Ecef) -> [[f64; 3]; 3] {
        let (r, theta, _lambda) = spherical(position);
        let p = self.synthesize(r, theta, _lambda);
        let (u, t) = (theta.sin(), theta.cos());
        let (r2, cot) = (r * r, t / u);
        // Physical second-derivative components in the spherical (r, θ, λ) frame.
        let g_rr = p.vrr;
        let g_tt = p.vtt / r2 + p.vr / r;
        let g_ll = p.vll / (r2 * u * u) + cot / r2 * p.vt + p.vr / r;
        let g_rt = p.vrt / r - p.vt / r2;
        let g_rl = p.vrl / (r * u) - p.vl / (r2 * u);
        let g_tl = (p.vtl - cot * p.vl) / (r2 * u);
        // Rotate into East-North-Up: E = λ̂, N = −θ̂, U = r̂.
        let (ee, nn, uu) = (g_ll, g_tt, g_rr);
        let (en, eu, nu) = (-g_tl, g_rl, -g_rt);
        [[ee, en, eu], [en, nn, nu], [eu, nu, uu]]
    }

    /// Synthesis: the potential `V` and all its first and second spherical partial derivatives at a
    /// geocentric radius `r`, colatitude `θ` and longitude `λ`. Holmes–Featherstone stable recursion
    /// with global scaling.
    #[allow(clippy::needless_range_loop)] // n also drives the column recursion, not just col1[n]
    fn synthesize(&self, r: f64, theta: f64, lambda: f64) -> Partials {
        let (t, u) = (theta.cos(), theta.sin()); // t = cosθ, u = sinθ
        let rho = self.a / r;
        let nmax = self.nmax;

        // Σ_m accumulators (after applying the trig for each order).
        let (mut v, mut vr, mut vrr) = (0.0, 0.0, 0.0);
        let (mut vt, mut vtt, mut vrt) = (0.0, 0.0, 0.0);
        let (mut vl, mut vll, mut vrl, mut vtl) = (0.0, 0.0, 0.0, 0.0);

        // The order-1 column and its θ-derivative: the zonal derivatives below are taken from it,
        // because the generic form divides a vanishing difference by sinθ near the poles.
        let mut col1 = vec![(0.0, 0.0); nmax + 1];
        if nmax >= 1 {
            let (mut p_n, mut p_nm1) = (self.sectorial(1, u), 0.0);
            for (n, slot) in col1.iter_mut().enumerate().skip(1) {
                *slot = (p_n, dp_dtheta(n, 1, t, u, p_n, p_nm1));
                if n < nmax {
                    let p_next = self.column_next(n + 1, 1, t, p_n, p_nm1);
                    p_nm1 = p_n;
                    p_n = p_next;
                }
            }
        }

        for m in 0..=nmax {
            let mut rho_n = rho.powi(m as i32);
            // Per-order C̄ / S̄ lumped sums, weighted by the radial/derivative factor each needs.
            let (mut am, mut bm) = (0.0, 0.0); // × P
            let (mut am_r, mut bm_r) = (0.0, 0.0); // × (n+1) P
            let (mut am_rr, mut bm_rr) = (0.0, 0.0); // × (n+1)(n+2) P
            let (mut am_t, mut bm_t) = (0.0, 0.0); // × dP/dθ
            let (mut am_tt, mut bm_tt) = (0.0, 0.0); // × d²P/dθ²
            let (mut am_rt, mut bm_rt) = (0.0, 0.0); // × (n+1) dP/dθ

            let mut p_n = self.sectorial(m, u); // n = m
            let mut p_nm1 = 0.0; // n = m-1 → 0
            for n in m..=nmax {
                let c = self.cnm[idx(n, m)];
                let s = self.snm[idx(n, m)];
                let (dp, ddp) = if n == 0 {
                    (0.0, 0.0)
                } else if m == 0 {
                    // dP̄_n0/dθ = −√(n(n+1)/2)·P̄_n1, exact at the poles; differentiated once more
                    // for the second derivative.
                    let k = -((n * (n + 1)) as f64 / 2.0).sqrt();
                    (k * col1[n].0, k * col1[n].1)
                } else {
                    let dp = dp_dtheta(n, m, t, u, p_n, p_nm1);
                    // d²P̄/dθ² = −(t/u)·dP̄/dθ − (n(n+1) − m²/u²)·P̄.
                    let ddp = -(t / u) * dp - ((n * (n + 1)) as f64 - (m * m) as f64 / (u * u)) * p_n;
                    (dp, ddp)
                };
                let w = rho_n;
                let nn1 = (n + 1) as f64;
                let nn2 = (n + 1) as f64 * (n + 2) as f64;
                am += w * c * p_n;
                bm += w * s * p_n;
                am_r += nn1 * w * c * p_n;
                bm_r += nn1 * w * s * p_n;
                am_rr += nn2 * w * c * p_n;
                bm_rr += nn2 * w * s * p_n;
                am_t += w * c * dp;
                bm_t += w * s * dp;
                am_tt += w * c * ddp;
                bm_tt += w * s * ddp;
                am_rt += nn1 * w * c * dp;
                bm_rt += nn1 * w * s * dp;

                if n < nmax {
                    let p_next = self.column_next(n + 1, m, t, p_n, p_nm1);
                    p_nm1 = p_n;
                    p_n = p_next;
                    rho_n *= rho;
                }
            }

            let mf = m as f64;
            let (cml, sml) = ((mf * lambda).cos(), (mf * lambda).sin());
            let even = |a: f64, b: f64| a * cml + b * sml; // combos that keep the trig
            let dl = |a: f64, b: f64| mf * (-a * sml + b * cml); // ∂/∂λ combos
            v += even(am, bm);
            vr += even(am_r, bm_r);
            vrr += even(am_rr, bm_rr);
            vt += even(am_t, bm_t);
            vtt += even(am_tt, bm_tt);
            vrt += even(am_rt, bm_rt);
            vl += dl(am, bm);
            vll += -mf * mf * even(am, bm);
            vrl += dl(am_r, bm_r);
            vtl += dl(am_t, bm_t);
        }

        // Un-scale and apply the GM/r prefactors.
        let s = HF_SCALE_INV;
        let (gmr, gmr2, gmr3) = (self.gm / r, self.gm / (r * r), self.gm / (r * r * r));
        Partials {
            v: gmr * v * s,
            vr: -gmr2 * vr * s,
            vrr: gmr3 * vrr * s,
            vt: gmr * vt * s,
            vtt: gmr * vtt * s,
            vrt: -gmr2 * vrt * s,
            vl: gmr * vl * s,
            vll: gmr * vll * s,
            vrl: -gmr2 * vrl * s,
            vtl: gmr * vtl * s,
        }
    }

    /// Scaled sectorial `P̄_mm` seeded per the Holmes–Featherstone convention.
    fn sectorial(&self, m: usize, u: f64) -> f64 {
        if m == 0 {
            return HF_SCALE; // P̄_00 = 1
        }
        // P̄_11 = √3·u, then P̄_mm = u·√((2m+1)/(2m))·P̄_{m-1,m-1}.
        let mut p = 3.0f64.sqrt() * u * HF_SCALE; // m = 1
        for k in 2..=m {
            p *= u * ((2 * k + 1) as f64 / (2 * k) as f64).sqrt();
        }
        p
    }

    /// One column-recursion step for `P̄_nm` (n ≥ m+1) from `P̄_{n-1,m}` and `P̄_{n-2,m}`.
    fn column_next(&self, n: usize, m: usize, t: f64, p_nm1: f64, p_nm2: f64) -> f64 {
        if n == m + 1 {
            // First off-sectorial: P̄_{m+1,m} = t·√(2m+3)·P̄_mm.
            t * ((2 * m + 3) as f64).sqrt() * p_nm1
        } else {
            let (nf, mf) = (n as f64, m as f64);
            let a = (((2.0 * nf - 1.0) * (2.0 * nf + 1.0)) / ((nf - mf) * (nf + mf))).sqrt();
            let b = (((2.0 * nf + 1.0) * (nf + mf - 1.0) * (nf - mf - 1.0))
                / ((nf - mf) * (nf + mf) * (2.0 * nf - 3.0)))
                .sqrt();
            a * t * p_nm1 - b * p_nm2
        }
    }
}

impl GravityModel for SphericalHarmonicModel {
    #[inline]
    fn gm(&self) -> f64 {
        self.gm
    }
    #[inline]
    fn reference_radius(&self) -> f64 {
        self.a
    }

    fn potential(&self, position: Ecef) -> f64 {
        let (r, theta, lambda) = spherical(position);
        self.synthesize(r, theta, lambda).v
    }

    fn acceleration(&self, position: Ecef) -> Vec3 {
        let (r, theta, lambda) = spherical(position);
        let p = self.synthesize(r, theta, lambda);
        let (st, ct) = theta.sin_cos();
        let (sl, cl) = lambda.sin_cos();
        // Spherical unit vectors in ECEF.
        let rhat = Vec3::new(st * cl, st * sl, ct);
        let thetahat = Vec3::new(ct * cl, ct * sl, -st);
        let lambdahat = Vec3::new(-sl, cl, 0.0);
        // ∇V = ∂V/∂r r̂ + (1/r)∂V/∂θ θ̂ + (1/(r sinθ))∂V/∂λ λ̂.
        rhat * p.vr + thetahat * (p.vt / r) + lambdahat * (p.vl / (r * st))
    }
}

/// `dP̄_nm/dθ = [n·t·P̄_nm − f·P̄_{n-1,m}]/u`, `f = √(((2n+1)/(2n-1))(n-m)(n+m))`, for `m ≥ 1`
/// (`P̄_nm` carries a factor `u^m`, so the division loses nothing).
fn dp_dtheta(n: usize, m: usize, t: f64, u: f64, p_n: f64, p_nm1: f64) -> f64 {
    let f = (((2 * n + 1) as f64 / (2 * n - 1) as f64) * ((n - m) * (n + m)) as f64).sqrt();
    (n as f64 * t * p_n - f * p_nm1) / u
}

fn parse_f64(s: &str) -> Option<f64> {
    // ICGEM files sometimes use Fortran 'D' exponents.
    s.replace(['D', 'd'], "e").parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{GRS80_A, GRS80_GM, GRS80_OMEGA};
    use geo4d::{GRS80, WGS84};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// A base normal field (degrees 0, 2, 4) plus an optional single-coefficient perturbation,
    /// all sharing GM/a/ω so their disturbing potential is exactly the perturbation.
    fn field(perturb: &[(usize, usize, f64, f64)]) -> SphericalHarmonicModel {
        let nmax = 4;
        let size = (nmax + 1) * (nmax + 2) / 2;
        let mut c = vec![0.0; size];
        let mut s = vec![0.0; size];
        c[idx(0, 0)] = 1.0;
        c[idx(2, 0)] = -0.484_166_774_985e-3;
        c[idx(4, 0)] = 0.790_303_733_511e-6;
        for &(n, m, dc, ds) in perturb {
            c[idx(n, m)] += dc;
            s[idx(n, m)] += ds;
        }
        SphericalHarmonicModel::new(GRS80_GM, GRS80_A, GRS80_OMEGA, nmax, c, s)
    }

    #[test]
    fn monopole_is_gm_over_r() {
        // Only C̄_00 = 1 → V = GM/r and |∇V| = GM/r².
        let gm = 3.986_004_418e14;
        let a = 6_378_137.0;
        let m = SphericalHarmonicModel::new(gm, a, 0.0, 0, vec![1.0], vec![0.0]);
        let p = Ecef::new(7_000_000.0, 0.0, 0.0);
        assert!(close(m.potential(p), gm / 7_000_000.0, 1e-3));
        assert!(close(m.acceleration(p).norm(), gm / 7_000_000.0_f64.powi(2), 1e-9));
        // The monopole acceleration points inward (toward the geocentre).
        assert!(m.acceleration(p).x < 0.0);
    }

    #[test]
    fn grs80_normal_field_reproduces_somigliana() {
        // The engine synthesized on the GRS80 even zonals must match closed-form normal gravity,
        // across latitudes and including the pole — a full end-to-end validation with no data.
        let model = SphericalHarmonicModel::grs80_normal_field();
        for lat in [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.9] {
            let on_ellipsoid = GRS80.to_ecef(Geodetic::new(lat, 20.0, 0.0));
            let synth = model.gravity_magnitude(on_ellipsoid);
            let closed = NormalGravity::GRS80.at_latitude(lat);
            assert!(close(synth, closed, 1e-5), "lat {lat}: {synth} vs {closed}");
        }
    }

    #[test]
    fn normal_field_has_zero_geoid_and_disturbance_against_itself() {
        // The disturbing potential of the normal field relative to itself is identically zero → with
        // w0 = U0 the geoid undulation (T/γ − (w0−U0)/γ) vanishes too.
        let model = SphericalHarmonicModel::grs80_normal_field();
        let normal = SphericalHarmonicModel::grs80_normal_field();
        let ng = NormalGravity::GRS80;
        for lat in [0.0, 45.0, 80.0] {
            let g = Geodetic::new(lat, 10.0, 0.0);
            assert!(close(model.geoid_undulation(g, &normal, &ng, ng.u0), 0.0, 1e-9));
            // A geoid potential 9.8 m²/s² below U0 raises the undulation by ~1 m (the zero-degree term).
            assert!(close(model.geoid_undulation(g, &normal, &ng, ng.u0 - 9.8), 1.0, 0.02));
        }
    }

    #[test]
    fn gfc_parser_round_trips_a_small_model() {
        let text = "\
begin_of_head
earth_gravity_constant 3.986004415e14
radius 6378136.3
max_degree 2
end_of_head
gfc 0 0 1.0 0.0
gfc 2 0 -4.84165e-4 0.0
gfc 2 2 2.439e-6 -1.400e-6
";
        let m = SphericalHarmonicModel::from_gfc(text).unwrap();
        assert_eq!(m.max_degree(), 2);
        assert!(close(m.gm(), 3.986004415e14, 1e3));
        assert!(close(m.reference_radius(), 6_378_136.3, 1e-3));
        // A missing header is an error.
        assert!(SphericalHarmonicModel::from_gfc("gfc 0 0 1 0").is_err());
    }

    #[test]
    fn gfc_parser_refuses_what_it_cannot_represent() {
        let head = "earth_gravity_constant 3.986004415e14\nradius 6378136.3\nmax_degree 2\nend_of_head\n";
        for body in [
            "gfct 2 0 -4.84e-4 0.0 20100101\ntrnd 2 0 -2e-11 0.0\n",
            "acos 2 0 8e-11 0.0 1.0\n",
        ] {
            let err = SphericalHarmonicModel::from_gfc(&format!("{head}{body}")).unwrap_err();
            assert!(matches!(err, GfcError::Unsupported(_)), "{err}");
        }
        let unnormalized = format!("norm unnormalized\n{head}gfc 2 0 -1.08e-3 0.0\n");
        assert!(matches!(
            SphericalHarmonicModel::from_gfc(&unnormalized),
            Err(GfcError::Unsupported(_))
        ));
    }

    #[test]
    fn acceleration_is_exact_on_the_rotation_axis() {
        // Near the poles the C̄21/S̄21 potential is ±GM·a²·√15·(C̄21·x + S̄21·y)/r⁴ (P̄21 ∝ sinθ·cosθ),
        // so on the axis its horizontal gravity is exactly ±GM·a²·√15·(C̄21, S̄21)/r⁴. This used to
        // be NaN.
        let dc = 1e-6;
        let model = field(&[(2, 1, dc, 2.0 * dc)]);
        let r: f64 = 6_400_000.0;
        let k = GRS80_GM * GRS80_A * GRS80_A * 15f64.sqrt() / r.powi(4);
        for z in [r, -r] {
            let k = k * z.signum();
            let g = model.acceleration(Ecef::new(0.0, 0.0, z));
            assert!(close(g.x, k * dc, 1e-9 * (k * dc).abs()), "{z}: {}", g.x);
            assert!(close(g.y, k * 2.0 * dc, 1e-9 * (k * dc).abs()), "{z}: {}", g.y);
            // 1 µm off the axis (where acos(z/r) used to round to 0) agrees, up to the ~1e-12 m/s²
            // the field really changes over 1 µm.
            let near = model.acceleration(Ecef::new(1e-6, 0.0, z));
            assert!(close(near.x, g.x, 1e-10) && close(near.y, g.y, 1e-10) && close(near.z, g.z, 1e-10));
        }
    }

    #[test]
    fn zonal_derivatives_stay_accurate_near_the_pole() {
        // The generic (n·t·P̄ − f·P̄_{n−1})/sinθ form lost ~2 digits per decade of θ here.
        let model = field(&[]);
        let (r, theta): (f64, f64) = (7.0e6, 1e-7);
        let p = model.synthesize(r, theta, 0.3);
        let (t, u) = (theta.cos(), theta.sin());
        let q = GRS80_A / r;
        let dp20 = -3.0 * 5f64.sqrt() * t * u;
        let dp40 = -3.0 / 8.0 * (140.0 * t.powi(3) - 60.0 * t) * u;
        let expect =
            GRS80_GM / r * (q.powi(2) * -0.484_166_774_985e-3 * dp20 + q.powi(4) * 0.790_303_733_511e-6 * dp40);
        assert!(close(p.vt, expect, 1e-9 * expect.abs()), "{} vs {expect}", p.vt);
    }

    #[test]
    fn gradient_tensor_is_trace_free_and_symmetric() {
        // Laplace's equation (trace = 0) outside the masses is an exact identity for every harmonic,
        // so a vanishing trace validates the synthesized *second* derivatives.
        let model = SphericalHarmonicModel::grs80_normal_field();
        let p = GRS80.to_ecef(Geodetic::new(35.0, 47.0, 250_000.0)); // above the surface
        let g = model.gradient_tensor_enu(p);
        let trace = g[0][0] + g[1][1] + g[2][2];
        assert!(trace.abs() < 1e-15, "Laplace trace {trace}");
        // Symmetric, and the components are the expected ~1e-6 s⁻² scale.
        assert_eq!((g[0][1], g[0][2], g[1][2]), (g[1][0], g[2][0], g[2][1]));
        assert!(g[2][2].abs() > 1e-7 && g[2][2].abs() < 1e-5);
    }

    #[test]
    fn anomaly_and_deflections_from_a_known_perturbation() {
        let normal = field(&[]);
        let pos = Geodetic::new(30.0, 40.0, 0.0);

        // The disturbing quantities of the field against itself are exactly zero.
        assert!(normal.gravity_anomaly(pos, &normal).abs() < 1e-12);
        let (xi0, eta0) = normal.deflections(pos, &normal, &NormalGravity::GRS80);
        assert!(xi0.abs() < 1e-18 && eta0.abs() < 1e-18);

        // A degree-3 zonal perturbation: Δg matches the analytic (n−1)·(GM/r²)(a/r)³·ΔC·P̄₃₀, and the
        // east deflection η is zero (no longitude dependence), while ξ is not.
        let dc = 1e-6;
        let zonal = field(&[(3, 0, dc, 0.0)]);
        let (r, theta, _) = spherical(WGS84.to_ecef(pos));
        let tc = theta.cos();
        let p30 = 7.0f64.sqrt() * 0.5 * (5.0 * tc.powi(3) - 3.0 * tc); // fully-normalized P̄₃₀
        let expected = 2.0 * (GRS80_GM / (r * r)) * (GRS80_A / r).powi(3) * dc * p30;
        assert!(close(zonal.gravity_anomaly(pos, &normal), expected, 1e-11));
        let (xi, eta) = zonal.deflections(pos, &normal, &NormalGravity::GRS80);
        assert!(eta.abs() < 1e-18 && xi.abs() > 1e-12, "xi {xi} eta {eta}");

        // A tesseral (order-1) perturbation produces a non-zero east deflection.
        let tess = field(&[(2, 1, dc, dc)]);
        let (_x, eta_t) = tess.deflections(pos, &normal, &NormalGravity::GRS80);
        assert!(eta_t.abs() > 1e-12, "{eta_t}");
    }
}
