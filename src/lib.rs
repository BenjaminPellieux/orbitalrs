use std::f64::consts::PI;

pub struct Tle {
    pub line1: String,
    pub line2: String,
}

pub struct StateVector {
    pub position: [f64; 3], // [X, Y, Z] in km
    pub velocity: [f64; 3], // [XDOT, YDOT, ZDOT] in km/s
}

pub struct OrbitalElements {
    pub inclination: f64,   // radians
    pub raan: f64,          // radians
    pub eccentricity: f64,  // unitless
    pub arg_perigee: f64,   // radians
    pub mean_anomaly: f64,  // radians
    pub mean_motion: f64,   // radians/min
    pub bstar: f64,         // 1/earth radii
    pub deep_space: bool,
}

// Constants
const TWOPI: f64 = 2.0 * std::f64::consts::PI;
const XKE: f64 = 0.0743669161;
const CK2: f64 = 5.413080e-4;
const XMNPDA: f64 = 1440.0;
const TOTHIRD: f64 = 2.0 / 3.0;
const XKMPER: f64 = 6378.135;

fn radians(deg: f64) -> f64 {
    deg * PI / 180.0
}

fn parse_real(line: &str, start: usize, len: usize) -> f64 {
    line.get(start - 1..start - 1 + len)
        .unwrap_or("0")
        .trim()
        .replace("-", "e-")
        .replace("+", "e+")
        .parse::<f64>()
        .unwrap_or(0.0)
}

pub fn convert_satellite_data(tle: &Tle) -> OrbitalElements {
    let line1: &String = &tle.line1;
    let line2: &String = &tle.line2;

    let _epoch: f64 = parse_real(line1, 19, 14);
    let _xndt2o: f64 = parse_real(line1, 34, 10);
    let _xndd6o: f64 = parse_real(line1, 45, 6);
    let _iexp: f64 = parse_real(line1, 51, 2);
    let bstar: f64 = parse_real(line1, 54, 6) * 1e-5 * 10f64.powf(parse_real(line1, 60, 2));

    let inclination: f64 = radians(parse_real(line2, 9, 8));
    let raan: f64 = radians(parse_real(line2, 18, 8));
    let eccentricity: f64 = parse_real(line2, 27, 7) * 1e-7;
    let arg_perigee: f64 = radians(parse_real(line2, 35, 8));
    let mean_anomaly: f64 = radians(parse_real(line2, 44, 8));
    let mean_motion: f64 = parse_real(line2, 53, 11);

    // convert mean motion to radians/min
    let xno: f64 = mean_motion * TWOPI / XMNPDA;
    let a1: f64 = (XKE / xno).powf(TOTHIRD);
    let temp: f64 = 1.5 * CK2 * (3.0 * inclination.cos().powi(2) - 1.0) / (1.0 - eccentricity.powi(2)).powf(1.5);
    let del1: f64 = temp / (a1 * a1);
    let ao: f64 = a1 * (1.0 - del1 * (0.5 * TOTHIRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo: f64 = temp / (ao * ao);
    let xnodp: f64 = xno / (1.0 + delo);

    let deep_space: bool = (TWOPI / xnodp) >= 225.0;

    OrbitalElements {
        inclination,
        raan,
        eccentricity,
        arg_perigee,
        mean_anomaly,
        mean_motion: xnodp,
        bstar,
        deep_space,
    }
}


/// Résout l’équation de Kepler : M = E - e * sin(E)
/// avec la méthode de Newton-Raphson
fn solve_kepler(mean_anomaly: f64, eccentricity: f64, tol: f64) -> f64 {
    let mut e: f64 = mean_anomaly; // Première estimation
    let mut delta:f64 = 1.0;

    while delta.abs() > tol {
        let f: f64 = e - eccentricity * e.sin() - mean_anomaly;
        let f_prime: f64 = 1.0 - eccentricity * e.cos();
        delta = f / f_prime;
        e -= delta;
    }

    e
}


pub fn sgp4(tsince: f64, elements: &OrbitalElements) -> StateVector {
    let a: f64 = (XKE / elements.mean_motion).powf(2.0 / 3.0); // demi-grand axe (rayons terrestres)
    let e: f64 = elements.eccentricity;
    let i: f64 = elements.inclination;
    let omega: f64 = elements.arg_perigee;
    let raan: f64 = elements.raan;

    // Anomalie moyenne + résolution de Kepler
    let m: f64 = (elements.mean_anomaly + elements.mean_motion * tsince) % TWOPI;
    let e_anomaly = solve_kepler(m, e, 1e-8);

    // Anomalie vraie
    let v: f64 = 2.0 * ((1.0 + e).sqrt() * (e_anomaly / 2.0).sin())
        .atan2((1.0 - e).sqrt() * (e_anomaly / 2.0).cos());

    // Distance (en rayons terrestres)
    let r: f64 = a * (1.0 - e * e_anomaly.cos());

    // Coordonnées dans le plan orbital
    let x_orb: f64 = r * v.cos();
    let y_orb: f64 = r * v.sin();

    // Vitesse dans le plan orbital
    let p: f64 = a * (1.0 - e * e); // paramètre orbital
    let r_dot: f64 = XKE * a.sqrt() * e * e_anomaly.sin() / r;
    let r_fi_dot: f64 = XKE * (p).sqrt() / (r * r);

    let vx_orb: f64 = r_dot * v.cos() - r * r_fi_dot * v.sin();
    let vy_orb: f64 = r_dot * v.sin() + r * r_fi_dot * v.cos();

    // Pré-calculs pour transformation inertielle
    let cos_omega: f64 = omega.cos();
    let sin_omega: f64 = omega.sin();
    let cos_raan: f64 = raan.cos();
    let sin_raan: f64 = raan.sin();
    let cos_i: f64 = i.cos();
    let sin_i: f64 = i.sin();

    // Position inertielle ECI (rayons terrestres)
    let x: f64 = x_orb * (cos_raan * cos_omega - sin_raan * sin_omega * cos_i)
          - y_orb * (cos_raan * sin_omega + sin_raan * cos_omega * cos_i);
    let y: f64 = x_orb * (sin_raan * cos_omega + cos_raan * sin_omega * cos_i)
          - y_orb * (sin_raan * sin_omega - cos_raan * cos_omega * cos_i);
    let z: f64 = x_orb * sin_omega * sin_i + y_orb * cos_omega * sin_i;

    // Vitesse inertielle ECI (rayons terrestres / min)
    let vx: f64 = vx_orb * (cos_raan * cos_omega - sin_raan * sin_omega * cos_i)
           - vy_orb * (cos_raan * sin_omega + sin_raan * cos_omega * cos_i);
    let vy: f64 = vx_orb * (sin_raan * cos_omega + cos_raan * sin_omega * cos_i)
           - vy_orb * (sin_raan * sin_omega - cos_raan * cos_omega * cos_i);
    let vz: f64 = vx_orb * sin_omega * sin_i + vy_orb * cos_omega * sin_i;

    StateVector {
        position: [x * XKMPER, y * XKMPER, z * XKMPER],                 // km
        velocity: [vx * XKMPER / 60.0, vy * XKMPER / 60.0, vz * XKMPER / 60.0], // km/s
    }
}