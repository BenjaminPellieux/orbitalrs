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
const TWOPI: f64 = 2.0 * PI;
const XKE: f64 = 0.0743669161;
const CK2: f64 = 5.413080e-4;
const AE: f64 = 1.0;
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
    let line1 = &tle.line1;
    let line2 = &tle.line2;

    let _epoch = parse_real(line1, 19, 14);
    let _xndt2o = parse_real(line1, 34, 10);
    let _xndd6o = parse_real(line1, 45, 6);
    let _iexp = parse_real(line1, 51, 2);
    let bstar = parse_real(line1, 54, 6) * 1e-5 * 10f64.powf(parse_real(line1, 60, 2));

    let inclination = radians(parse_real(line2, 9, 8));
    let raan = radians(parse_real(line2, 18, 8));
    let eccentricity = parse_real(line2, 27, 7) * 1e-7;
    let arg_perigee = radians(parse_real(line2, 35, 8));
    let mean_anomaly = radians(parse_real(line2, 44, 8));
    let mean_motion = parse_real(line2, 53, 11);

    // convert mean motion to radians/min
    let xno = mean_motion * TWOPI / XMNPDA;
    let a1 = (XKE / xno).powf(TOTHIRD);
    let temp = 1.5 * CK2 * (3.0 * inclination.cos().powi(2) - 1.0) / (1.0 - eccentricity.powi(2)).powf(1.5);
    let del1 = temp / (a1 * a1);
    let ao = a1 * (1.0 - del1 * (0.5 * TOTHIRD + del1 * (1.0 + 134.0 / 81.0 * del1)));
    let delo = temp / (ao * ao);
    let xnodp = xno / (1.0 + delo);

    let deep_space = (TWOPI / xnodp) >= 225.0;

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
    let mut e = mean_anomaly; // Première estimation
    let mut delta:f64 = 1.0;

    while delta.abs() > tol {
        let f = e - eccentricity * e.sin() - mean_anomaly;
        let f_prime = 1.0 - eccentricity * e.cos();
        delta = f / f_prime;
        e -= delta;
    }

    e
}


pub fn sgp4(tsince: f64, elements: &OrbitalElements) -> StateVector {
    // Constantes orbitales
    let a = (XKE / elements.mean_motion).powf(2.0 / 3.0); // demi-grand axe en rayons terrestres
    let e = elements.eccentricity;
    let i = elements.inclination;
    let omega = elements.arg_perigee;
    let raan = elements.raan;

    // Anomalie moyenne à t = tsince
    let m = elements.mean_anomaly + elements.mean_motion * tsince;

    // Résolution de l'équation de Kepler
    let e_anomaly = solve_kepler(m % (2.0 * PI), e, 1e-8);

    // Anomalie vraie
    let v = 2.0 * ((1.0 + e).sqrt() * (e_anomaly / 2.0).sin()).atan2((1.0 - e).sqrt() * (e_anomaly / 2.0).cos());

    // Distance au foyer
    let r = a * (1.0 - e * e_anomaly.cos());

    // Coordonnées dans le plan orbital
    let x_orb = r * v.cos();
    let y_orb = r * v.sin();

    // Passage au référentiel inertiel géocentrique (ECI)
    let cos_omega = omega.cos();
    let sin_omega = omega.sin();
    let cos_raan = raan.cos();
    let sin_raan = raan.sin();
    let cos_i = i.cos();
    let sin_i = i.sin();

    let x = x_orb * (cos_raan * cos_omega - sin_raan * sin_omega * cos_i)
          - y_orb * (cos_raan * sin_omega + sin_raan * cos_omega * cos_i);
    let y = x_orb * (sin_raan * cos_omega + cos_raan * sin_omega * cos_i)
          - y_orb * (sin_raan * sin_omega - cos_raan * cos_omega * cos_i);
    let z = x_orb * sin_omega * sin_i + y_orb * cos_omega * sin_i;

    StateVector {
        position: [x * XKMPER, y * XKMPER, z * XKMPER],         // km
        velocity: [0.0, 0.0, 0.0], // à compléter (dérivée orbitale complète si besoin)
    }
}


pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
}
