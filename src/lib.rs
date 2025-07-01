//! A library for performing orbital mechanics calculations using the SGP4 model.
//! This library provides functionality to convert Two-Line Element (TLE) data into orbital elements
//! and compute the state vector (position and velocity) of a satellite at a given time.

use std::f64::consts::PI;

/// Represents a Two-Line Element set (TLE) for a satellite.
/// TLEs are used to describe the orbits of Earth-orbiting objects.
pub struct Tle {
    /// First line of the TLE data.
    pub line1: String,
    /// Second line of the TLE data.
    pub line2: String,
}

/// Represents the state vector of a satellite, including its position and velocity.
pub struct StateVector {
    /// Position of the satellite in kilometers (X, Y, Z).
    pub position: [f64; 3],
    /// Velocity of the satellite in kilometers per second (XDOT, YDOT, ZDOT).
    pub velocity: [f64; 3],
}

/// Represents the orbital elements of a satellite.
pub struct OrbitalElements {
    /// Inclination of the orbit in radians.
    pub inclination: f64,
    /// Right Ascension of the Ascending Node (RAAN) in radians.
    pub raan: f64,
    /// Eccentricity of the orbit, unitless.
    pub eccentricity: f64,
    /// Argument of perigee in radians.
    pub arg_perigee: f64,
    /// Mean anomaly in radians.
    pub mean_anomaly: f64,
    /// Mean motion in radians per minute.
    pub mean_motion: f64,
    /// Bstar drag term in 1/earth radii.
    pub bstar: f64,
    /// Flag indicating if the orbit is in deep space.
    pub deep_space: bool,
}

/// Constant representing 2 * PI.
const TWOPI: f64 = 2.0 * std::f64::consts::PI;
/// Earth's gravitational constant.
const XKE: f64 = 0.0743669161;
/// Second zonal harmonic coefficient for Earth.
const CK2: f64 = 5.413080e-4;
/// Minutes per day.
const XMNPDA: f64 = 1440.0;
/// Two-thirds constant.
const TOTHIRD: f64 = 2.0 / 3.0;
/// Earth's radius in kilometers.
const XKMPER: f64 = 6378.135;

/// Converts degrees to radians.
///
/// # Arguments
/// * `deg` - Angle in degrees.
///
/// # Returns
/// * Angle in radians.
fn radians(deg: f64) -> f64 {
    deg * PI / 180.0
}

/// Parses a substring from a TLE line and converts it to a real number.
///
/// # Arguments
/// * `line` - The TLE line to parse.
/// * `start` - The starting index of the substring (1-based).
/// * `len` - The length of the substring.
///
/// # Returns
/// * The parsed real number.
fn parse_real(line: &str, start: usize, len: usize) -> f64 {
    line.get(start - 1..start - 1 + len)
        .unwrap_or("0")
        .trim()
        .replace("-", "e-")
        .replace("+", "e+")
        .parse::<f64>()
        .unwrap_or(0.0)
}

/// Converts satellite TLE data into orbital elements.
///
/// # Arguments
/// * `tle` - The Two-Line Element set for the satellite.
///
/// # Returns
/// * Orbital elements derived from the TLE data.
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

    // Convert mean motion to radians per minute
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

/// Solves Kepler's equation: M = E - e * sin(E) using the Newton-Raphson method.
///
/// # Arguments
/// * `mean_anomaly` - Mean anomaly in radians.
/// * `eccentricity` - Eccentricity of the orbit.
/// * `tol` - Tolerance for the solution.
///
/// # Returns
/// * Eccentric anomaly in radians.
fn solve_kepler(mean_anomaly: f64, eccentricity: f64, tol: f64) -> f64 {
    let mut e: f64 = mean_anomaly; // Initial estimate
    let mut delta: f64 = 1.0;

    while delta.abs() > tol {
        let f: f64 = e - eccentricity * e.sin() - mean_anomaly;
        let f_prime: f64 = 1.0 - eccentricity * e.cos();
        delta = f / f_prime;
        e -= delta;
    }

    e
}

/// Computes the state vector (position and velocity) of a satellite using the SGP4 model.
///
/// # Arguments
/// * `tsince` - Time since epoch in minutes.
/// * `elements` - Orbital elements of the satellite.
///
/// # Returns
/// * State vector containing the position and velocity of the satellite.
pub fn sgp4(tsince: f64, elements: &OrbitalElements) -> StateVector {
    let a: f64 = (XKE / elements.mean_motion).powf(2.0 / 3.0); // Semi-major axis (earth radii)
    let e: f64 = elements.eccentricity;
    let i: f64 = elements.inclination;
    let omega: f64 = elements.arg_perigee;
    let raan: f64 = elements.raan;

    // Mean anomaly + Kepler's equation solution
    let m: f64 = (elements.mean_anomaly + elements.mean_motion * tsince) % TWOPI;
    let e_anomaly = solve_kepler(m, e, 1e-8);

    // True anomaly
    let v: f64 = 2.0 * ((1.0 + e).sqrt() * (e_anomaly / 2.0).sin())
        .atan2((1.0 - e).sqrt() * (e_anomaly / 2.0).cos());

    // Distance (earth radii)
    let r: f64 = a * (1.0 - e * e_anomaly.cos());

    // Coordinates in the orbital plane
    let x_orb: f64 = r * v.cos();
    let y_orb: f64 = r * v.sin();

    // Velocity in the orbital plane
    let p: f64 = a * (1.0 - e * e); // Semi-latus rectum
    let r_dot: f64 = XKE * a.sqrt() * e * e_anomaly.sin() / r;
    let r_fi_dot: f64 = XKE * (p).sqrt() / (r * r);

    let vx_orb: f64 = r_dot * v.cos() - r * r_fi_dot * v.sin();
    let vy_orb: f64 = r_dot * v.sin() + r * r_fi_dot * v.cos();

    // Pre-calculations for inertial transformation
    let cos_omega: f64 = omega.cos();
    let sin_omega: f64 = omega.sin();
    let cos_raan: f64 = raan.cos();
    let sin_raan: f64 = raan.sin();
    let cos_i: f64 = i.cos();
    let sin_i: f64 = i.sin();

    // Inertial position ECI (earth radii)
    let x: f64 = x_orb * (cos_raan * cos_omega - sin_raan * sin_omega * cos_i)
        - y_orb * (cos_raan * sin_omega + sin_raan * cos_omega * cos_i);
    let y: f64 = x_orb * (sin_raan * cos_omega + cos_raan * sin_omega * cos_i)
        - y_orb * (sin_raan * sin_omega - cos_raan * cos_omega * cos_i);
    let z: f64 = x_orb * sin_omega * sin_i + y_orb * cos_omega * sin_i;

    // Inertial velocity ECI (earth radii per minute)
    let vx: f64 = vx_orb * (cos_raan * cos_omega - sin_raan * sin_omega * cos_i)
        - vy_orb * (cos_raan * sin_omega + sin_raan * cos_omega * cos_i);
    let vy: f64 = vx_orb * (sin_raan * cos_omega + cos_raan * sin_omega * cos_i)
        - vy_orb * (sin_raan * sin_omega - cos_raan * cos_omega * cos_i);
    let vz: f64 = vx_orb * sin_omega * sin_i + vy_orb * cos_omega * sin_i;

    StateVector {
        position: [x * XKMPER, y * XKMPER, z * XKMPER], // km
        velocity: [vx * XKMPER / 60.0, vy * XKMPER / 60.0, vz * XKMPER / 60.0], // km/s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests the conversion of TLE data to orbital elements.
    #[test]
    fn test_convert_satellite_data() {
        let tle = Tle {
            line1: "1 25544U 98067A   21135.57634567  .00002418  00000-0  50843-4 0  9993".to_string(),
            line2: "2 25544  51.6443 126.6639 0006738  34.7758 325.3542 15.48913328283873".to_string(),
        };

        let elements = convert_satellite_data(&tle);

        assert!(elements.inclination > 0.0);
        assert!(elements.raan > 0.0);
        assert!(elements.eccentricity >= 0.0 && elements.eccentricity < 1.0);
        assert!(elements.arg_perigee > 0.0);
        assert!(elements.mean_anomaly > 0.0);
        assert!(elements.mean_motion > 0.0);
    }

    /// Tests the solution of Kepler's equation.
    #[test]
    fn test_solve_kepler() {
        let mean_anomaly: f64 = 1.0;
        let eccentricity: f64 = 0.1;
        let tol: f64 = 1e-8;
        let e: f64 = solve_kepler(mean_anomaly, eccentricity, tol);
        let expected_e: f64 = 1.0885977523978936;
        assert!((e - expected_e).abs() < tol, "Kepler's equation solution is not within the expected tolerance");
    }

    /// Tests the SGP4 model computation of the state vector.
    #[test]
    fn test_sgp4() {
        let tle = Tle {
            line1: "1 25544U 98067A   21135.57634567  .00002418  00000-0  50843-4 0  9993".to_string(),
            line2: "2 25544  51.6443 126.6639 0006738  34.7758 325.3542 15.48913328283873".to_string(),
        };

        let elements = convert_satellite_data(&tle);
        let tsince = 0.0; // minutes since epoch
        let state = sgp4(tsince, &elements);

        assert!(state.position.iter().all(|&x| x.abs() < 10000.0)); // Check if position values are reasonable
        assert!(state.velocity.iter().all(|&x| x.abs() < 10.0)); // Check if velocity values are reasonable
    }
}
