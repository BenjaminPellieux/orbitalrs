//! Main application for comparing satellite positions and velocities using SGP4 model.
//! This application reads Two-Line Element (TLE) data and compares the computed satellite state vectors
//! with reference data.

use sgp4_rust::{Tle, convert_satellite_data, sgp4, OrbitalElements};
use std::fs::read_to_string;
use colored::*;

/// Displays a comparison of satellite positions.
///
/// This function prints a table comparing the reference positions with the simulated positions
/// at various times since epoch (TSINCE). The differences are color-coded based on their magnitude.
///
/// # Arguments
/// * `tsince_values` - A slice of time values since epoch in minutes.
/// * `positions` - A slice of reference positions corresponding to the time values.
/// * `elements` - The orbital elements of the satellite.
fn afficher_positions(tsince_values: &[f64], positions: &[[f64; 3]], elements: &OrbitalElements) {
    println!("\nComparaison des positions :");
    println!("┌── TSINCE ──┬───────────── REF POSITION [km] ────────┬───────────── SIMULATED [km] ───────────┬──────────── DELTA [km] ──────────┬─ Total Δ [km] ─┐");

    for (i, &tsince) in tsince_values.iter().enumerate() {
        let expected = positions[i];
        let state = sgp4(tsince, elements);

        let dx = (expected[0] - state.position[0]).abs();
        let dy = (expected[1] - state.position[1]).abs();
        let dz = (expected[2] - state.position[2]).abs();

        let fmt = |delta: f64| -> ColoredString {
            if delta > 150.0 && delta < 300.0 {
                format!("{:>7.2}", delta).yellow()
            } else if delta > 300.0 {
                format!("{:>7.2}", delta).red()
            } else {
                format!("{:>7.2}", delta).green()
            }
        };

        let fmt_total = |delta: f64| -> ColoredString {
            if delta > 200.0 && delta < 700.0 {
                format!("{:>7.2}", delta).yellow()
            } else if delta > 700.0 {
                format!("{:>7.2}", delta).red()
            } else {
                format!("{:>7.2}", delta).green()
            }
        };

        println!(
            "│ {:>10.3} │ X={:>10.2} Y={:>10.2} Z={:>10.2} │ X={:>10.2} Y={:>10.2} Z={:>10.2} │ ΔX={} ΔY={} ΔZ={} │   ΔR={}   │",
            tsince,
            expected[0], expected[1], expected[2],
            state.position[0], state.position[1], state.position[2],
            fmt(dx), fmt(dy), fmt(dz),
            fmt_total(dx + dy + dz),
        );
    }

    println!("└────────────┴────────────────────────────────────────┴────────────────────────────────────────┴──────────────────────────────────┴────────────────┘");
}

/// Displays a comparison of satellite velocities.
///
/// This function prints a table comparing the reference velocities with the simulated velocities
/// at various times since epoch (TSINCE). The differences are color-coded based on their magnitude.
///
/// # Arguments
/// * `tsince_values` - A slice of time values since epoch in minutes.
/// * `velocities` - A slice of reference velocities corresponding to the time values.
/// * `elements` - The orbital elements of the satellite.
fn afficher_vitesses(tsince_values: &[f64], velocities: &[[f64; 3]], elements: &OrbitalElements) {
    println!("\nComparaison des vitesses :");
    println!("┌─── TSINCE ────┬──────── REF VELOCITY [km/s] ──────┬──────────── SIMULATED [km/s] ───────────┬───────── DELTA [km/s] ───────┬─ Total Δ [km] ─┐");

    for (i, &tsince) in tsince_values.iter().enumerate() {
        let reference = velocities[i];
        let state = sgp4(tsince, elements);

        let dvx = (reference[0] - state.velocity[0]).abs();
        let dvy = (reference[1] - state.velocity[1]).abs();
        let dvz = (reference[2] - state.velocity[2]).abs();

        let fmt = |delta: f64| -> ColoredString {
            if delta > 0.03 && delta < 0.05 {
                format!("{:>6.4}", delta).yellow()
            } else if delta > 0.05 {
                format!("{:>6.4}", delta).red()
            } else {
                format!("{:>6.4}", delta).green()
            }
        };

        let fmt_total = |delta: f64| -> ColoredString {
            if delta > 0.050 && delta < 0.100 {
                format!("{:>6.4}", delta).yellow()
            } else if delta > 0.100 {
                format!("{:>6.4}", delta).red()
            } else {
                format!("{:>6.4}", delta).green()
            }
        };

        println!(
            "│ t = {:>5} min │ Ref = [{:>7.4}, {:>7.4}, {:>7.4}] │ Simulated = [{:>7.4}, {:>7.4}, {:>7.4}] │ Δ = [{}, {}, {}] │    ΔR={}   │",
            tsince,
            reference[0], reference[1], reference[2],
            state.velocity[0], state.velocity[1], state.velocity[2],
            fmt(dvx), fmt(dvy), fmt(dvz),
            fmt_total(dvx + dvy + dvz),
        );
    }

    println!("└───────────────┴───────────────────────────────────┴─────────────────────────────────────────┴──────────────────────────────┴────────────────┘");
}

/// Main function to read TLE data, compute satellite state vectors, and display comparisons.
///
/// This function reads the TLE data from a file, converts it to orbital elements, and then
/// extracts the reference positions and velocities. It then displays the comparisons of these
/// positions and velocities with the computed values.
fn main() {
    let input = read_to_string("data/sample.txt").expect("Could not read file");
    let lines: Vec<&str> = input.lines().collect();

    let tle = Tle {
        line1: lines[0].to_string(),
        line2: lines[1].to_string(),
    };

    let elements = convert_satellite_data(&tle);

    // Extract TSINCE, positions, and velocities
    let mut tsince_values: Vec<f64> = Vec::new();
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut velocities: Vec<[f64; 3]> = Vec::new();

    let mut mode = "";
    for line in &lines[3..] {
        let line = line.trim();
        if line.is_empty() {
            continue;
        } else if line.starts_with("SDP4") {
            mode = "positions";
            continue;
        } else if line.starts_with("XDOT") {
            mode = "velocities";
            continue;
        }

        let parts: Vec<f64> = line
            .split_whitespace()
            .filter_map(|x| x.parse::<f64>().ok())
            .collect();

        match mode {
            "positions" if parts.len() == 4 => {
                tsince_values.push(parts[0]);
                positions.push([parts[1], parts[2], parts[3]]);
            }
            "velocities" if parts.len() == 3 => {
                velocities.push([parts[0], parts[1], parts[2]]);
            }
            _ => (),
        }
    }

    afficher_positions(&tsince_values, &positions, &elements);
    afficher_vitesses(&tsince_values, &velocities, &elements);
}
