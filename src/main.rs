use sgp4_rust::{Tle, convert_satellite_data, sgp4};
use std::fs::read_to_string;

use colored::*;

fn afficher_positions(tsince_values: &[f64], positions: &[[f64; 3]], elements: &sgp4_rust::OrbitalElements) {

    println!("\nComparaison des positions :");
    println!("┌── TSINCE ──┬───────────── REF POSITION [km] ────────┬───────────── SIMULATED [km] ───────────┬──────────── DELTA [km] ──────────┬─ Total Δ [km] ─┐");
    for (i, &tsince) in tsince_values.iter().enumerate() {
        let expected: [f64; 3] = positions[i];
        let state: sgp4_rust::StateVector = sgp4(tsince, elements);

        let dx: f64 = (expected[0] - state.position[0]).abs();
        let dy: f64 = (expected[1] - state.position[1]).abs();
        let dz: f64 = (expected[2] - state.position[2]).abs();

        let fmt = |delta: f64| {
            if delta > 150.0 && delta < 300.0 {
                format!("{:>7.2}", delta).yellow()
            } else if delta > 300.0 {
                  format!("{:>7.2}", delta).red()
            } else {
                format!("{:>7.2}", delta).green()
            }
        };
        let fmt_total = |delta: f64| {
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



fn afficher_vitesses(tsince_values: &[f64], velocities: &[[f64; 3]], elements: &sgp4_rust::OrbitalElements) {

    println!("\nComparaison des vitesses :");
    println!("┌─── TSINCE ────┬──────── REF VELOCITY [km/s] ──────┬──────────── SIMULATED [km/s] ───────────┬───────── DELTA [km/s] ───────┬─ Total Δ [km] ─┐");

    for (i, &tsince) in tsince_values.iter().enumerate() {
        let reference = velocities[i];
        let state = sgp4(tsince, elements);

        let dvx = (reference[0] - state.velocity[0]).abs();
        let dvy = (reference[1] - state.velocity[1]).abs();
        let dvz = (reference[2] - state.velocity[2]).abs();

        let fmt = |delta: f64| {
            if delta > 0.03 && delta <  0.05 {
                format!("{:>6.4}", delta).yellow()
            } else if delta > 0.05 {
                  format!("{:>6.4}", delta).red()
            } else {
                format!("{:>6.4}", delta).green()
            }
        };


        let fmt_total = |delta: f64| {
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



fn main() {
    let input: String = read_to_string("data/sample.txt").expect("Could not read file");

    let lines: Vec<&str> = input.lines().collect();
    let tle: Tle = Tle {
        line1: lines[0].to_string(),
        line2: lines[1].to_string(),
    };
    let elements = convert_satellite_data(&tle);

    // Extraction des TSINCE
    let mut tsince_values: Vec<f64> = Vec::new();
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut velocities: Vec<[f64; 3]> = Vec::new();

    let mut mode: &'static str = "";
    for line in &lines[3..] {
        let line: &str = line.trim();
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
            .filter_map(|x: &str| x.parse::<f64>().ok())
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
