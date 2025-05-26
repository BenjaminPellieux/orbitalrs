use sgp4_rust::{Tle, convert_satellite_data, sgp4};

fn main() {
    // Ex. ISS (ZARYA) TLE
    let tle = Tle {
        line1: "1 25544U 98067A   24145.62488018  .00014250  00000+0  25510-3 0  9993".to_string(),
        line2: "2 25544  51.6404  27.6709 0004698  35.2185  72.1356 15.50013818398921".to_string(),
    };

    let elements = convert_satellite_data(&tle);
    let state = sgp4(0.0, &elements); // tsince = 0 minutes after epoch

    println!("Position (km): X: {:.3}, Y: {:.3}, Z: {:.3}", state.position[0], state.position[1], state.position[2]);
    println!("Velocity (km/s): XDOT: {:.6}, YDOT: {:.6}, ZDOT: {:.6}", state.velocity[0], state.velocity[1], state.velocity[2]);
}
