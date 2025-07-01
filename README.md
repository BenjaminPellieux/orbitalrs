# SGP4 Satellite Orbit Model

A Rust library for performing orbital mechanics calculations using the Simplified General Perturbations Satellite Orbit Model (SGP4). This library provides functionality to convert Two-Line Element (TLE) data into orbital elements and compute the state vector (position and velocity) of a satellite at a given time.

## Features

- Parse Two-Line Element (TLE) data.
- Convert TLE data into orbital elements.
- Compute satellite state vectors (position and velocity) using the SGP4 model.
- Compare computed satellite positions and velocities with reference data.

## Installation
To use this library and application, you need to have Rust installed on your machine. If you don't have Rust installed, follow the instructions [here](https://www.rust-lang.org/tools/install).

1. Clone the repository:

```bash
git clone <repository-url>
cd sgp4-rust
```    

2. Build the project:
```bash
cargo build --release
```

3. Run the application:

```bash
cargo run --release
```

## Usage

Prepare your TLE data file and place it in the data directory of the project. The file should have the following format:

```
1 LINE1_DATA
2 LINE2_DATA
SDP4
tsince1 posX1 posY1 posZ1
tsince2 posX2 posY2 posZ2
...
XDOT
velX1 velY1 velZ1
velX2 velY2 velZ2
...

```

    
Modify the path in the main.rs file if your data file is named differently than sample.txt.

Run the application to see the comparison between reference and computed satellite positions and velocities.

```bash
cargo run --release
```

## Comparing Positions and Velocities

To compare computed satellite positions and velocities with reference data, you can use the provided main.rs example:

```rs
use sgp4_rust::{Tle, convert_satellite_data, sgp4};
use std::fs::read_to_string;

fn main() {
    let input = read_to_string("data/sample.txt").expect("Could not read file");
    let lines: Vec<&str> = input.lines().collect();

    let tle = Tle {
        line1: lines[0].to_string(),
        line2: lines[1].to_string(),
    };

    let elements = convert_satellite_data(&tle);

    // Extract TSINCE, positions, and velocities from the file
    // ...

    afficher_positions(&tsince_values, &positions, &elements);
    afficher_vitesses(&tsince_values, &velocities, &elements);
}
```


## Documentation

The code is documented with comments and docstrings to help understand its functionality. You can generate the documentation locally using:

```bash
cargo doc --open

```
This will open the documentation in your default web browser.

## Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes. Ensure that your code adheres to the existing style and includes appropriate tests.
License

This project is licensed under the MIT License. See the LICENSE file for details.
Acknowledgements

Thanks to the Rust community for their extensive documentation and helpful resources.
Inspiration and reference data taken from various SGP4 implementations and orbital mechanics resources.
