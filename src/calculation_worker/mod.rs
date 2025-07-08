//
// calculation_worker/main.rs: Run calculation code
//
// Copyright (c) 2025 Pablo Álvarez Martín <pablo.alvmar12@gmail.com>
// Copyright (c) 2025 Gonzalo J. Carracedo <BatchDrake@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//

use ordered_float::NotNan;
use crate::simulation::Simulation;
use crate::simulation::SimulationParams;
use crate::spectrum::*;
use crate::curve::*;
use crate::sky_model::*;
use crate::instrument_model::*;
use crate::helpers::*;
use crate::detector::*;
use crate::curve::CurveAxis::*;
use crate::data_file_manager::*;
use crate::spectrum::SpectrumAxisOperations;
use std::collections::BTreeMap;

use std::time::UNIX_EPOCH;
use std::time::SystemTime;
use rand::prelude::*;
use rand_distr::{Poisson, Normal, Distribution};

#[derive(Clone)]
pub struct SNRcurve {
    initialized: bool,
    pub wl_to_pixel: Curve,
    pub wavelength: Vec<f64>,
    pub signal: Curve, 
    pub noise: Curve,
    pub counts: Curve
}

impl SNRcurve {
    pub fn new() -> SNRcurve{
        let mut snr_curve = SNRcurve {
            initialized: false,
            wl_to_pixel: Curve::default(),
            wavelength: Vec::new(),
            signal: Curve::default(),
            noise: Curve::default(),
            counts: Curve::default(),
        };

        snr_curve.wavelength.resize(DETECTOR_PIXELS as usize, 0.0); // can i put it here?

        return snr_curve;
    }

    pub fn snr_curve(&mut self) {
        self.wavelength.resize(DETECTOR_PIXELS as usize, 0.0);
    }
}

pub struct CalculationProduct {
    red_arm: SNRcurve,
    blue_arm: SNRcurve,
}
pub struct CalculationWorker {
    pub simulation: Simulation,
    sim_params: SimulationParams,
    new_spectrum: bool,
    pub input_spectrum: Spectrum,
}

impl CalculationWorker {
    pub fn init() -> CalculationWorker {
        let mut calculation_worker = CalculationWorker {
            simulation: Simulation::new(),
            new_spectrum: false,
            sim_params: SimulationParams::new(),
            input_spectrum: Spectrum::default(),
        };

        return calculation_worker;
    }

    fn simulate_arm(&mut self, arm: InstrumentArm, curve: &mut SNRcurve, input_simulation: &mut Simulation) { // WORKS WITH AN ALREADY DONE SIMULATION AS INPUT, AND ONLY MODIFIES IT
        let ron: f64;
        let inv_gain: f64;

        // generator seed -> now time
        let now = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
        let seed = now.as_secs() * 1_000_000 + now.subsec_micros() as u64;
        let mut rng = StdRng::seed_from_u64(seed);

        inv_gain = 1. / input_simulation.gain();
        ron = input_simulation.read_out_noise(); // In counts
        let binding = input_simulation.clone();
        let slice = binding.get_params().get_slice();

        curve.wl_to_pixel = input_simulation.wl_to_pixel_curve(arm);
        let electrons = input_simulation.clone().get_det().electrons_crv().get_curve();
        //println!("electrons: {:?}", electrons);
        let signal = input_simulation.clone().get_det().signal_crv().get_curve();
        let noise = input_simulation.clone().get_det().noise_crv().get_curve();

        for i in 0..DETECTOR_PIXELS as i32 {
            let wl = input_simulation.get_tarsis_model_mut().px_to_wavelength_val(arm, *slice, NotNan::new(i as f64).expect("x should not be NaN")).unwrap();

            curve.wavelength[i as usize] = input_simulation.px_to_wavelength(arm, NotNan::new(i as f64).expect("x should not be NaN"));
            
            // Poisson-distributed shot noise
            let poisson = Poisson::new(electrons.get_point(NotNan::new(wl as f64).expect("x should not be NaN"))).unwrap();
            let shot_electrons = poisson.sample(&mut rng);
            // Normally-distributed readout noise
            let normal = Normal::new(0.0, 1.0).unwrap();
            let readout_noise = ron * normal.sample(&mut rng);

            // Final simulated count
            let counts = (inv_gain * shot_electrons as f64 + readout_noise).round() as i64;
            curve.counts.get_map_mut().insert(NotNan::new(wl as f64).expect("x should not be NaN"), counts as f64);
        }

        curve.signal = signal;
        curve.noise = noise;

        curve.initialized = true;

    }

    fn single_shot(&mut self, curve: &mut SNRcurve, input_simulation: &mut Simulation, arm: InstrumentArm) -> SNRcurve {

        self.simulate_arm(arm, curve, input_simulation);// &mut curve);

        return curve.clone();

    }

    fn all_slices(&mut self, curve: &mut SNRcurve, input_simulation: &mut Simulation, arm: InstrumentArm) -> SNRcurve {

        for i in 0..TARSIS_SLICES {
            *input_simulation.get_params_mut().get_slice_mut() = i;

            self.simulate_arm(arm, curve, input_simulation);//, &mut curve);
        }

        return curve.clone();
    }

    pub fn simulate(&mut self, curve: &mut SNRcurve, input_simulation: &mut Simulation, arm: InstrumentArm) -> SNRcurve{

        if *input_simulation.get_params().get_slice() < 0 {
            self.all_slices(curve, input_simulation, arm);
        } else {
            self.single_shot(curve, input_simulation, arm);
        }

        return curve.clone();

    }

}