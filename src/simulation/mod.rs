//
// simulation/mod.cpp: Run instrument simulations
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
use crate::spectrum::*;
use crate::curve::*;
use crate::sky_model::*;
use crate::instrument_model::*;
use crate::helpers::*;
use crate::detector::*;
use crate::curve::CurveAxis::*;
use crate::spectrum::SpectrumAxisOperations;

#[derive(Clone)]
pub struct SimulationParams {
    prog_name: String,
    blue_detector: String,
    red_detector: String, 
    airmass: f64,
    moon: f64, 
    exposure: f64,
    r_ab_mag: f64,
    slice: usize,
}
impl SimulationParams {
    // Getter functions (for unit tests)
    pub fn get_prog_name(&self) -> &String {
        &self.prog_name
    }
    pub fn get_blue_detector(&self) -> &String {
        &self.blue_detector
    }
    pub fn get_red_detector(&self) -> &String {
        &self.red_detector
    }
    pub fn get_airmass(&self) -> &f64 {
        &self.airmass
    }
    pub fn get_moon(&self) -> &f64 {
        &self.moon
    }
    pub fn get_exposure(&self) -> &f64 {
        &self.exposure
    }
    pub fn get_r_ab_mag(&self) -> &f64 {
        &self.r_ab_mag
    }
    pub fn get_slice(&self) -> &usize {
        &self.slice
    }
    pub fn new() -> SimulationParams {
        SimulationParams {
            prog_name: String::from(""),
            blue_detector: String::from("ML15"),
            red_detector: String::from("ML15"),
            airmass: 1.0,
            moon: 0.0,
            exposure: 3600.0,
            r_ab_mag: 18.0,
            slice: 20,
        }
    }
}

pub struct Simulation {
    input: Spectrum,
    cousins_r: Curve,
    cousins_re_quiv_bw: f64,
    sky: Spectrum,
    sky_model: SkyModel,
    tarsis_model: InstrumentModel,
    det: Detector,
    params: SimulationParams,
}

impl Simulation {
    // Getter functions (for unit tests)
    pub fn get_input(&self) -> &Spectrum {
        &self.input
    }
    pub fn get_cousins_r(&self) -> &Curve {
        &self.cousins_r
    }
    pub fn get_cousins_re_quiv_bw(&self) -> &f64 {
        &self.cousins_re_quiv_bw
    }
    pub fn get_sky(&self) -> &Spectrum {
        &self.sky
    }
    pub fn get_sky_model(&self) -> &SkyModel {
        &self.sky_model
    }
    pub fn get_tarsis_model(&self) -> &InstrumentModel {
        &self.tarsis_model
    }
    pub fn get_det(&self) -> &Detector {
        &self.det
    }
    pub fn get_params(&self) -> &SimulationParams {
        &self.params
    }

    pub fn new() -> Simulation {
        let mut simulation = Simulation {
            input: Spectrum::default(),
            cousins_r: Curve::default(),
            cousins_re_quiv_bw: 0.0,
            sky: Spectrum::default(),
            sky_model: SkyModel::new(),
            tarsis_model: InstrumentModel::new(),
            det: Detector::new(),
            params: SimulationParams::new(),
        };

        let _ = simulation.cousins_r.load_curve("src/simulation/Generic_Cousins_R.txt");
        simulation.cousins_r.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.cousins_r.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency

        simulation.cousins_re_quiv_bw = simulation.cousins_r.integral();

        return simulation;
    }

    pub fn set_input(&mut self, spec: Spectrum) {
        self.input = spec;
    }

    pub fn normalize_to_r_mag(&mut self, mag_r: f64) {
        let mut filtered = Spectrum::default();
        filtered.get_curve().from_existing(&self.input.get_curve(), 1.0);
        //filtered.invert_axis(XAxis, SPEED_OF_LIGHT);
        SpectrumAxisOperations::invert_axis(&mut filtered, XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());
        filtered.multiply(&self.cousins_r);

        let desired_sb = surface_brightness_ab2freq_radiance(mag_r);
        let mean_sb = filtered.integral() / self.cousins_re_quiv_bw;

        self.input.scale_axis(YAxis, desired_sb / mean_sb);
    }

    pub fn set_params(&mut self, params: SimulationParams) {
        self.params = params.clone();

        self.sky_model.set_airmass(params.airmass);
        self.sky_model.set_moon(params.moon);

        if !self.sky.get_curve().get_curve().is_empty() {
            self.sky.get_curve().clear();
            self.sky = Spectrum::default();
        }

        // Update sky spectrum
        self.sky = self.sky_model.make_sky_spectrum(&self.input);

        // Update detector config
        self.det.set_exposure_time(params.exposure)
 
    }

    pub fn simulate_arm(&mut self, arm: InstrumentArm) { 
        let mut tarsis_prop: &mut InstrumentProperties = self.tarsis_model.properties_mut();
        let mut det_name: &str;

        match arm {
            InstrumentArm::BlueArm => {
                det_name = &self.params.blue_detector;
            }
            InstrumentArm::RedArm => {
                det_name = &self.params.red_detector;
            }
        }

        // Set coating
        *tarsis_prop.get_coating_mut() = self.det.get_spec().get_coating().clone();
        self.tarsis_model.set_input(arm, self.sky.clone());
        let flux = self.tarsis_model.make_pixel_photon_flux(self.params.slice);
        self.det.set_pixel_photon_flux(flux);        
    }

    pub fn signal(&mut self, px: NotNan<f64>) -> f64 {
        return self.det.signal_px(px);
    }

    pub fn noise(&mut self, px: NotNan<f64>) -> f64 {
        return self.det.noise(px);
    }

    pub fn electrons(&mut self, px: NotNan<f64>) -> f64 {
        return self.det.electrons_px(px);
    }

    pub fn read_out_noise(&mut self) -> f64 {
        return self.det.read_out_noise();
    }

    pub fn gain(&mut self) -> f64 {
        let spec = self.det.get_spec();
        
        if spec == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return *spec.get_gain();
    }
    pub fn px_to_wavelength(&mut self, arm: InstrumentArm, px: NotNan<f64>) -> f64 {
        return self.tarsis_model.px_to_wavelength_val(arm, self.params.slice, px).unwrap();

    }
    pub fn wl_to_pixel_curve(&mut self, arm: InstrumentArm) -> Curve {
        return self.tarsis_model.wavelength_to_pix_crv(arm, self.params.slice).unwrap();

    }
}

