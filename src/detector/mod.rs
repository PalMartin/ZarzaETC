//
// detector/mod.cpp: Detector model
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
use std::fs::File;
use std::io::Read;
use serde_yaml; 
use ordered_float::Pow;
use serde::{Serialize, Deserialize};

pub const DETECTOR_PIXELS: f64 = 2048.0;     
pub const DETECTOR_TEMPERATURE: f64 = 193.0;   // K

#[derive(Default, Debug, Clone, Serialize, Deserialize, PartialEq)]  
pub struct DetectorSpec {
    coating: String,
    pixel_side: f64,
    read_out_noise: f64,
    gain: f64,
    q_e: f64,
}
impl DetectorSpec {
    // Getter functions
    pub fn get_coating(&self) -> &String {
        &self.coating
    }
    pub fn get_coating_mut(&mut self) -> &mut String {
        &mut self.coating
    }
    pub fn get_pixel_side(&self) -> &f64 {
        &self.pixel_side
    }
    pub fn get_read_out_noise(&self) -> &f64 {
        &self.read_out_noise
    }
    pub fn get_gain(&self) -> &f64 {
        &self.gain
    }
    pub fn get_q_e(&self) -> &f64 {
        &self.q_e
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]  
pub struct DetectorProperties {
    ccd231_84_0_s77: DetectorSpec,
    ccd231_84_0_h69: DetectorSpec,
}
impl DetectorProperties {
    pub fn new() -> Self {
        let mut file = File::open("src/detector/detector.yaml").expect("Failed to open detector.yaml");
        let mut yaml_content = String::new();
        file.read_to_string(&mut yaml_content).expect("Failed to read file");
        let detector_properties: DetectorProperties = serde_yaml::from_str(&yaml_content).expect("Failed to parse YAML");
        return detector_properties;
    }
    // Getter functions
    pub fn get_ccd231_84_0_s77(&self) -> &DetectorSpec {
        &self.ccd231_84_0_s77
    }
    pub fn get_ccd231_84_0_h69(&self) -> &DetectorSpec {
        &self.ccd231_84_0_h69
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]  
pub struct Detector {
    properties: DetectorProperties,
    detector: DetectorSpec,
    exposure_time: f64,
    photon_flux_per_pixel: Spectrum,  // ph / (px m^2 s)
    photons_per_pixel: Spectrum,      // ph/px
    electrons_per_pixel: Spectrum,    // e/px
    signal: Spectrum,                 // c
}
impl Detector {
    pub fn get_properties(&self) -> &DetectorProperties {
        &self.properties
    }
    pub fn get_detector(&self) -> &DetectorSpec {
        &self.detector
    }
    pub fn get_exposure_time(&self) -> &f64 {
        &self.exposure_time
    }
    pub fn get_photon_flux_per_pixel(&self) -> &Spectrum {
        &self.photon_flux_per_pixel
    }
    pub fn get_photons_per_pixel(&self) -> &Spectrum {
        &self.photons_per_pixel
    }
    pub fn get_electrons_per_pixel(&self) -> &Spectrum {
        &self.electrons_per_pixel
    }
    pub fn get_signal(&self) -> &Spectrum {
        &self.signal
    }
    pub fn new() -> Detector {
        Detector {
            properties: DetectorProperties::new(),
            detector: Default::default(),
            exposure_time: 1.0,
            photon_flux_per_pixel: Default::default(),
            photons_per_pixel: Default::default(),
            electrons_per_pixel: Default::default(),
            signal: Default::default(),
        }
    }
    pub fn properties(&self) -> &DetectorProperties {
        &self.properties
    }

    pub fn set_pixel_photon_flux(&mut self, flux: Spectrum) {
        self.photon_flux_per_pixel.assign(&flux.get_curve());
        self.recalculate();
    }

    pub fn set_exposure_time(&mut self, t: f64) {
        self.exposure_time = t
    }

    pub fn set_detector(&mut self, name: &str) {
        if name == "CCD231-84-0-S77" {
            self.detector = self.properties.ccd231_84_0_s77.clone();
        } else if name == "CCD231-84-0-H69" {
            self.detector = self.properties.ccd231_84_0_h69.clone();
        } else {
            panic!("Unknown detector {}; expected 'CCD231-84-0-S77' or 'CCD231-84-0-H69'.", name);
        }
    }    

    pub fn get_spec(&self) -> DetectorSpec {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.detector.clone();
    }

    pub fn dark_electrons(&self, t: f64) -> f64 {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set");
        }

        let area  = self.detector.pixel_side * self.detector.pixel_side;
        let qd0   = 6.2415091e+13 * area; // e/s/m^2 = 1 nA / cm^2;
        let t_beta = 6400.0;
        let slope = 122.0;

        let qd = self.exposure_time * qd0 * slope * t.pow(3) * (-t_beta / t).exp();

        return qd;
    }

    pub fn recalculate(&mut self) {  
        if self.detector == DetectorSpec::default() {
            println!("Detector is not set.");
        }

        let inv_gain = 1.0 / self.detector.gain;

        // Compute electrons in each pixel by means of the exposure time
        self.photons_per_pixel.from_existing(&self.photon_flux_per_pixel.get_curve(), self.exposure_time * self.detector.pixel_side * self.detector.pixel_side);
        // Compute number of electrons by means of the quantum efficiency
        self.electrons_per_pixel.from_existing(&self.photons_per_pixel.get_curve(), self.detector.q_e);
        // Turn this into counts
        self.signal.from_existing(&self.electrons_per_pixel.get_curve(), inv_gain);
        // Add dark electrons to the electron-per-pixel curve
        self.electrons_per_pixel.add(self.dark_electrons(DETECTOR_TEMPERATURE));
    }

    pub fn signal_px(&self, px: NotNan<f64>) -> f64 {
        return self.signal.get_curve().get_point(px);
    }

    pub fn electrons_px(&self, px: NotNan<f64>) -> f64 {
        return self.electrons_per_pixel.get_curve().get_point(px);
    }

    pub fn signal_crv(&self) -> Spectrum {
        return self.signal.clone();
    }

    pub fn electrons_crv(&self) -> Spectrum {
        return self.electrons_per_pixel.clone();
    }

    pub fn read_out_noise(&self) -> f64 { 
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }
        return self.detector.read_out_noise / self.detector.gain;
    }

    pub fn noise(&self, px: NotNan<f64>) -> f64 {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        let mut ron2 = self.detector.read_out_noise;
        let inv_gain_2 = 1.0 / (self.detector.gain * self.detector.gain);

        ron2 *= ron2;

        return (inv_gain_2 * self.electrons_px(px) + ron2).sqrt();
    }

    pub fn snr(&self, px: NotNan<f64>) -> f64 {
        return self.signal.get_curve().get_point(px) / self.noise(px);
    }
}

