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
use serde_yaml; 
use serde_yaml::Value; // this is for the YAML nodes as a generic type
use ordered_float::Pow;
use std::collections::HashMap;

use crate::spectrum::*;
use crate::curve::*;
use crate::config_manager::*;
use crate::data_file_manager::*;
use crate::curve::CurveAxis::XAxis;
use crate::curve::CurveAxis::YAxis;

pub const DETECTOR_PIXELS: f64 = 2048.0;     
pub const DETECTOR_TEMPERATURE: f64 = 193.0;   // K

#[derive(Default, Debug, Clone, PartialEq)]  
pub struct DetectorSpec {

    pub config: Config,

    pub coating: String,
    pub pixel_side: f64,
    pub read_out_noise: f64,
    pub gain: f64,
    //pub q_e: f64, -> is now a curve in Detector

}

impl DetectorSpec {

    pub fn new(name: String) -> Self {
        let mut config = Config::new(name.clone());
        config.load();

        let mut spec = DetectorSpec {
            config,
            coating: "ML15".to_string(),
            pixel_side: 15e-6,
            read_out_noise: 0.0,
            gain: 1.0,
        };

        spec.deserialize(); // override with values from YAML if they exist
        return spec;
    } 

    pub fn serialize(&mut self) -> bool { // CONSIDER WHAT HAPPENS IF YAML_CONFIG IS NOT A MAPPING -> ERROR HANDLING?

        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "coating".into(), Config::serialize_field(&self.coating)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "pixel_side".into(), Config::serialize_field(&self.pixel_side)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "read_out_noise".into(), Config::serialize_field(&self.read_out_noise)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "gain".into(), Config::serialize_field(&self.gain)
        );

        return true
    }

    pub fn deserialize(&mut self) -> bool {

        self.config.deserialize_field(&mut self.coating, "coating");
        self.config.deserialize_field(&mut self.pixel_side, "pixel_side");
        self.config.deserialize_field(&mut self.read_out_noise, "read_out_noise");
        self.config.deserialize_field(&mut self.gain, "gain");

        return true; 
    }

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

}

#[derive(Default, Debug, Clone)]
pub struct DetectorProperties {

    pub config: Config,

    pub detectors: HashMap<String, DetectorSpec> // not sure if the pointer in "std::map<std::string, DetectorSpec *> detectors;" is accounted for.

}
impl DetectorProperties {

    pub fn new(name: String) -> DetectorProperties {
        let mut properties = DetectorProperties {
            config: Config::new(name),
            detectors: HashMap::new(),
        };

        properties.deserialize();
        return properties;
    }

    pub fn serialize(&mut self) -> bool {

        let mut yaml_detectors: HashMap<String, Value> = HashMap::new();

        for (k, v) in self.detectors.iter_mut() {
            v.serialize();
            yaml_detectors.insert(k.into(), v.config.yaml_node());
        }

        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "detectors".into(), Config::serialize_field(&yaml_detectors)
        );

        return true;
    }

    pub fn deserialize(&mut self) -> bool {

        let mut yaml_detectors: HashMap<String, Value> = HashMap::new();

        self.clear_detectors();

        if self.config.deserialize_field(&mut yaml_detectors, "detectors") {
            for (k, v) in yaml_detectors.iter_mut() {
                let mut spec = DetectorSpec::new("detectors.".to_string() + &k);
                if !spec.config.deserialize_yaml_node(&v) {
                    eprintln!("{}: failed to deserialize detector.", k);
                    // spec dropped
                }
                self.detectors.insert(k.to_string(), spec);
            }
        }
    
        return true;
    }

    pub fn clear_detectors(&mut self) {
        self.detectors.clear();
    }

    // getter functions
    pub fn get_detectors(&self) -> &HashMap<String, DetectorSpec> {
        &self.detectors
    }
}

#[derive(Default, Debug, Clone)]  
pub struct Detector {

    properties: DetectorProperties,
    detector: DetectorSpec,

    //q_e: Curve, I THINK THE QUANTUM EFFICIENCY IS INCLUDED IN THE TRANSMISSION CURVE, REVIEW AND UNCOMMENT IF NECESSARY
    exposure_time: f64,
    obj_photon_flux_per_pixel: Spectrum,  // ph / (px m^2 s)
    sky_photon_flux_per_pixel: Spectrum,  // ph / (px m^2 s)
    obj_photons_per_pixel: Spectrum,      // ph/px
    sky_photons_per_pixel: Spectrum,      // ph/px
    obj_electrons_per_pixel: Spectrum,    // e/px
    sky_electrons_per_pixel: Spectrum,    // e/px
    electrons_per_pixel: Spectrum,    // e/px
    signal: Spectrum,                 // c

}
impl Detector {

    pub fn new() -> Detector {

        let mut properties = DetectorProperties::new("detectors".to_string());
        properties.config.load();
        properties.deserialize();

        let detector = Detector {
            properties,
            detector: Default::default(),

            //q_e: Default::default(),
            exposure_time: 3600.0,
            obj_photon_flux_per_pixel: Default::default(),
            sky_photon_flux_per_pixel: Default::default(),
            obj_photons_per_pixel: Default::default(),
            sky_photons_per_pixel: Default::default(),
            obj_electrons_per_pixel: Default::default(),
            sky_electrons_per_pixel: Default::default(),
            electrons_per_pixel: Default::default(),
            signal: Default::default(),
        };

        return detector;
    }

    pub fn properties(&self) -> &DetectorProperties {
        &self.properties
    }

    pub fn get_spec(&self) -> DetectorSpec {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.detector.clone();
    }

    pub fn set_pixel_photon_flux(&mut self, flux_obj: Spectrum, flux_sky: Spectrum) {
        self.obj_photon_flux_per_pixel.assign(&flux_obj.get_curve());
        self.sky_photon_flux_per_pixel.assign(&flux_sky.get_curve());
        self.recalculate();
    }

    pub fn set_exposure_time(&mut self, t: f64) {
        self.exposure_time = t
    }

    pub fn set_detector(&mut self, det: String) -> bool {  // IMPORTANT: CHANGE IF THE Q_E CURVE IS SEPARATED FROM THE TRANSMISSION CURVE

        //if det == "ccd231_84_0_s77".to_string() {
        //    self.q_e.load_curve(&DataFileManager::data_file(&"NBB_QE_2.csv".to_string())).unwrap();
        //    self.q_e.scale_axis(XAxis, 1e-9); // X axis was in nm
        //} else if det == "ccd231_84_0_h69".to_string() {
        //    self.q_e.load_curve(&DataFileManager::data_file(&"ML15_QE_2.csv".to_string())).unwrap();
        //    self.q_e.scale_axis(XAxis, 1e-9); // X axis was in nm
        //} else if det == "ccd231_84_0_h69+dd".to_string() {
        //    self.q_e.load_curve(&DataFileManager::data_file(&"ML15+DD_QE".to_string())).unwrap(); // For HR-R
        //    self.q_e.scale_axis(XAxis, 1e-9); // X axis was in nm
        //} else {
        //    eprintln!("Unknown detector {}; expected 'ccd231_84_0_s77', 'ccd231_84_0_h69' or 'ccd231_84_0_h69'.", det);
        //}
        
        let it = self.properties.detectors.get(&det);

        match it {
            Some(spec) => {
                self.detector = spec.clone();
                self.detector.deserialize();
                return true;
            }
            None => {
                self.detector = Default::default();
                return false;
            }
        }
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
        self.obj_photons_per_pixel.from_existing(&self.obj_photon_flux_per_pixel.get_curve(), self.exposure_time * self.detector.pixel_side * self.detector.pixel_side);
        self.sky_photons_per_pixel.from_existing(&self.sky_photon_flux_per_pixel.get_curve(), self.exposure_time * self.detector.pixel_side * self.detector.pixel_side);
        // Compute number of electrons by means of the quantum efficiency
        self.obj_electrons_per_pixel.from_existing(&self.obj_photons_per_pixel.get_curve(), 1.0); // y_units = gain?
        self.sky_electrons_per_pixel.from_existing(&self.sky_photons_per_pixel.get_curve(), 1.0); // y_units = gain?
        // Turn this into counts
        self.signal.from_existing(&self.obj_electrons_per_pixel.get_curve(), inv_gain);
        // Add dark electrons to the electron-per-pixel curve
        self.electrons_per_pixel.from_existing(&self.obj_electrons_per_pixel.get_curve(), 1.0); 
        self.electrons_per_pixel.add(&self.sky_electrons_per_pixel.get_curve());
        self.electrons_per_pixel.add(self.dark_electrons(DETECTOR_TEMPERATURE));
    }

    pub fn signal_value(&self, wl: NotNan<f64>) -> f64 {
        return self.signal.get_point(wl);
    }

    pub fn signal_crv(&self) -> Spectrum {
        return self.signal.clone();
    }

    pub fn electrons_value(&self, wl: NotNan<f64>) -> f64 {
        return self.electrons_per_pixel.get_point(wl);
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

    pub fn noise_value(&self, wl: NotNan<f64>) -> f64 {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }
    
        let mut ron2 = self.detector.read_out_noise;
        let inv_gain_2 = 1.0 / (self.detector.gain * self.detector.gain);
        ron2 *= ron2;
    
        return (inv_gain_2 * self.electrons_value(wl) + ron2).sqrt(); // MISSING THE RANDOM NOISE
    }

    pub fn noise_crv(&self) -> Spectrum { 
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }
        
        let mut ron2 = self.detector.read_out_noise;
        let inv_gain_2 = 1.0 / (self.detector.gain * self.detector.gain);
        ron2 *= ron2;
    
        let mut elec = self.electrons_crv();
        elec.multiply(inv_gain_2);
        elec.add(ron2);
        // MISSING THE RANDOM NOISE
        for (_x, y) in elec.get_curve_mut().get_map_mut().iter_mut() {
            *y = y.sqrt();  // Modify the value directly
        }
    
        return elec;
    }

    pub fn snr_value(&self, wl: NotNan<f64>) -> f64 { 
        return self.signal.get_point(wl) / self.noise_value(wl);
    }

    pub fn snr_crv(&mut self) -> &mut Curve { 
        let noise_curve = self.noise_crv(); // binding
        let mut noise_inv = noise_curve;
        noise_inv.invert_axis(YAxis, NotNan::new(1.0).unwrap());
        
        let snr = self.signal.get_curve_mut(); 
        snr.multiply(&*noise_inv.get_curve_mut());
    
        return snr;
    }

    

    // getter functions
    pub fn get_properties(&self) -> &DetectorProperties {
        &self.properties
    }
    pub fn get_detector(&self) -> &DetectorSpec {
        &self.detector
    }
    //pub fn get_q_e(&self) -> &Curve {
    //    &self.q_e
    //}
    pub fn get_exposure_time(&self) -> &f64 {
        &self.exposure_time
    }
    pub fn get_obj_photon_flux_per_pixel(&self) -> &Spectrum {
        &self.obj_photon_flux_per_pixel
    }
    pub fn get_sky_photon_flux_per_pixel(&self) -> &Spectrum {
        &self.sky_photon_flux_per_pixel
    }
    pub fn get_obj_photons_per_pixel(&self) -> &Spectrum {
        &self.obj_photons_per_pixel
    }
    pub fn get_sky_photons_per_pixel(&self) -> &Spectrum {
        &self.sky_photons_per_pixel
    }
    pub fn get_obj_electrons_per_pixel(&self) -> &Spectrum {
        &self.obj_electrons_per_pixel
    }
    pub fn get_sky_electrons_per_pixel(&self) -> &Spectrum {
        &self.sky_electrons_per_pixel
    }
    pub fn get_electrons_per_pixel(&self) -> &Spectrum {
        &self.electrons_per_pixel
    }
    pub fn get_signal(&self) -> &Spectrum {
        &self.signal
    }
}

