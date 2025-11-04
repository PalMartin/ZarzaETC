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
use serde_yaml::Value;
use ordered_float::Pow;
use std::collections::HashMap;
use crate::spectrum::*;
use crate::curve::*;
use crate::config_manager::*;
use crate::curve::CurveAxis::YAxis;

// Constants
/// Number of pixels in the detector
pub const DETECTOR_PIXELS: f64 = 2048.0;    
/// Temperature of the detector in Kelvin 
pub const DETECTOR_TEMPERATURE: f64 = 193.0;   // K

/// Structure to hold the specifications of a detector
#[derive(Default, Debug, Clone, PartialEq)]  
pub struct DetectorSpec {
    /// Configuration manager instance for loading and saving specifications.
    pub config: Config,

    /// Coating type of the detector.
    pub coating: String,
    /// Side length of a pixel in meters.
    pub pixel_side: f64,
    /// Read-out noise of the detector in electrons.
    pub read_out_noise: f64,
    /// Gain of the detector (electrons per ADU).
    pub gain: f64,

}

/// Implementations of methods for DetectorSpec.
impl DetectorSpec {

    /// Constructor to create a new DetectorSpec with default values and load configuration from a file.
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

    /// Serialize the current state of the DetectorSpec to the configuration manager for saving to file.
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

    /// Deserialize the DetectorSpec from the configuration manager, loading values from the YAML file.
    pub fn deserialize(&mut self) -> bool {

        self.config.deserialize_field(&mut self.coating, "coating");
        self.config.deserialize_field(&mut self.pixel_side, "pixel_side");
        self.config.deserialize_field(&mut self.read_out_noise, "read_out_noise");
        self.config.deserialize_field(&mut self.gain, "gain");

        return true; 
    }

    /// Inmutable acces to the coating of the detector.
    pub fn get_coating(&self) -> &String {
        &self.coating
    }
    /// Mutable acces to the coating of the detector.
    pub fn get_coating_mut(&mut self) -> &mut String {
        &mut self.coating
    }
    /// Inmutable acces to the pixel side of the detector.
    pub fn get_pixel_side(&self) -> &f64 {
        &self.pixel_side
    }
    /// Inmutable acces to the read out noise of the detector.
    pub fn get_read_out_noise(&self) -> &f64 {
        &self.read_out_noise
    }
    /// Inmutable acces to the gain of the detector.
    pub fn get_gain(&self) -> &f64 {
        &self.gain
    }

}

/// Structure to hold the properties of all detectors
#[derive(Default, Debug, Clone)]
pub struct DetectorProperties {
    /// Configuration manager instance for loading and saving properties.
    pub config: Config,
    /// Map of detector names to their specifications.
    pub detectors: HashMap<String, DetectorSpec> // not sure if the pointer in "std::map<std::string, DetectorSpec *> detectors;" is accounted for.

}

/// Implementations of methods for DetectorProperties.
impl DetectorProperties {

    /// Create a new DetectorProperties instance with default values and load configuration from a file.
    pub fn new(name: String) -> DetectorProperties {
        /// Initialize configuration manager with the given name.
        let mut properties = DetectorProperties {
            config: Config::new(name),
            detectors: HashMap::new(),
        };

        properties.deserialize();
        return properties;
    }

    /// Serialize the current properties to the configuration manager for saving to file.
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

    /// Deserialize properties from the configuration manager, loading values from file.
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

    /// Clear all detectors from the properties.
    pub fn clear_detectors(&mut self) {
        self.detectors.clear();
    }

    /// Inmutable access to the map of detectors.
    pub fn get_detectors(&self) -> &HashMap<String, DetectorSpec> {
        &self.detectors
    }
}

/// Main Detector structure
#[derive(Default, Debug, Clone)]  
pub struct Detector {

    /// Detector properties.
    properties: DetectorProperties,
    /// Detector specifications.
    detector: DetectorSpec,

    /// Exposure time in seconds.
    exposure_time: f64,
    /// Object photon flux per pixel (photons per pixel per square meter per second).
    obj_photon_flux_per_pixel: Spectrum,  // ph / (px m^2 s)
    /// Sky photon flux per pixel (photons per pixel per square meter per second).
    sky_photon_flux_per_pixel: Spectrum,  // ph / (px m^2 s)
    /// Object photons per pixel (photons per pixel).
    obj_photons_per_pixel: Spectrum,      // ph/px
    /// Sky photons per pixel (photons per pixel).
    sky_photons_per_pixel: Spectrum,      // ph/px
    /// Object electrons per pixel (electrons per pixel).
    obj_electrons_per_pixel: Spectrum,    // e/px
    /// Sky electrons per pixel (electrons per pixel).
    sky_electrons_per_pixel: Spectrum,    // e/px
    /// Total electrons per pixel (electrons per pixel).
    electrons_per_pixel: Spectrum,    // e/px
    /// Signal (ADU per pixel).
    signal: Spectrum,                 // c

}

/// Implementations of methods for Detector.
impl Detector {

    /// Create a new Detector instance with default values and load properties from a file.
    pub fn new() -> Detector {

        // Load properties from file
        let mut properties = DetectorProperties::new("detectors".to_string());
        properties.config.load();
        properties.deserialize();

        // Initialize the Detector with default values
        let detector = Detector {
            properties,
            detector: Default::default(),

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

    /// Inmutable access to the detector properties.
    pub fn properties(&self) -> &DetectorProperties {
        &self.properties
    }

    /// Get a clone of the current detector specification.
    pub fn get_spec(&self) -> DetectorSpec {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.detector.clone();
    }

    /// Set the photon flux per pixel for both the object and sky, and recalculate dependent values.
    pub fn set_pixel_photon_flux(&mut self, flux_obj: Spectrum, flux_sky: Spectrum) {
        self.obj_photon_flux_per_pixel.assign(&flux_obj.get_curve());
        self.sky_photon_flux_per_pixel.assign(&flux_sky.get_curve());
        self.recalculate();
    }

    /// Set the exposure time and recalculate dependent values.
    pub fn set_exposure_time(&mut self, t: f64) {
        self.exposure_time = t
    }

    /// Set the detector by name from the properties. Returns true if successful, false if the detector name is not found.
    pub fn set_detector(&mut self, det: String) -> bool {  
        
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

    /// Calculate the number of dark electrons based on the detector temperature and exposure time.
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

    /// Recalculate all dependent values based on the current detector specifications, exposure time, and photon fluxes. 
    pub fn recalculate(&mut self) { 
        if self.detector == DetectorSpec::default() {
            println!("Detector is not set.");
        }

        let inv_gain = 1.0 / self.detector.gain;
        // Compute electrons in each pixel by means of the exposure time.
        self.obj_photons_per_pixel.from_existing(&self.obj_photon_flux_per_pixel.get_curve(), self.exposure_time * self.detector.pixel_side * self.detector.pixel_side);
        self.sky_photons_per_pixel.from_existing(&self.sky_photon_flux_per_pixel.get_curve(), self.exposure_time * self.detector.pixel_side * self.detector.pixel_side);
        // Compute number of electrons by means of the quantum efficiency (thq auqntum efficiency is accounted for in the transmission curve of the provided data, therefore we use a value of 1.0).
        self.obj_electrons_per_pixel.from_existing(&self.obj_photons_per_pixel.get_curve(), 1.0); // y_units = gain?
        self.sky_electrons_per_pixel.from_existing(&self.sky_photons_per_pixel.get_curve(), 1.0); // y_units = gain?
        // Convert electrons to ADU (signal).
        self.signal.from_existing(&self.obj_electrons_per_pixel.get_curve(), inv_gain);
        // Add dark electrons to the electron-per-pixel curve.
        self.electrons_per_pixel.from_existing(&self.obj_electrons_per_pixel.get_curve(), 1.0); 
        self.electrons_per_pixel.add(&self.sky_electrons_per_pixel.get_curve());
        self.electrons_per_pixel.add(self.dark_electrons(DETECTOR_TEMPERATURE));
    }

    /// Get the signal value at a specific wavelength.
    pub fn signal_value(&self, wl: NotNan<f64>) -> f64 {
        return self.signal.get_point(wl);
    }

    /// Get a clone of the signal curve.
    pub fn signal_crv(&self) -> Spectrum {
        return self.signal.clone();
    }

    /// Get the electrons value at a specific wavelength.
    pub fn electrons_value(&self, wl: NotNan<f64>) -> f64 {
        return self.electrons_per_pixel.get_point(wl);
    }

    /// Get a clone of the electrons curve.
    pub fn electrons_crv(&self) -> Spectrum {
        return self.electrons_per_pixel.clone();
    }

    /// Get the read out noise in ADU.
    pub fn read_out_noise(&self) -> f64 { 
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }
        return self.detector.read_out_noise / self.detector.gain;
    }

    /// Get the noise value at a specific wavelength.
    pub fn noise_value(&self, wl: NotNan<f64>) -> f64 {
        if self.detector == DetectorSpec::default() {
            panic!("Detector is not set.");
        }
    
        let mut ron2 = self.detector.read_out_noise;
        let inv_gain_2 = 1.0 / (self.detector.gain * self.detector.gain);
        ron2 *= ron2;
    
        return (inv_gain_2 * self.electrons_value(wl) + ron2).sqrt(); // MISSING THE RANDOM NOISE
    }

    /// Get a clone of the noise curve.
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

    /// Get the signal-to-noise ratio value at a specific wavelength.
    pub fn snr_value(&self, wl: NotNan<f64>) -> f64 { 
        return self.signal.get_point(wl) / self.noise_value(wl);
    }

    /// Get a mutable reference to the signal-to-noise ratio curve.
    pub fn snr_crv(&mut self) -> &mut Curve { 
        let noise_curve = self.noise_crv(); // binding
        let mut noise_inv = noise_curve;
        noise_inv.invert_axis(YAxis, NotNan::new(1.0).unwrap());
        
        let snr = self.signal.get_curve_mut(); 
        snr.multiply(&*noise_inv.get_curve_mut());
    
        return snr;
    }

    

    // Inmutable access to the detector properties.
    pub fn get_properties(&self) -> &DetectorProperties {
        &self.properties
    }
    // Inmutable access to the current detector specification.
    pub fn get_detector(&self) -> &DetectorSpec {
        &self.detector
    }
    // Inmutable access to the exposure time.
    pub fn get_exposure_time(&self) -> &f64 {
        &self.exposure_time
    }
    // Inmutable access to the object photon flux per pixel.
    pub fn get_obj_photon_flux_per_pixel(&self) -> &Spectrum {
        &self.obj_photon_flux_per_pixel
    }
    // Inmutable access to the sky photon flux per pixel.
    pub fn get_sky_photon_flux_per_pixel(&self) -> &Spectrum {
        &self.sky_photon_flux_per_pixel
    }
    // Inmutable access to the object photons per pixel.
    pub fn get_obj_photons_per_pixel(&self) -> &Spectrum {
        &self.obj_photons_per_pixel
    }
    // Inmutable access to the sky photons per pixel.
    pub fn get_sky_photons_per_pixel(&self) -> &Spectrum {
        &self.sky_photons_per_pixel
    }
    // Inmutable access to the object electrons per pixel.
    pub fn get_obj_electrons_per_pixel(&self) -> &Spectrum {
        &self.obj_electrons_per_pixel
    }
    // Inmutable access to the sky electrons per pixel.
    pub fn get_sky_electrons_per_pixel(&self) -> &Spectrum {
        &self.sky_electrons_per_pixel
    }
    // Inmutable access to the total electrons per pixel.
    pub fn get_electrons_per_pixel(&self) -> &Spectrum {
        &self.electrons_per_pixel
    }
    // Inmutable access to the signal.
    pub fn get_signal(&self) -> &Spectrum {
        &self.signal
    }
}

