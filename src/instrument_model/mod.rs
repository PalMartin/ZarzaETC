//
// instrument_model/mod.cpp: Instrument model
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


use crate::curve;
use crate::spectrum::*;
use crate::curve::*;
use curve::CurveAxis::*;
use crate::helpers::*;
use crate::config_manager::*;
use crate::data_file_manager::*;
use std::str::FromStr;

// Constants
/// Diameter of the CAHA telescope aperture in meters.
pub const CAHA_APERTURE_DIAMETER: f64 = 3.5;      // m
/// Focal length of the CAHA telescope in meters.
pub const CAHA_FOCAL_LENGTH: f64 = 12.195;        // m
/// Aperture area of the CAHA telescope in square meters.
pub const CAHA_APERTURE_AREA: f64 = 0.25 * std::f64::consts::PI * CAHA_APERTURE_DIAMETER * CAHA_APERTURE_DIAMETER;
/// Effective area of the CAHA telescope in square meters.
pub const CAHA_EFFECTIVE_AREA: f64 = 9.093;       // m^2
/// Number of slices in TARSIS instrument.
pub const TARSIS_SLICES: usize = 40;              // Number of TARSIS slices per FOV
/// Spectral pixel length (number of pixels in the spectral direction).
pub const SPECTRAL_PIXEL_LENGTH: f64 = 2048.0;    // Guessed from the resolution elements and ranges


/// Structure to hold instrument properties loaded from configuration files and used in the instrument model.
#[derive(Default, Clone)]
pub struct InstrumentProperties {
    /// Configuration manager instance for loading and saving properties.
    pub config: Config,

    /// Focal ratio of the telescope.
    f_num: f64,
    /// Efficiency of the telescope's aperture.
    ap_efficiency: f64,
    /// Coating type.
    coating: String,
}

/// Enum to represent the two arms of the instrument.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum InstrumentArm {
    BlueArm,
    RedArm,
}

/// Implement string parsing for InstrumentArm to allow easy conversion from string inputs.
impl FromStr for InstrumentArm {
    type Err = String;

    /// Converts a string to an InstrumentArm enum variant.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "red" => Ok(InstrumentArm::RedArm),
            "blue" => Ok(InstrumentArm::BlueArm),
            _ => Err(format!("Invalid arm: {}", s)),
        }
    }
}

/// Structure representing the instrument model, including properties and various curves for dispersion and resolution.
#[derive(Clone)]
pub struct InstrumentModel {
    /// Instrument properties loaded from configuration.
    properties: InstrumentProperties, 
    /// Transmission curve for the blue arm with ML15 coating.
    blue_ml15: Curve, 
    /// Transmission curve for the blue arm with NBB coating.
    blue_nbb: Curve, 
    /// Transmission curve for the red arm with ML15 coating.
    red_ml15: Curve, 

    ///////////////////////////// Blue arm ///////////////////////////
    /// Spectral dispersion curves for each slice in the blue arm.
    blue_disp: Vec<Curve>,    // Spectral dispersion
    /// Resolution element curves for each slice in the blue arm.
    blue_repx: Vec<Curve>,    // Resolution element (in pixels)
    /// Wavelength to pixel curves for each slice in the blue arm.
    blue_w2px: Vec<Curve>,    // int of inv of sd
    /// Pixel to wavelength curves for each slice in the blue arm.
    blue_px2w: Vec<Curve>,    // Inverse of above

    ////////////////////////////// Red arm ///////////////////////////
    /// Spectral dispersion curves for each slice in the red arm.
    red_disp: Vec<Curve>,     // Spectral dispersion
    /// Resolution element curves for each slice in the red arm.
    red_repx: Vec<Curve>,     // Resolution element (in pixels)
    /// Wavelength to pixel curves for each slice in the red arm.
    red_w2px: Vec<Curve>,     // int of inv of sd
    /// Pixel to wavelength curves for each slice in the red arm.
    red_px2w: Vec<Curve>,     // Inverse of above

    /// Attenuated input spectrum before the dispersor.
    atten_spectrum: Spectrum,       // Attenuated before dispersor
    /// Current arm/path of the attenuated spectrum.
    current_path: InstrumentArm,    // Path of the attenuated spectrum
}

/// Implementation of methods for InstrumentProperties, including loading and saving configuration.
impl InstrumentProperties {
    /// Create a new InstrumentProperties instance with default values and load configuration from a file.
    pub fn new(name: String) -> InstrumentProperties {
        // Initialize configuration manager with the given name.
        let mut config = Config::new(name.clone());
        config.load();

        // Set default values for the instrument properties.
        let mut instrument_properties = InstrumentProperties {
            // Initialize with the configuration manager.
            config,
            // Default focal ratio based on CAHA telescope parameters.
            f_num: CAHA_FOCAL_LENGTH / CAHA_APERTURE_DIAMETER,
            // Default aperture efficiency based on CAHA telescope parameters.
            ap_efficiency: CAHA_EFFECTIVE_AREA / CAHA_APERTURE_AREA,
            // Default coating type.
            coating: "ML15".to_string(),
        };

        // Load values from the configuration file, overriding defaults from cotrol YAML if present.
        instrument_properties.deserialize();
        return instrument_properties;
    }

    /// Serialize the current properties to the configuration manager for saving to file.
    pub fn serialize(&mut self) -> bool {

        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "f_num".into(), Config::serialize_field(&self.f_num)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "ap_efficiency".into(), Config::serialize_field(&self.ap_efficiency)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "coating".into(), Config::serialize_field(&self.coating)
        );

        return true;
    }

    /// Deserialize properties from the configuration manager, loading values from file.
    pub fn deserialize(&mut self) -> bool {

        self.config.deserialize_field(&mut self.f_num, "f_num");
        self.config.deserialize_field(&mut self.ap_efficiency, "ap_efficiency");
        self.config.deserialize_field(&mut self.coating, "coating");

        return true;
    }

    /// Inmutable access to the focal ratio.
    pub fn get_f_num(&self) -> &f64 {
        &self.f_num
    }
    /// Inmutable access to the aperture efficiency.
    pub fn get_ap_efficiency(&self) -> &f64 {
        &self.ap_efficiency
    }
    /// Inmutable access to the coating type.
    pub fn get_coating(&self) -> &String {
        &self.coating
    }
    /// Mutable access to the coating type.
    pub fn get_coating_mut(&mut self) -> &mut String {
        &mut self.coating
    }
}

/// Implementation of methods for InstrumentModel, including initialization and various utility functions.
impl InstrumentModel {
    /// Create a new InstrumentModel instance, loading necessary curves and initializing properties.
    pub fn new() -> InstrumentModel {
        // Load instrument properties from configuration file.
        let mut properties = InstrumentProperties::new("instrument".to_string());
        properties.config.load();
        properties.deserialize();

        // Initialize the InstrumentModel with default values and loaded properties.
        let mut instrument_model = InstrumentModel {
            // Set instrument properties.
            properties,

            // Initialize curves and vectors with default values or capacities.
            blue_ml15: Default::default(),
            blue_nbb: Default::default(),
            red_ml15: Default::default(),
            blue_disp: Vec::with_capacity(TARSIS_SLICES),  
            blue_repx: Vec::with_capacity(TARSIS_SLICES),  
            blue_w2px: Vec::with_capacity(TARSIS_SLICES),  
            blue_px2w: Vec::with_capacity(TARSIS_SLICES),  
            red_disp: Vec::with_capacity(TARSIS_SLICES),  
            red_repx: Vec::with_capacity(TARSIS_SLICES),  
            red_w2px: Vec::with_capacity(TARSIS_SLICES),  
            red_px2w: Vec::with_capacity(TARSIS_SLICES), 
            atten_spectrum: Default::default(),
            current_path: InstrumentArm::BlueArm,  // Blue arm is the default value
        };
        
        // Load transmission curves from data files and scale the X axis from nm to m
        let _ = instrument_model.blue_ml15.load_curve(&DataFileManager::data_file(&"blueTransmission.csv".to_string()));
        instrument_model.blue_ml15.scale_axis(XAxis, 1e-9);

        let _ = instrument_model.blue_nbb.load_curve(&DataFileManager::data_file(&"blueTransmission2.csv".to_string()));
        instrument_model.blue_nbb.scale_axis(XAxis, 1e-9);

        let _ = instrument_model.red_ml15.load_curve(&DataFileManager::data_file(&"redTransmission.csv".to_string()));
        instrument_model.red_ml15.scale_axis(XAxis, 1e-9);

        // Load dispersion and resolution curves for each slice, scaling axes and preparing for use
        // in wavelength to pixel and pixel to wavelength conversions.
        // We also compute the integrals of the dispersion curves to get the wavelength to pixel curves
        // and flip them to get the pixel to wavelength curves.
        for i in 0..TARSIS_SLICES {
            ///////////////////////////// Blue arm ///////////////////////////
            // Spectral dispersion.
            let mut blue_disp_curve: Curve = Default::default();
            let _ = blue_disp_curve.load_slice_dat(&DataFileManager::data_file(&"dispersionBlue.csv".to_string()), i); 
            let _ = blue_disp_curve.extend_right(); // Extend the curve to the right to cover the full range.
            let _ = blue_disp_curve.extend_left();  // Extend the curve to the left to cover the full range.
            let _ = blue_disp_curve.scale_axis(XAxis, 1e-9); // Scale X axis from nm to m.
            let _ = blue_disp_curve.scale_axis(YAxis, 1e-9); // Scale Y axis from nm/px to m/px.
            let _ = blue_disp_curve.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN")); // Invert Y axis to get (m -> px/m).    
            let _ = instrument_model.blue_disp.push(blue_disp_curve.clone());
            // Resolution element (in pixels).
            let mut blue_repx_curve: Curve = Default::default();
            let _ = blue_repx_curve.load_slice_dat(&DataFileManager::data_file(&"pxResolutionBlue.csv".to_string()), i);
            let _ = blue_repx_curve.extend_right(); // Extend the curve to the right to cover the full range.
            let _ = blue_repx_curve.extend_left(); // Extend the curve to the left to cover the full range.
            let _ = blue_repx_curve.scale_axis(XAxis, 1e-9); // Scale X axis from nm to m.
            let _ = blue_repx_curve.invert_axis(YAxis, NotNan::new(STD2FWHM).expect("x should not be NaN")); // Invert Y axis and convert from std to fwhm.
            let _ = instrument_model.blue_repx.push(blue_repx_curve.clone());
            // Wavelength to pixel and pixel to wavelength curves.
            let mut blue_w2px_curve: Curve = Default::default();
            let _ = blue_w2px_curve.assign(&blue_disp_curve.clone()); // Start from the dispersion curve.
            let _ = blue_w2px_curve.integrate(0.0);  // Integrate to get wavelength to pixel.
            let _ = instrument_model.blue_w2px.push(blue_w2px_curve.clone()); 
            // Pixel to wavelength is the inverse of the above.
            let mut blue_px2w_curve: Curve = Default::default();
            let _ = blue_px2w_curve.assign(&blue_w2px_curve.clone()); // Start from the wavelength to pixel curve.
            let _ = blue_px2w_curve.flip(); // Flip to get pixel to wavelength.
            let _ = instrument_model.blue_px2w.push(blue_px2w_curve.clone());
            ////////////////////////////// Red arm ///////////////////////////
            // Spectral dispersion.
            let mut red_disp_curve: Curve = Default::default();
            let _ = red_disp_curve.load_slice_dat(&DataFileManager::data_file(&"dispersionRed.csv".to_string()), i);
            let _ = red_disp_curve.extend_right(); // Extend the curve to the right to cover the full range.
            let _ = red_disp_curve.extend_left(); // Extend the curve to the left to cover the full range.
            let _ = red_disp_curve.scale_axis(XAxis, 1e-9); // Scale X axis from nm to m.
            let _ = red_disp_curve.scale_axis(YAxis, 1e-9); // Scale Y axis from nm/px to m/px.
            let _ = red_disp_curve.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN")); // Invert Y axis to get (m -> px/m).
            let _ = instrument_model.red_disp.push(red_disp_curve.clone());
            // Resolution element (in pixels).
            let mut red_repx_curve: Curve = Default::default();
            let _ = red_repx_curve.load_slice_dat(&DataFileManager::data_file(&"pxResolutionRed.csv".to_string()), i);
            let _ = red_repx_curve.extend_right(); // Extend the curve to the right to cover the full range.
            let _ = red_repx_curve.extend_left(); // Extend the curve to the left to cover the full range.
            let _ = red_repx_curve.scale_axis(XAxis, 1e-9); // Scale X axis from nm to m.
            let _ = red_repx_curve.invert_axis(YAxis, NotNan::new(STD2FWHM).expect("x should not be NaN")); // Invert Y axis and convert from std to fwhm.
            let _ = instrument_model.red_repx.push(red_repx_curve.clone());
            // Wavelength to pixel and pixel to wavelength curves.
            let mut red_w2px_curve: Curve = Default::default();
            let _ = red_w2px_curve.assign(&red_disp_curve.clone()); // Start from the dispersion curve.
            let _ = red_w2px_curve.integrate(0.0); // Integrate to get wavelength to pixel.
            let _ = instrument_model.red_w2px.push(red_w2px_curve.clone());
            // Pixel to wavelength is the inverse of the above.
            let mut red_px2w_curve: Curve = Default::default();
            let _ = red_px2w_curve.assign(&red_w2px_curve.clone()); // Start from the wavelength to pixel curve.
            let _ = red_px2w_curve.flip(); // Flip to get pixel to wavelength.
            let _ = instrument_model.red_px2w.push(red_px2w_curve.clone());
        }
        // Return the fully initialized instrument model.
        return instrument_model;
    }
    /// Inmutable acces to the blue ML15 transmission curve.
    pub fn get_blue_ml15(&self) -> &Curve {
        &self.blue_ml15
    }
    /// Inmutable acces to the blue NBB transmission curve.
    pub fn get_blue_nbb(&self) -> &Curve {
        &self.blue_nbb
    }
    /// Inmutable acces to the red ML15 transmission curve.
    pub fn get_red_ml15(&self) -> &Curve {
        &self.red_ml15
    }
    /// Inmutable access to the dispersion curve for the blue arm.
    pub fn get_blue_disp(&self) -> &Vec<Curve> {
        &self.blue_disp
    }
    /// Inmutable access to the resolution element curve for the blue arm.
    pub fn get_blue_repx(&self) -> &Vec<Curve> {
        &self.blue_repx
    }
    /// Inmutable access to the wavelength to pixel curve for the blue arm.
    pub fn get_blue_w2px(&self) -> &Vec<Curve> {
        &self.blue_w2px
    }
    /// Inmutable access to the pixel to wavelength curve for the blue arm.
    pub fn get_blue_px2w(&self) -> &Vec<Curve> {
        &self.blue_px2w
    }
    /// Inmutable access to the dispersion curve for the red arm.
    pub fn get_red_disp(&self) -> &Vec<Curve> {
        &self.red_disp
    }
    /// Inmutable access to the resolution element curve for the red arm.
    pub fn get_red_repx(&self) -> &Vec<Curve> {
        &self.red_repx
    }
    /// Inmutable access to the wavelength to pixel curve for the red arm.
    pub fn get_red_w2px(&self) -> &Vec<Curve> {
        &self.red_w2px
    }
    /// Inmutable access to the pixel to wavelength curve for the red arm.
    pub fn get_red_px2w(&self) -> &Vec<Curve> {
        &self.red_px2w
    }
    /// Inmutable access to the attenuated input spectrum.
    pub fn get_atten_spectrum(&self) -> &Spectrum {
        &self.atten_spectrum
    }
    /// Inmutable access to the current arm/path of the attenuated spectrum.
    pub fn get_current_path(&self) -> &InstrumentArm {
        &self.current_path
    }
    /// Inmutable acces to the instrument properties.
    pub fn properties(&self) -> &InstrumentProperties {
        &self.properties
    }
    /// Mutable acces to the instrument properties.
    pub fn properties_mut(&mut self) -> &mut InstrumentProperties {
        &mut self.properties
    }

    /// Turns a pixel into wavelength in a given pixel.
    pub fn px_to_wavelength_val(&mut self, arm: InstrumentArm, slice: usize, pixel: NotNan<f64>) -> Result<f64, String> {
        // Check slice bounds.
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        // Select the appropriate arm and get the wavelength for the given pixel.
        match arm {
            InstrumentArm::BlueArm => {
                // Get the pixel to wavelength value of the curve for the specified slice.
                let curve = &self.blue_px2w[slice];
                Ok(curve.get_point(pixel))  
            }
            InstrumentArm::RedArm => {
                // Get the pixel to wavelength value of the curve for the specified slice.
                let curve = &self.red_px2w[slice];
                Ok(curve.get_point(pixel)) 
            }
        }
    }

    /// Turns a pixel into wavelength for all pixels.
    pub fn px_to_wavelength_crv(&mut self, arm: InstrumentArm, slice: usize) -> Result<Curve, String> {
        // Check slice bounds.
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        // Select the appropriate arm and return the full pixel to wavelength curve.
        match arm {
            InstrumentArm::BlueArm => {
                // Get the pixel to wavelength curve for the specified slice.
                let curve = &self.blue_px2w[slice];
                Ok(curve.clone())  
            }
            InstrumentArm::RedArm => {
                // Get the pixel to wavelength curve for the specified slice.
                let curve = &self.blue_px2w[slice];
                Ok(curve.clone()) 
            }
        }
    }

    /// Turns a wavelength into pixel in a given wavelength.
    pub fn wavelength_to_px_val(&mut self, arm: InstrumentArm, slice: usize, lambda: NotNan<f64>) -> Result<f64, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                // Get the wavelength to pixel value of the curve for the specified slice.
                let curve = &self.blue_w2px[slice];
                Ok(curve.get_point(lambda))  
            }
            InstrumentArm::RedArm => {
                // Get the wavelength to pixel value of the curve for the specified slice.
                let curve = &self.red_w2px[slice];
                Ok(curve.get_point(lambda)) 
            }
        }
    }

    /// Turns a wavelength into pixel for all wavelengths.
    pub fn wavelength_to_pix_crv(&mut self, arm: InstrumentArm, slice: usize) -> Result<Curve, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                // Get the wavelength to pixel curve for the specified slice.
                let curve = &self.blue_w2px[slice];
                Ok(curve.clone())  
            }
            InstrumentArm::RedArm => {
                // Get the wavelength to pixel curve for the specified slice.º
                let curve = &self.red_w2px[slice];
                Ok(curve.clone()) 
            }
        }
    }

    /// Set the input (object) spectrum. The input spectrum must be in radiance units,
    /// with a *wavelength* spectral axis* i.e. J / (s * m^2 * sr * m)
    pub fn set_input(&mut self, arm: InstrumentArm, input: Spectrum) -> () {
        // Initialize variables for calculations.
        let light_cone_sr: f64;
        let aperture_ang_radius: f64;
        let total_scale: f64;
        let transmission: Curve;
        let coating = self.properties.coating.clone();

        // Select the appropriate transmission curve based on the arm and coating type.
        match arm {
            InstrumentArm::BlueArm => {
                if coating == "ML15" {
                    transmission = self.blue_ml15.clone();
                } else if coating == "NBB" {
                    transmission = self.blue_nbb.clone();
                } else {
                    println!("{} is not a supported AR coating for the blue arm", coating);
                    return;
                }
            }
            InstrumentArm::RedArm => {
                if coating == "ML15" {
                    transmission = self.red_ml15.clone();
                } else {
                    println!("{} is not a supported AR coating for the blue arm", coating);
                    return;
                }
            }
        }
        // Calculate the light cone solid angle, total scale factor, and apply to the input spectrum.
        aperture_ang_radius = (0.5 / self.properties.f_num).atan();
        light_cone_sr = std::f64::consts::PI * aperture_ang_radius * aperture_ang_radius;
        total_scale = light_cone_sr * self.properties.ap_efficiency;
        self.current_path = arm;
        self.atten_spectrum.from_existing(&input.get_curve(), 1.0);   // Set input radiance
        self.atten_spectrum.scale_axis_factor(YAxis, total_scale);           // To irradiance
        self.atten_spectrum.multiply(&transmission);                  // Attenuate by transmission
    }

    /// Convolve a curve around a given x0 with a Gaussian of given inv_sigma (1/std) and oversampling factor.
    pub fn convolve_around(curve: &Curve, x0: f64, inv_sigma: f64, oversample: u32) -> f64 {
        // Initialize variables for convolution.
        let mut scale = 0.0;
        let half_prec = 0.5 * inv_sigma * inv_sigma;
        let dx = STD2FWHM / (inv_sigma * oversample as f64);
        let mut y = 0.0;

        // Ensure oversample is odd to have a central point.
        let half_width = oversample / 2;

        // Perform convolution by sampling points around x0.
        for i in 0..oversample {
            let offset = (i as i32 - half_width as i32) as f64 * dx;
            let x = x0 + offset;
            let weight = (-half_prec * offset * offset).exp();
            let not_nan_x = NotNan::new(x).expect("x should not be NaN");
            let dy = curve.get_point(not_nan_x) * weight;
            y += dy;
            scale += weight;
        }

        y /= scale;

        return y;
    }

    /// Switch the current arm/path of the attenuated spectrum.
    pub fn switch_current_arm(&mut self, arm: InstrumentArm) {
        self.current_path = arm;
    }

    /// Returns the per-pixel photon flux,in units of in ph / (s m^2)
    pub fn make_pixel_photon_flux(&self, slice: usize) -> Spectrum {
        // Initialize variables for calculations.
        let mut disp_spectrum: Spectrum = Default::default();
        let mut w2px_ptr: Curve = Default::default();
        let mut px2w_ptr: Curve = Default::default();
        let mut res_el_ptr: Curve = Default::default();
        let mut disp_ptr: Curve = Default::default();

        // Assign the appropriate curves based on the current arm/path and slice.
        match self.current_path {
            InstrumentArm::BlueArm => {
                disp_ptr = self.blue_disp[slice].clone();
                w2px_ptr = self.blue_w2px[slice].clone();
                px2w_ptr = self.blue_px2w[slice].clone();
                res_el_ptr = self.blue_repx[slice].clone();
            }
            InstrumentArm::RedArm => {
                disp_ptr = self.red_disp[slice].clone();
                w2px_ptr = self.red_w2px[slice].clone();
                px2w_ptr = self.red_px2w[slice].clone();
                res_el_ptr = self.red_repx[slice].clone();
            }
        }

        // Create references to the curves for easier access.
        let w2px = w2px_ptr;
        let px2w = px2w_ptr;
        let res_el = res_el_ptr;
        let disp = disp_ptr;

        // Prepare the dispersion spectrum by copying the attenuated spectrum and scaling the X axis.   
        disp_spectrum.from_existing(&self.atten_spectrum.get_curve(), 1.0);
        disp_spectrum.scale_axis_curve_diff(XAxis, &w2px, &disp);

        // Initialize the pixel flux spectrum.
        let mut pixel_flux: Spectrum = Default::default();
        
        // Loop over each pixel to compute the photon flux.
        for i in 0..SPECTRAL_PIXEL_LENGTH.round() as i32 {
            // Get the wavelength corresponding to the current pixel.
            let wl = px2w.get_point(i.into());
            // Convert energy flux to photon flux.
            let to_photons = wl / (PLANCK_CONSTANT * SPEED_OF_LIGHT);

            // Convolve the dispersion spectrum around the current wavelength with a Gaussian
            // defined by the resolution element at that wavelength. Add the result to the pixel flux spectrum. 
            if !wl.is_nan() {
                pixel_flux.get_curve_mut().get_map_mut().insert(NotNan::new(wl as f64).expect("x should not be NaN"), 
                Self::convolve_around(&disp_spectrum.get_curve(), i as f64,
                    res_el.get_point(NotNan::new(wl).expect("x should not be NaN")),
                    11) * to_photons);
            }
        }
        
        return pixel_flux;

    }
}


