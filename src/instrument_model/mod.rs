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

use crate::curve;

use ordered_float::NotNan;
use crate::spectrum::*;
use crate::curve::*;
use curve::CurveAxis::*;
use crate::helpers::*;


pub const CAHA_APERTURE_DIAMETER: f64 = 3.5;      // m
pub const CAHA_FOCAL_LENGTH: f64 = 12.195;        // m
pub const CAHA_APERTURE_AREA: f64 = 0.25 * std::f64::consts::PI * CAHA_APERTURE_DIAMETER * CAHA_APERTURE_DIAMETER;
pub const CAHA_EFFECTIVE_AREA: f64 = 9.093;       // m^2
pub const TARSIS_SLICES: usize = 40;              // Number of TARSIS slices per FOV
pub const SPECTRAL_PIXEL_LENGTH: f64 = 2048.0;    // Guessed from the resolution elements and ranges

#[derive(Clone)]
pub struct InstrumentProperties {
    f_num: f64,
    ap_efficiency: f64,
    coating: String,
}

#[derive(Debug, PartialEq)]
pub enum InstrumentArm {
    BlueArm,
    RedArm,
}

pub struct InstrumentModel {
    properties: InstrumentProperties, 
    blue_ml15: Curve, 
    blue_nbb: Curve, 
    red_ml15: Curve, 

    ///////////////////////////// Blue arm ///////////////////////////
    blue_disp: Vec<Curve>,    // Spectral dispersion
    blue_repx: Vec<Curve>,    // Resolution element (in pixels)
    blue_w2px: Vec<Curve>,    // int of inv of sd
    blue_px2w: Vec<Curve>,    // Inverse of above

    ////////////////////////////// Red arm ///////////////////////////
    red_disp: Vec<Curve>,     // Spectral dispersion
    red_repx: Vec<Curve>,     // Resolution element (in pixels)
    red_w2px: Vec<Curve>,     // int of inv of sd
    red_px2w: Vec<Curve>,     // Inverse of above

    atten_spectrum: Spectrum,           // Attenuated before dispersor
    current_path: InstrumentArm,    // Path of the attenuated spectrum
}

impl InstrumentProperties {
    pub fn new() -> InstrumentProperties {
        InstrumentProperties {
            f_num: CAHA_FOCAL_LENGTH / CAHA_APERTURE_DIAMETER,
            ap_efficiency: CAHA_EFFECTIVE_AREA / CAHA_APERTURE_AREA,
            coating: String::from("ML15"),
        }
    }
    pub fn get_f_num(&self) -> &f64 {
        &self.f_num
    }
    pub fn get_ap_efficiency(&self) -> &f64 {
        &self.ap_efficiency
    }
    pub fn get_coating(&self) -> &String {
        &self.coating
    }
    pub fn get_coating_mut(&mut self) -> &mut String {
        &mut self.coating
    }
}

impl InstrumentModel {
    // Creates an Instrument Model
    pub fn new() -> InstrumentModel {
        let mut instrument_model = InstrumentModel {
            properties: InstrumentProperties::new(),
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
        
        let _ = instrument_model.blue_ml15.load_curve("src/instrument_model/blueTransmission.csv");
        instrument_model.blue_ml15.scale_axis(XAxis, 1e-9);

        let _ = instrument_model.blue_nbb.load_curve("src/instrument_model/blueTransmission2.csv");
        instrument_model.blue_nbb.scale_axis(XAxis, 1e-9);

        let _ = instrument_model.red_ml15.load_curve("src/instrument_model/redTransmission.csv");
        instrument_model.red_ml15.scale_axis(XAxis, 1e-9);

        // Iteration over tarsis slices
        for i in 0..TARSIS_SLICES {
            let mut blue_disp_curve: Curve = Default::default();
            let _ = blue_disp_curve.load_slice_dat("src/instrument_model/dispersionBlue.csv", i);
            let _ = blue_disp_curve.extend_right();
            let _ = blue_disp_curve.extend_left();
            let _ = blue_disp_curve.scale_axis(XAxis, 1e-9);
            let _ = blue_disp_curve.scale_axis(YAxis, 1e-9);
            let _ = blue_disp_curve.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN"));
            let _ = instrument_model.blue_disp.push(blue_disp_curve.clone());

            let mut blue_repx_curve: Curve = Default::default();
            let _ = blue_repx_curve.load_slice_dat("src/instrument_model/pxResolutionBlue.csv", i);
            let _ = blue_repx_curve.extend_right();
            let _ = blue_repx_curve.extend_left();
            let _ = blue_repx_curve.scale_axis(XAxis, 1e-9);
            let _ = blue_repx_curve.invert_axis(YAxis, NotNan::new(STD2FWHM).expect("x should not be NaN"));
            let _ = instrument_model.blue_repx.push(blue_repx_curve.clone());

            let mut blue_w2px_curve: Curve = Default::default();
            let _ = blue_w2px_curve.from_existing(&blue_disp_curve.clone(), 1.0);
            let _ = blue_w2px_curve.integrate(0.0); //?
            let _ = instrument_model.blue_w2px.push(blue_w2px_curve.clone());

            let mut blue_px2w_curve: Curve = Default::default();
            let _ = blue_px2w_curve.from_existing(&blue_w2px_curve.clone(), 1.0);
            let _ = blue_px2w_curve.flip();
            let _ = instrument_model.blue_px2w.push(blue_px2w_curve.clone());


            let mut red_disp_curve: Curve = Default::default();
            let _ = red_disp_curve.load_slice_dat("src/instrument_model/dispersionRed.csv", i);
            let _ = red_disp_curve.extend_right();
            let _ = red_disp_curve.extend_left();
            let _ = red_disp_curve.scale_axis(XAxis, 1e-9);
            let _ = red_disp_curve.scale_axis(YAxis, 1e-9);
            let _ = red_disp_curve.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN"));
            let _ = instrument_model.red_disp.push(red_disp_curve.clone());

            let mut red_repx_curve: Curve = Default::default();
            let _ = red_repx_curve.load_slice_dat("src/instrument_model/pxResolutionRed.csv", i);
            let _ = red_repx_curve.extend_right();
            let _ = red_repx_curve.extend_left();
            let _ = red_repx_curve.scale_axis(XAxis, 1e-9);
            let _ = red_repx_curve.invert_axis(YAxis, NotNan::new(STD2FWHM).expect("x should not be NaN"));
            let _ = instrument_model.red_repx.push(red_repx_curve.clone());

            let mut red_w2px_curve: Curve = Default::default();
            let _ = red_w2px_curve.from_existing(&red_disp_curve.clone(), 1.0);
            let _ = red_w2px_curve.integrate(0.0);
            let _ = instrument_model.red_w2px.push(red_w2px_curve.clone());

            let mut red_px2w_curve: Curve = Default::default();
            let _ = red_px2w_curve.from_existing(&red_w2px_curve.clone(), 1.0);
            let _ = instrument_model.red_px2w.push(red_px2w_curve.clone());
        }
        return instrument_model;
    }
    // Getter functions
    pub fn get_blue_ml15(&self) -> &Curve {
        &self.blue_ml15
    }
    pub fn get_blue_nbb(&self) -> &Curve {
        &self.blue_nbb
    }
    pub fn get_red_ml15(&self) -> &Curve {
        &self.red_ml15
    }
    pub fn get_blue_disp(&self) -> &Vec<Curve> {
        &self.blue_disp
    }
    pub fn get_blue_repx(&self) -> &Vec<Curve> {
        &self.blue_repx
    }
    pub fn get_blue_w2px(&self) -> &Vec<Curve> {
        &self.blue_w2px
    }
    pub fn get_blue_px2w(&self) -> &Vec<Curve> {
        &self.blue_px2w
    }
    pub fn get_red_disp(&self) -> &Vec<Curve> {
        &self.red_disp
    }
    pub fn get_red_repx(&self) -> &Vec<Curve> {
        &self.red_repx
    }
    pub fn get_red_w2px(&self) -> &Vec<Curve> {
        &self.red_w2px
    }
    pub fn get_red_px2w(&self) -> &Vec<Curve> {
        &self.red_px2w
    }
    pub fn get_atten_spectrum(&self) -> &Spectrum {
        &self.atten_spectrum
    }
    pub fn get_current_path(&self) -> &InstrumentArm {
        &self.current_path
    }

    pub fn properties(&self) -> &InstrumentProperties {
        &self.properties
    }
    pub fn properties_mut(&mut self) -> &mut InstrumentProperties {
        &mut self.properties
    }
    // Turns a pixel into lambda in a given pixel 
    pub fn px_to_wavelength_val(&mut self, arm: InstrumentArm, slice: usize, pixel: NotNan<f64>) -> Result<f64, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                let curve = &self.blue_px2w[slice];
                Ok(curve.get_point(pixel))  
            }
            InstrumentArm::RedArm => {
                let curve = &self.blue_px2w[slice];
                Ok(curve.get_point(pixel)) 
            }
        }
    }
    // Turns a pixel into lambda for all pixels
    pub fn px_to_wavelength_crv(&mut self, arm: InstrumentArm, slice: usize) -> Result<Curve, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                let curve = &self.blue_px2w[slice];
                Ok(curve.clone())  
            }
            InstrumentArm::RedArm => {
                let curve = &self.blue_px2w[slice];
                Ok(curve.clone()) 
            }
        }
    }
    // Turns a lambda into pixel in a given pixel 
    pub fn wavelength_to_px_val(&mut self, arm: InstrumentArm, slice: usize, lambda: NotNan<f64>) -> Result<f64, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                let curve = &self.blue_w2px[slice];
                Ok(curve.get_point(lambda))  
            }
            InstrumentArm::RedArm => {
                let curve = &self.red_w2px[slice];
                Ok(curve.get_point(lambda)) 
            }
        }
    }
    // Turns a lambda into pixel for all pixels
    pub fn wavelength_to_pix_crv(&mut self, arm: InstrumentArm, slice: usize) -> Result<Curve, String> {
        if slice >= TARSIS_SLICES {
            return Err(format!("Slice {} out of bounds", slice + 1));
        }
        match arm {
            InstrumentArm::BlueArm => {
                let curve = &self.blue_w2px[slice];
                Ok(curve.clone())  
            }
            InstrumentArm::RedArm => {
                let curve = &self.red_w2px[slice];
                Ok(curve.clone()) 
            }
        }
    }
    // Set input spectrum. The input spectrum must be in radiance units,
    // with a *wavelength* spectral axis* i.e. J / (s * m^2 * sr * m)
    pub fn set_input(&mut self, arm: InstrumentArm, input: Spectrum) -> () {
        let light_cone_sr: f64;
        let aperture_ang_radius: f64;
        let total_scale: f64;
        let transmission: Curve;
        let coating = self.properties.coating.clone();

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
        aperture_ang_radius = (0.5 / self.properties.f_num).atan();
        light_cone_sr = std::f64::consts::PI * aperture_ang_radius * aperture_ang_radius;
        total_scale = light_cone_sr * self.properties.ap_efficiency;
        self.current_path = arm;
        self.atten_spectrum.from_existing(&input.get_curve(), 1.0);           // Set input radiance
        self.atten_spectrum.scale_axis(YAxis, total_scale); // To irradiance
        self.atten_spectrum.multiply(&transmission);         // Attenuate by transmission
    }

    pub fn convolve_around(curve: &Curve, x0: f64, inv_sigma: f64, oversample: u32) -> f64 {
        let mut scale = 0.0;
        let half_prec = 0.5 * inv_sigma * inv_sigma;
        let dx = STD2FWHM / (inv_sigma * oversample as f64);
        let mut y = 0.0;

        let half_width = (oversample as f64 / 2.0).round() as i32;

        for i in -half_width..=half_width {
            let offset = i as f64 * dx;
            let x = x0 + offset;
            let weight = (-half_prec * offset * offset).exp();
            y += curve.get_point(NotNan::new(x).expect("x should not be NaN")) * weight;
            scale += weight;
        }

        y /= scale;

        return y;
    }

    pub fn switch_current_arm(&mut self, arm: InstrumentArm) {
        self.current_path = arm;
    }

    // Returns the per-pixel photon flux,in units of in ph / (s m^2)
    pub fn make_pixel_photon_flux(&self, slice: usize) -> Spectrum {
        let mut disp_spectrum: Spectrum = Default::default();
        let mut w2px_ptr: Curve = Default::default();
        let mut px2w_ptr: Curve = Default::default();
        let mut res_el_ptr: Curve = Default::default();
        let mut disp_ptr: Curve = Default::default();

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

        let w2px = w2px_ptr;
        let px2w = px2w_ptr;
        let res_el = res_el_ptr;
        let disp = disp_ptr;

        disp_spectrum.from_existing(&self.atten_spectrum.get_curve(), 1.0);
        disp_spectrum.scale_axis(XAxis, &w2px);
        disp_spectrum.scale_axis(XAxis, &disp);

        let mut pixel_flux: Spectrum = Default::default();
        
        for i in 0..SPECTRAL_PIXEL_LENGTH.round() as i32 {
            let wl = px2w.get_point(i.into());
            let to_photons = wl / (PLANCK_CONSTANT * SPEED_OF_LIGHT);

            if !wl.is_nan() {
                pixel_flux.get_curve_mut().get_map_mut().insert(NotNan::new(i as f64).expect("x should not be NaN"), 
                Self::convolve_around(&disp_spectrum.get_curve(), i as f64,
                    res_el.get_point(NotNan::new(wl).expect("x should not be NaN")),
                    11) * to_photons);
            }
        }

        return pixel_flux;

    }
}


