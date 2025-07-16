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
use crate::data_file_manager::*;
use crate::spectrum::SpectrumAxisOperations;
//use std::arch::aarch64::float64x1x3_t;
use std::collections::BTreeMap;
//use ordered_float::Pow;
use std::cmp;

pub const GRID_SIZE: usize = 40;      // m


#[derive(Clone, Debug)]
pub struct SimulationParams {
    prog_name: String,
    blue_detector: String,
    red_detector: String, 
    airmass: f64,
    moon: f64, 
    dit: f64,
    ndit: i32,
    //exposure: f64,
    snr: f64,
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
    pub fn get_airmass_mut(&mut self) -> &mut f64 {
        &mut self.airmass
    }
    pub fn get_moon(&self) -> &f64 {
        &self.moon
    }
    pub fn get_moon_mut(&mut self) -> &mut f64 {
        &mut self.moon
    }
    pub fn get_dit(&self) -> &f64 {
        &self.dit
    }
    pub fn get_dit_mut(&mut self) -> &mut f64 {
        &mut self.dit
    }
    pub fn get_ndit(&self) -> &i32 {
        &self.ndit
    }
    pub fn get_ndit_mut(&mut self) -> &mut i32 {
        &mut self.ndit
    }
    //pub fn get_exposure_time(&self) -> &f64 {
    //    &(self.dit * self.ndit as f64)
    //}
    pub fn get_r_ab_mag(&self) -> &f64 {
        &self.r_ab_mag
    }
    pub fn get_r_ab_mag_mut(&mut self) -> &mut f64 {
        &mut self.r_ab_mag
    }
    pub fn get_slice(&self) -> &usize {
        &self.slice
    }
    pub fn get_slice_mut(&mut self) -> &mut usize {
        &mut self.slice
    }
    pub fn new() -> SimulationParams {
        SimulationParams {
            prog_name: String::from(""),
            blue_detector: String::from("ML15"),
            red_detector: String::from("ML15"),
            airmass: 1.0,
            moon: 0.0,
            dit: 3600.0, 
            ndit: 1, 
            //exposure: 0.0,
            snr: 10.0,
            r_ab_mag: 18.0,
            slice: 20,
        }
    }
}

#[derive(Clone)]
pub struct Simulation {
    input: Spectrum,
    filter: Curve,
    filter_equiv_bw: f64,
    sky: Spectrum,
    obj: Spectrum,
    sky_model: SkyModel,
    tarsis_model: InstrumentModel,
    det: Detector,
    params: SimulationParams,
}

//#[derive(Clone, Debug)]
pub enum SpatialDistribution {
    Infinite,
    Uniform { center: f64, radius: f64 },
    Sersic { center: f64, r_e: f64, n: f64 },
    Exponential { center: f64, r_e: f64 },
    Gaussian { center: f64, sigma: f64 },
    Point { x0: usize, y0: usize },
}


impl Simulation {
    // Getter functions (for unit tests)
    pub fn get_input(&self) -> &Spectrum {
        &self.input
    }
    pub fn get_input_mut(&mut self) -> &mut Spectrum {
        &mut self.input
    }
    pub fn get_filter(&self) -> &Curve {
        &self.filter
    }
    pub fn get_filter_equiv_bw(&self) -> &f64 {
        &self.filter_equiv_bw
    }
    //pub fn get_cousins_r(&self) -> &Curve {
    //     &self.cousins_r
    // }
    // pub fn get_cousins_r_equiv_bw(&self) -> &f64 {
    //     &self.cousins_r_equiv_bw
    // }
    pub fn get_sky(&self) -> &Spectrum {
        &self.sky
    }
    pub fn get_sky_model(&self) -> &SkyModel {
        &self.sky_model
    }
    pub fn get_tarsis_model(&self) -> &InstrumentModel {
        &self.tarsis_model
    }
    pub fn get_tarsis_model_mut(&mut self) -> &mut InstrumentModel {
        &mut self.tarsis_model
    }
    pub fn get_det(&self) -> &Detector {
        &self.det
    }
    pub fn get_det_mut(&mut self) -> &mut Detector {
        &mut self.det
    }
    pub fn get_params(&self) -> &SimulationParams {
        &self.params
    }
    pub fn get_params_mut(&mut self) -> &mut SimulationParams {
        &mut self.params
    }

    pub fn new() -> Simulation {
        let simulation = Simulation {
            input: Spectrum::default(),
            filter: Curve::default(),
            filter_equiv_bw: 0.0,
            sky: Spectrum::default(),
            obj: Spectrum::default(),
            sky_model: SkyModel::new(),
            tarsis_model: InstrumentModel::new(),
            det: Detector::new(),
            params: SimulationParams::new(),
        };

        return simulation;
    }

    pub fn set_filter(&mut self, filter_path: &String) {
        // FILTER PATH: the path of the desired filter (cousins_r, gjc or sdss):
        // Generic_Cousins_R.csv; Generic_Johnson_UBVRIJHKL.U.csv; Generic_Johnson_UBVRIJHKL.B.csv; 
        // Generic_Johnson_UBVRIJHKL.V.csv; Generic_Johnson_UBVRIJHKL.R.csv; Generic_Johnson_UBVRIJHKL.I.csv; 
        // SLOAN_SDSS.u.csv; SLOAN_SDSS.g.csv; SLOAN_SDSS.r.csv; SLOAN_SDSS.i.csv

        let _ = self.filter.load_curve(&DataFileManager::data_file(&filter_path.to_string())).unwrap();
        self.filter.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        self.filter.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        self.filter_equiv_bw = self.filter.integral();
    }

    pub fn set_input_user(&mut self, spec_file_path: &String) { // USER MUST CHECK ID THEIR DATA IS NORMALIZED AND IN THE CORRECT UNITS
        let mut spec = Spectrum::default();
        let _ = spec.load_curve(spec_file_path).unwrap();
        self.input = spec;
    }

    pub fn set_input_template(&mut self, template_path: &String) { // CHECK UNIT CONVERSION!!
        // TEMPLATE PATH: the path of the desired spectral template (galactic (swire), active galaxy (AGN atlas) or stellar (pickles)):
        // Ell2_template.csv; ll5_template.csv; Ell13_template.csv; S0_template.csv; Sa_template.csv; Sb_template.csv;
        // Sc_template.csv; Sd_template.csv; Sdm_template.csv; Spi4_template.csv; liner_template.csv; eyfert1_template.csv;
        // seyfert2_template.csv; qso_template.csv; pickles_uk_O9V.csv; pickles_uk_BOV.csv; pickles_uk_B3V.csv; pickles_uk_A0V.csv;
        // pickles_uk_A2V.csv; pickles_uk_A5V.csv; pickles_uk_F0V.csv; pickles_uk_F2V.csv; pickles_uk_F5V.csv; pickles_uk_G0V.csv;
        // pickles_uk_G2V.csv; pickles_uk_G5V.csv; pickles_uk_K0V.csv; pickles_uk_K2V.csv; pickles_uk_K5V.csv; pickles_uk_K7V.csv;
        // pickles_uk_M0V.csv; pickles_uk_M2V.csv; pickles_uk_M4V.csv; pickles_uk_M5V.csv; 

        let mut spec = Spectrum::default();
        // X axis is Angstrom
        // Y Axis is 1e0 Erg / (s cm^2 A) per 2.7 arcsec diam fiber. I.e.,
        //            1.7466e-1 erg / (s cm^2 A arcsec^2), i.e. 
        //             7.4309394e7 W / (m^2 A sr)  
        let _ = spec.load_curve(&DataFileManager::data_file(&template_path.to_string())).unwrap(); // erg s-1 cm-2 A-1 -> CHECK IF ITS REAL FOR EVERY TEMPLATE
        let _ = spec.scale_axis_factor(YAxis, 7.4309394e7); // To SI units  -> FIX THE ANGULAR UNITS: THIS IS CONSIDERING ALL THE EMISSION OF THE GALAXY IS PUNTUAL
        let _ = spec.scale_axis_factor(XAxis, 1e-10); // Convert angstrom to meters
        self.input = spec; // spectra in radiance units (J s-1 m-2 m-1 sr-1)

    }

    // FIX UNITS IN SET_INPUT_LINE
    // OJO, ARREGLAR: HAY QUE PONER LA OPCIÓN DE A QUE LAMBDA CORTAR LA EMISIÓN DEL CONTINUO (ROYO CONSIDERAR UNA LÍNEA ACOMPAÑADA DEL CONTINUO DE UNA LAMBDA ESPECIFICA HASTA OTRA)
    pub fn set_input_line(&mut self, central_wl: NotNan<f64>, int_flux: f64, fwhm: f64, cont_flux: f64, line_type: &str, res: &str, arm: InstrumentArm, wl_min: f64, wl_max: f64) {
        // wl in angstrom and input fluxes must be given in phot/s/cm2/angstrom, same units as the galaxy templates

        let mut line = BTreeMap::new();
        let mut sigma = 0.0;
        if res == String::from("resolved") {
            sigma = fwhm / STD2FWHM; // FWHM in Ángstrom​​
        } 
        else if res == String::from("unresolved") {
            let mut fwhm_ins = Curve::default();
            match arm {
                InstrumentArm::BlueArm => {
                    fwhm_ins = self.tarsis_model.get_blue_repx()[self.params.slice].clone();
                }
                InstrumentArm::RedArm => { // IMPLEMENT HIGH RESOLUTION DETECTOR
                    fwhm_ins = self.tarsis_model.get_red_repx()[self.params.slice].clone();
                }
            }
            sigma = fwhm_ins.get_point(central_wl) / STD2FWHM;
        } else {
            panic!("Unexpected type format for {}. Expected 'resolved' or 'unresolved'.", res)
        }

        let num_points = 500;
        for i in 0..num_points {
            let wl_val = wl_min + (i as f64) * (wl_max - wl_min) / (num_points as f64 - 1.0);
            let wl = NotNan::new(wl_val).expect("wl should not be NaN");
            line.insert(wl, cont_flux);
        }
    
        match line_type {
            "emission_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = wl_min + (i as f64) * ((wl_max - wl_min) / (num_points as f64 - 1.0));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = int_flux * (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux + flux);
                    line.insert(central_wl, int_flux + cont_flux); 
                }
            }
            "absorption_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = wl_min + (i as f64) * ((wl_max - wl_min) / (num_points as f64 - 1.0));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = -int_flux * (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux - flux);
                    line.insert(central_wl, int_flux + cont_flux);
                }
            }
            _ => {
                eprintln!("Warning: Unknown line type '{}'. Expected: 'emission_line' or 'absorption_line'.", line_type);
            }
        }

        let mut spec = Spectrum::default();
        *spec.get_curve_mut().get_map_mut() = line.clone();
        let _ = spec.scale_axis_factor(YAxis, 7.4309394e7); // To SI units  -> FIX THE ANGULAR UNITS: THIS IS CONSIDERING ALL THE EMISSION OF THE GALAXY IS PUNTUAL
        let _ = spec.scale_axis_factor(XAxis, 1e-10); // Convert angstrom to meters
        self.input = spec; // spectra in radiance units (J s-1 m-2 m-1 sr-1)


    }

    pub fn set_spatial_distribution(&mut self, spatial_distr: SpatialDistribution, pixel_position: usize) {

        let grid: Vec<Vec<f64>> = match spatial_distr {
            SpatialDistribution::Infinite => self.clone().infinite_profile(),
            SpatialDistribution::Uniform { center, radius } => {
                self.clone().uniform_profile(center, radius)
            }
            SpatialDistribution::Sersic { center, r_e, n } => {
                self.clone().sersic_profile(center, r_e, n)
            }
            SpatialDistribution::Exponential { center, r_e } => {
                self.clone().exponential_profile(center, r_e)
            }
            SpatialDistribution::Gaussian { center, sigma } => {
                self.clone().gaussian_profile(center, sigma)
            }
            SpatialDistribution::Point { x0, y0 } => {
                self.clone().point_source(x0, y0)
            }
        };

        self.input.scale_axis_factor(YAxis, grid[pixel_position][*self.get_params().get_slice()]); // NOT SHURE IF THIS IS THE CORRECT ORDER OF COORDENATES (SLICE-PIXEL_POSITION)
        //println!("Spatial Distr Factor {}", grid[pixel_position][*self.get_params().get_slice()])
    }

    // pub fn normalize_to_r_mag(&mut self, mag_r: f64) { // ONLY FOR WHEN FILTER -> COUSINS_R IS SELECTED
    //     let mut filtered = Spectrum::default();
    //     filtered.from_existing(&self.input.get_curve(), 1.0);
    //     filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());
    //     filtered.multiply(&self.filter); 
        
    //     let desired_sb = surface_brightness_ab2freq_radiance(mag_r);

    //     let mean_sb = filtered.integral() / self.filter_equiv_bw;

    //     self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
    // }

    // pub fn normalize_to_filter_mag(&mut self, mag: f64, filter: &str) {
    //     // Normalize to Johnson-Cousins U,B,V,R,I or SDSS u,g,r,i bands
    //     let mut filtered = Spectrum::default();
    //     filtered.from_existing(&self.input.get_curve(), 1.0);
    //     filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); 
    //     filtered.multiply(by_what);

    //     let desired_sb = surface_brightness_ab2freq_radiance(mag);
    //     let mean_sb = filtered.integral() / self.filter_equiv_bw;
    //     self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
    // }

    // IT IS NECESSARY TO PREVIOUSLY USE "SET_FILTER"
    pub fn normalize_to_filter_mag(&mut self, mag: f64) { // ONLY FOR WHEN FILTER -> COUSINS_R IS SELECTED
        let mut filtered = Spectrum::default();
        filtered.from_existing(&self.input.get_curve(), 1.0);
        filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());
        filtered.multiply(&self.filter); 
        
        let desired_sb = surface_brightness_ab2freq_radiance(mag);

        let mean_sb = filtered.integral() / self.filter_equiv_bw;

        self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
    }

    pub fn set_params(&mut self, params: SimulationParams) {
        self.params = params.clone();

        self.sky_model.set_airmass(params.airmass);
        self.sky_model.set_moon(params.moon);

        if !self.sky.get_curve().get_map().is_empty() {
            self.sky.clear();
            self.sky = Spectrum::default();
        }

        // Update sky spectrum
        self.sky = self.sky_model.make_sky_spectrum(&self.input).0;
        self.obj = self.sky_model.make_sky_spectrum(&self.input).1;

        // Update detector config
        self.det.set_exposure_time(params.dit * params.ndit as f64)
    }

    pub fn simulate_arm(&mut self, arm: InstrumentArm) { 
        let tarsis_prop: &mut InstrumentProperties = self.tarsis_model.properties_mut();
        // let mut det_name: &str;

        // match arm {
        //     InstrumentArm::BlueArm => {
        //         det_name = &self.params.blue_detector;
        //     }
        //     InstrumentArm::RedArm => {
        //         det_name = &self.params.red_detector;
        //     }
        // }

        // Set coating
        *tarsis_prop.get_coating_mut() = self.det.get_spec().get_coating().clone();

        self.tarsis_model.set_input(arm.clone(), self.obj.clone());
        let flux_obj = self.tarsis_model.make_pixel_photon_flux(self.params.slice);       
        self.tarsis_model.set_input(arm.clone(), self.sky.clone());
        let flux_sky = self.tarsis_model.make_pixel_photon_flux(self.params.slice);     
        self.det.set_pixel_photon_flux(flux_obj, flux_sky);      
    }

    pub fn read_out_noise(&self) -> f64 {
        let spec = self.det.get_spec();
        if spec == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.det.read_out_noise();
    }

    pub fn gain(&self) -> f64 {
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

    pub fn signal(self, wl: NotNan<f64>) -> f64 {
        return self.det.signal_value(wl);
    }

    pub fn noise(self, wl: NotNan<f64>) -> f64 {
        return self.det.noise_value(wl);
    }

    pub fn electrons(self, wl: NotNan<f64>) -> f64 {
        return self.det.electrons_value(wl);

    }

    // Returns DIT or NDIT depending in what the user providad (first element of the touple), as well as INT (=DIT*NDIT) (second element of the touple)
    pub fn texp_from_snr(&mut self, snr: f64, lambda_ref: NotNan<f64>) -> f64 {//, input: f64) -> (f64, f64) { // CHECK, ITS WRONG

        let inv_gain = 1.0 / self.det.get_detector().get_gain();

        let obj = self.det.get_obj_photon_flux_per_pixel().get_curve().get_point(lambda_ref) * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_pixel_side() * 1.0 * inv_gain;
        let sky = self.det.get_sky_photon_flux_per_pixel().get_curve().get_point(lambda_ref) * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_pixel_side() * 1.0 * inv_gain;

        let a = (obj).powf(2.0);
        let b = snr * snr * ((obj) + (sky) + (self.det.dark_electrons(DETECTOR_TEMPERATURE) / self.det.get_exposure_time()));
        let c = snr * snr * self.read_out_noise() * self.read_out_noise();

        let t = (b + ((b * b) + (4.0 * a * c)).sqrt()) / (2.0 * a);

        return t;
        //return (t / input, t);
    }

    pub fn lim_mag(&mut self, lambda_ref: NotNan<f64>, arm: InstrumentArm) -> f64 { 
        //let px = self.tarsis_model.wavelength_to_px_val(arm, self.params.slice, lambda_ref).unwrap();
        // Total throughput: sky extinction + instrument transmission (telescope + instrument + QE) (CONSIDERED QE ADDED INTO THE TRANSMISSION CURVE)
        let ext_frac = mag2frac(self.get_sky_model().get_sky_ext().get_point(lambda_ref) * self.get_sky_model().get_airmass());
        let trans_tot: f64;
        match arm { // THIS IS CONSIDERING THE QE IS ADDED INTO THE TRANSMISSION CURVE
            InstrumentArm::RedArm => {
                trans_tot = self.get_tarsis_model().get_red_ml15().get_point(lambda_ref);
            },
            InstrumentArm::BlueArm => {
                trans_tot = self.get_tarsis_model().get_blue_nbb().get_point(lambda_ref); // USE ONLY NBB OR CONSIDER ML15?
            },
        }

        let res_in_wl: f64;
        match arm {
            InstrumentArm::BlueArm => {
                let res_in_px = (self.get_tarsis_model().get_blue_repx()[self.params.slice].get_point(NotNan::new(lambda_ref.into_inner()).expect("x should not be NaN"))).ceil();
                res_in_wl =  self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
            InstrumentArm::RedArm => {
                let res_in_px = (self.get_tarsis_model().get_red_repx()[self.params.slice].get_point(NotNan::new(lambda_ref.into_inner()).expect("x should not be NaN"))).ceil();
                res_in_wl = self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
        }

        let area = CAHA_APERTURE_AREA;
        let n = ext_frac * trans_tot;

        let h = PLANCK_CONSTANT;
        let c = SPEED_OF_LIGHT;
        let lambda = lambda_ref.into_inner(); // in meters
        let d_lambda = res_in_wl; // in meters
        let f_0 = 3.631e-23; // AB zeropoint in W·m⁻²·Hz⁻¹
        let f_lambda = f_0 * c / (lambda * lambda); // Convert to W·m⁻²·m⁻¹

        let energy_per_photon = h * c / lambda;
        let total_time = self.params.ndit as f64 * self.params.dit;

        let k = f_lambda * d_lambda * area * total_time * n / energy_per_photon; // in counts (Signal and noise are in couns so this must also be in counts: in n is the QE so by converting to photons we automatically get counts)

        let mag_lim = -2.5 * ((self.params.snr * (self.params.snr * (self.params.snr * self.params.snr + (4.0 * (self.det.noise_value(lambda_ref) * self.det.noise_value(lambda_ref) - self.det.signal_value(lambda_ref)) / k)).sqrt())) / (2.0 * k)).log10();

        return mag_lim;
    }

    pub fn lim_flux(&mut self, lambda_ref: NotNan<f64>, arm: InstrumentArm) -> f64 {  // IN  W·m⁻²·m⁻¹
        // Total throughput: sky extinction + instrument transmission (telescope + instrument + QE) (CONSIDERED QE ADDED INTO THE TRANSMISSION CURVE)
        let ext_frac = mag2frac(self.get_sky_model().get_sky_ext().get_point(lambda_ref) * self.get_sky_model().get_airmass());
        let trans_tot: f64;
        match arm { // THIS IS CONSIDERING THE QE IS ADDED INTO THE TRANSMISSION CURVE
            InstrumentArm::RedArm => {
                trans_tot = self.get_tarsis_model().get_red_ml15().get_point(lambda_ref);
            },
            InstrumentArm::BlueArm => {
                trans_tot = self.get_tarsis_model().get_blue_nbb().get_point(lambda_ref); // USE ONLY NBB OR CONSIDER ML15?
            },
        }

        let res_in_wl: f64;
        match arm {
            InstrumentArm::BlueArm => {
                let res_in_px = (self.get_tarsis_model().get_blue_repx()[self.params.slice].get_point(NotNan::new(lambda_ref.into_inner()).expect("x should not be NaN"))).ceil();
                res_in_wl =  self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
            InstrumentArm::RedArm => {
                let res_in_px = (self.get_tarsis_model().get_red_repx()[self.params.slice].get_point(NotNan::new(lambda_ref.into_inner()).expect("x should not be NaN"))).ceil();
                res_in_wl = self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
        }

        let area = CAHA_APERTURE_AREA;
        let n = ext_frac * trans_tot;

        let h = PLANCK_CONSTANT; // Planck constant in J·s
        let c = SPEED_OF_LIGHT; // Speed of light in m/s
        let lambda = lambda_ref.into_inner(); // in meters
        let d_lambda = res_in_wl; // in meters

        let energy_per_photon = h * c / lambda;
        let total_time = self.params.ndit as f64 * self.params.dit;
        let k = d_lambda * area * n * total_time / energy_per_photon; // counts per W·m⁻²·m⁻¹ (Signal and noise are in couns so this must also be in counts: in n is the QE so by converting to photons we automatically get counts)

        //let k = area * n * self.params.ndit as f64 * self.params.dit;

        let flux_lim = (self.params.snr * (self.params.snr + (self.params.snr * self.params.snr + (4.0 * (self.det.noise_value(lambda_ref) * self.det.noise_value(lambda_ref) - self.det.signal_value(lambda_ref)) / k)).sqrt())) / (2.0 * k);

        return flux_lim;
    }

    pub fn rad_vel_unc(&mut self, arm: InstrumentArm, obj_size: f64) -> f64 { // change the object size specification when the spatial distribution is implemented
        // only works if the input is an emission/absorption line generated with "set_input_line"
        let tot: f64 = self.det.signal_crv().get_curve().get_map().values().sum();
        let line_cntr = self.det.signal_crv().get_curve().get_map().iter().map(|(x, v)| x.into_inner() * (v / tot)).sum();
        //println!("line cnt: {}", line_cntr);
        let res_in_wl: f64;
        match arm {
            InstrumentArm::BlueArm => {
                let res_in_px = (self.get_tarsis_model().get_blue_repx()[self.params.slice].get_point(NotNan::new(line_cntr).expect("x should not be NaN"))).ceil();
                res_in_wl =  self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_blue_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
            InstrumentArm::RedArm => {
                let res_in_px = (self.get_tarsis_model().get_red_repx()[self.params.slice].get_point(NotNan::new(line_cntr).expect("x should not be NaN"))).ceil();
                res_in_wl = self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(res_in_px).expect("x should not be NaN")) - self.get_tarsis_model().get_red_px2w()[self.params.slice].get_point(NotNan::new(0.0).expect("x should not be NaN"));
            }
        }
        //println!("RES: {}", res_in_wl);

        let mut w : Vec<f64> = Vec::new();
        let mut a_0 : Vec<f64> = Vec::new();
        let mut noise : Vec<f64> = Vec::new();
        for i in 0..SPECTRAL_PIXEL_LENGTH.round() as i32 {
            let wl = self.tarsis_model.px_to_wavelength_val(arm, self.params.slice, NotNan::new(i as f64).expect("x should not be NaN")).unwrap();
            w.push(wl * wl * self.det.signal_crv().get_curve().get_diff(NotNan::new(wl).expect("x should not be NaN")) * self.det.signal_crv().get_curve().get_diff(NotNan::new(wl).expect("x should not be NaN")) / (self.det.signal_value(NotNan::new(wl).expect("x should not be NaN")) + self.det.noise_value(NotNan::new(wl).expect("x should not be NaN"))));
            a_0.push(self.det.signal_value(NotNan::new(wl).expect("x should not be NaN")));
            noise.push(self.det.noise_value(NotNan::new(wl).expect("x should not be NaN")));
        }
        let W: f64 = w.into_iter().sum();
        let A_0: f64 = a_0.into_iter().sum();
        let N: f64 = noise.into_iter().sum();
        let Q: f64 = W.sqrt() / A_0.sqrt();

        let dv_rms: f64 = 0.5 * (res_in_wl + cmp::min(obj_size.ceil() as i32 , 4) as f64 + (SPEED_OF_LIGHT / (Q * (A_0 + N).sqrt()))); // 4 is the size of the spaxel side; need the size of the object as input (obj_size)

        return dv_rms;
    }


    // SPATIAL FLUX DISTRIBUTION

    // ¡¡¡¡¡¡ Do a correct scaling so that the grid multipy the flux at each pixel so that is the correct percentage of the galaxy total flux (integrate the spectral template?)!!!
    // I found that i must multiply each value of the mask by the inverse of the summation of all values of the mask

    fn infinite_profile(self) -> Vec<Vec<f64>> {
        let value = 1.0 / (GRID_SIZE * GRID_SIZE) as f64;
        let grid = vec![vec![value; GRID_SIZE]; GRID_SIZE];
        return grid;
    }

    fn uniform_profile(self, center: f64, radius: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; GRID_SIZE]; GRID_SIZE]; 
        let mut total = 0.0;
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                let r = ((x as f64 - center).powi(2) + (y as f64 - center).powi(2)).sqrt();
                if r <= radius {
                    grid[y][x] = 1.0;
                    total += 1.0;
                } else {
                    grid[y][x] = 0.0;
                }
            }
        }
        // normalization
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                grid[y][x] /= total;
            }
        }

        return grid;
    }

    fn sersic_profile(self, center: f64, r_e: f64, n: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; GRID_SIZE]; GRID_SIZE]; 
        let b_n = 2.0 * n - (1.0 / 3.0) + (4.0 / (405.0 * n)) + (46.0 / (25515.0 * n.powf(2.0))); // Approximation
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                let r = ((x as f64 - center).powi(2) + (y as f64 - center).powi(2)).sqrt();
                //grid[y][x] = i_e * ((-b_n * ((r / r_e).powf(1.0 / n) - 1.0)).exp());
                grid[y][x] = (-b_n * ((r / r_e).powf(1.0 / n))).exp(); // NORMALIZED
            }
        }
        // Normalization
        let total: f64 = grid.iter().flatten().sum();
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                grid[y][x] /= total;
            }
        }

        return grid;
    }

    fn exponential_profile(self, center: f64, r_e: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; GRID_SIZE]; GRID_SIZE]; 
        //let b_n = 2.0 * n - (1.0 / 3.0); // Approximation for b_n only for n > 8.
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                let r = ((x as f64 - center).powf(2.0) + (y as f64 - center).powf(2.0)).sqrt();
                //grid[y][x] = i_e * (((r / r_e).powf(1.0 / n) - 1.0).exp());
                grid[y][x] = (-r / r_e).exp(); // NORMALIZED
            }
        }
        // Normalization
        let total: f64 = grid.iter().flatten().sum();
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                grid[y][x] /= total;
            }
        }

        return grid;
    }

    fn gaussian_profile(self, center: f64, sigma: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; GRID_SIZE]; GRID_SIZE];
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                let r = ((x as f64 - center).powf(2.0) + (y as f64 - center).powf(2.0)).sqrt();
                //grid[y][x] = i_e * ((- r.powf(2.0) / (2.0 * sigma.powf(2.0))).exp());
                grid[y][x] = (-(r.powi(2)) / (2.0 * sigma.powi(2))).exp();// NORMALIZED
            }
        }
        // Normalization
        let total: f64 = grid.iter().flatten().sum();
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                grid[y][x] /= total;
            }
        }

        return grid;
    }

    fn point_source(self, x0: usize, y0: usize) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; GRID_SIZE]; GRID_SIZE]; 
        for y in 0..GRID_SIZE {
            for x in 0..GRID_SIZE {
                if x == x0 && y == y0 {
                    grid[y][x] = 1.0;
                } else {
                    grid[y][x] = 0.0;
                }
            }
        }
        return grid;
    }
}

