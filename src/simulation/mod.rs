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
use std::collections::BTreeMap;
//use ordered_float::Pow;

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

pub struct Simulation {
    input: Spectrum,
    filter: Curve,
    filter_equiv_bw: f64,
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

    pub fn new(filter_path: &String) -> Simulation {
        // FILTER PATH: the path of the desired filter (cousins_r, gjc or sdss):
        // Generic_Cousins_R.csv; Generic_Johnson_UBVRIJHKL.U.csv; Generic_Johnson_UBVRIJHKL.B.csv; 
        // Generic_Johnson_UBVRIJHKL.V.csv; Generic_Johnson_UBVRIJHKL.R.csv; Generic_Johnson_UBVRIJHKL.I.csv; 
        // SLOAN_SDSS.u.csv; SLOAN_SDSS.g.csv; SLOAN_SDSS.r.csv; SLOAN_SDSS.i.csv
        let mut simulation = Simulation {
            input: Spectrum::default(),
            filter: Curve::default(),
            filter_equiv_bw: 0.0,
            sky: Spectrum::default(),
            sky_model: SkyModel::new(),
            tarsis_model: InstrumentModel::new(),
            det: Detector::new(),
            params: SimulationParams::new(),
        };

        let _ = simulation.filter.load_curve(&DataFileManager::data_file(&filter_path.to_string())).unwrap();
        simulation.filter.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.filter.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.filter_equiv_bw = simulation.filter.integral();

        return simulation;
    }

    pub fn set_input_user(&mut self, spec: Spectrum) {
        self.input = spec;
        // multiply by sky transmission

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

        //match spatial_distr {
        //    "infinite" => {
        //        let intensity_grid = infinite_profile(grid_size);
        //    }
        //    "uniform" => {
        //        let intensity_grid = uniform_profile(grid_size, center, radius);
        //    }
        //    "sersic" => {
        //        let intensity_grid = sersic_profile(grid_size, center, r_e, n);
        //    }
        //    "exponential" => {
        //        let intensity_grid = exponential_profile(grid_size, center, r_e);
        //    }
        //    "gaussian" => {
        //        let intensity_grid = gaussian_profile(grid_size, center, sigma);
        //    }
        //    "point" => {
        //        let intensity_grid = point_source(grid_size, x0, y0);
        //    }
        //}
        //for (y, row) in intensity_grid.iter().enumerate() {
        //    for (x, &val) in row.iter().enumerate() {
        //        self.input.get_point(x).multiply(val) ; // esto es lo de la distribución espacial
        //        self.input[x][y] = 
        //    }
        //}
    }

    // FIX UNITS IN SET_INPUT_LINE
    pub fn set_input_line(&mut self, central_wl: NotNan<f64>, int_flux: f64, fwhm: f64, cont_flux: f64, line_type: &str, res: &str, arm: InstrumentArm) {
        // wl in angstrom and input fluxes must be given in phot/s/cm2/angstrom
        let mut line = BTreeMap::new();
        let mut sigma = 0.0;
        if res == String::from("Resolved") {
            sigma = fwhm / STD2FWHM; // FWHM in Ángstrom​​
        } 
        else if res == String::from("Unresolved") {
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
            panic!("Unexpected type format for {}. Expected 'Resolved' or 'Unresolved'.", res)
        }
        let wl_min = central_wl - 5.0 * sigma;
        let wl_max = central_wl + 5.0 * sigma;
        let num_points = 1000; // The unit test dont work if there are less points
    
        match line_type {
            "emission_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = *(wl_min + (i as f64) * *(((wl_max - wl_min) / (num_points as f64 - 1.0))));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = int_flux * (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux + flux);
                    line.insert(central_wl, int_flux + cont_flux); 
                }
                *self.input.get_curve_mut().get_map_mut() = line;
            }
            "absorption_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = *(wl_min + (i as f64) * *(((wl_max - wl_min) / (num_points as f64 - 1.0))));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = -int_flux * (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux - flux);
                    line.insert(central_wl, int_flux + cont_flux);
                }
                *self.input.get_curve_mut().get_map_mut() = line;
            }
            _ => {
                eprintln!("Warning: Unknown line type '{}'. Expected: 'emission_line' or 'absorption_line'.", line_type);
            }
        }
    }

    pub fn normalize_to_r_mag(&mut self, mag_r: f64) { // ONLY FOR WHEN FILTER -> COUSINS_R IS SELECTED
        let mut filtered = Spectrum::default();
        filtered.from_existing(&self.input.get_curve(), 1.0);
        filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());
        //println!("FILTERED: {:?}", filtered.get_curve().get_curve());
        filtered.multiply(&self.filter); 
        
        let desired_sb = surface_brightness_ab2freq_radiance(mag_r);
        //println!("desired_sb {}:", desired_sb);

        let mean_sb = filtered.integral() / self.filter_equiv_bw;
        //println!("cousins_re_quiv_bw {}:", self.cousins_re_quiv_bw);
        //println!("mean_sb {}:", mean_sb);

        self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
    }

    pub fn normalize_to_filter_mag(&mut self, mag: f64, filter: &str) {
        // Normalize to Johnson-Cousins U,B,V,R,I or SDSS u,g,r,i bands
        let mut filtered = Spectrum::default();
        filtered.from_existing(&self.input.get_curve(), 1.0);
        filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); 

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
        self.sky = self.sky_model.make_sky_spectrum(&self.input);
        //println!("INPUT: {:?}", &self.input);

        // Update detector config
        self.det.set_exposure_time(params.dit * params.ndit as f64)
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
        self.tarsis_model.set_input(arm.clone(), self.sky.clone());
        //println!("sky: {:?}", self.sky);
        let flux = self.tarsis_model.make_pixel_photon_flux(self.params.slice);
        //println!("flux: {:?}", flux);
        //self.det.set_pixel_photon_flux(flux, self.params.slice, arm);        
        self.det.set_pixel_photon_flux(flux);      
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

    pub fn texp_from_snr_px(&self, px: NotNan<f64>, arm: InstrumentArm) -> f64 { // CHECK EXPRESSION
        let n = self.det.get_photon_flux_per_pixel().get_point(px) * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_pixel_side() * self.det.get_q_e().get_point(px) / self.det.get_detector().get_gain();
        return self.det.noise_px(px) * self.params.snr / n;
    }

    pub fn texp_from_snr_crv(&mut self, arm: InstrumentArm) -> f64 { // CHECK EXPRESSION
        let mut texp_crv = BTreeMap::new();
        let mut n = 0.0;
        let mut texp_px = 0.0;
        for (x, y) in self.det.noise_crv().get_curve().get_map().iter() {
            n = self.det.get_photon_flux_per_pixel().get_point(*x) * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_pixel_side() * self.det.get_q_e().get_point(*x) / self.det.get_detector().get_gain();
            texp_crv.insert(*x, *y * self.params.snr / n);
        }
        let texp = texp_crv.values().cloned().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(0.0);
        return texp; 
    }

    pub fn lim_flux(&self, px: NotNan<f64>, fwhm: f64) -> f64 { // CHECK FWHM OF THE LINE
        let mut ron2 = self.det.read_out_noise() * self.det.read_out_noise();
        let noise = (ron2 + self.sky.get_point(px) * (self.params.dit * self.params.ndit as f64)).sqrt(); // ASSUMPTION: CONTRIBUTION OF THE EMISSION LINE TO THE NOISE IS NEGLIGIBLE
        let flux = self.params.snr * noise * fwhm / (self.params.dit * self.params.ndit as f64); // in ADU/s
        return flux * self.det.get_detector().get_gain() / self.det.get_q_e().get_point(px) // in fotons/s
    }

    pub fn lim_mag(&self, px: NotNan<f64>, fwhm: f64) -> f64 { // CHECK -> I calculated this from the limiting (line) flux
        return -2.5 * self.lim_flux(px, fwhm).log10() - 48.6; // AB MAGNITUDE
    }

    pub fn rad_vel_unc(&mut self, central_wl: NotNan<f64>, fwhm: f64, px: NotNan<f64>, res: &str, arm: InstrumentArm) -> f64 {
        let mut sigma = 0.0;
        if res == String::from("Resolved") { 
            sigma = STD2FWHM / fwhm; // FWHM in Ángstrom​​
        } 
        else if res == String::from("Unresolved") {
            let mut fwhm_ins = Curve::default();
            match arm {
                InstrumentArm::RedArm => {
                    fwhm_ins = self.tarsis_model.get_red_repx()[self.params.slice].clone();
                }
                InstrumentArm::BlueArm => {
                    fwhm_ins = self.tarsis_model.get_blue_repx()[self.params.slice].clone();
                }
            }
            sigma = STD2FWHM / fwhm_ins.get_point(central_wl); 
        } else {
            panic!("Unexpected type format for {}. Expected 'Resolved' or 'Unresolved'.", res)
        }
        return SPEED_OF_LIGHT * fwhm / *((central_wl * self.det.snr_px(px))); // CHECK IF THIS EXPRESSION IS OK 
    }

    // SPATIAL FLUX DISTRIBUTION

    // ¡¡¡¡¡¡ Do a correct scaling so that the grid multipy the flux at each pixel so that is the correct percentage of the galaxy total flux (integrate the spectral template?)!!!
    // I found that i must multiply each value of the mask by the inverse of the summation of all values of the mask

    fn infinite_profile(grid_size: usize) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![1.0; grid_size]; grid_size];
        return grid;
    }

    fn uniform_profile(grid_size: usize, center: f64, radius: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; grid_size]; grid_size]; 
        for y in 0..grid_size {
            for x in 0..grid_size {
                let r = ((x as f64 - center).powi(2) + (y as f64 - center).powi(2)).sqrt();
                if r <= radius {
                    grid[y][x] = 1.0;
                } else {
                    grid[y][x] = 0.0;
                }
            }
        }
        return grid;
    }

    fn sersic_profile(grid_size: usize, center: f64, r_e: f64, n: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; grid_size]; grid_size]; 
        let b_n = 2.0 * n - (1.0 / 3.0) + (4.0 / (405.0 * n)) + (46.0 / (25515.0 * n.powf(2.0))); // Approximation
        for y in 0..grid_size {
            for x in 0..grid_size {
                let r = ((x as f64 - center).powi(2) + (y as f64 - center).powi(2)).sqrt();
                //grid[y][x] = i_e * ((-b_n * ((r / r_e).powf(1.0 / n) - 1.0)).exp());
                grid[y][x] = (-b_n * ((r / r_e).powf(1.0 / n))).exp(); // NORMALIZED
            }
        }
        return grid;
    }

    fn exponential_profile(grid_size: usize, center: f64, r_e: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; grid_size]; grid_size]; 
        //let b_n = 2.0 * n - (1.0 / 3.0); // Approximation for b_n only for n > 8.
        for y in 0..grid_size {
            for x in 0..grid_size {
                let r = ((x as f64 - center).powf(2.0) + (y as f64 - center).powf(2.0)).sqrt();
                //grid[y][x] = i_e * (((r / r_e).powf(1.0 / n) - 1.0).exp());
                grid[y][x] = (-r / r_e).exp(); // NORMALIZED
            }
        }
        return grid;
    }

    fn gaussian_profile(grid_size: usize, center: f64, sigma: f64) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; grid_size]; grid_size];
        for y in 0..grid_size {
            for x in 0..grid_size {
                let r = ((x as f64 - center).powf(2.0) + (y as f64 - center).powf(2.0)).sqrt();
                //grid[y][x] = i_e * ((- r.powf(2.0) / (2.0 * sigma.powf(2.0))).exp());
                grid[y][x] = - r.powf(2.0) / (2.0 * sigma.powf(2.0)).exp(); // NORMALIZED
            }
        }
        return grid;
    }

    fn point_source(grid_size: usize, x0: usize, y0: usize) -> Vec<Vec<f64>> {
        let mut grid = vec![vec![0.0; grid_size]; grid_size]; 
        for y in 0..grid_size {
            for x in 0..grid_size {
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

