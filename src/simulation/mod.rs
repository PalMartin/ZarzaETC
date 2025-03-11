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
use std::collections::BTreeMap;
use ordered_float::Pow;

#[derive(Clone, Debug)]
pub struct SimulationParams {
    prog_name: String,
    blue_detector: String,
    red_detector: String, 
    airmass: f64,
    moon: f64, 
    exposure: f64,
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
    pub fn get_exposure(&self) -> &f64 {
        &self.exposure
    }
    pub fn get_exposure_mut(&mut self) -> &mut f64 {
        &mut self.exposure
    }
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
            exposure: 3600.0,
            snr: 10.0,
            r_ab_mag: 18.0,
            slice: 20,
        }
    }
}

pub struct Simulation {
    input: Spectrum,
    cousins_r: Curve,
    cousins_r_equiv_bw: f64,
    gjc_u: Curve,
    gjc_u_equiv_bw: f64,
    gjc_b: Curve,
    gjc_b_equiv_bw: f64,
    gjc_v: Curve,
    gjc_v_equiv_bw: f64,
    gjc_r: Curve,
    gjc_r_equiv_bw: f64,
    gjc_i: Curve,
    gjc_i_equiv_bw: f64,
    sdss_u: Curve,
    sdss_u_equiv_bw: f64,
    sdss_g: Curve,
    sdss_g_equiv_bw: f64,
    sdss_r: Curve,
    sdss_r_equiv_bw: f64,
    sdss_i: Curve,
    sdss_i_equiv_bw: f64,
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
    pub fn get_cousins_r(&self) -> &Curve {
        &self.cousins_r
    }
    pub fn get_cousins_r_equiv_bw(&self) -> &f64 {
        &self.cousins_r_equiv_bw
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

    pub fn new() -> Simulation {
        let mut simulation = Simulation {
            input: Spectrum::default(),
            cousins_r: Curve::default(),
            cousins_r_equiv_bw: 0.0,
            gjc_u: Curve::default(),
            gjc_u_equiv_bw: 0.0,
            gjc_b: Curve::default(),
            gjc_b_equiv_bw: 0.0,
            gjc_v: Curve::default(),
            gjc_v_equiv_bw: 0.0,
            gjc_r: Curve::default(),
            gjc_r_equiv_bw: 0.0,
            gjc_i: Curve::default(),
            gjc_i_equiv_bw: 0.0,
            sdss_u: Curve::default(),
            sdss_u_equiv_bw: 0.0,
            sdss_g: Curve::default(),
            sdss_g_equiv_bw: 0.0,
            sdss_r: Curve::default(),
            sdss_r_equiv_bw: 0.0,
            sdss_i: Curve::default(),
            sdss_i_equiv_bw: 0.0,
            sky: Spectrum::default(),
            sky_model: SkyModel::new(),
            tarsis_model: InstrumentModel::new(),
            det: Detector::new(),
            params: SimulationParams::new(),
        };

        let _ = simulation.cousins_r.load_curve("src/simulation/filters/Generic_Cousins_R.csv").unwrap();
        simulation.cousins_r.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.cousins_r.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.cousins_r_equiv_bw = simulation.cousins_r.integral();

        let _ = simulation.gjc_u.load_curve("src/simulation/filters/Generic_Johnson_UBVRIJHKL.U.csv").unwrap();
        simulation.gjc_u.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.gjc_u.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.gjc_u_equiv_bw = simulation.gjc_u.integral();

        let _ = simulation.gjc_b.load_curve("src/simulation/filters/Generic_Johnson_UBVRIJHKL.B.csv").unwrap();
        simulation.gjc_b.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.gjc_b.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.gjc_b_equiv_bw = simulation.gjc_b.integral();

        let _ = simulation.gjc_v.load_curve("src/simulation/filters/Generic_Johnson_UBVRIJHKL.V.csv").unwrap();
        simulation.gjc_v.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.gjc_v.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.gjc_v_equiv_bw = simulation.gjc_v.integral();

        let _ = simulation.gjc_r.load_curve("src/simulation/filters/Generic_Johnson_UBVRIJHKL.R.csv").unwrap();
        simulation.gjc_r.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.gjc_r.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.gjc_r_equiv_bw = simulation.cousins_r.integral();

        let _ = simulation.gjc_i.load_curve("src/simulation/filters/Generic_Johnson_UBVRIJHKL.I.csv").unwrap();
        simulation.gjc_i.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.gjc_i.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.gjc_i_equiv_bw = simulation.gjc_i.integral();

        let _ = simulation.sdss_u.load_curve("src/simulation/filters/SLOAN_SDSS.u.csv").unwrap();
        simulation.sdss_u.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.sdss_u.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.sdss_u_equiv_bw = simulation.sdss_u.integral();

        let _ = simulation.sdss_g.load_curve("src/simulation/filters/SLOAN_SDSS.g.csv").unwrap();
        simulation.sdss_g.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.sdss_g.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.sdss_g_equiv_bw = simulation.sdss_g.integral();

        let _ = simulation.sdss_r.load_curve("src/simulation/filters/SLOAN_SDSS.r.csv").unwrap();
        simulation.sdss_r.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.sdss_r.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.sdss_r_equiv_bw = simulation.sdss_r.integral();

        let _ = simulation.sdss_i.load_curve("src/simulation/filters/SLOAN_SDSS.i.csv").unwrap();
        simulation.sdss_i.scale_axis(XAxis, 1e-10); // X axis was in angstrom
        simulation.sdss_i.invert_axis(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap()); // To frequency
        simulation.sdss_i_equiv_bw = simulation.sdss_i.integral();

        return simulation;
    }

    pub fn set_input_user(&mut self, spec: Spectrum) {
        self.input = spec;
    }

    pub fn set_input_template(&mut self, template: &str) { // CHECK UNITS
        match template {
            "Ell2" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Ell2_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Ell5" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Ell5_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Ell13" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Ell13_template_norm.csv").unwrap();
                self.input = spec;
            }
            "S0" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/S0_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Sa" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Sa_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Sb" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Sb_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Sc" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Sc_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Sd" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Sd_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Sdm" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Sdm_template_norm.csv").unwrap();
                self.input = spec;
            }
            "Spi4" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/Spi4_template_norm.csv").unwrap();
                self.input = spec;
            }
            "LINER" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/liner_template.csv").unwrap();
                self.input = spec;
            }
            "Sy1" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/seyfert1_template.csv").unwrap();
                self.input = spec;
            }
            "Sy2" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/seyfert2_template.csv").unwrap();
                self.input = spec;
            }
            "QSO" => {
                let mut spec = Spectrum::default();
                let _ = spec.load_curve("src/simulation/SEDs/qso_template.csv").unwrap();
                self.input = spec;
            }
            _ => {
                eprintln!("Warning: Unknown template '{}'. Expected: 'Ell2', 'Ell5', 'Ell13', 'S0', 'Sa', 'Sb', 'Sc', 'Sd', 'Sdm', 'Spi4', 'LINER', 'Sy1', 'Sy2', or 'QSO'.", template);
            }
            // ADD THE STELLAR LIBRARIES
        }
    }


    pub fn set_input_line(&mut self, central_wl: NotNan<f64>, int_flux: f64, fwhm: f64, cont_flux: f64, line_type: &str, res: &str, arm: InstrumentArm) {
        // wl in  nm anf lux in ???
        let mut line = BTreeMap::new();
        let mut sigma = 0.0;
        if res == String::from("Resolved") {
            sigma = fwhm / STD2FWHM; // FWHM in Ángstrom​​
        } 
        else if res == String::from("Unresolved") {
            let mut fwhm_ins = Curve::default();
            match arm {
                InstrumentArm::BlueArm => {
                    let _ = fwhm_ins.load_curve("src/simulation/fwhm_blue.csv").unwrap(); 
                }
                InstrumentArm::RedArm => { // IMPLEMENT HIGH RESOLUTION DETECTOR
                    let _ = fwhm_ins.load_curve("src/simulation/fwhm_red.csv").unwrap();
                }
            }
            sigma = fwhm_ins.get_point(central_wl) / STD2FWHM;
        } else {
            panic!("Unexpected type format for {}. Expected 'Resolved' or 'Unresolved'.", res)
        }
        let wl_min = central_wl - 5.0 * sigma;
        let wl_max = central_wl + 5.0 * sigma;
        let num_points = 100;
    
        match line_type {
            "emission_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = *(wl_min + (i as f64) * *(((wl_max - wl_min) / (num_points as f64 - 1.0))));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = (int_flux / (sigma * (2.0 * std::f64::consts::PI).sqrt())) *
                        (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux + flux);
                }
                *self.input.get_curve_mut().get_map_mut() = line;
            }
            "absorption_line" => {
                for i in 0..num_points {
                    let wl_val: f64 = *(wl_min + (i as f64) * *(((wl_max - wl_min) / (num_points as f64 - 1.0))));  
                    let wl = NotNan::new(wl_val).expect("x should not be NaN");
    
                    let flux = -(int_flux / (sigma * (2.0 * std::f64::consts::PI).sqrt())) *
                        (-(f64::from(wl) - f64::from(central_wl)).powf(2.0) / (2.0 * sigma.powf(2.0))).exp();
                    line.insert(wl, cont_flux - flux);
                }
                *self.input.get_curve_mut().get_map_mut() = line;
            }
            _ => {
                eprintln!("Warning: Unknown line type '{}'. Expected: 'emission_line' or 'absorption_line'.", line_type);
            }
        }
    }

    pub fn normalize_to_r_mag(&mut self, mag_r: f64) {
        let mut filtered = Spectrum::default();
        filtered.get_curve_mut().from_existing(&self.input.get_curve(), 1.0);
        filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());
        //println!("FILTERED: {:?}", filtered.get_curve().get_curve());
        filtered.get_curve_mut().multiply(&self.cousins_r); 
        
        let desired_sb = surface_brightness_ab2freq_radiance(mag_r);
        //println!("desired_sb {}:", desired_sb);

        let mean_sb = filtered.integral() / self.cousins_r_equiv_bw;
        //println!("cousins_re_quiv_bw {}:", self.cousins_re_quiv_bw);
        //println!("mean_sb {}:", mean_sb);

        self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
    }

    pub fn normalize_to_filter_mag(&mut self, mag: f64, filter: &str) {
        // Normalize to Johnson-Cousins U,B,V,R,I or SDSS u,g,r,i bands
        let mut filtered = Spectrum::default();
        filtered.get_curve_mut().from_existing(&self.input.get_curve(), 1.0);
        filtered.invert_axis_spec(XAxis, NotNan::new(SPEED_OF_LIGHT).unwrap());  
        match filter {
            "J-C_U" => {  
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.gjc_u_equiv_bw;
            }
            "J-C_B" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.gjc_b_equiv_bw;
            }

            "J-C_V" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.gjc_v_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "J-C_R" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.gjc_r_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "J-C_I" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.gjc_i_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "SDSS_u" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.sdss_u_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "SDSS_g" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.sdss_g_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "SDSS_r" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.sdss_r_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            }

            "SDSS_i" => {
                let desired_sb = surface_brightness_ab2freq_radiance(mag);
                let mean_sb = filtered.integral() / self.sdss_i_equiv_bw;
                self.input.scale_axis_factor(YAxis, desired_sb / mean_sb);
            },
            _ => {
            eprintln!("Warning: Unknown filter '{}'. Expected: 'J-C_U', 'J-C_B', 'J-C_V', 'J-C_R', 'J-C_I', 'SDSS_u', 'SDSS_g', 'SDSS_r' or 'SDSS_i'.", filter)
            }
        }
    }

    pub fn set_params(&mut self, params: SimulationParams) {
        self.params = params.clone();

        self.sky_model.set_airmass(params.airmass);
        self.sky_model.set_moon(params.moon);

        if !self.sky.get_curve().get_map().is_empty() {
            self.sky.get_curve_mut().clear();
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

    pub fn signal(&self, px: NotNan<f64>) -> f64 {
        let spec = self.det.get_spec();
        if spec == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.det.signal_px(px);
    }

    pub fn noise(&self, px: NotNan<f64>) -> f64 {
        let spec = self.det.get_spec();
        if spec == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.det.noise(px);
    }

    pub fn electrons(&self, px: NotNan<f64>) -> f64 {
        let spec = self.det.get_spec();
        if spec == DetectorSpec::default() {
            panic!("Detector is not set.");
        }

        return self.det.electrons_px(px);
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

    pub fn snr_from_texp_px(&self, px: NotNan<f64>) -> f64 {
        return self.signal(px) / self.noise(px);
    }

    pub fn texp_from_snr_px(&self, px: NotNan<f64>) -> f64 { // CHECK EXPRESSION
        let n = self.det.get_photon_flux_per_pixel().get_curve().get_point(px) * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_pixel_side() * self.det.get_detector().get_q_e() / self.det.get_detector().get_gain();
        return self.noise(px) * self.params.snr / n;
    }

    pub fn lim_flux(&self, px: NotNan<f64>, fwhm: f64) -> f64 { // CHECK FWHM OF THE LINE
        let mut ron2 = self.det.read_out_noise() * self.det.read_out_noise();
        let noise = (ron2 + self.sky.get_curve().get_point(px) * self.params.exposure).sqrt(); // ASSUMPTION: CONTRIBUTION OF THE EMISSION LINE TO THE NOISE IS NEGLIGIBLE
        let flux = self.params.snr * noise * fwhm / self.params.exposure; // in ADU/s
        return flux * self.det.get_detector().get_gain() / self.det.get_detector().get_q_e() // in fotons/s
    }

    pub fn lim_mag(&self, px: NotNan<f64>, fwhm: f64) -> f64 { // CHECK -> I calculated this from the limiting (line) flux
        return -2.5 * self.lim_flux(px, fwhm).log10() - 48.6; // AB MAGNITUDE
    }

    pub fn rad_vel_unc(&mut self, central_wl: NotNan<f64>, fwhm: f64, px: NotNan<f64>, res: &str, arm: InstrumentArm) -> f64 {
        let mut sigma = 0.0;
        if res == String::from("Resolved") { 
            sigma = fwhm / STD2FWHM; // FWHM in Ángstrom​​
        } 
        else if res == String::from("Unresolved") {
            let mut fwhm_ins = Curve::default();
            match arm {
                InstrumentArm::BlueArm => {
                    let _ = fwhm_ins.load_curve("src/simulation/fwhm_blue.csv").unwrap();
                }
                InstrumentArm::RedArm => { // IMPLEMENT HIGH RESOLUTION DETECTOR
                    let _ = fwhm_ins.load_curve("src/simulation/fwhm_red.csv").unwrap();
                }
            }
            sigma = fwhm_ins.get_point(central_wl) / STD2FWHM; 
        } else {
            panic!("Unexpected type format for {}. Expected 'Resolved' or 'Unresolved'.", res)
        }
        return SPEED_OF_LIGHT * fwhm / *((central_wl * self.snr_from_texp_px(px))); // CHECK IF THIS EXPRESSION IS REAL
    }
}

