use ordered_float::NotNan;
use crate::detector::*;
use crate::curve::*;
use crate::curve::CurveAxis::*;
use crate::spectrum::*;
use crate::simulation::*;
use crate::instrument_model::*;
use crate::helpers::*;
use crate::instrument_model::InstrumentArm::RedArm;
use crate::instrument_model::InstrumentArm::BlueArm;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_simulation_params() {
        let simul_params = SimulationParams::new();
        //assert_eq!(simul_params.get_blue_detector(), "ML15");
    }
    
    #[test]
    fn test_new_simulation() {
        let simul = Simulation::new();
        //assert!(simul.get_input().get_curve().get_curve().is_empty());
        //assert!(simul.get_cousins_r().clone().get_curve().is_empty());
        //assert_eq!(*simul.get_cousins_re_quiv_bw(), 0.0);
        //assert!(simul.get_sky().get_curve().get_curve().is_empty());
        //assert_eq!(simul.get_params().get_prog_name(), "");
        //assert_eq!(simul.get_params().get_blue_detector(), "ML15");
        //assert_eq!(simul.get_params().get_red_detector(), "ML15");
        //assert_eq!(*simul.get_params().get_airmass(), 1.0);
        //assert_eq!(*simul.get_params().get_moon(), 0.0);
        //assert_eq!(*simul.get_params().get_exposure(), 3600.0);
        //assert_eq!(*simul.get_params().get_r_ab_mag(), 18.0);
        //assert_eq!(*simul.get_params().get_slice(), 20);
    }

    #[test]
    fn test_set_input() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //assert!(simul.get_input().get_curve().get_point(NotNan::new(2.0).unwrap()) - spectrum.get_curve().get_point(NotNan::new(2.0).unwrap()) > -spectrum.get_curve().get_point(NotNan::new(2.0).unwrap()));
        //assert!(simul.get_input().get_curve().get_point(NotNan::new(2.0).unwrap()) - spectrum.get_curve().get_point(NotNan::new(2.0).unwrap()) < 10e-4);
    }

    #[test]
    fn test_set_input_line() {
        //Halpha emission
        //let mut simul = Simulation::new();
        //simul.set_input_line(NotNan::new(656.3).unwrap(), 10e-13, 0.1, 1e-16, &String::from("emission_line"), &String::from("Resolved"), RedArm); 
        //let _ = simul.get_input().get_curve().save_curve("test-data/generated_halpha_emission_line.csv");
        //assert!(!simul.get_input().get_curve().get_map().is_empty());

        //NaDI absorption
        //let mut simul = Simulation::new();
        //simul.set_input_line(NotNan::new(589.59).unwrap(), 10e-15, 0.0196, 1e-16, &String::from("absorption_line"), &String::from("Resolved"), RedArm); 
        //let _ = simul.get_input().get_curve().save_curve("test-data/generated_NaD1_absorption_line.csv");
        //assert!(!simul.get_input().get_curve().get_map().is_empty());
    }

    #[test]
    fn test_normalize_to_r_mag() { // FUNCIONA MUY RARO -> FUNCIONA UNA VEZ Y LUEJA DE IR, VUELVE A FUNCIONAR SI ALADO O QUITO get_curve_mut() EN LA LINEA 115 DE SIMULATION
        let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //spectrum.scale_axis(XAxis, 1e-9); // X axis was in nm
        //simul.set_input_line(NotNan::new(656.3).unwrap(), 10e-13, 0.1, 1e-16, &String::from("emission_line"), &String::from("Resolved"), RedArm); 
        //simul.get_input_mut().scale_axis(XAxis, 1e-9); // X axis was in nm
        //println!("Debug spectrum:");
        //spectrum.get_curve().debug();

        //simul.set_input(spectrum.clone());
        //println!("Debug input:");
        //simul.get_input().get_curve().debug();

       //simul.normalize_to_r_mag(16.0);
        //println!("Debug input after norm:");
        //simul.get_input().get_curve().debug();

        //assert!(!simul.get_input().get_curve().get_map().is_empty());
    }

    #[test]
    fn test_set_params() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //let params = SimulationParams::new();
        //simul.set_params(params.clone());

        //assert_eq!(simul.get_params().get_prog_name(), params.get_prog_name());
        //assert_eq!(simul.get_params().get_blue_detector(), params.get_blue_detector());
        //assert_eq!(simul.get_params().get_red_detector(), params.get_red_detector());
        //assert_eq!(simul.get_params().get_airmass(), params.get_airmass());
        //assert_eq!(simul.get_params().get_moon(), params.get_moon());
        //assert_eq!(simul.get_params().get_exposure(), params.get_exposure());
        //assert_eq!(simul.get_params().get_r_ab_mag(), params.get_r_ab_mag());
        //assert_eq!(simul.get_params().get_slice(), params.get_slice());
    }

    #[test]
    fn test_simulate_arm() { // FIX 
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        //simul.set_input(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum.clone());
        //simul.get_tarsis_model_mut().set_input(InstrumentArm::BlueArm, spectrum);
        //simul.simulate_arm(InstrumentArm::BlueArm);
        //assert!(!simul.get_det().get_photon_flux_per_pixel().get_curve().get_curve().is_empty());
    }

    #[test]
    fn test_signal() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum);
        //let px = NotNan::new(10.0).unwrap();
        //let signal = simul.signal(px);
        //assert!(signal >= 0.0);
    }

    #[test]
    fn test_noise() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum);
        //let px = NotNan::new(10.0).unwrap();
        //let noise = simul.noise(px);
        //assert!(noise >= 0.0);
    }

    #[test]
    fn test_electrons() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum);
        //let px = NotNan::new(10.0).unwrap();
        //let electrons = simul.electrons(px);
        //assert!(electrons >= 0.0);
    }

    #[test]
    fn test_read_out_noise() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum);
        //let ron = simul.read_out_noise();
        //assert!(ron >= 0.0);
    }

    #[test]
    fn test_gain() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //simul.get_det_mut().set_detector("CCD231-84-0-S77");
        //simul.get_det_mut().set_pixel_photon_flux(spectrum);
        //let gain = simul.gain();
        //assert!(gain >= 0.0);
    }

    #[test]
    fn test_px_to_wavelength() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //let px = NotNan::new(10.0).unwrap();
        //let wl = simul.px_to_wavelength(InstrumentArm::BlueArm, px);
        //assert!(wl > 0.0);
    }

    #[test]
    fn test_wl_to_pixel_curve() {
        //let mut simul = Simulation::new();
        //let mut spectrum = Spectrum::default();
        //let _ = spectrum.get_curve_mut().load_curve("test-data/flat.csv").unwrap();
        //simul.set_input_user(spectrum.clone());
        //let crv = simul.wl_to_pixel_curve(InstrumentArm::BlueArm);
        //assert!(!crv.get_map().is_empty());
    }
 }

