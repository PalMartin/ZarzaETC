mod cli;

use clap::Parser;
use ordered_float::NotNan;
use rand_distr::uniform;
use rand_distr::Uniform;
use zarza_etc::calculation_worker::CalculationProduct;
use zarza_etc::detector::*;
use zarza_etc::curve::*;
use zarza_etc::curve::CurveAxis::*;
use zarza_etc::spectrum::*;
use zarza_etc::simulation::*;
use zarza_etc::instrument_model::*;
use zarza_etc::helpers::*;
use zarza_etc::instrument_model::InstrumentArm::RedArm;
use zarza_etc::instrument_model::InstrumentArm::BlueArm;
use zarza_etc::calculation_worker::*;
use std::sync::Once;
use std::env;


use cli::{Cli, Commands};


// environments
static INIT: Once = Once::new();

pub fn init_env() {
    INIT.call_once(|| {
        env::set_var("ZARZA_ETC_DATA_DIR", "data"); // así
    });
}

fn spectral_flux_distr(simul: &mut Simulation, arm: String) {

    let cli = Cli::parse();

    match cli.input.spectral.as_str() {
                "template" => {
                    let temp = cli.input.temp.as_ref().expect("Missing --temp");
                    simul.set_input_template(&format!("SEDs/{}", temp));
                }
                "line" => {
                    simul.set_input_line(
                        NotNan::new(cli.input.lambda.expect("Missing --lambda")).unwrap(),
                        cli.input.flux.expect("Missing --flux"),
                        cli.input.fwhm.expect("Missing --fwhm"),
                        cli.input.continuum.expect("Missing --continuum"),
                        &cli.input.line_type.clone().expect("Missing --line-type"),
                        &cli.input.resolution.clone().expect("Missing --resolution"),
                        if arm == "red" {
                            InstrumentArm::RedArm
                        } else {
                            InstrumentArm::BlueArm
                        },
                        cli.input.lambda_min.expect("Missing --lambda-min"),
                        cli.input.lambda_max.expect("Missing --lambda-max"),
                    );
                }

                "user" => { // TODO
                    let file = cli.input.file.as_ref().expect("Missing --file");
                    simul.set_input_user(file);
                }
                _ => panic!("Unknown input type. Must be: 'template', 'line' or 'user'."),
            }
}

fn spatial_flux_distr(simul: &mut Simulation, pixel_position: usize) {

    let cli = Cli::parse();

    match cli.spdtr.spatial.as_str() {
                "infinite" => {
                    simul.set_spatial_distribution(SpatialDistribution::Infinite, pixel_position);
                }
                "uniform" => {
                    simul.set_spatial_distribution(
                        SpatialDistribution::Uniform { 
                            center: cli.spdtr.center.expect("Missing --center"), 
                            radius: cli.spdtr.radius.expect("Missing --radius") 
                        }, pixel_position);
                }
                "sersic" => {
                    simul.set_spatial_distribution(
                        SpatialDistribution::Sersic { 
                            center: cli.spdtr.center.expect("Missing --center"), 
                            r_e: cli.spdtr.r_e.expect("Missing --r_e"), 
                            n: cli.spdtr.n.expect("Missing --n") 
                        }, pixel_position);
                }
                "exponential" => {
                    simul.set_spatial_distribution(
                        SpatialDistribution::Exponential { 
                            center: cli.spdtr.center.expect("Missing --center"), 
                            r_e: cli.spdtr.r_e.expect("Missing --r_e") 
                        }, pixel_position);
                }
                "gaussian" => {
                    simul.set_spatial_distribution(SpatialDistribution::Gaussian { 
                        center: cli.spdtr.center.expect("Missing --center"), 
                        sigma: cli.spdtr.sigma.expect("Missing --sigma")
                    }, pixel_position);
                }
                "point" => {
                    simul.set_spatial_distribution(
                        SpatialDistribution::Point { 
                            x0: cli.spdtr.x0.expect("Missing --x0"), 
                            y0: cli.spdtr.y0 .expect("Missing --y0")
                        }, pixel_position);
                }
                _ => panic!("Unknown input type. Must be 'infinite', 'uniform', 'sersic', 'exponential', 'gaussian' or 'point'."),
            }
}

fn filter_name(filter: &str) -> &str {

    match filter {
                "Rc" => {
                    return "Generic_Cousins_R.csv"
                }
                "U" => {
                    return "Generic_Johnson_UBVRIJHKL.B.csv"
                }
                "B" => {
                    return "Generic_Johnson_UBVRIJHKL.I.csv" 
                }
                "V" => {
                    return "Generic_Johnson_UBVRIJHKL.R.csv"
                }
                "R" => {
                    return "Generic_Johnson_UBVRIJHKL.U.csv"
                }
                "I" => {
                    return "Generic_Johnson_UBVRIJHKL.V.csv"
                }
                "u" => {
                    return "SLOAN_SDSS.g.csv"
                }
                "g" => {
                    return "SLOAN_SDSS.i.csv"
                }
                "r" => {
                    return "SLOAN_SDSS.r.csv"
                }
                "i" => {
                    return "SLOAN_SDSS.u.csv"
                }
                _ => panic!("Unknown filter. Must be 'Rc', 'U', 'B', 'V', 'R', 'u', 'g', 'r', 'i'."),
            }
}

fn make_simulation(slice: usize, pixpos: usize, airmass: f64, moon: f64, ndit:i32, dit: f64, det: String, arm: String, filter: &str, mag: f64) -> Simulation {
    // Initialize simulation parameters
    let mut simul_params = SimulationParams::new();
    // Modify/set slice
    *simul_params.get_slice_mut() = slice;
    // Modify airmass
    *simul_params.get_airmass_mut() = airmass;
    // Modify moon fraction
    *simul_params.get_moon_mut() = moon;
    // Modify DIT
    *simul_params.get_dit_mut() = dit;
    // Modify NDIT
    *simul_params.get_ndit_mut() = ndit as i32;
    
    // Initialize simulation
    let mut simul = Simulation::new();
    // Set detector
    simul.get_det_mut().set_detector(det.to_string());     
    // Set spectral flux distribution
    spectral_flux_distr(&mut simul, arm.clone());
    // Set spatial flux distribution
    spatial_flux_distr(&mut simul, pixpos);
    // Set simulation parameters into the simulation
    simul.set_params(simul_params.clone());
    // Set filter for normalization
    simul.set_filter(&format!("filters/{}", filter_name(filter)));
    // Normalize to filter magnitude
    simul.normalize_to_filter_mag(mag);

    return simul
}
fn main() {

    init_env();

    let cli = Cli::parse();

    match cli.modes {
        Commands::Snr { slice, pixpos, airmass, moon, ndit, dit, det, arm, filter, mag } => {
            
            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm
            if arm == "red" {
                simul.simulate_arm(RedArm);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.")
            }

            // Computation outputs
            let _ = simul.get_det_mut().signal_crv().save_curve("signal_curve.csv");
            let _ = simul.get_det_mut().noise_crv().save_curve("noise_curve.csv");
            let _ = simul.get_det_mut().snr_crv().save_curve("snr_curve.csv");
            println!("The ressults of the computation of the SNR (the signal, noise and SNR curves), have been saved in the archives of names 'signal_curve.csv', 'noise_curve.csv' and 'snr_curve.csv' respectively.")
        }
        Commands::TExp { slice, pixpos, airmass, moon, ndit, dit, det, arm, snr, lambda_ref, filter, mag } => {

            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm
            if arm == "red" {
                simul.simulate_arm(RedArm);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.")
            }

            let exposure_time = simul.texp_from_snr(snr, NotNan::new(lambda_ref).expect("x should not be NaN"));
            println!("The exposure time necessary to achieve a SNR of {:.3} at a wavelength of {:.3}[nm] is: {:.3}[s].", snr, lambda_ref * 1e9, exposure_time);
        }
        Commands::MagLim { slice, pixpos, airmass, moon, ndit, dit, det, arm, lambda_ref, filter, mag } => {

            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm
            let lim_mag: f64;
            if arm == "red" {
                simul.simulate_arm(RedArm);
                lim_mag = simul.lim_mag(NotNan::new(lambda_ref).expect("x should not be NaN"), RedArm);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
                lim_mag = simul.lim_mag(NotNan::new(lambda_ref).expect("x should not be NaN"), BlueArm);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.");
                lim_mag = 0.0;
            }

            println!("The limiting magnitude at a wavelength of {:.3}[nm] for an exposure time of {:.3}[s] is: {:.3}", lambda_ref * 1e9, *simul.get_params().get_ndit() as f64 * simul.get_params().get_dit(), lim_mag);
        } 
        Commands::FluxLim { slice, pixpos, airmass, moon, ndit, dit, det, arm, lambda_ref, filter, mag } => {

            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm
            let lim_flux: f64;
            if arm == "red" {
                simul.simulate_arm(RedArm);
                lim_flux = simul.lim_flux(NotNan::new(lambda_ref).expect("x should not be NaN"), RedArm);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
                lim_flux = simul.lim_flux(NotNan::new(lambda_ref).expect("x should not be NaN"), BlueArm);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.");
                lim_flux = 0.0;
            }

            println!("The limiting flux at a wavelength of {:.3}[nm] for an exposure time of {:.3}[s] is: {:.3}", lambda_ref * 1e9, *simul.get_params().get_ndit() as f64 * simul.get_params().get_dit(), lim_flux);
        }
        Commands::RadVelUnc { slice, pixpos, airmass, moon, ndit, dit, det, arm, obj_size, filter, mag } => {

            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm
            let rad_vel_unc: f64;
            if arm == "red" {
                simul.simulate_arm(RedArm);
                rad_vel_unc = simul.rad_vel_unc(RedArm, obj_size);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
                rad_vel_unc = simul.rad_vel_unc(BlueArm, obj_size);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.");
                rad_vel_unc = 0.0;
            }

            println!("The radial velocity uncertainty for an object that fills {:.3}[px] on the detector is: {:.3}", obj_size, rad_vel_unc);
        }
        Commands::CountSimul { slice, pixpos, airmass, moon, ndit, dit, det, arm, filter, mag } => {

            // Make simulation
            let mut simul = make_simulation(slice, pixpos, airmass, moon, ndit, dit, det, arm.clone(), &filter, mag);
            // Simulate instrument arm and calculation worker
            let mut snr_curve = SNRcurve::new();
            let mut calculation_worker = CalculationWorker::init();
            if arm == "red" {
                simul.simulate_arm(RedArm);
                calculation_worker.simulate(&mut snr_curve, &mut simul, RedArm);
            } else if arm == "blue" {
                simul.simulate_arm(BlueArm);
                calculation_worker.simulate(&mut snr_curve, &mut simul, BlueArm);
            } else {
                eprint!("Unknown arm. Must be 'red' or 'blue'.")
            }
            //let _ = snr_curve.signal.save_curve("calc_work_simulated_signal_curve.csv");
            //let _ = snr_curve.noise.save_curve("calc_work_simulated_noise_curve.csv");
            let _ = snr_curve.counts.save_curve("calc_work_simulated_counts_curve.csv");
            println!("The ressults of the simulated counts with noise has been saved in the archive of name 'calc_work_simulated_counts_curve.csv'.")
        }
    }
}
