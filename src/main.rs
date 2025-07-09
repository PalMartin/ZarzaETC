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

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {

    /// Slice number
    #[arg(short, long, default_value_t = 20)]
    slice: usize,

    /// Number of exposures
    #[arg(short, long, default_value_t = 3)]
    ndit: i32,

    /// Time per exposure
    #[arg(short, long, default_value_t = 1200.0)]
    dit: f64,

    /// Select detector for the simulation
    #[arg(short, long)]
    choose_det: String,

    /// Select instrument arm for the simulation
    #[arg(short, long, default_value_t = 1)]
    arm: i32,

    /// Object spectral template
    #[arg(short, long)]
    template: String,
    
}

// environments
static INIT: Once = Once::new();

pub fn init_env() {
    INIT.call_once(|| {
        env::set_var("ZARZA_ETC_DATA_DIR", "data"); // así
    });
}

fn main() {

    init_env();

    let args = Args::parse();

    // Initialize simulation parameters
    let mut simul_params = SimulationParams::new();
    // Modify/set slice
    *simul_params.get_slice_mut() = args.slice;
    // Modify DIT
    *simul_params.get_dit_mut() = args.dit;
    // Modify NDIT
    *simul_params.get_ndit_mut() = args.ndit as i32;

    // Initialize simulation
    let mut simul = Simulation::new();
    // Set detector
    simul.get_det_mut().set_detector(args.choose_det.to_string());
    // Set input object spectrum
    simul.set_input_template(&String::from("SEDs/".to_string() + &args.template));
    //simul.set_input_line(NotNan::new(6563.0).unwrap(), 9.8e-14, 0.1, 2.8e-14, &String::from("emission_line"), &String::from("Resolved"), RedArm, 6500.0, 6610.0);
    // Set spatial distribution
    simul.set_spatial_distribution(SpatialDistribution::Infinite, 20);
    // for other with different inputs it works like this: simul.set_spatial_distribution(SpatialDistribution::Uniform{center: 40.0, radius: 15.0}, pixel_position);

    // Set simulation parameters into the simulation
    simul.set_params(simul_params.clone());
    // Simulate instrument arm
    if args.arm == 1 {
        simul.simulate_arm(RedArm);
    } else {
        simul.simulate_arm(BlueArm);
    }

    let _ = simul.get_det_mut().signal_crv().save_curve("simulated_signal_curve.csv");
    let _ = simul.get_det_mut().noise_crv().save_curve("simulated_noise_curve.csv");
    let _ = simul.get_det_mut().snr_crv().save_curve("simulated_snr_curve.csv");

    println!("Exposure time for NDIT = 3: {:?}", simul.texp_from_snr(10.0, NotNan::new(6.563e-7).expect("x should not be NaN"), 3.0));
    //println!("Limiting magnitude: {:.3}", simul.lim_mag(NotNan::new(6.5e-7).expect("x should not be NaN"), RedArm));
    println!("Limiting magnitude: {:.3}", simul.lim_mag(NotNan::new(6.563e-7).expect("x should not be NaN"), RedArm));
    println!("Limiting flux: {:.3e} [W m-2 m-1]", simul.lim_flux(NotNan::new(6.563e-7).expect("x should not be NaN"), RedArm));
    println!("Radial velocity uncertainty: {:.3}, [km s] (check units)", simul.rad_vel_unc(RedArm, 150.0) / 1000.0);

    //// SIMULATION WITH CALCULATION WORKER TO ADD RANDOM NOISE TO THE COUNTS ////

    // Initialize simulation parameters
    let mut simul_params = SimulationParams::new();
    // Modify/set slice
    *simul_params.get_slice_mut() = args.slice;
    // Modify DIT
    *simul_params.get_dit_mut() = args.dit;
    // Modify NDIT
    *simul_params.get_ndit_mut() = args.ndit as i32;

    // Initialize simulation
    let mut simul_2 = Simulation::new();
    // Set detector
    simul_2.get_det_mut().set_detector(args.choose_det.to_string());
    // Set input object spectrum
    simul_2.set_input_template(&String::from("SEDs/".to_string() + &args.template));
    //simul_2.set_input_line(NotNan::new(6563.0).unwrap(), 9.8e-14, 0.1, 2.8e-14, &String::from("emission_line"), &String::from("Resolved"), RedArm, 6500.0, 6610.0);
    // Set simulation parameters into the simulation
    simul_2.set_params(simul_params.clone());
    // Simulate instrument arm
    if args.arm == 1 {
        simul_2.simulate_arm(RedArm);
    } else {
        simul_2.simulate_arm(BlueArm);
    }

    let mut snr_curve = SNRcurve::new();
    let mut calculation_worker = CalculationWorker::init();
    if args.arm == 1 {
        calculation_worker.simulate(&mut snr_curve, &mut simul_2, RedArm);
    } else {
        calculation_worker.simulate(&mut snr_curve, &mut simul_2, BlueArm);
    }
    let _ = snr_curve.signal.save_curve("calc_work_simulated_signal_curve.csv");
    let _ = snr_curve.noise.save_curve("calc_work_simulated_noise_curve.csv");
    let _ = snr_curve.counts.save_curve("calc_work_simulated_counts_curve.csv");


}

// THIS THING GETS RUNNED LIKE THIS: 
// cargo run -- --help
// cargo run -- --arg1 ARG1 ...