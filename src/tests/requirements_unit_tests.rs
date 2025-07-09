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
use crate::data_file_manager::*;
use crate::config_manager::*;
use std::sync::Once;
use std::env;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_f64_eq;
    use crate::assert_f32_eq;

    static INIT: Once = Once::new();

    pub fn init_env() {
        INIT.call_once(|| {
            env::set_var("ZARZA_ETC_DATA_DIR", "data"); // así
        });
    }


    #[test]
    fn transmission_requirement_blueml15() {

        init_env();

        // SIMULATION:
        //let simul = Simulation::new(&"filters/Generic_Cousins_R.csv".to_string());
        let simul = Simulation::new();

        // TRANSMISSION:
        // Extract the transmission:
        let trans = simul.get_tarsis_model().get_blue_ml15().clone().get_map();
        let count = trans.len();
        let avg = trans.values().copied().sum::<f64>() / count as f64;
        //println!("AVG TRANSMISSION DETECTOR: {}", avg/0.787);
        println!("AVG: {}", avg);

        //assert!(trans.clone().values().all(|&v| v > 0.338) , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (ML15)\x1b[0m");
        //assert!(avg > 0.266 , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (ML15)\x1b[0m");
        assert!(0.327 < avg && avg < 0.349 , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (ML15)\x1b[0m");

    }

    #[test]
    fn transmission_requirement_bluenbb() {

        init_env();
        // SIMULATION:
        //let simul = Simulation::new(&"filters/Generic_Cousins_R.csv".to_string());
        let simul = Simulation::new();


        // TRANSMISSION:
        // Extract the transmission:
        let trans = simul.get_tarsis_model().get_blue_nbb().clone().get_map();
        let count = trans.len();
        let avg = trans.values().copied().sum::<f64>() / count as f64;
        println!("AVG: {}", avg);
        
        //assert!(trans.clone().values().all(|&v| v > 0.338) , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (NBB)\x1b[0m");
        //assert!(avg > 0.266 , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (NBB)\x1b[0m");
        assert!(0.327 < avg && avg < 0.349 , "\x1b[1;5;31mMinimum transmission requirements not met on the blue arm (NBB)\x1b[0m");
    }

    #[test]
    fn transmission_requirement_redml15() {
        // SIMULATION:
        let simul = Simulation::new();

        // TRANSMISSION:
        // Extract the transmission:
        let trans = simul.get_tarsis_model().get_red_ml15().clone().get_map();
        let count = trans.len();
        let avg = trans.values().copied().sum::<f64>() / count as f64;
        println!("AVG: {}", avg);

        //assert!(trans.clone().values().all(|&v| v > 0.402) , "\x1b[1;5;31mMinimum transmission requirements not met on the red arm\x1b[0m");
        //assert!(avg > 0.319 , "\x1b[1;5;31mMinimum transmission requirements not met on the red arm\x1b[0m");
        assert!(0.389 < avg && avg < 0.415 , "\x1b[1;5;31mMinimum transmission requirements not met on the red arm\x1b[0m");
    }

    #[test]
    fn blue_disperssion_requirement() {
        // SIMULATION:
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        for i in 0..TARSIS_SLICES {
            // Modify/set slice
            *simul_params.get_slice_mut() = i; 
            // Modify NDIT
            //*simul_params.get_ndit_mut() = 3;
            // Initialize simulator
            let simul = Simulation::new();

            let mut disp = simul.get_tarsis_model().get_blue_disp()[i].clone();
            // undo operations made in code
            disp.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN"));
            disp.scale_axis(YAxis, 1e9);
            
            assert!(disp.get_point(NotNan::new(3.2e-7).unwrap()) > 0.108 && disp.get_point(NotNan::new(3.2e-7).unwrap()) < 0.120, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
            assert!(disp.get_point(NotNan::new(4.2e-7).unwrap()) > 0.101 && disp.get_point(NotNan::new(4.2e-7).unwrap()) < 0.113, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
            assert!(disp.get_point(NotNan::new(5.2e-7).unwrap()) > 0.100 && disp.get_point(NotNan::new(5.2e-7).unwrap()) < 0.112, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
        }
    }

    #[test]
    fn red_disperssion_requirement() {
        // SIMULATION:
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        for i in 0..TARSIS_SLICES {
            // Modify/set slice
            *simul_params.get_slice_mut() = i; 
            // Modify NDIT
            //*simul_params.get_ndit_mut() = 3;
            // Initialize simulator
            let simul = Simulation::new();

            let mut disp = simul.get_tarsis_model().get_red_disp()[i].clone();
            // undo operations made in code
            disp.invert_axis(YAxis, NotNan::new(1.0).expect("x should not be NaN"));
            disp.scale_axis(YAxis, 1e9);
            
            assert!(disp.get_point(NotNan::new(5.1e-7).unwrap()) > 0.170 && disp.get_point(NotNan::new(5.1e-7).unwrap()) < 0.188, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
            assert!(disp.get_point(NotNan::new(6.6e-7).unwrap()) > 0.168 && disp.get_point(NotNan::new(6.6e-7).unwrap()) < 0.186, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
            assert!(disp.get_point(NotNan::new(8.1e-7).unwrap()) > 0.167 && disp.get_point(NotNan::new(8.1e-7).unwrap()) < 0.185, "\x1b[1;5;31mMinimum dispersion requirements not met on the blue arm\x1b[0m");
        }
    }



 }

