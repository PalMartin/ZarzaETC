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


//pub mod test_helpers;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_f64_eq;
    use crate::assert_f32_eq;

    #[test]
    fn input_line_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        // Modify NDIT
        *simul_params.get_ndit_mut() = 3;    
        let mut simul = Simulation::new();
        // Set input object spectrum
        simul.set_input_line(NotNan::new(6563.0).unwrap(), 1e-4, 0.1, 1e-5, &String::from("emission_line"), &String::from("Resolved"), RedArm);
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());

        assert_f64_eq!(simul.get_input().get_curve().get_point(NotNan::new(6563.0).unwrap()), 0.00011) // fix to do the same as for galaxy templates -> A to m
    }

    #[test]
    fn input_template_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        // Modify NDIT
        *simul_params.get_ndit_mut() = 3;
        let mut simul = Simulation::new();
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        println!("PUNTO: {}", simul.get_input().get_curve().get_point(NotNan::new(6.57e-7).unwrap()));

        assert_f64_eq!(simul.get_input().get_point(NotNan::new(6.57e-7).unwrap()), 2.8268339999999998e-14 * 1e10 * 1e-3)
    }

    #[test]
    fn para_params() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Modify NDIT
        *simul_params.get_ndit_mut() = 3;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);


        // CALCULAR A MANO LA SEÑAL ESPERADA PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC

        //assert_f64_eq!()
        println!("QE: {:?}", simul.get_det().get_q_e().get_point(NotNan::new(6.55e-7).unwrap()));
        println!("TRANS: {}", simul.get_tarsis_model().get_red_ml15().get_point(NotNan::new(6.55e-7).unwrap()));
        println!("DISP: {}", simul.get_tarsis_model().get_red_disp()[20].get_point(NotNan::new(6.55e-7).unwrap()));
        println!("W2PX: {}", simul.get_tarsis_model().get_red_w2px()[20].get_point(NotNan::new(6.55e-7).unwrap()));
        println!("PX2W: {}", simul.get_tarsis_model().get_red_px2w()[20].get_point(NotNan::new(6.55e-7).unwrap()));
        println!("RES: {}", simul.get_tarsis_model().get_red_repx()[20].get_point(NotNan::new(6.55e-7).unwrap()));
        println!("EXT: {}", simul.get_sky_model().get_sky_ext().get_point(NotNan::new(6.55e-7).unwrap()));

        println!("SIGNAL: {:?}", simul.get_det_mut().signal_crv().get_point(NotNan::new(6.55e-7).unwrap()));
        println!("NOISE: {:?}", simul.get_det_mut().noise_crv().get_point(NotNan::new(6.55e-7).unwrap()));
        println!("SNR: {:?}", simul.get_det_mut().snr_crv().get_point(NotNan::new(6.55e-7).unwrap()));
        //println!("NOISE: {:?}", simul.get_det_mut().noise_crv());
        //println!("SNR: {:?}", simul.get_det_mut().snr_crv());

        assert_f64_eq!(simul.get_input().get_curve().get_point(NotNan::new(6.57e-7).unwrap()), 1.0 as f64)

    }

    #[test]
    fn signal_px_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO LA SEÑAL ESPERADA PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC


        assert_f64_eq!(simul.get_det_mut().signal_px(NotNan::new(6.57e-7).unwrap()), 2.3729193744677297e-7)

    }

    #[test]
    fn signal_crv_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO LA SEÑAL ESPERADA PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC, PERO EXTRAYENDO EL PUNTO DEL CALCULO CON CURVAS

        //assert_f64_eq!()

    }

    #[test]
    fn noise_px_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO EL RUIDO ESPERADO PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC

        //assert_f64_eq!()

    }

    #[test]
    fn noise_crv_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO EL RUIDO ESPERADO PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC, PERO EXTRAYENDO EL PUNTO DEL CALCULO CON CURVAS

        //assert_f64_eq!()

    }

    #[test]
    fn snr_px_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO LA SNR ESPERADO PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC

        //assert_f64_eq!()

    }

    #[test]
    fn snr_crv_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        // CALCULAR A MANO LA SNR ESPERADA PARA EL MISMO PUNTO QUE ANTES (6.57e-7), TENIENDO EN CUENTA LA EMISIÓN DEL CIELO ETC, PERO EXTRAYENDO EL PUNTO DEL CALCULO CON CURVAS

        //assert_f64_eq!()

    }

    #[test]
    fn t_exp_test() {
        // Initialize simulation parameters
        let mut simul_params = SimulationParams::new();
        // Modify/set slice
        *simul_params.get_slice_mut() = 20;
        // Initialize simulator
        let mut simul = Simulation::new();
        // Set detector
        simul.get_det_mut().set_detector("CCD231-84-0-H69");
        // Set input object spectrum
        simul.set_input_template(&String::from("Sa"));
        // Set simulation parameters into the simulation
        simul.set_params(simul_params.clone());
        // Simulate instrument arm
        simul.simulate_arm(RedArm);

        //assert_f64_eq!()

    }
 }

