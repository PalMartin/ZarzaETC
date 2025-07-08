//use ordered_float::NotNan;
use crate::detector::*;
use crate::curve::*;
//use crate::curve::CurveAxis::*;
// use crate::spectrum::*;
// use crate::simulation::*;
use crate::instrument_model::*;
use crate::sky_model::*;
// use crate::helpers::*;
//use crate::instrument_model::InstrumentArm::RedArm;
//use crate::instrument_model::InstrumentArm::BlueArm;
use crate::data_file_manager::*;
use crate::config_manager::*;
use std::sync::Once;
use std::env;

//pub mod test_helpers;


#[cfg(test)]
mod tests {
    //use crate::instrument_model;

    use super::*;
    // use crate::assert_f64_eq;
    // use crate::assert_f32_eq;

    static INIT: Once = Once::new();

    pub fn init_env() {
        INIT.call_once(|| {
            env::set_var("ZARZA_ETC_DATA_DIR", "data"); // así
        });
    }

    #[test]
    fn data_loading_test() {
        
        init_env();

        //println!("SEARCH PATHS: {:?}", DataFileManager::instance().search_paths());

        let data_file_path = DataFileManager::data_file(&"blueTransmission.csv".to_string());
        //println!("PATH FOUND: {}", data_file_path);

        let mut my_curve = Curve::default();
        let _ = my_curve.load_curve(&DataFileManager::data_file(&"blueTransmission.csv".to_string())).unwrap();
        //println!("LOADED CURVE: {:?}", my_curve);

    }

    #[test]
    fn config_test(){

        let mut config = Config::new("detectors".to_string());
        config.load();
        //println!("CONFIG TEST: {:?}", config);

    }

    #[test]
    fn new_detector_test(){

        let detector = Detector::new();
        
        //println!("DETECTOR TEST: {:?}", detector);

    }

    #[test]
    fn set_detector_test(){

        let mut detector = Detector::new();
        detector.set_detector("ccd231_84_0_s77".to_string());

        //println!("DETECTOR PROPERTIES: {:?}", detector.get_properties());
        //println!("DETECTOR SELECTED: {:?}", detector.get_detector());

    }

    #[test]
    fn new_instrument_model() {

        let instrument = InstrumentModel::new();
        println!("INSTRUMENT MODEL PROPERTIES -> f_num: {:?}", instrument.properties().get_f_num()); 

    }

    #[test]
    fn new_sky_model() {

        let sky: SkyModel = SkyModel::new();
        println!("SKY MODEL PROPERTIES -> airmass: {:?}", sky.properties.sky_emission); 
        
    }

 }

