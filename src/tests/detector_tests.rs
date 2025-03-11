use ordered_float::NotNan;
use crate::detector::*;
use crate::curve::*;
use crate::spectrum::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_detector_properties() {
        let properties = DetectorProperties::new();
        assert_eq!(properties.get_ccd231_84_0_s77().get_coating(), "NBB");
    }

    #[test]
    fn test_new_detector() {
        let detector = Detector::new();
        assert_eq!(*detector.get_exposure_time(), 1.0);
        assert_eq!(detector.get_properties().get_ccd231_84_0_s77().get_coating(), "NBB");
    }

    #[test]
    fn test_set_detector_and_get_spec() {
        let mut detector = Detector::new();
        detector.set_detector(&"CCD231-84-0-S77".to_string());
        assert_eq!(detector.get_spec(), detector.get_properties().get_ccd231_84_0_s77().clone());
        assert_eq!(detector.get_spec().get_coating(), detector.get_properties().get_ccd231_84_0_s77().get_coating());
    }

    #[test]
    fn test_dark_electrons() {
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        let dark_e = detector.dark_electrons(DETECTOR_TEMPERATURE);
        assert!(dark_e >= 0.0);
    }

    #[test]
    fn test_read_out_noise() {
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        let ron = detector.read_out_noise();
        assert!(ron >= 0.0);
    }

    #[test]
    fn test_set_pixel_photon_flux() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let _ = detector.get_photon_flux_per_pixel().get_curve().save_curve("test-data/curve_constant_pixel_photon_flux.csv").unwrap();
        assert!(!detector.get_photon_flux_per_pixel().get_curve().get_map().is_empty())
    }

    #[test]
    fn test_set_exposure_time() {
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_exposure_time(1.0);
        assert!(detector.get_exposure_time() - 1.0 > -0.9);
        assert!(detector.get_exposure_time() - 1.0 < 10e-4);
    }

    #[test]
    fn test_signal_px() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let _ = detector.get_signal().get_curve().save_curve("test-data/curve_constant_signal.csv").unwrap();
        assert!(detector.signal_px(NotNan::new(2.0).unwrap()) - 0.000000001125 < 10e-4);
        assert!(detector.signal_px(NotNan::new(2.0).unwrap()) - 0.000000001125 > -0.000000001125);
    }
    
    #[test]
    fn test_electrons_px() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let _ = detector.get_electrons_per_pixel().get_curve().save_curve("test-data/curve_constant_electron.csv").unwrap();
        assert!(detector.signal_px(NotNan::new(2.0).unwrap()) - 0.04886848172065906 < 10e-4);
        assert!(detector.signal_px(NotNan::new(2.0).unwrap()) - 0.04886848172065906 > -0.04886848172065906);
    }

    #[test]
    fn test_signal_crv() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let _ = detector.get_signal().get_curve().save_curve("test-data/curve_constant_signal.csv").unwrap();
        assert!(detector.get_signal().get_curve().get_point(NotNan::new(2.0).unwrap()) - 0.000000001125  > -0.000000001125);
        assert!(detector.get_signal().get_curve().get_point(NotNan::new(2.0).unwrap()) - 0.000000001125 < 1e-4);
        
    }

    #[test]
    fn test_electrons_crv() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let _ = detector.get_electrons_per_pixel().get_curve().save_curve("test-data/curve_constant_signal.csv").unwrap();
        assert!(detector.get_electrons_per_pixel().get_curve().get_point(NotNan::new(2.0).unwrap()) - 0.04886848172065906  > -0.04886848172065906);
        assert!(detector.get_electrons_per_pixel().get_curve().get_point(NotNan::new(2.0).unwrap()) - 0.04886848172065906 < 1e-4);
    }

    #[test]
    fn test_noise() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let noise = detector.noise(NotNan::new(2.0).unwrap());
        
    }

    #[test]
    fn test_snr() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        let mut detector = Detector::new();
        detector.set_detector("CCD231-84-0-S77");
        detector.set_pixel_photon_flux(spectrum);
        let noise = detector.snr(NotNan::new(2.0).unwrap());

    }
}

