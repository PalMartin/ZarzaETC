use ordered_float::NotNan;
use crate::instrument_model::*;
use crate::curve::*;
use crate::spectrum::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_instrument_model_initialization() {
        let model = InstrumentModel::new();
        let _ = model.get_blue_ml15().save_curve("test-data/blue_ml15_test.csv");
        let _ = model.get_blue_nbb().save_curve("test-data/blue_nbb_test.csv");
        let _ = model.get_red_ml15().save_curve("test-data/red_ml15_test.csv");
        let _ = model.get_blue_w2px()[3].save_curve("test-data/blue_disp_test.csv");
    }

    #[test]
    fn test_properties() {
        let model = InstrumentModel::new();
        let props = model.properties();
        assert!(*props.get_f_num() > 0.0);
        assert!(*props.get_ap_efficiency() > 0.0);
        assert!(!props.get_coating().is_empty());
    }

    #[test]
    fn test_px_to_wavelength_val() {
        let mut model = InstrumentModel::new();
        let pixel = NotNan::new(100.0).unwrap();
        let result = model.px_to_wavelength_val(InstrumentArm::BlueArm, 0, pixel);
        assert!(result.is_ok());
    }

    #[test]
    fn test_px_to_wavelength_crv() {
        let mut model = InstrumentModel::new();
        let result = model.px_to_wavelength_crv(InstrumentArm::BlueArm, 0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_wavelength_to_px_val() {
        let mut model = InstrumentModel::new();
        let wavelength = NotNan::new(500e-9).unwrap();
        let result = model.wavelength_to_px_val(InstrumentArm::RedArm, 0, wavelength);
        assert!(result.is_ok());
    }

    #[test]
    fn test_wavelength_to_pix_crv() {
        let mut model = InstrumentModel::new();
        let result = model.wavelength_to_pix_crv(InstrumentArm::RedArm, 0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_set_input() {
        let mut model = InstrumentModel::new();
        let mut spectrum = Spectrum::default();
        let _ = spectrum.load_curve("test-data/curve_test.csv").unwrap();
        model.set_input(InstrumentArm::BlueArm, spectrum);
        assert_eq!(*model.get_current_path(), InstrumentArm::BlueArm);
    }

    #[test]
    fn test_convolve_around() {
        let mut curve = Curve::default();
        let _ = curve.load_curve("test-data/curve_test.csv").unwrap();
        let x0 = 500.0;
        let inv_sigma = 1.0;
        let oversample = 10;
        let result = InstrumentModel::convolve_around(&curve, x0, inv_sigma, oversample);
        assert!(result.is_finite());
    }

    #[test]
    fn test_switch_current_arm() {
        let mut model = InstrumentModel::new();
        model.switch_current_arm(InstrumentArm::RedArm);
        assert_eq!(*model.get_current_path(), InstrumentArm::RedArm);
    }

    #[test]
    fn test_make_pixel_photon_flux() {
        
    }
}
