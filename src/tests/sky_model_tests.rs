use ordered_float::NotNan;

use crate::sky_model::*;
use crate::curve::*;
use crate::spectrum::*;
use crate::helpers::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_sky_properties_creation() {
        let props = SkyProperties::new();
        assert_eq!(props.sky_emission, "src/sky_model/CAHASKy.csv");
        assert_eq!(props.sky_emission_ref_airmass, 1.0);
        assert_eq!(props.sky_extinction, "src/sky_model/CAHASkyExt.csv");
    }
    #[test]
    fn test_sky_model_initialization() {
        let mut model = SkyModel::new();
        assert!(!model.get_sky_spectrum().get_curve().get_curve().is_empty());
        assert!(!model.get_sky_ext().clone().get_curve().is_empty());
        assert!(!model.get_moon_to_mag().clone().get_curve().is_empty());
        assert_eq!(*model.get_airmass(), 1.0);
        assert_eq!(*model.get_moon_fraction(), 0.0);
    }
    #[test]
    fn test_set_moon() {
        let mut model = SkyModel::new();
        model.set_moon(50.0);
        assert_eq!(*model.get_moon_fraction(), 50.0);
    }
    #[test]
    fn test_set_airmass() {
        let mut model = SkyModel::new();
        model.set_airmass(2.0);
        assert_eq!(*model.get_airmass(), 2.0);
    }
    #[test]
    fn test_set_zenith_distance() {
        let mut model = SkyModel::new();
        model.set_zenith_distance(60.0);
        assert!((model.get_airmass() - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_make_sky_spectrum() {
        let model = SkyModel::new();
        let mut object_spectrum: Spectrum = Default::default();
        let _ = object_spectrum.load_curve("test-data/curve_test.csv").unwrap();
        let mut result_spectrum = model.make_sky_spectrum(&object_spectrum);
        //Debugging
        println!("Sky Model: {:?}", result_spectrum.get_curve_mut().get_curve_mut());
        let _ = result_spectrum.get_curve_mut().save_curve("test-data/sky_model_test.csv");

        assert!(!result_spectrum.x_points().is_empty());
    }
}

