use crate::spectrum::*;
use crate::curve::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_spectrum_initialization() {
        let mut spectrum = Spectrum::default();
        assert!(spectrum.get_curve().get_map_mut().is_empty());
        assert_eq!(spectrum.clone().get_curve().get_oob_left(), 0.0);
        assert_eq!(spectrum.clone().get_curve().get_oob_right(), 0.0);
    }
    #[test]
    fn test_spectrum_load() {
        let mut spectrum = Spectrum::default();
        let _ = spectrum.get_curve_mut().load_curve("test-data/curve_constant.csv").unwrap();
        assert!(!spectrum.get_curve().get_map().is_empty());
        assert!(spectrum.get_curve().integral() - 25.0 < 1e-4);
    }
    #[test]
    fn test_scale_axis_factor_spec() {

    }
    #[test]
    fn test_scale_axis_curve_spec() {

    }
    #[test]
    fn test_scale_axis_curve_diff_spec() {

    }
    #[test]
    fn test_invert_axis_spec() {

    }
}

