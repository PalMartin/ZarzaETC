use ordered_float::NotNan;

use crate::curve::*;
use crate::curve::CurveAxis::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_curve_initialization() {
        let mut curve = <Curve as CurveDefault>::default();
        assert!(curve.get_curve_mut().is_empty());
        assert_eq!(curve.clone().get_oob_left(), 0.0);
        assert_eq!(curve.clone().get_oob_right(), 0.0);
    }
    #[test]
    fn test_is_oob() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        assert_eq!(my_curve.is_oob(NotNan::new(-1.0).unwrap()), true);
        assert_eq!(my_curve.is_oob(NotNan::new(10.0).unwrap()), true);
        assert_eq!(my_curve.is_oob(NotNan::new(4.0).unwrap()), false);
    }
    #[test]
    fn test_set_units() {

    }
    #[test]
    fn test_integral() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        assert!(my_curve.integral() - 25.0 < 1e-4)
    }
    #[test]
    fn test_dist_mean() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        assert!(my_curve.dist_mean() - 2.5 < 1e-4)
    }
    #[test]
    fn test_integrate() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_2 = <Curve as CurveDefault>::default();
        let _ = my_curve_2.load_curve("test-data/curve_cum_int.csv").unwrap();
        my_curve.integrate(2.0);
        assert_eq!(my_curve.get_point(NotNan::new(5.0).unwrap()), my_curve_2.get_point(NotNan::new(5.0).unwrap()))
    }
    #[test]
    fn test_flip() {

    }
    #[test]
    fn test_invert_axis() {
        
    }
    #[test]
    fn test_x_points() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        assert_eq!(my_curve.x_points(), vec![0.0,1.0,2.0,3.0,4.0,5.0]);
    }
    #[test]
    fn test_set_and_get_point() {
        let mut curve = <Curve as CurveDefault>::default();
        curve.set(NotNan::new(1.0).unwrap(), 2.0);
        assert_eq!(curve.get_point(NotNan::new(1.0).unwrap()), 2.0);
    }
    #[test]
    fn test_get_diff() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        assert!(my_curve.get_diff(NotNan::new(2.5).unwrap()) < 1e-4) // slope of 0.0, uniform curve
    }
    #[test]
    fn test_assign() {

    }
    #[test]
    fn test_from_existing() {

    }
    #[test]
    fn test_scale_axis_1() {

    }
    #[test]
    fn test_multiply_1() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_2 = <Curve as CurveDefault>::default();
        let _ = my_curve_2.load_curve("test-data/curve_constant_2.csv").unwrap();
        my_curve.multiply(5.0);
        let _ = my_curve.save_curve("test-data/mult_test.csv");
        assert!(my_curve.integral() - my_curve_2.integral() < 1e-4);
    }
    #[test]
    fn test_add_1() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_2 = <Curve as CurveDefault>::default();
        let _ = my_curve_2.load_curve("test-data/curve_constant_2.csv").unwrap();
        my_curve.add(20.0);
        let _ = my_curve.save_curve("test-data/mult_test.csv");
        println!("{}", my_curve.integral());
        println!("{}", my_curve_2.integral());
        assert!(my_curve.get_point(NotNan::new(5.0).unwrap()) - my_curve_2.get_point(NotNan::new(5.0).unwrap()) < 1e-4);
    }
    #[test]
    fn test_scale_axis_2() {

    }
    #[test]
    fn test_multiply_2() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_mult = <Curve as CurveDefault>::default();
        let _ = my_curve_mult.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_2 = <Curve as CurveDefault>::default();
        let _ = my_curve_2.load_curve("test-data/curve_constant_2.csv").unwrap();
        my_curve.multiply(my_curve_mult);
        assert!(my_curve.integral() - my_curve_2.integral() < 1e-4);
    }
    #[test]
    fn test_add_2() {
        let mut my_curve = <Curve as CurveDefault>::default();
        let _ = my_curve.load_curve("test-data/curve_constant.csv").unwrap();
        let mut my_curve_3 = <Curve as CurveDefault>::default();
        let _ = my_curve_3.load_curve("test-data/curve_constant_3.csv").unwrap();
        let mut my_curve_2 = <Curve as CurveDefault>::default();
        let _ = my_curve_2.load_curve("test-data/curve_constant_2.csv").unwrap();
        my_curve.add(my_curve_3);
        assert!(my_curve.integral() - my_curve_2.integral() < 1e-4);
    }
}