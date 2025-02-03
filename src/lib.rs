pub mod curve;

#[cfg(test)]
pub mod tests {

    #[test]
    fn default_load_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");
    }
    #[test]
    fn integral_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.integral();
    }
    #[test]
    fn dist_mean_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.dist_mean();
    }
    #[test]
    fn integrate_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.integrate(5.0)
    }
    #[test]
    fn flip_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.flip()
    }
    #[test]
    fn extend_left_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.extend_left()
    }
    #[test]
    fn extend_right_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_tests.csv");

        my_curve.extend_right()
    }
    #[test]
    fn multiply_by_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let mut my_curve_2 = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");
        let _ = my_curve_2.load_curve("./curve_tes_2.csv");

        my_curve.multiply_by(&mut my_curve_2);
    }
    #[test]
    fn add_curve_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let mut my_curve_2 = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");
        let _ = my_curve_2.load_curve("./curve_test_2.csv");

        my_curve.add_curve(&mut my_curve_2);
    }
    #[test]
    fn add_value_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv");

        my_curve.add_value(5.0);
    }

}
