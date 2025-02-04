pub mod curve;

#[cfg(test)]
pub mod tests {

    #[test]
    fn default_load_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap(); //.unwrap() is necessary to fail the test if the loading of the curve has an error
    }
    #[test]
    fn integral_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.integral();
    }
    #[test]
    fn dist_mean_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.dist_mean();
    }
    #[test]
    fn integrate_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.integrate(5.0)
    }
    #[test]
    fn flip_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.flip()
    }
    #[test]
    fn extend_left_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.extend_left()
    }
    #[test]
    fn extend_right_works() {
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.extend_right()
    }
    #[test]
    fn multiply_by_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let mut my_curve_2 = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();
        let _ = my_curve_2.load_curve("./curve_test_2.csv").unwrap();

        my_curve.multiply_by(&mut my_curve_2);

        let _ = my_curve.save_curve("./prueba.csv").unwrap();
    }
    #[test]
    fn add_curve_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let mut my_curve_2 = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();
        let _ = my_curve_2.load_curve("./curve_test_2.csv").unwrap();

        my_curve.add_curve(&mut my_curve_2);
    }
    #[test]
    fn add_value_works(){
        use crate::curve::{Curve, CurveOperations};

        let mut my_curve = <Curve as CurveOperations>::default();
        let _ = my_curve.load_curve("./curve_test.csv").unwrap();

        my_curve.add_value(5.0);
    }

}
