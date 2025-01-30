pub mod curve;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let val1 = 0.0;
        let val2 = vec![0.0, 1.0, 2.0];
        let val3 = "test";
        let mut my_curve = super::curve::Curve{m_oob_right: val1, m_oob_left: val1, m_units_x: val3.to_string(), m_units_y: val3.to_string(), m_curve_x: val2.clone(), m_curve_y: val2.clone()};
        super::curve::CurveOperations::flip_2(&mut my_curve);
    }
}
