use std::error::Error;
use csv::ReaderBuilder;
use csv::WriterBuilder;
use std::collections::BTreeMap;
use ordered_float::NotNan;

pub struct Curve {
    oob_right: f64,
    oob_left: f64,
    units_x: String,
    units_y: String,
    curve: BTreeMap<NotNan<f64>,f64>,
}

pub enum CurveAxis {
    XAxis,
    YAxis,
}

pub trait CurveOperations { 
    fn default() -> Self;
    fn is_oob(&self, x: NotNan<f64>) -> bool;
    fn set_units(&mut self, axis: CurveAxis, units: String);
    fn integral(&mut self) -> f64;
    fn dist_mean(&mut self) -> f64;
    fn integrate(&mut self, k: f64) -> ();
    fn flip(&mut self) -> ();
    fn extend_left(&mut self) -> ();
    fn extend_right(&mut self) -> ();
    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: NotNan<f64>) -> ();
    fn scale_axis_curve(&mut self, axis: CurveAxis, other: &mut Curve) -> ();
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) -> ();
    fn x_points(&self) -> Vec<NotNan<f64>>;
    fn set(&mut self, x: NotNan<f64>, y: f64) -> ();
    fn get_point(&self, x: NotNan<f64>) -> f64;
    fn get_diff(&self, x: NotNan<f64>) -> f64;
    fn multiply_by(&mut self, other: &mut Curve) -> ();
    fn add_curve(&mut self, other: &mut Curve) -> ();
    fn add_value(&mut self, val: f64) -> ();
    fn assign(&mut self, other: &Curve) -> ();
    fn from_existing(&mut self, other: &Curve, y_units: f64) -> ();
    fn clear(&mut self) -> ();
    fn debug(&self) -> ();
    fn load_curve(&mut self, path: &str) -> Result<BTreeMap<NotNan<f64>,f64>, Box<dyn Error>>;
    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>>;
}

impl CurveOperations for Curve {

    fn default() -> Curve {
        Curve {
            oob_right: 0.0,
            oob_left: 0.0,
            units_x: String::new(),
            units_y: String::new(),
            curve: BTreeMap::new(),
        }
    }

    fn is_oob(&self, x: NotNan<f64>) -> bool {
        if self.curve.is_empty() { 
            return true;
        }
        if x < *self.curve.keys().next().unwrap() { 
            return true;
        }
        if *self.curve.keys().next_back().unwrap() < x { 
            return true;
        }
        else {
            return false;
        }
    }

    fn set_units(&mut self, axis: CurveAxis, units: String) -> () {
        match axis {
            CurveAxis::XAxis => self.units_x = Some(units).unwrap(),
            CurveAxis::YAxis => self.units_y = Some(units).unwrap(),
        }
    }

    fn integral(&mut self) -> f64 {  
        // Initialize variables
        let mut accum = 0.0;
        let mut err = 0.0;

        // Check if out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return f64::INFINITY * (self.oob_left + self.oob_right);
        } 

        if self.curve.is_empty() {
            return 0.0;
        }

        else{
            // Get first point and remove it for the loop
            let (x0, y0) = self.curve.pop_first().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            for (x, y) in self.curve.iter() {
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);

                // Kahan summation
                let term_real = my * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;

                err = term_err - term_real;
                accum = tmp;

                x_prev = *x;
                y_prev = *y;
            }
            // Put the first point back
            self.curve.insert(x0, y0);
            return accum;
        }
    }

    fn dist_mean(&mut self) -> f64 {
        // Initialize variables
        let mut accum = 0.0;
        let mut err = 0.0;

        // Curve is out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return 0.5 * (self.oob_left + self.oob_right);
        } 

        if self.curve.is_empty() {
            return 0.0;
        }

        else{
            // Get first point and remove it for the loop
            let (x0, y0) = self.curve.pop_first().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            for (x, y) in self.curve.iter() {
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);
                let x_0 = 0.5 * *((x + x_prev));

                // Kahan summation
                let term_real = x_0 * my * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;

                err = term_err - term_real;
                accum = tmp;

                x_prev = *x;
                y_prev = *y;
            }
            // Put the first point back
            self.curve.insert(x0, y0);
            return accum/self.integral();
        }
    }

    fn integrate(&mut self, k: f64) -> () { 
        // Initialize variables
        let mut accum = 0.0;
        let mut err = 0.0;

        self.oob_left = k;

        if self.curve.is_empty() {
            self.oob_right = k;
            return;
        }

        if self.curve.len() == 1 {
            let x = *self.curve.keys().next().unwrap();
            self.curve.clear();
            self.curve.insert(x, k);
            self.oob_right = k;
            return;
        }

        else {
            // Get first point and remove it for the loop
            let (x0, y0) = self.curve.pop_first().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            let mut dx: NotNan<f64> = NotNan::new(0.0).unwrap();

            for (x, y) in self.curve.iter() {
                dx = x - x_prev;

                // Kahan summation
                let term_real = 0.5 * (y + y_prev) * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp-accum;
    
                err = term_err - term_real;
                accum = tmp;
    
                x_prev = *x;
                y_prev = *y;
            }
            self.curve.insert(x0, k);
            let x_next = *self.curve.keys().last().unwrap() + dx;
            self.curve.insert(x_next, *self.curve.values().last().unwrap() * *((*self.curve.keys().last().unwrap() - x_prev + dx)) + self.get_point(x_prev));

            self.oob_right = accum;
        }
    }

    fn flip(&mut self) {
        // Unit flip
        std::mem::swap(&mut self.units_x, &mut self.units_y);
        // x & y coordinates flip
        let mut flip_crv = BTreeMap::new();
        for (&x, &y) in &self.curve {
            flip_crv.insert(x,y);
        }
        self.curve = flip_crv;
    }

    fn extend_left(&mut self) -> () {
        if self.curve.is_empty() {
            return;
        }
        self.oob_left = *self.curve.values().next().unwrap();
    }

    fn extend_right(&mut self) -> () {
        if self.curve.is_empty() {
            return;
        }
        self.oob_right = *self.curve.values().next_back().unwrap(); 
    }

    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: NotNan<f64>) -> () { 
        match axis {
            CurveAxis::XAxis => {
                for &(mut x) in self.curve.keys() {
                    x *= factor;
                }
            }
            CurveAxis::YAxis => {
                for y in self.curve.values_mut() {
                    *y *= *factor;
                }
            }
        }
        self.oob_left *= *factor;
        self.oob_right *= *factor;
    }

    fn scale_axis_curve(&mut self, axis: CurveAxis, other: &mut Curve) -> () { 
        match axis {
            CurveAxis::XAxis => {
                for (ref mut x_1, x_2) in self.curve.keys().zip(other.curve.keys()) {
                    *x_1 = x_2;
                }
            }
            CurveAxis::YAxis => {
                for (ref mut y_1, y_2) in self.curve.values_mut().zip(other.curve.values_mut()) {
                    *y_1 = y_2;
                }
            }
        }
        self.oob_left = other.oob_left;
        self.oob_right = other.oob_right;
    }

    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) -> () {
        match axis {
            CurveAxis::XAxis => { 
                for ref mut x in self.curve.keys() {
                    *x = &(factor / *x);
                }
            }
            CurveAxis::YAxis => {
                for ref mut y in self.curve.values_mut() {
                    **y *= &(*factor / **y);
                }
            }
        }
        self.oob_left = 1.0 / self.oob_left;
        self.oob_right = 1.0 / self.oob_right;
    }

    fn x_points(&self) -> Vec<NotNan<f64>> {
        return self.curve.clone().into_keys().collect();
    } 

    fn set(&mut self, x: NotNan<f64>, y: f64) -> () {
        self.curve.insert(x, y);
    }

    fn get_point(&self, x: NotNan<f64>) -> f64 { // NOT WORKING WELL
        if self.curve.contains_key(&x) {
            return *self.curve.get(&x).unwrap();
        }
        else {
            let (&x0, &y0) = self.curve.range(..x).next_back().unwrap_or_else(|| panic!("No previous point found for interpolation"));
            let (&x1, &y1) = self.curve.range(x..).next().unwrap_or_else(|| panic!("No next point found for interpolation"));
            let y = y0 + *((x - x0) * (y1 - y0) / (x1 - x0));
            return y;
        }
    }

    fn get_diff(&self, x: NotNan<f64>) -> f64 {
        let (&prev,_) = self.curve.range(..x).next_back().unwrap_or_else(|| panic!("No previous point found for the provided value"));
        let (&next,_) = self.curve.range(x..).next().unwrap_or_else(|| panic!("No next point found for the provided value"));
        let Some((&x0,&y0)) = self.curve.get_key_value(&prev) else { panic!("Key not found") };
        let Some((&x1,&y1)) = self.curve.get_key_value(&next) else { panic!("Key not found") };
        if x1 == *self.curve.keys().next().unwrap() {
            return 0.0;
        }
        if x1 == *self.curve.keys().next_back().unwrap() {
            return 0.0;
        }
        else {
            // Two cases
            if x1 != x {
                // Middle of segment
                let diff = (y1 - y0) / *((x1 - x0));
                return diff;
            }
            else {
                // Edge of a segment
                let (&x2, &y2) = self.curve.range(next..).nth(1).unwrap_or_else(|| panic!("No next-next point found for the provided value"));
                if x2 == *self.curve.keys().next_back().unwrap() {
                    return 0.0
                }
                else {
                    let diff = (y2 - y0) / *((x2 - x0));
                    return diff;
                }
            }
        }
    }

    fn multiply_by(&mut self, other: &mut Curve) { 
        self.curve.append(&mut other.curve);
        for (x, ref mut y) in self.curve.iter() {
            *y = &(self.get_point(*x) * other.get_point(*x));
        }
        self.oob_left *= other.oob_left;
        self.oob_right *= other.oob_right;
    }

    fn add_curve(&mut self, other: &mut Curve) { 
        self.curve.append(&mut other.curve);
        for (x, ref mut y) in self.curve.iter() {
            *y = &(self.get_point(*x) + other.get_point(*x));
        }
        self.oob_left += other.oob_left;
        self.oob_right += other.oob_right;
    }

    fn add_value(&mut self, val: f64) -> () {
        for (x, ref mut y) in self.curve.iter() {
            *y = &(self.get_point(*x) + val);
        }
        self.oob_left += val;
        self.oob_right += val;
    }

    fn assign(&mut self, other: &Curve) -> () { 
        let Some((&own_first_x, &_)) = self.curve.first_key_value() else { todo!() };
        let Some((&crv_first_x, &_)) = other.curve.first_key_value() else { todo!() };
        let Some((&own_last_x, &_)) = self.curve.last_key_value() else { todo!() };
        let Some((&crv_last_x, &_)) = other.curve.last_key_value() else { todo!() };
        // Nothing to add
        if crv_first_x == crv_last_x {
            return;
        }
        // No curve
        if own_first_x == own_last_x {
            return;
        }
        // Assign middle part
        for (x, ref mut y) in self.curve.iter() {
            if other.is_oob(*x) == false {
                *y = &other.get_point(*x);
            }
        }
        // Curve longer to the left
        if crv_first_x < own_first_x {
            self.oob_left = other.oob_left;
        }
        // Curve longer to the right
        if own_last_x < crv_last_x {
            self.oob_right = other.oob_right;
        }
        // Assign whole curve
        let mut new_pairs = BTreeMap::new();
        for (x, y) in self.curve.iter_mut() {
            new_pairs.insert(*x, *y);
        }
        self.curve.clear();
        self.curve.append(&mut new_pairs)
    }

    fn from_existing(&mut self, other: &Curve, y_units: f64) -> () {
        self.curve = other.curve.clone();

        if y_units != 1.0 {
            for (_, y) in self.curve.iter_mut() {
                *y *= y_units;
            }
        }
    }

    fn clear(&mut self) -> () {
        self.curve.clear();
        self.oob_left = 0.0;
        self.oob_right = 0.0;
    }

    fn debug(&self) -> () {
        for (x, y) in self.curve.iter() {
            print!("{} = {}, ", x, y)
        }
        println!();
    }

    fn load_curve(&mut self, path: &str) -> Result<BTreeMap<NotNan<f64>,f64>, Box<dyn Error>> {

        let mut rdr = ReaderBuilder::new().has_headers(false).from_path(path)?;
        for result in rdr.records() {
            let record = result?;

            let x: NotNan<f64> = record.get(0).unwrap().parse()?;
            let y: f64 = record.get(1).unwrap().parse()?;

            self.curve.insert(x, y);
        }
        Ok(self.curve.clone())
    }

    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>> {
        
        let mut wtr = WriterBuilder::new().has_headers(false).from_path(path).unwrap();
        //wtr.write_record(&["X", "Y"])?; 
        for (x, y) in self.curve.iter() {
            wtr.write_record(&[x.to_string(), y.to_string()]).unwrap();
        }
        wtr.flush().unwrap(); 
        Ok(())
    }
}
