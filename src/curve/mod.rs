use std::error::Error;
use csv::ReaderBuilder;
use csv::WriterBuilder;
use std::collections::BTreeMap;
use ordered_float::NotNan;

#[derive(Default, Clone, Debug)]
pub struct Curve {
    oob_right: f64,
    oob_left: f64,
    units_x: String,
    units_y: String,
    curve: BTreeMap<NotNan<f64>,f64>,
}

impl Curve {
    pub fn get_curve_mut(&mut self) -> &mut BTreeMap<NotNan<f64>, f64> {
        &mut self.curve
    }
    pub fn get_oob_left_mut(&mut self) -> &mut f64 {
        &mut self.oob_left
    }
    pub fn get_oob_right_mut(&mut self) -> &mut f64 {
        &mut self.oob_right
    }
    pub fn get_curve(self) -> BTreeMap<NotNan<f64>, f64> {
        self.curve
    }
    pub fn get_oob_left(self) -> f64 {
        self.oob_left.clone()
    }
    pub fn get_oob_right(self) -> f64 {
        self.oob_right.clone()
    }
}

pub enum CurveAxis {
    XAxis,
    YAxis,
}

//pub trait CurveDefault {
//    fn default() -> Self;
//}

pub trait CurveOperations {
    fn is_oob(&self, x: NotNan<f64>) -> bool;
    fn set_units(&mut self, axis: CurveAxis, units: String);
    fn integral(&self) -> f64;
    fn dist_mean(&self) -> f64;
    fn integrate(&mut self, k: f64);
    fn flip(&mut self);
    fn extend_left(&mut self);
    fn extend_right(&mut self);
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>);
    fn x_points(&self) -> Vec<NotNan<f64>>;
    fn set(&mut self, x: NotNan<f64>, y: f64);
    fn get_point(&self, x: NotNan<f64>) -> f64;
    fn get_diff(&self, x: NotNan<f64>) -> f64;
    fn assign(&mut self, other: &Curve);
    fn from_existing(&mut self, other: &Curve, y_units: f64);
    fn clear(&mut self);
    fn debug(&self);
    fn load_curve(&mut self, path: &str) -> Result<(), Box<dyn Error>>;
    fn load_slice_dat(&mut self, path: &str, row_index: usize) -> Result<(), Box<dyn Error>>;
    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>>;
}

pub trait BinaryCurveOperations<T> {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: T);
    fn multiply(&mut self, by_what: T);
    fn add(&mut self, what: T);
}

impl CurveOperations for Curve {

    fn is_oob(&self, x: NotNan<f64>) -> bool {
        if self.curve.is_empty() { 
            return true;
        }
        let Some((&first_x, &_)) = self.curve.first_key_value() else { todo!() };
        if x < first_x { 
            return true;
        }
        let Some((&last_x, &_)) = self.curve.last_key_value() else { todo!() };
        if last_x < x { 
            return true;
        } else {
            return false;
        }
    }

    fn set_units(&mut self, axis: CurveAxis, units: String) {
        match axis {
            CurveAxis::XAxis => self.units_x = units,
            CurveAxis::YAxis => self.units_y = units,
        }
    }

    fn integral(&self) -> f64 {  
        // Initialize variables
        let mut accum = 0.0;
        let mut err = 0.0;

        // Check if out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return f64::INFINITY * (self.oob_left + self.oob_right);
        } 

        if self.curve.is_empty() {
            return 0.0;
        } else {
            // Initialize iterator
            let mut iter = self.curve.iter();
            // Get first point 
            let (&x0, &y0) = iter.next().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            for (&x, &y) in iter {
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);

                // Kahan summation
                let term_real = my * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;

                err = term_err - term_real;
                accum = tmp;

                x_prev = x;
                y_prev = y;
            }
            return accum;
        }
    }

    fn dist_mean(&self) -> f64 {
        // Initialize variables
        let mut accum = 0.0;
        let mut err = 0.0;

        // Curve is out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return 0.5 * (self.oob_left + self.oob_right);
        } 

        if self.curve.is_empty() {
            return 0.0;
        } else {
            // Initialize iterator
            let mut iter = self.curve.iter();
            // Get first point 
            let (&x0, &y0) = iter.next().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            for (&x, &y) in iter {
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);
                let x_0 = 0.5 * *((x + x_prev));

                // Kahan summation
                let term_real = x_0 * my * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;

                err = term_err - term_real;
                accum = tmp;

                x_prev = x;
                y_prev = y;
            }
            return accum/self.integral();
        }
    }

    fn integrate(&mut self, k: f64) {
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
        // Initialize iterator
        let mut iter = self.curve.iter();
        // Get first point 
        let (&x0, &y0) = iter.next().unwrap();
        let mut x_prev = x0;
        let mut y_prev = y0;

        let mut new_pairs = BTreeMap::new();
        let mut dx: NotNan<f64> = NotNan::new(0.0).unwrap();

        for (&x, &y) in iter {
            dx = x - x_prev;

            // Kahan summation
            let term_real = 0.5 * (y + y_prev) * *dx - err;
            let tmp = accum + term_real;
            let term_err = tmp - accum;
    
            err = term_err - term_real;
            accum = tmp;
            
            new_pairs.insert(x, accum);

            x_prev = x;
            y_prev = y;
        }            
        new_pairs.insert(x0, k);
        let x_next = *self.curve.keys().last().unwrap() + dx;
        new_pairs.insert(x_next, *self.curve.values().last().unwrap() * *((*self.curve.keys().last().unwrap() - x_prev + dx)) + self.get_point(x_prev));
        
        self.curve = new_pairs;
        self.oob_right = accum;
    }

    fn flip(&mut self) {
        // Unit flip
        std::mem::swap(&mut self.units_x, &mut self.units_y);
        // x & y coordinates flip
        let mut flip_crv = BTreeMap::new();
        for (&x, &y) in &self.curve {
            let y_notnan = NotNan::new(y).expect("y should not be NaN");
            flip_crv.insert(y_notnan, *x);
        }
        self.curve = flip_crv;
    }

    fn extend_left(&mut self) {
        if self.curve.is_empty() {
            return;
        }
        self.oob_left = *self.curve.values().next().unwrap();
    }

    fn extend_right(&mut self) {
        if self.curve.is_empty() {
            return;
        }
        self.oob_right = *self.curve.values().next_back().unwrap();
    }

    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.iter() {
                    let new_x = NotNan::new(*factor / x.into_inner());
                    new_pairs.insert(new_x.expect("Attempted to insert NaN"), *y);
                }
                self.curve = new_pairs;
            }
            CurveAxis::YAxis => {
                for (_, y) in self.curve.iter_mut() {
                    *y = *factor / *y;
                }
            }
        }

        self.oob_left = 1.0 / self.oob_left;
        self.oob_right = 1.0 / self.oob_right;
    }

    fn x_points(&self) -> Vec<NotNan<f64>> {
        return self.curve.clone().into_keys().collect();
    } 

    fn set(&mut self, x: NotNan<f64>, y: f64) {
        self.curve.insert(x, y);
    }

    fn get_point(&self, x: NotNan<f64>) -> f64 {
        if self.curve.is_empty() {
            panic!("The curve is empty.")
        }
        let Some((&first_x, &_)) = self.curve.first_key_value() else { todo!() };
        if x < first_x { 
            return self.oob_left;
        }
        let Some((&last_x, &_)) = self.curve.last_key_value() else { todo!() };
        if last_x < x { 
            return self.oob_right;
        }
        if self.curve.contains_key(&x) {
            return *self.curve.get(&x).unwrap();
        } else {
            let (&x0, &y0) = self.curve.range(..x).next_back().unwrap();
            let (&x1, &y1) = self.curve.range(x..).next().unwrap();
            let y = y0 + *((x - x0) * (y1 - y0) / (x1 - x0));
            return y;
        }
    }

    fn get_diff(&self, x: NotNan<f64>) -> f64 {
        let (&prev,_) = self.curve.range(..x).next_back().unwrap();
        let (&next,_) = self.curve.range(x..).next().unwrap();
        let Some((&x0,&y0)) = self.curve.get_key_value(&prev) else { panic!("Key not found") };
        let Some((&x1,&y1)) = self.curve.get_key_value(&next) else { panic!("Key not found") };
        if x1 == *self.curve.keys().next().unwrap() {
            return 0.0;
        }
        if x1 == *self.curve.keys().next_back().unwrap() {
            return 0.0;
        } else {
            // Two cases
            if x1 != x {
                // Middle of segment
                let diff = (y1 - y0) / *((x1 - x0));
                return diff;
            }
            else {
                // Edge of a segment
                let (&x2, &y2) = self.curve.range(next..).nth(1).unwrap();
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

    fn assign(&mut self, other: &Curve) { 
        let Some((&own_first_x, &_)) = self.curve.first_key_value() else { return; };
        let Some((&crv_first_x, &_)) = other.curve.first_key_value() else { return; };
        let Some((&own_last_x, &_)) = self.curve.last_key_value() else { return; };
        let Some((&crv_last_x, &_)) = other.curve.last_key_value() else { return; };
        // Nothing to add
        if other.curve.is_empty() {
            return;
        }
        // No curve
        if self.curve.is_empty() {
            self.curve = other.curve.clone();
        }
        // Assign middle part
        let mut new_pairs = BTreeMap::new();
        for (x, _) in self.curve.iter() {
            if !other.is_oob(*x) {
                new_pairs.insert(*x, other.get_point(*x));
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
        for (x, y) in other.curve.iter() {
            new_pairs.insert(*x, *y);
        }
        self.curve = new_pairs;
    }

    fn from_existing(&mut self, other: &Curve, y_units: f64) {
        self.curve = other.curve.clone();

        if y_units != 1.0 {
            for (_, y) in self.curve.iter_mut() {
                *y *= y_units;
            }
        }
    }

    fn clear(&mut self) {
        self.curve.clear();
        self.oob_left = 0.0;
        self.oob_right = 0.0;
    }

    fn debug(&self) {
        for (x, y) in self.curve.iter() {
            print!("{} = {}, ", x, y)
        }
        println!();
    }

    fn load_curve(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b',').from_path(path)?;
        let mut btree: BTreeMap<NotNan<f64>, f64> = BTreeMap::new();
        for result in rdr.records() {
            let record = result?;

            // Enserue CSV has 2 columns
            if record.len() != 2 { 
                return Err(format!("Invalid row format").into());
            }

            //  Show an error if a value is not a valid number
            let x: f64 = record[0].parse::<f64>().map_err(|_| format!("Invalid number"))?;
            let y: f64 = record[1].parse::<f64>().map_err(|_| format!("Invalid number"))?;

            let x_notnan = NotNan::new(x)?;//.map_err(|_| "NaN value encountered")?; -> not necessary with the previous error
            btree.insert(x_notnan, y);
        }
        self.curve = btree.clone();

        Ok(())
    }  

    fn load_slice_dat(&mut self, path: &str, row_index: usize) -> Result<(), Box<dyn Error>> {
        let mut rdr = ReaderBuilder::new().has_headers(true).delimiter(b',').from_path(path)?;
        let headers = rdr.headers()?.clone(); 
        let mut btree: BTreeMap<NotNan<f64>, f64> = BTreeMap::new();
        
        for (idx, result) in rdr.records().enumerate() {
            let record = result?;
            if idx == row_index { 
                for (i, header) in headers.iter().enumerate() { 
                    let x: f64 = header.parse().map_err(|_| format!("Invalid x value: {}", header))?;
                    let y: f64 = record[i].parse().map_err(|_| format!("Invalid y value at column {}", i))?;
    
                    let x_notnan = NotNan::new(x)?;
                    btree.insert(x_notnan, y);
                }
                self.curve = btree;
                return Ok(()); 
            }
        }
        Err(format!("Row index {} out of bounds", row_index).into())
    }

    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>> {
        let mut wtr = WriterBuilder::new().has_headers(false).from_path(path).unwrap();
        //wtr.write_record(&["X", "Y"])?; -> uncomment if headers are necessary
        for (x, y) in self.curve.iter() {
            wtr.write_record(&[x.to_string(), y.to_string()]).unwrap();
        }
        wtr.flush().unwrap(); 
        Ok(())
    }
}

impl BinaryCurveOperations<f64> for Curve {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: f64) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.iter() {
                    new_pairs.insert(*x * by_what, *y);
                }
                self.curve = new_pairs;
            }
            CurveAxis::YAxis => {
                for (_, y) in self.curve.iter_mut() {
                    *y *= by_what;
                }
            } 
        }

        self.oob_left *= by_what;
        self.oob_right *= by_what;
    } 

    fn multiply(&mut self, by_what: f64) {
        let mut new_pairs = BTreeMap::new();
        for (x, _) in self.curve.iter() {
            new_pairs.insert(*x, self.get_point(*x) * by_what);
        }
        self.curve = new_pairs;

        self.oob_left *= by_what;
        self.oob_right *= by_what;
    }

    fn add(&mut self, what: f64) {
        let mut new_pairs = BTreeMap::new();
        for (x, _) in self.curve.iter() {
            new_pairs.insert(*x, self.get_point(*x) + what);
        }
        self.curve = new_pairs;

        self.oob_left += what;
        self.oob_right += what;
    }
}

impl BinaryCurveOperations<Curve> for Curve {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: Curve) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.iter() {
                    new_pairs.insert(*x * by_what.get_point(*x), *y);
                }
                self.curve = new_pairs;
            }
            CurveAxis::YAxis => {
                for (x, y) in self.curve.iter_mut() {
                    *y *= by_what.get_point(*x);
                }
            }
        }

        self.oob_left = by_what.oob_left;
        self.oob_right = by_what.oob_right;
    } 

    fn multiply(&mut self, by_what: Curve) {
        // Unite the keys of both curves
        let mut x_vals: Vec<_> = self.curve.keys().cloned().collect();
        x_vals.extend(by_what.curve.keys().cloned());
        x_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        x_vals.dedup();
        // Iterate on them
        let mut new_pairs = BTreeMap::new();
        for x in x_vals.iter() {
            new_pairs.insert(*x, self.get_point(*x) * by_what.get_point(*x));
        }
        self.curve = new_pairs;

        self.oob_left *= by_what.oob_left;
        self.oob_right *= by_what.oob_right;
    }

    fn add(&mut self, what: Curve) {
        // Unite the keys of both curves
        let mut x_vals: Vec<_> = self.curve.keys().cloned().collect();
        x_vals.extend(what.curve.keys().cloned());
        x_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        x_vals.dedup();
        // Iterate on them
        let mut new_pairs = BTreeMap::new();
        for x in x_vals.iter() {
            new_pairs.insert(*x, self.get_point(*x) + what.get_point(*x));
        }
        self.curve = new_pairs;

        self.oob_left += what.oob_left;
        self.oob_right += what.oob_right;
    }
}

