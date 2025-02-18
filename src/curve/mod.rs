//
// curve/mod.cpp: Generic physical curve, with units
//
// Copyright (c) 2025 Pablo Álvarez Martín <pablo.alvmar12@gmail.com>
// Copyright (c) 2023 Gonzalo J. Carracedo <BatchDrake@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//

use std::error::Error;
use csv::ReaderBuilder;
use csv::WriterBuilder;
use std::collections::BTreeMap;
use std::collections::HashSet;
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

pub trait CurveOperations {
    fn is_oob(&self, x: NotNan<f64>) -> bool;
    fn bounds(&self) -> Option<(NotNan<f64>, NotNan<f64>)>;
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
    fn bounds(&self) -> Option<(NotNan<f64>, NotNan<f64>)> {
        if self.curve.is_empty() {
            return None;
        }

        let first = *self.curve.first_key_value().unwrap().0;
        let last  = *self.curve.last_key_value().unwrap().0;

        return Some((first, last));
    }

    fn is_oob(&self, x: NotNan<f64>) -> bool {
        (match self.curve.first_key_value() {
            Some(first_x) => x < *first_x.0,
            None          => true,
        }) || (match self.curve.last_key_value() {
            Some(last_x)  => *last_x.0 < x,
            None          => true,
        })
    }

    fn set_units(&mut self, axis: CurveAxis, units: String) {
        match axis {
            CurveAxis::XAxis => self.units_x = units,
            CurveAxis::YAxis => self.units_y = units,
        };
    }

    fn integral(&self) -> f64 {  
        // Check if out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return f64::INFINITY * (self.oob_left + self.oob_right);
        } else if self.curve.is_empty() {
            return 0.0;
        } else {
            let mut accum = 0.0;
            let mut err = 0.0;

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
                let tmp       = accum + term_real;
                let term_err  = tmp - accum;

                err = term_err - term_real;
                accum = tmp;

                x_prev = x;
                y_prev = y;
            }
            return accum;
        }
    }

    fn dist_mean(&self) -> f64 {
        // Curve is out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            return 0.5 * (self.oob_left + self.oob_right);
        } else if self.curve.is_empty() {
            return 0.0;
        } else {
            let mut accum = 0.0;
            let mut err = 0.0;

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
        } else if self.curve.len() == 1 {
            let x = *self.curve.keys().next().unwrap();
            self.curve.clear();
            self.curve.insert(x, k);
            self.oob_right = k;
        } else {
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
    }

    fn flip(&mut self) {
        // Unit flip
        std::mem::swap(&mut self.units_x, &mut self.units_y);

        // x & y coordinates flip
        let mut flip_crv = BTreeMap::new();
        for (&x, &y) in &self.curve {
            let y_notnan = NotNan::new(y).expect("flip(): NaN values in the X axis are forbidden");
            flip_crv.insert(y_notnan, *x);
        }
        self.curve = flip_crv;
    }

    fn extend_left(&mut self) {
        match self.bounds() {
            None               => panic!("extend_left(): curve is empty"),
            Some((first_x, _)) => self.oob_left = *first_x
        };
    }

    fn extend_right(&mut self) {
        match self.bounds() {
            None              => panic!("extend_right(): curve is empty"),
            Some((_, last_x)) => self.oob_right = *last_x
        };
    }

    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        //
        // While we can alter the values of the Y axis in-place, this is not
        // the case for the X axis, as this would break the ordering of the
        // underlaying tree. We have to instantiate a new BTreeMap, transferring
        // and converting the old values as new values.
        //
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.iter() {
                    let new_x = NotNan::new(*factor / x.into_inner());
                    new_pairs.insert(
                        new_x.expect(
                        "invert_axis(): NaN values in the X axis are forbidden"),
                        *y);
                }

                self.curve = new_pairs;
            }
            CurveAxis::YAxis => {
                for (_, y) in self.curve.iter_mut() {
                    *y = *factor / *y;
                }
            }
        }

        self.oob_left  = 1.0 / self.oob_left;
        self.oob_right = 1.0 / self.oob_right;
    }

    fn x_points(&self) -> Vec<NotNan<f64>> {
        self.curve.clone().into_keys().collect()
    } 

    fn set(&mut self, x: NotNan<f64>, y: f64) {
        self.curve.insert(x, y);
    }

    fn get_point(&self, x: NotNan<f64>) -> f64 {
        let Some((first_x, last_x)) = self.bounds() else {
            panic!("get_point(x): curve is empty");
        };

        if x < first_x {
            self.oob_left
        } else if x > last_x {
            self.oob_right
        } else {
            match self.curve.get(&x) {
                Some(y) => *y,
                None    => {
                    let (&x0, &y0) = self.curve.range(..x).next_back().unwrap();
                    let (&x1, &y1) = self.curve.range(x..).next().unwrap();
                    y0 + *((x - x0) * (y1 - y0) / (x1 - x0))
                }
            }
        }
    }

    fn get_diff(&self, x: NotNan<f64>) -> f64 {
        // There are several cases here:
        // 1. Diffing out of bounds returns zero
        // 2. Diffing in bounds depends whether we are in the edge between
        //    two segments or not.

        let Some((first_x, last_x)) = self.bounds() else {
            panic!("get_point(x): curve is empty");
        };

        if x < first_x || x >= last_x {
            0.0
        } else {
            // These are guaranteed to exist
            let mut left   = self.curve.range(..x);
            let mut right  = self.curve.range(x..);

            let (&x0, &y0) = left.next_back().unwrap();
            let (&x1, &y1) = right.next().unwrap();

            if x != x1 || x1 == last_x {
                // Middle of the segment, or edge of the last segment
                (y1 - y0) / *((x1 - x0))
            } else {
                // Edge of a segment. This is also guaranteed to exist, as
                // the previous condition is also applied for the last segment.
                let (&x2, &y2) = right.next().unwrap();
                (y2 - y0) / *((x2 - x0))
            }
        }
    }

    // TODO: Check units!!!
    fn assign(&mut self, other: &Curve) {
        let self_bounds  = self.bounds();
        let other_bounds = other.bounds();

        // Check if there's something to assign
        if !other_bounds.is_none() {
            if self_bounds.is_none() {
                // Self is empty, assign directly
                self.curve     = other.curve.clone();
                self.oob_left  = other.oob_left;
                self.oob_right = other.oob_right;
            } else {
                // Self is not empty, check overlaps

                // Assign points that already exist in our curve
                for (x, y) in self.curve.iter_mut() {
                    if !other.is_oob(*x) {
                        *y = other.get_point(*x);
                    }
                }

                let (self_first_x, self_last_x)   = self_bounds.unwrap();
                let (other_first_x, other_last_x) = other_bounds.unwrap();

                // Curve longer to the left
                if other_first_x < self_first_x {
                    self.oob_left = other.oob_left;
                }

                // Curve longer to the right
                if self_last_x < other_last_x {
                    self.oob_right = other.oob_right;
                }

                // Assign everything else
                for (x, y) in other.curve.iter() {
                    self.curve.insert(*x, *y);
                }
            }
        }
    }

    fn from_existing(&mut self, other: &Curve, y_units: f64) {
        self.curve = other.curve.clone();
        self.oob_left  = other.oob_left;
        self.oob_right = other.oob_right;

        if y_units != 1.0 {
            self.oob_left  *= y_units;
            self.oob_right *= y_units;

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
                self.oob_left  *= by_what;
                self.oob_right *= by_what;
            }
            CurveAxis::YAxis => self.multiply(by_what),
        }

        
    } 

    fn multiply(&mut self, by_what: f64) {
        for (_, y) in self.curve.iter_mut() {
            *y *= by_what;
        }

        self.oob_left  *= by_what;
        self.oob_right *= by_what;
    }

    fn add(&mut self, what: f64) {
        for (_, y) in self.curve.iter_mut() {
            *y += what;
        }

        self.oob_left  += what;
        self.oob_right += what;
    }
}

fn x_union(curve_a: &Curve, curve_b: &Curve) -> HashSet<NotNan<f64>> {
    let mut set: HashSet<NotNan<f64>> = HashSet::from_iter(
        curve_a.curve.keys().cloned());

    for (x, _) in curve_b.curve.iter() {
        set.insert(*x);
    }
    set
}

impl BinaryCurveOperations<&Curve> for Curve {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: &Curve) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.iter() {
                    new_pairs.insert(*x * by_what.get_point(*x), *y);
                }
                self.curve = new_pairs;
                self.oob_left = by_what.oob_left;
                self.oob_right = by_what.oob_right;
            }
            CurveAxis::YAxis => {
                for (x, y) in self.curve.iter_mut() {
                    *y *= by_what.get_point(*x);
                }
            }
        }
    } 

    fn multiply(&mut self, by_what: &Curve) {
        let x_vals = x_union(self, by_what);
        for x in x_vals.iter() {
            self.curve.insert(*x, self.get_point(*x) * by_what.get_point(*x));
        }

        self.oob_left  *= by_what.oob_left;
        self.oob_right *= by_what.oob_right;
    }

    fn add(&mut self, what: &Curve) {
        let x_vals = x_union(self, what);

        for x in x_vals.iter() {
            self.curve.insert(*x, self.get_point(*x) + what.get_point(*x));
        }

        self.oob_left += what.oob_left;
        self.oob_right += what.oob_right;
    }
}

