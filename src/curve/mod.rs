//
// curve/mod.cpp: Generic physical curve, with units
//
// Copyright (c) 2025 Pablo Álvarez Martín <pablo.alvmar12@gmail.com>
// Copyright (c) 2025 Gonzalo J. Carracedo <BatchDrake@gmail.com>
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
use serde::{Serialize, Deserialize};
use std::mem;


/// A generic curve structure, defined as a set of (x,y) points, with units for
/// each axis, and out-of-bounds (OOB) values for the left and right sides.
/// The curve is stored as a BTreeMap, which keeps the points sorted by X values.
/// NaN values are not allowed in the X axis, as they would break the ordering
/// of the tree. Out-of-bounds values are used when the curve is evaluated outside 
/// the defined range. The units are stored as strings.  The curve supports linear 
/// interpolation between points, and flat extrapolation using the OOB values.
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct Curve {
    /// Out-of-bounds values for the right side.
    oob_right: f64,
    /// Out-of-bounds values for the left side.
    oob_left: f64,
    /// Units for the X axis.
    units_x: String,
    /// Units for the Y axis.
    units_y: String,
    /// The actual curve data as a sorted map, stored as a BTreeMap of (x,y) points. 
    curve: BTreeMap<NotNan<f64>,f64>,
}

/// Accessor methods for the Curve struct.
impl Curve {
    /// Mutable access to the curve map.
    pub fn get_map_mut(&mut self) -> &mut BTreeMap<NotNan<f64>, f64> {
        &mut self.curve
    }
    /// Mutable access to OOB values.
    pub fn get_oob_left_mut(&mut self) -> &mut f64 {
        &mut self.oob_left
    }
    /// Mutable access to OOB values.
    pub fn get_oob_right_mut(&mut self) -> &mut f64 {
        &mut self.oob_right
    }
    /// Immutable access to the curve map.
    pub fn get_map(self) -> BTreeMap<NotNan<f64>, f64> {
        self.curve
    }
    /// Immutable access to OOB values.
    pub fn get_oob_left(&self) -> f64 {
        self.oob_left.clone()
    }
    /// Immutable access to OOB values.
    pub fn get_oob_right(&self) -> f64 {
        self.oob_right.clone()
    }
}

/// Axis selector for curve operations.
pub enum CurveAxis {
    /// Select X axis.
    XAxis,
    /// Select Y axis.
    YAxis,
}

/// Trait defining operations that can be performed on a Curve.
pub trait CurveOperations {
    /// Check if a given x value is out of bounds.
    fn is_oob(&self, x: NotNan<f64>) -> bool;
    /// Returns the minimum and maximum X values, if any.
    fn bounds(&self) -> Option<(NotNan<f64>, NotNan<f64>)>;
    /// Set the units for a given axis.
    fn set_units(&mut self, axis: CurveAxis, units: String);
    /// Calculate the integral of the curve.
    fn integral(&self) -> f64;
    /// Computes the mean value of the curve's distribution.
    fn dist_mean(&self) -> f64;
    /// Integrates the curve, starting from a given value f, and modifying it in place.
    fn integrate(&mut self, k: f64);
    /// Flips the curve, swapping X and Y axes and their units.
    fn flip(&mut self);
    /// Extends the curve to the left using the first defined point.
    fn extend_left(&mut self) -> Result<(), String>;
    /// Extends the curve to the right using the last defined point.
    fn extend_right(&mut self) -> Result<(), String>;
    /// Inverts the values of a given axis, multiplying by a factor.
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>);
    /// Returns a vector of all X points in the curve.
    fn x_points(&self) -> Vec<NotNan<f64>>;
    /// Sets a point in the curve, adding it if it doesn't exist.
    fn set(&mut self, x: NotNan<f64>, y: f64);
    /// Gets the Y value for a given X, interpolating if necessary.
    fn get_point(&self, x: NotNan<f64>) -> f64;
    /// Gets the derivative (slope) of the curve at a given X.
    fn get_diff(&self, x: NotNan<f64>) -> f64;
    /// Assigns values from another curve, interpolating as necessary.
    fn assign(&mut self, other: &Curve);
    /// Copies values from another curve, scaling Y values by a factor.
    fn from_existing(&mut self, other: &Curve, y_units: f64);
    /// Clears the curve, removing all points and resetting OOB values.
    fn clear(&mut self);
    /// Prints debug information about the curve.
    fn debug(&self);
    /// Loads curve data from a CSV file.
    fn load_curve(&mut self, path: &str) -> Result<(), Box<dyn Error>>;
    /// Loads curve data from a specific row in a CSV file with headers. Only works for CSV files with slice format.
    fn load_slice_dat(&mut self, path: &str, row_index: usize) -> Result<(), Box<dyn Error>>;
    /// Saves the curve data to a CSV file.
    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>>;
}

/// Trait defining binary operations that can be performed on a Curve: 
/// Operations that involve two curves or a curve and a scalar. These
/// operations modify the curve in place.                   
pub trait BinaryCurveOperations<T> {
    /// Scales a given axis by a scalar or another curve.
    fn scale_axis(&mut self, axis: CurveAxis, by_what: T);
    /// Multiplies the Y values of the curve by a scalar or another curve.
    fn multiply(&mut self, by_what: T);
    /// Adds a scalar or another curve to the Y values of the curve.
    fn add(&mut self, what: T);
}

/// Implementation of CurveOperations trait for the Curve struct.
impl CurveOperations for Curve {
    /// Returns the minimum and maximum X values, if any.
    fn bounds(&self) -> Option<(NotNan<f64>, NotNan<f64>)> {
        // Check if the curve is empty.
        if self.curve.is_empty() {
            return None;
        }

        // Get the first and last keys (X values) from the BTreeMap.
        let first = *self.curve.first_key_value().unwrap().0;
        let last  = *self.curve.last_key_value().unwrap().0;

        // Return the bounds as a tuple. 
        return Some((first, last));
    }

    /// Check if a given x value is out of bounds.
    fn is_oob(&self, x: NotNan<f64>) -> bool {
        (match self.curve.first_key_value() {
            Some(first_x) => x < *first_x.0,
            None          => true,
        }) || (match self.curve.last_key_value() {
            Some(last_x)  => *last_x.0 < x,
            None          => true,
        })
    }

    /// Set the units for a given axis.
    fn set_units(&mut self, axis: CurveAxis, units: String) {
        // Set the units for the specified axis.
        match axis {
            CurveAxis::XAxis => self.units_x = units,
            CurveAxis::YAxis => self.units_y = units,
        };
    }

    /// Calculate the integral of the curve.
    fn integral(&self) -> f64 {  
        // Check if urve is out of bounds   
        if self.oob_right + self.oob_left != 0.0 {
            // If so, return infinity
            return f64::INFINITY * (self.oob_left + self.oob_right);
        // If the curve is empty, return 0
        } else if self.curve.is_empty() {
            return 0.0;
        } else {
            // Initialize variables for Kahan summation.
            let mut accum = 0.0;
            let mut err = 0.0;

            // Initialize iterator.
            let mut iter = self.curve.iter();
            // Get first point.
            let (&x0, &y0) = iter.next().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            // Iterate over the curve points.
            for (&x, &y) in iter {
                // Calculate the width and average height of the trapezoid.
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);

                // Kahan summation.
                let term_real = my * *dx - err;
                let tmp       = accum + term_real;
                let term_err  = tmp - accum;

                // Update error and accumulator.
                err = term_err - term_real;
                accum = tmp;

                // Update previous point.
                x_prev = x;
                y_prev = y;
            }
            // Return the accumulated integral value.
            return accum;
        }
    }

    /// Computes the mean value of the curve's distribution.
    fn dist_mean(&self) -> f64 {
        // Check if curve is out of bounds
        if self.oob_right + self.oob_left != 0.0 {
            // If so, return the mean of the OOB values
            return 0.5 * (self.oob_left + self.oob_right);
        // If the curve is empty, return 0
        } else if self.curve.is_empty() {
            return 0.0;
        } else {
            // Initialize variables for Kahan summation.
            let mut accum = 0.0;
            let mut err = 0.0;

            // Initialize iterator.
            let mut iter = self.curve.iter();
            // Get first point.
            let (&x0, &y0) = iter.next().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            // Iterate over the curve points.
            for (&x, &y) in iter {
                // Calculate the width and average height of the trapezoid.
                let dx = x - x_prev;
                let my = 0.5 * (y + y_prev);
                // Average x value in the trapezoid (for mean calculation).
                let x_0 = 0.5 * *((x + x_prev));

                // Kahan summation.
                let term_real = x_0 * my * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;

                // Update error and accumulator.
                err = term_err - term_real;
                accum = tmp;

                // Update previous point.
                x_prev = x;
                y_prev = y;
            }
            // Return the mean value (integral divided by total integral).
            return accum/self.integral();
        }
    }

    /// Integrates the curve, starting from a given value k, and modifying it in place.
    fn integrate(&mut self, k: f64) { 
        // Initialize variables for Kahan summation.
        let mut accum = 0.0;
        let mut err = 0.0;

        // Set left OOB value to k.
        self.oob_left = k;

        // Handle special cases for empty or single-point curves.
        if self.curve.is_empty() {
            // If the curve is empty, set right OOB to k.
            self.oob_right = k;
        } else if self.curve.len() == 1 {
            // If the curve has only one point, set that point's Y to k and right OOB to k.
            let x = *self.curve.keys().next().unwrap();
            self.curve.clear();
            self.curve.insert(x, k);
            self.oob_right = k;
        } else {
            // Initialize iterator.
            let mut iter = self.curve.iter();
            
            // Get first point. 
            let (&x0, &y0) = iter.next().unwrap();
            let mut x_prev = x0;
            let mut y_prev = y0;

            // Create a new BTreeMap to store the integrated points.
            let mut new_pairs = BTreeMap::new();
            // Initialize dx variable.
            let mut dx: NotNan<f64> = NotNan::new(0.0).unwrap();

            // Iterate over the curve points.
            for (&x, &y) in iter {
                // Calculate the width of the trapezoid.
                dx = x - x_prev;

                // Kahan summation.
                let term_real = 0.5 * (y + y_prev) * *dx - err;
                let tmp = accum + term_real;
                let term_err = tmp - accum;
        
                // Update error and accumulator.
                err = term_err - term_real;
                accum = tmp;
                
                // Insert the new point into the new BTreeMap.
                new_pairs.insert(x, accum);

                // Update previous point.
                x_prev = x;
                y_prev = y;
            }       
            // Handle the last point separately to ensure it's included.     
            new_pairs.insert(x0, k);
            
            // Update the curve with the new integrated points and right OOB value with the final accumulated value.
            self.curve = new_pairs;
            self.oob_right = accum;
        }
    }

    /// Flips the curve, swapping X and Y axes and their units.
    fn flip(&mut self) {
        // Swap units of X and Y axes.
        std::mem::swap(&mut self.units_x, &mut self.units_y);

        // Create a new BTreeMap to store the flipped points.
        let mut flip_crv = BTreeMap::new();
        // Iterate over the current curve and insert flipped points into the new BTreeMap.
        for (&x, &y) in &self.curve {
            let y_notnan = NotNan::new(y).expect("flip(): NaN values in the X axis are forbidden");
            flip_crv.insert(y_notnan, *x);
        }
        // Update the curve with the flipped points.
        self.curve = flip_crv;

        // Swap out-of-bounds values.
        self.oob_left = *self.curve.first_key_value().unwrap().1;
        self.oob_right = *self.curve.last_key_value().unwrap().1;
    }

    /// Extend the curve to the left using the first defined point.
    fn extend_left(&mut self) -> Result<(), String> {
        match self.bounds() {
            // If the curve is empty, return an error.
            None               => return Err(format!("extend_left(): curve is empty")),
            // Set the left OOB value to the Y value of the first point.
            Some((first_x, _)) => self.oob_left = *self.curve.first_key_value().unwrap().1,
        };
        Ok(())
    }

    /// Extend the curve to the roght using the last defined point.
    fn extend_right(&mut self) -> Result<(), String> {
        match self.bounds() {
            // If the curve is empty, return an error.
            None              => return Err(format!("extend_right(): curve is empty")),
            // Set the right OOB value to the Y value of the last point.
            Some((_, last_x)) => self.oob_right = *self.curve.last_key_value().unwrap().1,
        };
        Ok(())
    }

    /// Inverts the values of a given axis, multiplying by a factor.
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        //
        // While we can alter the values of the Y axis in-place, this is not
        // the case for the X axis, as this would break the ordering of the
        // underlaying tree. We have to instantiate a new BTreeMap, transferring
        // and converting the old values as new values.
        //
        match axis {
            // Invert X axis values.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to store the inverted points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over the current curve and insert inverted points into the new BTreeMap.
                for (x, y) in self.curve.iter() {
                    // Calculate the new X value by inverting the current X value.
                    let new_x = NotNan::new(*factor / x.into_inner());
                    // Insert the new point into the new BTreeMap.
                    new_pairs.insert(
                        new_x.expect(
                        "invert_axis(): NaN values in the X axis are forbidden"),
                        *y);
                }
                // Update the curve with the new inverted points.
                self.curve = new_pairs;

                // Swap out-of-bounds values.
                mem::swap(&mut self.oob_left, &mut self.oob_right);
            }
            // Invert Y axis values.
            CurveAxis::YAxis => {
                for (_, y) in self.curve.iter_mut() {
                    *y = *factor / *y;
                }
                // Invert out-of-bounds values.
                self.oob_left  = 1.0 / self.oob_left;
                self.oob_right = 1.0 / self.oob_right;
            }
        }
    }

    /// Returns a vector of all X points in the curve.
    fn x_points(&self) -> Vec<NotNan<f64>> {
        // Collect all X values from the curve into a vector and return it.
        self.curve.clone().into_keys().collect()
    } 

    /// Sets a point in the curve, adding it if it doesn't exist.
    fn set(&mut self, x: NotNan<f64>, y: f64) {
        // Insert or update the point (x, y) in the curve.
        self.curve.insert(x, y);
    }

    /// Gets the Y value for a given X, interpolating if necessary.
    fn get_point(&self, x: NotNan<f64>) -> f64 {
        // Check if the curve is empty and panic if so.
        let Some((first_x, last_x)) = self.bounds() else {
            panic!("get_point(x): curve is empty");
        };

        // Handle out-of-bounds cases.
        if x < first_x {
            self.oob_left
        } else if x > last_x {
            self.oob_right
        } else {
            // Interpolate or return exact value if it exists.
            match self.curve.get(&x) {
                // Exact match found, return the corresponding Y value.
                Some(y) => *y,
                // No exact match, perform linear interpolation.
                None => {
                    let (&x0, &y0) = self.curve.range(..x).next_back().unwrap();
                    let (&x1, &y1) = self.curve.range(x..).next().unwrap();
                    if x1 == x0 {
                    }
                    y0 + *((x - x0) * (y1 - y0) / (x1 - x0))
                }
            }
        }
    }

    /// Gets the derivative (slope) of the curve at a given X.
    fn get_diff(&self, x: NotNan<f64>) -> f64 {
        // There are several cases here:
        // 1. Diffing out of bounds returns zero
        // 2. Diffing in bounds depends whether we are in the edge between
        //    two segments or not.

        // Check if the curve is empty and panic if so.
        let Some((first_x, last_x)) = self.bounds() else {
            panic!("get_point(x): curve is empty");
        };

        // Handle out-of-bounds cases.
        if x < first_x || x >= last_x {
            return 0.0;
        } else {
            // This is guaranteed to exist, as we are in bounds.
            // Interpolate or return exact value if it exists.
            // We need to check if we are at the edge of a segment.
            // If so, we need to use the next segment for the diff.
            // If not, we can use the current segment. 
            // If we are at the edge of the last segment, we return zero.
            let mut left   = self.curve.range(..x);
            let mut right  = self.curve.range(x..);

            // Get the last point on the left side and the first point on the right side.
            let Some((&x0, &y0)) = left.next_back() else {
                return 0.0; // No valid left point
            };
            let Some((&x1, &y1)) = right.next() else {
                return 0.0; // No valid right point
            };

            if x != x1 || x1 == last_x {
                // Check if we are at the edge of a segment.
                // Not at the edge of a segment, or at the edge of the last segment
                if x1 != x0 {
                    // Normal case, return the diff between the two points.
                    return (y1 - y0) / *((x1 - x0));
                }
            } else {
                // Edge of a segment. This is also guaranteed to exist, as
                // the previous condition is also applied for the last segment.
                if let Some((&x2, &y2)) = right.next() {
                    // We are at the edge of a segment, use the next segment for the diff.
                    if x2 != x0 {
                        return (y2 - y0) / *((x2 - x0));
                    }
                }
            }
        }
        // Fallback return value to avoid panic
        return 0.0;
    }

    /// Assigns values from another curve, interpolating as necessary.
    fn assign(&mut self, other: &Curve) {
        // Get bounds of both curves
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

                // Get bounds again, they must exist now as both curves are non-empty.
                let (self_first_x, self_last_x)   = self_bounds.unwrap();
                let (other_first_x, other_last_x) = other_bounds.unwrap();

                // Curve longer to the left.
                if other_first_x < self_first_x {
                    // assign left OOB value.
                    self.oob_left = other.oob_left;
                }

                // Curve longer to the right.
                if self_last_x < other_last_x {
                    // assign right OOB value.
                    self.oob_right = other.oob_right;
                }

                // Assign everything else.
                for (x, y) in other.curve.iter() {
                    self.curve.insert(*x, *y);
                }
            }
        }
    }

    /// Copies values from another curve, scaling Y values by a factor.
    fn from_existing(&mut self, other: &Curve, y_units: f64) {
        // Directly copy the curve and OOB values.
        self.curve = other.curve.clone();
        self.oob_left  = other.oob_left;
        self.oob_right = other.oob_right;

        // Scale Y values if necessary.
        if y_units != 1.0 {
            // Scale OOB values.
            self.oob_left  *= y_units;
            self.oob_right *= y_units;

            // Scale Y values in the curve.
            for (_, y) in self.curve.iter_mut() {
                *y *= y_units;
            }
        }
    }

    /// Clears the curve, removing all points and resetting OOB values.
    fn clear(&mut self) {
        // Clear the curve and reset OOB values.
        self.curve.clear();
        self.oob_left = 0.0;
        self.oob_right = 0.0;
    }

    /// Prints debug information about the curve.
    fn debug(&self) {
        // Print debug information about the curve.
        println!("Min X value: {:?}; Max X value: {:?}", self.curve.first_key_value().map(|(k, _)| k), self.curve.last_key_value().map(|(k, _)| k));
        println!("Left bound: {}; Right bound: {}", self.oob_left, self.oob_right);
        println!();
    }

    /// Loads curve data from a CSV file.
    fn load_curve(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        // Create a CSV reader with no headers and comma delimiter.
        let mut rdr = ReaderBuilder::new().has_headers(false).delimiter(b',').from_path(path)?;
        // Create a new BTreeMap to store the loaded points.
        let mut btree: BTreeMap<NotNan<f64>, f64> = BTreeMap::new();
        // Iterate over the records in the CSV file.
        for result in rdr.records() {
            let record = result?;

            // Enserue CSV has 2 columns
            if record.len() != 2 { 
                return Err(format!("Invalid row format").into());
            }

            //  Show an error if a value is not a valid number
            let x: f64 = record[0].parse::<f64>().map_err(|_| format!("Invalid number"))?;
            let y: f64 = record[1].parse::<f64>().map_err(|_| format!("Invalid number"))?;

            // Insert the point into the BTreeMap, ensuring no NaN values in X axis.
            let x_notnan = NotNan::new(x)?;
            btree.insert(x_notnan, y);
        }
        // Update the curve with the loaded points.
        self.curve = btree.clone();

        Ok(())
    }  

    /// Loads curve data from a specific row in a CSV file with headers. Only works for CSV files with slice format.
    fn load_slice_dat(&mut self, path: &str, row_index: usize) -> Result<(), Box<dyn Error>> {
        // Create a CSV reader with headers and comma delimiter.
        let mut rdr = ReaderBuilder::new().has_headers(true).delimiter(b',').from_path(path)?;
        // Get headers to use as X values.
        let headers = rdr.headers()?.clone(); 
        // Create a new BTreeMap to store the loaded points.
        let mut btree: BTreeMap<NotNan<f64>, f64> = BTreeMap::new();
        
        // Iterate over the records in the CSV file.
        for (idx, result) in rdr.records().enumerate() {
            let record = result?;
            // Ensure the record has the same number of columns as headers.
            if idx == row_index { 
                for (i, header) in headers.iter().enumerate() { 
                    let x: f64 = header.parse().map_err(|_| format!("Invalid x value: {}", header))?;
                    let y: f64 = record[i].parse().map_err(|_| format!("Invalid y value at column {}", i))?;
    
                    // Insert the point into the BTreeMap, ensuring no NaN values in X axis.
                    let x_notnan = NotNan::new(x)?;
                    btree.insert(x_notnan, y);
                }
                // Update the curve with the loaded points.
                self.curve = btree;
                return Ok(()); 
            }
        }
        // If we reach here, the specified row index was not found.
        Err(format!("Row index {} out of bounds", row_index).into())
    }

    /// Saves the curve data to a CSV file.
    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>> {
        // Create a CSV writer with no headers.
        let mut wtr = WriterBuilder::new().has_headers(false).from_path(path).unwrap();
        //wtr.write_record(&["X", "Y"])?; -> UNCOMMENT IF HEADERS ARE NECESSARY 
        for (x, y) in self.curve.iter() {
            // Write each point as a record in the CSV file.
            wtr.write_record(&[x.to_string(), y.to_string()]).unwrap();
        }
        // Flush the writer to ensure all data is written to the file.
        wtr.flush().unwrap(); 
        Ok(())
    }
}

/// Implementation of BinaryCurveOperations trait for the Curve struct with f64 as the type.
impl BinaryCurveOperations<f64> for Curve {
    /// Scales a given axis by a scalar.
    fn scale_axis(&mut self, axis: CurveAxis, by_what: f64) {
        match axis {
            // Scale X axis values.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to store the scaled points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over the current curve and insert scaled points into the new BTreeMap.
                for (x, y) in self.curve.iter() {
                    // Calculate the new X value by scaling the current X value.
                    new_pairs.insert(*x * by_what, *y);
                }
                // Update the curve with the new scaled points.
                self.curve = new_pairs;
            }
            // Scale Y axis values by multiplying by the scalar.
            CurveAxis::YAxis => self.multiply(by_what),
        }

        
    } 

    /// Multiplies the Y values of the curve by a scalar.
    fn multiply(&mut self, by_what: f64) {
        // Scale Y values in the curve.
        for (_, y) in self.curve.iter_mut() {
            // Multiply each Y value by the given factor.
            *y *= by_what;
        }

        // Scale out-of-bounds values.
        self.oob_left  *= by_what;
        self.oob_right *= by_what;
    }

    /// Adds a scalar to the Y values of the curve.
    fn add(&mut self, what: f64) {
        // Add the given value to Y values in the curve.
        for (_, y) in self.curve.iter_mut() {
            // Add the given value to each Y value.
            *y += what;
        }

        // Add the given value to the out-of-bounds values.
        self.oob_left  += what;
        self.oob_right += what;
    }
}

/// Helper function to get the union of X values from two curves.
fn x_union(curve_a: &Curve, curve_b: &Curve) -> HashSet<NotNan<f64>> {
    /// Create a HashSet to store the union of X values.
    let mut set: HashSet<NotNan<f64>> = HashSet::from_iter(
        // Insert X values from the first curve.
        curve_a.curve.keys().cloned());

    /// Insert X values from the second curve.
    for (x, _) in curve_b.curve.iter() {
        set.insert(*x);
    }
    // Return the union of X values. 
    set
}

/// Implementation of BinaryCurveOperations trait for the Curve struct with &Curve as the type.
impl BinaryCurveOperations<&Curve> for Curve {
    /// Scales a given axis by another curve.
    fn scale_axis(&mut self, axis: CurveAxis, by_what: &Curve) {
        match axis {
            // Scale X axis values.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to store the scaled points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over the current curve and insert scaled points into the new BTreeMap.
                for (x, y) in self.curve.iter() {
                    // Calculate the new X value by scaling the current X value using the other curve.
                    new_pairs.insert(*x * by_what.get_point(*x), *y);
                }
                // Update the curve with the new scaled points.
                self.curve = new_pairs;
            }
            // Scale Y axis values.
            CurveAxis::YAxis => {
                for (x, y) in self.curve.iter_mut() {
                    // Multiply each Y value by the corresponding value from the other curve.
                    *y *= by_what.get_point(*x);
                }
            }
        }
    } 

    /// Multiplies the Y values of the curve by another curve.
    fn multiply(&mut self, by_what: &Curve) {
        // Get the union of X values from both curves.
        let x_vals = x_union(self, by_what);
        // Create a new BTreeMap to store the multiplied points.
        let mut new_pairs = BTreeMap::new();
        
        // Iterate over the union of X values and multiply corresponding Y values.
        for x in x_vals.iter() {
            // Multiply the Y values from both curves at the given X value.
            new_pairs.insert(*x, self.get_point(*x) * by_what.get_point(*x));
        }
        // Update the curve with the new multiplied points.
        self.curve = new_pairs;

        // Scale out-of-bounds values.
        self.oob_left  *= by_what.oob_left;
        self.oob_right *= by_what.oob_right;
    }

    /// Adds another curve to the Y values of the curve.
    fn add(&mut self, what: &Curve) {
        // Get the union of X values from both curves.
        let x_vals = x_union(self, what);
        // Create a new BTreeMap to store the added points.
        let mut new_pairs = BTreeMap::new();

        // Iterate over the union of X values and add corresponding Y values.
        for x in x_vals.iter() {
            // Add the Y values from both curves at the given X value.
            new_pairs.insert(*x, self.get_point(*x) + what.get_point(*x));
        }
        // Update the curve with the new added points.
        self.curve = new_pairs;

        // Add the out-of-bounds values.
        self.oob_left += what.oob_left;
        self.oob_right += what.oob_right;
    }
}

