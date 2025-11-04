//
// spectrum/mod.cpp: Generic spectrum curve
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

use crate::curve;

use std::error::Error;
use ordered_float::NotNan;
use std::collections::BTreeMap;
use crate::curve::Curve;
use crate::curve::CurveAxis;
use curve::CurveOperations;
use curve::BinaryCurveOperations;
use serde::{Serialize, Deserialize};

/// A generic spectrum structure built upon a Curve structure.
/// The "spectrum" encapsulates a "Curve" and provides additional operations for spectra, 
/// such as scaling by a factor or another curve, and axis inversion consistent with physical units.
//
/// Spectra are commonly used to represent distributions (e.g., intensity, energy)
/// across a domain such as wavelength or frequency. This struct enables both
/// general curve manipulation and domain-specific spectral operations.
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct Spectrum {
    /// The underlying curve representing the spectrum data.
    curve: Curve,
}

/// Accessor methods for the Spectrum struct.
impl Spectrum {
    /// Mutable access to the underlying curve.
    pub fn get_curve_mut(&mut self) -> &mut Curve {
        &mut self.curve
    }
    /// Immutable access to the underlying curve.
    pub fn get_curve(&self) -> Curve {
        self.curve.clone()
    }
}

/// Trait defining operations specific to spectrum axes.
pub trait SpectrumAxisOperations {
    /// Scales the specified axis by a scalar factor.
    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: f64);
    /// Scales the specified axis using another curve.
    fn scale_axis_curve(&mut self, axis: CurveAxis, other: Curve);
    /// Scales the specified axis using another curve and its derivative.
    fn scale_axis_curve_diff(&mut self, axis: CurveAxis, other: &Curve, diff: &Curve);
    /// Inverts the specified axis using a factor.
    fn invert_axis_spec(&mut self, axis: CurveAxis, factor: NotNan<f64>);
}

/// Implementing CurveOperations for Spectrum by delegating to the inner Curve.
impl CurveOperations for Spectrum {
    fn is_oob(&self, x: NotNan<f64>) -> bool {
        self.curve.is_oob(x)
    }
    fn bounds(&self) -> Option<(NotNan<f64>, NotNan<f64>)> {
        self.curve.bounds()
    }
    fn set_units(&mut self, axis: CurveAxis, units: String) {
        self.curve.set_units(axis, units);
    }
    fn integral(&self) -> f64 {
        self.curve.integral()
    }
    fn dist_mean(&self) -> f64 {
        self.curve.dist_mean()
    }
    fn integrate(&mut self, k: f64) {
        self.curve.integrate(k);
    }
    fn flip(&mut self) {
        self.curve.flip();
    }
    fn extend_left(&mut self) -> Result<(), String> {
        self.curve.extend_left()
    }
    fn extend_right(&mut self) -> Result<(), String> {
        self.curve.extend_right()
    }
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        self.curve.invert_axis(axis, factor);
    }
    fn x_points(&self) -> Vec<NotNan<f64>> {
        self.curve.x_points()
    }
    fn set(&mut self, x: NotNan<f64>, y: f64) {
        self.curve.set(x, y);
    }
    fn get_point(&self, x: NotNan<f64>) -> f64 {
        self.curve.get_point(x)
    }
    fn get_diff(&self, x: NotNan<f64>) -> f64 {
        self.curve.get_diff(x)
    }
    fn assign(&mut self, other: &Curve) {
        self.curve.assign(other);
    }
    fn from_existing(&mut self, other: &Curve, y_units: f64) {
        self.curve.from_existing(other, y_units);
    }
    fn clear(&mut self) {
        self.curve.clear();
    }
    fn debug(&self) {
        self.curve.debug();
    }
    fn load_curve(&mut self, path: &str) -> Result<(), Box<dyn Error>> {
        self.curve.load_curve(path)
    }
    fn load_slice_dat(&mut self, path: &str, row_index: usize) -> Result<(), Box<dyn Error>> {
        self.curve.load_slice_dat(path, row_index)
    }
    fn save_curve(&self, path: &str) -> Result<(), Box<dyn Error>> {
        self.curve.save_curve(path)
    }
}

/// Implementing BinaryCurveOperations for Spectrum by delegating to the inner Curve for the scalar case.
impl BinaryCurveOperations<f64> for Spectrum {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: f64) {
        self.curve.scale_axis(axis, by_what)
    }
    fn multiply(&mut self, by_what: f64){
        self.curve.multiply(by_what)
    }
    fn add(&mut self, what: f64) {
        self.curve.add(what)

    }
}

/// Implementing BinaryCurveOperations for Spectrum by delegating to the inner Curve for the curve case.
impl BinaryCurveOperations<&Curve> for Spectrum {
    fn scale_axis(&mut self, axis: CurveAxis, by_what: &Curve) {
        self.curve.scale_axis(axis, by_what)
    }
    fn multiply(&mut self, by_what: &Curve){
        self.curve.multiply(by_what)
    }
    fn add(&mut self, what: &Curve) {
        self.curve.add(what)

    }
}

/// Implementing SpectrumAxisOperations for Spectrum with specific logic for spectral data.
impl SpectrumAxisOperations for Spectrum {
    /// Scale the specified axis by a scalar factor.
    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: f64) {
        match axis {
            // Scaling the X axis by a factor involves multiplying X values and dividing Y values.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to hold the scaled points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over existing points, scaling X and Y accordingly.
                for (x, y) in self.curve.get_map_mut().iter() {
                    // Insert the scaled point into the new map.
                    new_pairs.insert(*x * factor, *y / factor);
                }
                // Replace the old map with the new scaled map.
                *self.curve.get_map_mut() = new_pairs;
                *self.curve.get_oob_left_mut() /= factor;
                *self.curve.get_oob_right_mut() /= factor;
            }
            // Scaling the Y axis by a factor involves multiplying Y values directly.
            CurveAxis::YAxis => self.curve.multiply(factor)
        }
    } 

    /// Scale the specified axis using another curve.
    fn scale_axis_curve(&mut self, axis: CurveAxis, mut other: Curve) {
        match axis {
            // Scaling the X axis by another curve involves adjusting Y values based on the other curve's values.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to hold the adjusted points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over existing points, adjusting Y based on the other curve.
                for (x, y) in self.curve.get_map_mut().iter_mut() {
                    // Calculate the derivative of the other curve at the current X value.
                    let diff = other.clone().get_diff(*x).abs();
                    // Only adjust if the derivative is non-zero to avoid division by zero.
                    if diff != 0.0 {
                        // Insert the adjusted point into the new map.
                        let x_notnan = NotNan::new(other.get_point(*x)).expect("x should not be NaN");
                        new_pairs.insert(x_notnan, *y / diff);
                    }
                }

                // Replace the old map with the new adjusted map.
                *self.curve.get_map_mut() = new_pairs;

                // Adjust out-of-bounds values based on the other curve's endpoints.
                let Some((&_, &crv_first_y)) = other.get_map_mut().first_key_value() else { todo!() };
                let Some((&_, &crv_last_y)) = other.get_map_mut().last_key_value() else { todo!() };
                if *self.curve.get_oob_left_mut() != 0.0 {
                    *self.curve.get_oob_left_mut() /= other.get_diff(NotNan::new(crv_first_y).expect("Attempted to insert NaN"))
                }
                if *self.curve.get_oob_right_mut() != 0.0 {
                    *self.curve.get_oob_right_mut() /= other.get_diff(NotNan::new(crv_last_y).expect("Attempted to insert NaN"))
                }
            }
            // Scaling the Y axis by another curve involves multiplying Y values directly.
            CurveAxis::YAxis => self.curve.multiply(&other)
        }
    }

    /// Scale the specified axis using another curve and its derivative.
    fn scale_axis_curve_diff(&mut self, axis: CurveAxis, other: &Curve, diff: &Curve) {
        match axis {
            // Scaling the X axis by another curve and its derivative involves adjusting Y values based on the derivative.
            CurveAxis::XAxis => {
                // Create a new BTreeMap to hold the adjusted points.
                let mut new_pairs = BTreeMap::new();
                // Iterate over existing points, adjusting Y based on the derivative of the other curve.
                for (x, y) in self.curve.get_map_mut().iter() {
                    let dfdx = diff.get_point(*x).abs();
                    // Only adjust if the derivative is non-zero to avoid division by zero.
                    if dfdx != 0.0 {
                        // Insert the adjusted point into the new map.
                        let x_notnan = NotNan::new(other.get_point(*x)).expect("x should not be NaN");
                        new_pairs.insert(x_notnan, y / dfdx);
                    }
                }
                // Replace the old map with the new adjusted map.
                *self.curve.get_map_mut() = new_pairs;

                // Adjust out-of-bounds values based on the derivative of the other curve's endpoints.
                let Some((&_, &crv_first_y)) = other.clone().get_map().first_key_value() else { todo!() };
                let Some((&_, &crv_last_y)) = other.clone().get_map().last_key_value() else { todo!() };
                if *self.curve.get_oob_left_mut() != 0.0 {
                    *self.curve.get_oob_left_mut() /= diff.get_point(NotNan::new(crv_first_y).expect("Attempted to insert NaN")).abs()
                }
                if *self.curve.get_oob_right_mut() != 0.0 {
                    *self.curve.get_oob_right_mut() /= diff.get_point(NotNan::new(crv_last_y).expect("Attempted to insert NaN")).abs()
                }
            }
            // Scaling the Y axis by another curve involves multiplying Y values directly.
            CurveAxis::YAxis => self.curve.multiply(other)
        }
    }

    /// Invert the specified axis using a factor.
    fn invert_axis_spec(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        // Check if the curve is empty.
        if self.curve.get_map_mut().is_empty() {
            return;
        } else {
            match axis {
                // Inverting the X axis involves transforming X values and adjusting Y values accordingly.
                CurveAxis::XAxis => {
                    // Ensure there are no negative X values, as inversion is not defined for them.
                    let Some((&first_x, &_)) = self.curve.get_map_mut().first_key_value() else { todo!() };
                    // If the first X value is negative, print a warning and skip inversion.
                    if first_x < NotNan::new(0.0).expect("Attempted to insert NaN") {
                        eprintln!("Inverting spectrums with negative values in the X axis not yet supported")
                    } else {
                        // Create a new BTreeMap to hold the inverted points.
                        let mut new_pairs = BTreeMap::new();
                        // Iterate over existing points, inverting X and adjusting Y accordingly.
                        for (x, y) in self.curve.get_map_mut().iter_mut() {
                            // Calculate the new inverted X value.
                            let new_x = NotNan::new(*factor / x.into_inner());
                            // Adjust Y based on the inversion formula.
                            *y *= *((*x * *x) / *factor);
                            // Insert the inverted point into the new map.
                            new_pairs.insert(new_x.expect("Attempted to insert NaN"), *y);
                        }
                        // Replace the old map with the new inverted map.
                        *self.curve.get_map_mut() = new_pairs;
                    }
                }
                // Inverting the Y axis involves transforming Y values directly.
                CurveAxis::YAxis => {
                    for (_, y) in self.curve.get_map_mut().iter_mut() {
                        *y = *factor / *y;
                    }
                }
            }
        }
    }
}