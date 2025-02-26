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

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct Spectrum {
    curve: Curve,
}

impl Spectrum {
    pub fn get_curve_mut(&mut self) -> &mut Curve {
        &mut self.curve
    }
    pub fn get_curve(&self) -> Curve {
        self.curve.clone()
    }
}

pub trait SpectrumAxisOperations {
    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: f64);
    fn scale_axis_curve(&mut self, axis: CurveAxis, other: Curve);
    fn scale_axis_curve_diff(&mut self, axis: CurveAxis, other: Curve, diff: Curve);
    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>);
}

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

impl SpectrumAxisOperations for Spectrum {
    fn scale_axis_factor(&mut self, axis: CurveAxis, factor: f64) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.get_curve_mut().iter() {
                    new_pairs.insert(*x * factor, *y / factor);
                }
                *self.curve.get_curve_mut() = new_pairs;
                *self.curve.get_oob_left_mut() /= factor;
                *self.curve.get_oob_right_mut() /= factor;
            }
            CurveAxis::YAxis => self.get_curve_mut().multiply(factor)
        }
    } 

    fn scale_axis_curve(&mut self, axis: CurveAxis, mut other: Curve) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.get_curve_mut().iter_mut() {
                    let diff = other.get_diff(*x).abs();
                    if diff != 0.0 {
                        let x_notnan = NotNan::new(other.get_point(*x)).expect("x should not be NaN");
                        new_pairs.insert(x_notnan, *y / diff);
                    }
                }

                *self.curve.get_curve_mut() = new_pairs;

                let Some((&_, &crv_first_y)) = other.get_curve_mut().first_key_value() else { todo!() };
                let Some((&_, &crv_last_y)) = other.get_curve_mut().last_key_value() else { todo!() };
                if *self.curve.get_oob_left_mut() != 0.0 {
                    *self.curve.get_oob_left_mut() /= other.get_diff(NotNan::new(crv_first_y).expect("Attempted to insert NaN"))
                }
                if *self.curve.get_oob_right_mut() != 0.0 {
                    *self.curve.get_oob_right_mut() /= other.get_diff(NotNan::new(crv_last_y).expect("Attempted to insert NaN"))
                }
            }
            CurveAxis::YAxis => self.get_curve_mut().multiply(&other)
        }
    }
    fn scale_axis_curve_diff(&mut self, axis: CurveAxis, mut other: Curve, diff: Curve) {
        match axis {
            CurveAxis::XAxis => {
                let mut new_pairs = BTreeMap::new();
                for (x, y) in self.curve.get_curve_mut().iter() {
                    let dfdx = diff.get_diff(*x).abs();
                    if dfdx != 0.0 {
                        let x_notnan = NotNan::new(other.get_point(*x)).expect("x should not be NaN");
                        new_pairs.insert(x_notnan, *y / dfdx);
                    }
                }
                *self.curve.get_curve_mut() = new_pairs;

                let Some((&_, &crv_first_y)) = other.get_curve_mut().first_key_value() else { todo!() };
                let Some((&_, &crv_last_y)) = other.get_curve_mut().last_key_value() else { todo!() };
                if *self.curve.get_oob_left_mut() != 0.0 {
                    *self.curve.get_oob_left_mut() /= diff.get_point(NotNan::new(crv_first_y).expect("Attempted to insert NaN")).abs()
                }
                if *self.curve.get_oob_right_mut() != 0.0 {
                    *self.curve.get_oob_right_mut() /= diff.get_point(NotNan::new(crv_last_y).expect("Attempted to insert NaN")).abs()
                }
            }
            CurveAxis::YAxis => self.get_curve_mut().multiply(&other)
        }
    }

    fn invert_axis(&mut self, axis: CurveAxis, factor: NotNan<f64>) {
        if self.curve.get_curve_mut().is_empty() {
            return;
        } else {
            match axis {
                CurveAxis::XAxis => {
                    let Some((&first_x, &_)) = self.curve.get_curve_mut().first_key_value() else { todo!() };
                    if first_x < NotNan::new(0.0).expect("Attempted to insert NaN") {
                        eprintln!("Inverting spectrums with negative values in the X axis not yet supported")
                    } else {
                        let mut new_pairs = BTreeMap::new();
                        for (x, y) in self.curve.get_curve_mut().iter_mut() {
                            let new_x = NotNan::new(*factor / x.into_inner());
                            *y *= *((*x * *x) / *factor);
                            new_pairs.insert(new_x.expect("Attempted to insert NaN"), *y);
                        }
                        *self.curve.get_curve_mut() = new_pairs;
                    }
                }
                CurveAxis::YAxis => {
                    for (_, y) in self.curve.get_curve_mut().iter_mut() {
                        *y = *factor / *y;
                    }
                }
            }
        }
    }
}