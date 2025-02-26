//
// helpers/mod.cpp: Helper functions.
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

pub const SPEED_OF_LIGHT: f64 = 299792458.0;      // m/s
pub const PLANCK_CONSTANT: f64 = 6.62607015e-34;  // J s
pub const JANSKY: f64 = 1e-26;                    // W / (m^2 Hz)
pub const AB_ZEROPOINT: f64 = 3631.0 * JANSKY;    // W / (m^2 Hz)
pub const ARCSEC: f64 = 4.8481368e-6;             // radian
pub const STD2FWHM: f64 = 2.35482004503095;       // sqrt(8 * ln(2))
pub const INVSQRT2PI: f64 = 0.39894228040143;     // 1 / sqrt(2 * pi)

// Converts magnitude to flux fraction
pub fn mag2frac(mag: f64) -> f64 {
    let base: f64 = 10.0;
    return base.powf(-0.4 * mag)
}

// Converts surface brightness in AB to frequency radiance (W / (m^2 Hz))
pub fn surface_brightness_ab2freq_radiance(mag: f64) -> f64 {
    let fnu = mag2frac(mag) * AB_ZEROPOINT / (ARCSEC * ARCSEC);
    return fnu
}

// Converts surface brightness in AB to radiance (W / (m^2 sr))
pub fn surface_brightness_ab2radiance(mag: f64, wl: f64) -> f64 {
    let fnu = surface_brightness_ab2freq_radiance(mag);
    let flambda = SPEED_OF_LIGHT / (wl * wl) * fnu;
    return flambda
}

// Converts surface brightness in Vega to radiance (W / (m^2 sr))
pub fn surface_brightness_vega2radiance(mag: f64) -> f64 {
    let flambda = mag2frac(mag) * 0.03631 / (ARCSEC * ARCSEC);
    return flambda
}