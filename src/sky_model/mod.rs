//
// sky_model/mod.cpp: Sky model
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

use ordered_float::NotNan;
use crate::spectrum::Spectrum;
use crate::curve::Curve;
use crate::curve::CurveAxis::YAxis;
use crate::curve::CurveAxis::XAxis;
use curve::CurveOperations;
use curve::BinaryCurveOperations;
use crate::helpers::*;

pub struct SkyProperties {
    pub sky_emission: String,
    pub sky_emission_ref_airmass: f64,
    pub sky_extinction: String,
}

pub struct SkyModel {
    properties: SkyProperties,
    sky_spectrum: Spectrum,
    sky_ext: Curve,
    moon_to_mag: Curve,
    airmass: f64,
    moon_fraction: f64,
}

impl SkyProperties {
    pub fn new() -> SkyProperties {
        SkyProperties {
            sky_emission: String::from("src/sky_model/CAHASKy.csv"),
            sky_emission_ref_airmass: 1.0,
            sky_extinction: String::from("src/sky_model/CAHASkyExt.csv"),
        }
    }
}

impl SkyModel {
    pub fn new() -> Self {
        let mut sky_model = SkyModel {
            properties: SkyProperties::new(),
            sky_spectrum: Default::default(),
            sky_ext: Default::default(),
            moon_to_mag: Default::default(),
            airmass: 1.0,
            moon_fraction: 0.0,
        };
        
        sky_model.load_data();
        return sky_model;
    }

    pub fn get_sky_spectrum_mut(&mut self) -> &mut Spectrum {
        &mut self.sky_spectrum
    }
    pub fn get_sky_ext_mut(&mut self) -> &mut Curve {
        &mut self.sky_ext
    }
    pub fn get_moon_to_mag_mut(&mut self) -> &mut Curve {
        &mut self.moon_to_mag
    }
    pub fn get_airmass_mut(&mut self) -> &mut f64 {
        &mut self.airmass
    }
    pub fn get_moon_fraction_mut(&mut self) -> &mut f64 {
        &mut self.moon_fraction
    }
    pub fn get_sky_spectrum(&self) -> &Spectrum {
        &self.sky_spectrum
    }
    pub fn get_sky_ext(&self) -> &Curve {
        &self.sky_ext
    }
    pub fn get_moon_to_mag(&self) -> &Curve {
        &self.moon_to_mag
    }
    pub fn get_airmass(&self) -> &f64 {
        &self.airmass
    }
    pub fn get_moon_fraction(&self) -> &f64 {
        &self.moon_fraction
    }

    pub fn load_data(&mut self) {
        //
        // http://www.caha.es/sanchez/sky/
        // X axis is Angstrom
        // Y Axis is 1e-16 Erg / (s cm^2 A) per 2.7 arcsec diam fiber. I.e.,
        //            1.7466e-17 erg / (s cm^2 A arcsec^2), i.e.
        //             7.4309394e-10 W / (m^2 A sr)  
        let _ = self.sky_spectrum.load_curve("src/sky_model/CAHASky.csv").unwrap();
        self.sky_spectrum.scale_axis(YAxis, 7.4309394e-10); // To SI units
        self.sky_spectrum.scale_axis(XAxis, 1e-10); // Convert angstrom to meters
        let _ = self.sky_ext.load_curve("src/sky_model/CAHASkyExt.csv").unwrap();
        let _ = self.moon_to_mag.load_curve("src/sky_model/moonBrightness.csv").unwrap();
    }

    pub fn set_moon(&mut self, moon: f64) {
        if moon < 0.0 || moon > 100.0 {
            println!("Moon percent out of bounds");
        }

        self.moon_fraction = moon;
    }

    pub fn set_airmass(&mut self, airmass: f64) {
        if airmass < 1.0 {
            println!("Airmass out of bounds");
        }

        self.airmass = airmass;
    }

    pub fn set_zenith_distance(&mut self, z: f64) {
        if z < 0.0 || z >= 90.0 {
            println!("Zenith distance out of bounds");
        }

        let z_rad = z * std::f64::consts::PI / 180.0;
        self.set_airmass(1.0 / z_rad.cos())
    }

    pub fn make_sky_spectrum(&self, object: &Spectrum) -> Spectrum {
        let mut spectrum: Spectrum = Default::default();
        let sky_ext = &self.sky_ext;
        let sky_bg = &self.sky_spectrum;
        let moon = &self.moon_to_mag;
 
        spectrum.from_existing(&sky_bg.get_curve(), 1.0);
        spectrum.scale_axis(YAxis, self.airmass);
        spectrum.get_curve_mut().add(&object.get_curve());

        let xp = spectrum.x_points();
        for p in xp {
            let ext_frac = mag2frac(sky_ext.get_point(p) * self.airmass);
            
            // Apply the model: I_sky = extinction(airmass) * (object + moon + background * airmass)
            spectrum.set(p, ext_frac * (spectrum.get_point(p) + surface_brightness_ab2radiance(moon.get_point(NotNan::new(self.moon_fraction).expect("x should not be NaN")), *p)));
        }

        return spectrum;
    }
}

