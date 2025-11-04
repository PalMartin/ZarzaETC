//
// sky_model/mod.rs: Sky model
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

use ordered_float::NotNan;

use crate::curve;
use crate::spectrum::*;
use crate::curve::*;
use crate::curve::CurveAxis::YAxis;
use crate::curve::CurveAxis::XAxis;
use curve::CurveOperations;
use curve::BinaryCurveOperations;
use crate::helpers::*;
use crate::config_manager::*;
use crate::data_file_manager::*;

/// Structure to hold sky properties loaded from configuration files and used in the sky model.
#[derive(Default, Clone)]
pub struct SkyProperties {
    /// Configuration manager instance for loading and saving properties.
    pub config: Config,

    /// Path to the sky emission data file.
    pub sky_emission: String,
    /// Reference airmass for the sky emission data.
    pub sky_emission_ref_airmass: f64,
    /// Path to the sky extinction data file.
    pub sky_extinction: String,
}

/// Structure representing the sky model, including properties and sky emission/extinction data.
#[derive(Clone)]
pub struct SkyModel {
    /// Sky properties loaded from configuration.
    pub properties: SkyProperties,
    /// Spectrum representing the sky emission.
    sky_spectrum: Spectrum,
    /// Spectrum representing the object's spectrum before sky effects.
    obj_spectrum: Spectrum,
    /// Curve representing the sky extinction as a function of wavelength.
    sky_ext: Curve,
    /// Curve representing the moon brightness as a function of moon phase.
    moon_to_mag: Curve,
    /// Current airmass value for the observation.
    airmass: f64,
    /// Current moon fraction (percentage of moon illumination).
    moon_fraction: f64,
}

/// Implementation of methods for managing sky properties, including loading and saving configuration.
impl SkyProperties {
    /// Create a new SkyProperties instance with default values and load configuration from file.
    pub fn new(name: String) -> SkyProperties {
        // Initialize configuration manager with the given name.
        let mut config = Config::new(name.clone());
        config.load();

        // Set default values for sky properties.
        let mut sky_properties = SkyProperties {
            // Initialize with the configuration manager.
            config,
            // Default sky emission data file.
            sky_emission: String::from("CAHASKy.csv"),
            // Default reference airmass for sky emission.
            sky_emission_ref_airmass: 1.0,
            // Default sky extinction data file.
            sky_extinction: String::from("CAHASkyExt.csv"),
        };

        // Load values from the configuration file, overriding defaults from cotrol YAML if present.
        sky_properties.deserialize(); 
        return sky_properties;
    }

    /// Serialize the current properties to the configuration manager for saving to file.
    pub fn serialize(&mut self) -> bool {

        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "sky_emission".into(), Config::serialize_field(&self.sky_emission)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "sky_emission_ref_airmass".into(), Config::serialize_field(&self.sky_emission_ref_airmass)
        );
        self.config.yaml_config.as_mapping_mut().unwrap().insert(
            "sky_extinction".into(), Config::serialize_field(&self.sky_extinction)
        );

        return true;
    }

    /// Deserialize properties from the configuration manager, loading values from file.
    pub fn deserialize(&mut self) -> bool {

        self.config.deserialize_field(&mut self.sky_emission, "sky_emission");
        self.config.deserialize_field(&mut self.sky_emission_ref_airmass, "sky_emission_ref_airmass");
        self.config.deserialize_field(&mut self.sky_extinction, "sky_extinction");

        return true;
    }
}

/// Implementation of methods for the SkyModel, including loading data and generating sky spectra.
impl SkyModel {
    /// Create a new SkyModel instance, loading properties and data from configuration and files.
    pub fn new() -> Self {
        // Load properties from configuration file.
        let mut properties = SkyProperties::new("sky".to_string());
        properties.config.load();
        properties.deserialize();

        // Initialize the SkyModel with loaded properties and default values for spectra and curves.
        let mut sky_model = SkyModel {
            // Set the loaded properties.
            properties,

            // Initialize spectra and curves with default values.
            sky_spectrum: Default::default(),
            obj_spectrum: Default::default(),
            sky_ext: Default::default(),
            moon_to_mag: Default::default(),
            // Set default airmass and moon fraction.
            airmass: 1.0,
            moon_fraction: 0.0,
        };

        // Load sky emission and extinction data from files.
        sky_model.load_data();
        return sky_model;
    }

    /// Mutable access to the sky emission spectrum.
    pub fn get_sky_spectrum_mut(&mut self) -> &mut Spectrum {
        &mut self.sky_spectrum
    }
    /// Mutable access to the object's spectrum before sky effects.
    pub fn get_obj_spectrum_mut(&mut self) -> &mut Spectrum {
        &mut self.obj_spectrum
    }
    /// Mutable access to the sky extinction curve.
    pub fn get_sky_ext_mut(&mut self) -> &mut Curve {
        &mut self.sky_ext
    }
    /// Mutable access to the moon brightness curve.
    pub fn get_moon_to_mag_mut(&mut self) -> &mut Curve {
        &mut self.moon_to_mag
    }
    /// Mutable access to the current airmass value.
    pub fn get_airmass_mut(&mut self) -> &mut f64 {
        &mut self.airmass
    }
    /// Mutable access to the current moon fraction value.
    pub fn get_moon_fraction_mut(&mut self) -> &mut f64 {
        &mut self.moon_fraction
    }
    /// Immutable access to the sky emission spectrum.
    pub fn get_sky_spectrum(&self) -> &Spectrum {
        &self.sky_spectrum
    }
    /// Immutable access to the object's spectrum before sky effects.
    pub fn get_sky_ext(&self) -> &Curve {
        &self.sky_ext
    }
    /// Immutable access to the moon brightness curve.
    pub fn get_moon_to_mag(&self) -> &Curve {
        &self.moon_to_mag
    }
    // Immutable access to the current airmass value.
    pub fn get_airmass(&self) -> &f64 {
        &self.airmass
    }
    /// Immutable access to the current moon fraction value.
    pub fn get_moon_fraction(&self) -> &f64 {
        &self.moon_fraction
    }

    /// Load sky emission and extinction data from the specified files.
    pub fn load_data(&mut self) {
        // Data obtained from CAHA website:
        // http://www.caha.es/sanchez/sky/
        // Units:
        // X axis is Angstrom
        // Y Axis is 1e-16 Erg / (s cm^2 A) per 2.7 arcsec diam fiber. I.e.,
        //            1.7466e-17 erg / (s cm^2 A arcsec^2), i.e.
        //             7.4309394e-10 W / (m^2 A sr)  
        let _ = self.sky_spectrum.load_curve(&DataFileManager::data_file(&"CAHASky.csv".to_string())).unwrap();
        self.sky_spectrum.scale_axis_factor(YAxis, 7.4309394e-10); // To SI units
        self.sky_spectrum.scale_axis_factor(XAxis, 1e-10); // Convert angstrom to meters
        let _ = self.sky_ext.load_curve(&DataFileManager::data_file(&"CAHASkyExt.csv".to_string())).unwrap();
        self.sky_ext.scale_axis(XAxis, 1e-9); // Convert nm to meters
        let _ = self.moon_to_mag.load_curve(&DataFileManager::data_file(&"moonBrightness.csv".to_string())).unwrap();
    }

    /// Set the moon fraction (percentage of moon illumination) for the sky model.
    pub fn set_moon(&mut self, moon: f64) {
        // Ensure the moon fraction is within valid bounds (0 to 100).
        if moon < 0.0 || moon > 100.0 {
            println!("Moon percent out of bounds");
        }

        self.moon_fraction = moon;
    }

    /// Set the airmass value for the sky model.
    pub fn set_airmass(&mut self, airmass: f64) {
        // Ensure the airmass is within valid bounds (greater than or equal to 1).
        if airmass < 1.0 {
            println!("Airmass out of bounds");
        }

        self.airmass = airmass;
    }

    // Set the zenith distance (angle from zenith) and compute the corresponding airmass.
    pub fn set_zenith_distance(&mut self, z: f64) {
        // Ensure the zenith distance is within valid bounds (0 to 90 degrees).
        if z < 0.0 || z >= 90.0 {
            println!("Zenith distance out of bounds");
        }

        // Convert zenith distance from degrees to radians and compute airmass.
        let z_rad = z * std::f64::consts::PI / 180.0;
        self.set_airmass(1.0 / z_rad.cos())
    }

    /// Generate the sky spectrum and the object's spectrum after applying sky effects.
    /// Returns a tuple containing the sky spectrum and the object's spectrum.
    pub fn make_sky_spectrum(&self, object: &Spectrum) -> (Spectrum, Spectrum) {
        // Create mutable copies of the sky and object spectra to apply transformations.
        let mut spectrum_sky: Spectrum = Default::default();
        let mut spectrum_obj: Spectrum = Default::default();
        // References to the sky extinction curve, sky emission spectrum and moon brightness curve.
        let sky_ext = &self.sky_ext;
        let sky_bg = &self.sky_spectrum;
        let moon = &self.moon_to_mag;
 
        // Rescale the sky background to the current airmass.
        spectrum_sky.from_existing(&sky_bg.get_curve(), 1.0);
        spectrum_sky.scale_axis_factor(YAxis, self.airmass);
        spectrum_obj.from_existing(&object.get_curve(), 1.0);

        // Apply extinction and moon brightness to the sky and object spectra.
        let xp = spectrum_sky.x_points();
        for p in xp {
            // Get the extinction factor for the current wavelength and airmass.
            let ext_frac = mag2frac(sky_ext.get_point(p) * self.airmass);
            // Apply extinction to the sky background and add moon contribution if applicable.
            spectrum_sky.set(p, ext_frac * (spectrum_sky.get_point(p) + surface_brightness_ab2radiance(moon.get_point(NotNan::new(self.moon_fraction).expect("x should not be NaN")), *p)));
        }
        let xp = spectrum_obj.x_points();
        for p in xp {
            let ext_frac = mag2frac(sky_ext.get_point(p) * self.airmass);
            // Apply extinction to the object's spectrum.
            spectrum_obj.set(p, ext_frac * spectrum_obj.get_point(p));
        }
        return (spectrum_sky, spectrum_obj);

    }
}

