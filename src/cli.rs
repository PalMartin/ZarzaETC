use clap::{Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "ZarzaETC CLI")]
#[command(about = "Command Line Interfaze of the TARSIS ETC 'ZarzaETC'.", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub modes: Commands,

    #[command(flatten)]
    pub input: InputArgs,

    #[command(flatten)]
    pub spdtr: SpatialDistrArgs,

}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Computation of the SNR
    Snr {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation.
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon (0.0–100.0).
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red" (default: "blue").
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computaton of the exposure time
    TExp {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation.
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon.
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red".
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Target SNR to be achieved for the calculation of the exposure time.
        #[arg(long)]
        snr: f64,
        /// Wavelength of reference to make the calculations of the expected exposure time.
        #[arg(long)]
        lambda_ref: f64,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
        
    },
    /// Computation of the limiting magnitude
    LimMag {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation.
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon.
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red" (default: "blue").
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Wavelength of reference to make the calculations of the limiting magnitude.
        #[arg(long)]
        lambda_ref: f64,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computation of the limiting flux
    LimFlux {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation.
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon.
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red".
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Wavelength of reference to make the calculations of the limiting flux.
        #[arg(long)]
        lambda_ref: f64,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computation of the uncertainty in the radial velocity
    RadVelUnc {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation.
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon.
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red" (default: "blue").
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Wavelength of reference to make the calculations of the uncertainty in the radial velocity.
        #[arg(long)]
        obj_size: f64,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Simulation of the detected counts
    CountSimul {
        /// Index of the slice used for the simulation.
        #[arg(long)]
        slice: usize,
        /// Pixel position within the selected slice.
        #[arg(long)]
        pixpos: usize,
        /// Airmass during the observation (default: 1.0).
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        /// Illuminated fraction of the Moon (0.0–100.0, default: 0.0).
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        /// Number of exposures.
        #[arg(long)]
        ndit: i32,
        /// Time per exposure.
        #[arg(long)]
        dit: f64,
        /// CCD detector model, must be "ccd231_84_0_s77" or "ccd231_84_0_h69".
        #[arg(long, default_value_t = String::from("ccd231_84_0_s77"))]
        det: String,
        /// Spectrograph arm, must be "blue" or "red" (default: "blue").
        #[arg(long, default_value_t = String::from("blue"))]
        arm: String,
        /// Photometric filter used for normalization.
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        /// Target magnitude in the selected filter for the normalization.
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
}

////////////////////////////////////////// SPEXTRAL AND SPATIAL FLUX DISTRIBUTION INPUTS //////////////////////////////////////////

#[derive(Args, Debug)]
pub struct InputArgs {
    /// Selection of the type of input for the comutation: "template", "line" or "user".
    #[arg(long)]
    pub spectral: String, 

    /// For "template" mode: Name of the template to be used.
    #[arg(long)]
    pub temp: Option<String>,
    /// For "line" mode: Central wavelength of the line. In Angstrom!!
    #[arg(long)]
    pub lambda: Option<f64>,
    /// For "line" mode: Peak flux of the line. In Angstrom!!
    #[arg(long)]
    pub flux: Option<f64>,
    /// For "line" mode: FWHM of the line. In Angstrom!!
    #[arg(long)]
    pub fwhm: Option<f64>,
    /// For "line" mode: Flux of the continuum.
    #[arg(long)]
    pub continuum: Option<f64>,
    /// For "line" mode: Type of line, 'emission' or 'absorption'.
    #[arg(long)]
    pub line_type: Option<String>,
    /// For "line" mode: Specify 'resolved' or 'unresolved' line.
    #[arg(long)]
    pub resolution: Option<String>,
    /// For "line" mode: Min wavelength of the wavelength window to generate de line. In Angstrom!!
    #[arg(long)]
    pub lambda_min: Option<f64>, 
    /// For "line" mode: Max wavelength of the wavelength window to generate de line. In Angstrom!!
    #[arg(long)]
    pub lambda_max: Option<f64>,
    /// For "user" mode: Peth of the CSV with the user spectrum.
    #[arg(long)]
    pub file: Option<String>,
}

#[derive(Args, Debug)]
pub struct SpatialDistrArgs {
    /// Selection of the spatial flux distribution of the object: "infinite", "uniform", "sersic", "exponential", "gaussian" or "point".
    #[arg(long)]
    pub spatial: String, 

    /// Center of the emission?
    #[arg(long)]
    pub center: Option<f64>,
    /// Radius of the emission?
    #[arg(long)]
    pub radius: Option<f64>,
    /// Effectivde radius?
    #[arg(long)]
    pub r_e: Option<f64>,
    /// Sersic index?
    #[arg(long)]
    pub n: Option<f64>,
    /// Sigma of the gaussian profile?
    #[arg(long)]
    pub sigma: Option<f64>,
    // /// X coordinates of the point emission
    // #[arg(long)]
    // pub x0: Option<usize>,
    // /// Y coordinates of the point emission
    // #[arg(long)]
    // pub y0: Option<usize>,
}
