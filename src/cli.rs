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
        /// 
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        ///
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computaton of the exposure time
    TExp {
        ///
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        ///
        #[arg(long)]
        snr: f64,
        ///
        #[arg(long)]
        lambda_ref: f64,
        ///
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
        #[arg(long, default_value_t = 18.0)]
        mag: f64
        
    },
    /// Computation of the limiting magnitude
    MagLim {
        ///
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        ///
        #[arg(long)]
        lambda_ref: f64,
        ///
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computation of the limiting flux
    FluxLim {
        ///
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        ///
        #[arg(long)]
        lambda_ref: f64,
        ///
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Computation of the uncertainty in the radial velocity
    RadVelUnc {
        ///
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        ///
        #[arg(long)]
        obj_size: f64,
        ///
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
        #[arg(long, default_value_t = 18.0)]
        mag: f64
    },
    /// Simulation of the detected counts
    CountSimul {
        ///
        #[arg(long)]
        slice: usize,
        ///
        #[arg(long)]
        pixpos: usize,
        ///
        #[arg(long, default_value_t = 1.0)]
        airmass: f64,
        ///
        #[arg(long, default_value_t = 0.0)]
        moon: f64,
        ///
        #[arg(long)]
        ndit: i32,
        ///
        #[arg(long)]
        dit: f64,
        ///
        #[arg(long)]
        det: String,
        ///
        #[arg(long)]
        arm: String,
        /// 
        #[arg(long, default_value_t = String::from("Rc"))]
        filter: String,
        ///
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
    /// X coordinates of the point emission
    #[arg(long)]
    pub x0: Option<usize>,
    /// Y coordinates of the point emission
    #[arg(long)]
    pub y0: Option<usize>,
}
