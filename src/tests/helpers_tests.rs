use crate::helpers::*;

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_mag2frac() {
        let mag = 0.0;
        assert_eq!(mag2frac(mag), 1.0);
        let mag = 1.0;
        assert!((mag2frac(mag) - 0.3981).abs() < 1e-4);
    }
    #[test]
    fn test_surface_brightness_ab2freq_radiance() {
        let mag = 0.0;
        let expected = AB_ZEROPOINT / (ARCSEC * ARCSEC);
        assert!((surface_brightness_ab2freq_radiance(mag) - expected).abs() < 1e-4);
    }
    #[test]
    fn test_surface_brightness_ab2radiance() {
        let mag = 0.0;
        let wl = 500e-9;
        let expected = (SPEED_OF_LIGHT / (wl * wl)) * surface_brightness_ab2freq_radiance(mag);
        assert!((surface_brightness_ab2radiance(mag, wl) - expected).abs() < 1e-4);
    }
    #[test]
    fn test_surface_brightness_vega2radiance() {
        let mag = 0.0;
        let expected = 0.03631 / (ARCSEC * ARCSEC);
        assert!((surface_brightness_vega2radiance(mag) - expected).abs() < 1e-4);
    }
}

