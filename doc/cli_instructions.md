# ZarzaETC - Command Line Interface (CLI)
`zarza_etc` is the Command Line Interface for the **TARSIS ETC (Exposure TIme Calculator)**

## Usage 

```bash
zarza_etc --spectral <SPECTRAL> --spatial <SPATIAL> <COMMAND>
```

## Required Options

These two options must be provided for every command:

| Option                  | Description                                                                                                 |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `--spectral <SPECTRAL>` | Defines the spectral input mode (`template`, `line`, or `user`)                                             |
| `--spatial <SPATIAL>`   | Defines the spatial flux distribution (`infinite`, `uniform`, `sersic`, `exponential`, `gaussian`, `point`) |

## Spectral Input Configuration (`--spectral`)

The available spectral distributions are the following, and depending on the spectral model selected, extra parameters may be required.

| Spatial Profile                  | Options                      |
| ------------------------------ | ---------------------------- |
| `template` | `--temp` |
| `line` |  `--lambda`, `--flux`,  `--fwhm`, `--continuum`, `--line-type`, `--resolution`, `--lambda-min`, `--lambda-max`|
| `user`                       | `--file`   |

Optional parameters:

| Parameter     | Description                          |
| ---------- | ------------------------------------ |
| `--temp` | For "template" mode: Name of the template to be used       |
| `--lambda` | For "line" mode: Central wavelength of the line. Must be in Angstrom |
| `--flux`    | For "line" mode: Peak flux of the line. Must be in cgs units   |
| `--fwhm`      | For "line" mode: FWHM of the line. Must be in Angstrom                        |
| `--continuum`  | For "line" mode: Flux of the continuum         |
| `--line-type`  | For "line" mode: Type of line, `emission` or `absorption                     |
| `--resolution`  | For "line" mode: Specify `resolved` or `unresolved line    |
| `--lambda-min`  | For "line" mode: Minimum wavelength of the wavelength window to generate de line. Must be in Angstrom                       |
| `--lambda-max`  | For "line" mode: Maximum wavelength of the wavelength window to generate de line. Must be in Angstrom                 |
| `--file`  | Gaussian sigma                       |


## Spatial Profile Configuration ('--spatial')

The available spatial profiles are the following, and depending on the spatial model selected, extra parameters may be required.

| Spatial Profile                  | Options                      |
| ------------------------------ | ---------------------------- |
| `infinite` | *(No additional parameters)* |
| `uniform` |  `--center`, `--radius` |
| `sersic`                       | `--center`, `--r-e`, `--n`   |
| `exponential`                  | `--center`, `--r-e`       |
| `gaussian`                     | `--center`, `--sigma`        |
|  `point` | *(No additional parameters)* |

Optional parameters:

| Parameter     | Description                          |
| ---------- | ------------------------------------ |
| `--center` | Center of emission                   |
| `--radius` | Radius of emission (for exponential) |
| `--r-e`    | Effective radius (Sérsic profile and Exponential profile)    |
| `--n`      | Sérsic index                         |
| `--sigma`  | Gaussian sigma                       |

## Available Commands

### `snr` — Signal-to-Noise Ratio

Computes the Signal-to-Noise Ratio expected for the observing conditions given by the user. This saves the signal, noise and SNR curves in .dsv

**Required parameters (besides `--spectral` and `--spatial`):**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |


**Optional parameters:**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |

**Example:**

```bash
zarza_etc --spectral "template" --spatial "infinite" --temp "Sa_template.csv" snr --slice 20 --pixpos 20 --ndit 3 --dit 3600.0 --det "ccd231_84_0_h69"
```

### `t-exp` — Exposure Time

Computes the exposure time required to achieve a given SNR. Returns the result in the console.

**Required parameters (besides `--spectral` and `--spatial`):**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--snr <FLOAT>` | Target SNR to be achieved for the calculation of the exposure time |
| `--lambda-ref <FLOAT>` | Wavelength of reference to make the calculations of the expected exposure time |

**Optional parameters:**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |

**Example:**

```bash
zarza_etc --spectral "template" --spatial "infinite" --temp "Sa_template.csv" t-exp --slice 20 --pixpos 20 --ndit 3 --dit 3600.0 --det "ccd231_84_0_h69" --arm "red" --mag 22.0 --filter "Rc" --airmass 1.2 --moon 50.0 --snr 10.0 --lambda-ref 6563
```

### `lim-mag` — Limiting Magnitude

Computation of the limiting magnitude detectable by the instrument for an observational setup.

**Required parameters (besides `--spectral` and `--spatial`):**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |
| `--lambda-ref <FLOAT>` | Wavelength of reference to make the calculations of the expected exposure time |

**Optional parameters:**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |-
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |

**Example:**

```bash
--spectral "line" --spatial "infinite" --lambda 6563 --flux 12.0e-14 --fwhm 0.10 --continuum 2.8e-14 --line-type "emission_line" --resolution "resolved"  --lambda-min 6000 --lambda-max 7000 lim-mag --slice 20 --pixpos 20 --ndit 3 --dit 3600.0 --det "ccd231_84_0_h69" --arm "red" --mag 22.0 --filter "Rc" --airmass 1.2 --moon 50.0 --lambda-ref 6563
```

### `lim-flux` — Limiting Flux

Computation of the limiting flux detectable by the instrument for an observational setup. Only works for the emission/absorption line spectral setup.

**Required parameters (besides `--spectral` and `--spatial`):**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |
| `--lambda-ref <FLOAT>` | Wavelength of reference to make the calculations of the expected exposure time |

**Optional parameters:**

| Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |-
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |

**Example:**

```bash

```

### `rad-vel-unc`  — Radial Velocity Uncertainty

Computation of the uncertainty in the radial velocity expected of the detector for an observational setup. Only works for the emission/absorption line spectral setup.
 
 **Required parameters (besides `--spectral` and `--spatial`):**
 
 | Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |
| `--obj-size <FLOAT>` | Size of the object on he detector to make the calculations of the uncertainty in the radial velocity |
 
 **Optional parameters:**
 
 | Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |-
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |

 
 **Example:**
 
 ```bash

```
 
### `count-simul`  — Simulated Counts 

Simulation of the detected counts by the instrument. The ressults of the simulation are saved in the archive with the name: "calc_work_simulated_counts_curve.csv".

**Required parameters (besides `--spectral` and `--spatial`):**

 | Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--slice <INT>`  | Index of the slice used for the simulation      |
| `--pixpos <INT>` | Pixel position within the selected slice           |
| `--ndit <INT>`   | Number of exposures                      |
| `--dit <FLOAT>`  | Time per exposure      |

**Optional parameters:**

 | Parameter           | Description                              |
| ---------------- | ---------------------------------------- |
| `--airmass <FLOAT>`  | Airmass during the observation (default: 1)     |
| `--moon <FLOAT>`  | Illuminated fraction of the Moon (0.0–100.0) (default: 0)    |-
| `--det <STRING>` | CCD detector model, must be `ccd231_84_0_s77` or `ccd231_84_0_h69` (default: `ccd231_84_0_s77`) |
| `--arm <STRING>`  | Spectrograph arm, must be "blue" or "red" (default: "blue") (default: blue)     |
| `--filter <STRING>`  | Photometric filter used for normalization (default: Rc) |
| `--mag <FLOAT>`  | Target magnitude in the selected filter for the normalization (default: 18)     |open 

**Example:**

```bash

```

### help         

Print this message or the help of the given subcommand(s).

Example of help for a command:

```bash
zarza_etc snr --help
```
or 

```bash
zarza_etc snr -h
```
