# O3DTM
> 3D thermal physical model for the lunar surface

## TYPICAL ROUTINES

### Generate new crater environment (e.g. create environment called 'abc')

1. Add lat/long limits to `generate_crater_environement` (i.e. case 'abc'...)
2. Run `generate_crater_environment`('abc', true) 

This will generate the crater environment and perform ray tracing.

### Run 3D model for crater environment 'abc' for scattering powers of 0, 0.5, 1, 1.5, 2 with maximum temperatures and seasonal variation

```matlab
simulate_crater_temperatures('max', 'seasons', 5, 0, 2)
```

### Generate map of H parameter

```matlab
generate_map('H', 'LOLA A0') % The 'LOLA A0' here makes the fitting routine use LOLA albedo data
```

### Set up 3D model to run for different time periods
	
1. Download new ephemeris data from JPL Horizons (see inputs/ephemerides/ for required output format)
2. Run `save_horizons_data`
3. Run with desired dtm parameters (either update defaults in define_parameters or use custom_parameters as function argument):
  - `diviner_start_dtm` & `diviner_end_dtm` define the time period to get results from
  - `phase_start_dtm` & `phase_end_dtm` define period where simulations begin (simulation begins at local midday which occurs between phase_start_dtm and phase_end_dtm)
  - > NOTE want ~ 5 years between phase_end_dtm and diviner_start_dtm so that model has time for temperatures to fully stabilise before being used to generate results


## ALL FUNCTIONS/SCRIPTS

### Main content in code/...

| NAME | DESCRIPTION |
|------|-------------|
|compare_1d_3d_models | Generates data file with 1D and 3D model temperatures for a given crater|
|convert_height_to_slope	|	Converts elevation matrix to slope matrix |
|create_static_path		|	Returns static file location (where big data files are stored)
|define_parameters	|		Defines model parameter values|
|divergent_colormap	|		Changes plot colormap to be divergent (+ve values = red, -ve values = blue, 0 = white)|
|find_diviner_extreme_temperatures|	Finds extreme diviner temperatures for given ppd|
|generate_compressed_crater_environment|	Compresses crater environment so that spatial resolution is same in NS and EW direction (used for high latitude craters)|
|generate_crater_environment|		Generates crater environment used for 3D modeling|
|generate_crater_ray_tracing|		Pre-calculates ray-tracing and solar angle calculations for crater environment for use in 3D model|
|generate_map		|		Generates parameter map using fitting routine and 1D model|
|generate_simulated_crater_environment|	Generates crater environment for simulated bowl shaped crater|
|git_commit			|	(Tool to make creating git commit and backing up to GitHub quicker)|
|hayne_1d_model	|			1D thermal model|
|plot_3d_surface	|			Used to plot 3D surfaces easily|
|plot_compare_model_diviner_temperatures|	Plots comparison of 1D model and Diviner temperatures|
|plot_diff_map	|			Plots map of difference between two different values (useful for comparing different results of generate_map for same variable)|
|plot_histogram	|			Plots histogram of data with bar colours consistent with map colourmap|
|plot_map		|		Plots map of value over lunar surface, generally used for results of generate_map|
|plot_ortho_map			|	Plots orthographic map of value over lunar surface|
|plot_report_figures		|	Plots figures used in project report|
|process_horizons_data		|	Converts downloaded JPL Horizons .txt data files to MATLAB matrices|
|read_diviner_temperatures|		Loads data file for max/min diviner temperature|
|read_horizons_data	|		Loads data file for ephemeris data|
|read_lola_A0_data		|	Loads LOLA albedo data|
|read_lola_height_data	|		Loads LOLA elevation data|
|save_diviner_extreme_temperatures|	Saves data files for diviner max/min temperature data|
|save_horizons_data	|		Saves processed JPL Horizons data into data file for use elsewhere|
|save_report_figure		|	Saves figure in format used in project report|
|scale_data		|		Scales a matrix so that the mean of the values is a specified value|
|scratchpad			|	(Script used to test snippets of code)|
|simulate_1d_temperatures|		Simulates 1D temperatures for whole lunar surface and saves result in data file|
|simulate_crater_temperatures|		Simulates crater temperatures using 3D model for different scattering powers and saves result in a data file|
|spectral	|			Colourmap varying from red > white > blue (used in divergent_colormap)|
|temperature_model_3d	|		3D thermal model|
|viewshed		|		[Copy of MATLAB viewshed function fixing roundoff bug which will be fixed in future MATLAB release]|
|viridis		|			Colourmp (used as default for most of my plots)|


### Additional content in code/useful_scripts/...

> may need to move these files to code/ for them to work properly

|NAME | DESCRIPTION |
|-----|-------------|
|b2r			|		Blue to red colourmap|
|call_generate_map|			Script to generate multiple parameter maps|
|find_albedo_dependence|			Finds temperature dependence on albedo|
|find_H_dependence|			Finds temperature dependence on H|
|find_temperature_dependence|		More general function to find temperature dependence of any given parameter|
|fit_scattering_power			|(Failed) attempt to fit scattering power from 3D model temperatures and Diviner temperatures|
|generate_lola_maps			|Generates various useful maps from LOLA data|
|generate_slope_temperature_maps|		Generates various temperature maps and animation of 1D model temperatures|
|graph_diviner_datetimes			|Plots graph of when extreme diviner temperatures occurred|
|inferno					|Alternative colourmap|
|interpolate_emissivity|			Interpolates emissivity to try to compare fitted values to Apollo measurements|
|magma					|Alternative colourmap|
|phasemap 			|	Alternative colourmap good for representing angles|
|plasma					|Alternative colourmap|
|plot_albedo_interpolation|		Graphically represents albedo interpolation process used in finding albedo for 128ppd craters|



### Various old files etc. are in code/old/ (if anything breaks, the function it needs has probably been archived to here)
