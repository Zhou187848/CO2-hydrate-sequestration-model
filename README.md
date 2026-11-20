# Near-well hydrate accumulation shapes injectivity during deepwater CO₂ storage

Version 1.0.1 contains the wellbore and reservoir solvers used in the study. The wellbore solver supplies the inlet conditions for a separate reservoir run. The [source workbook](data/Supplementary_Data_1_public.xlsx) indexes the shareable model results and published laboratory data by figure and table. Fig. 6 has 14 operating cases; Fig. 7 has 36 porosity–permeability combinations.

## Running the models

Start MRST 2025b in MATLAB R2025a. From the repository root:

```matlab
addpath('reservoir','wellbore','field')
[report, series] = compare_li_kinetics(8,false);
disp(report.RMSE_kg)
```

This example reproduces the laboratory comparison using the published Run 1 data (RMSE about 0.00678 kg). The reservoir solver also accepts user-defined geometry. For a short example with illustrative layer boundaries:

```matlab
c = modern_configuration(struct('sand_top',110,'sand_bottom',155, ...
    'duration',60,'snapshot_days',[],'save',false,'verbose',false));
out = axisymmetric_co2(c);
```

- `run_interface_reference`, `run_fig6`, `run_interface_grid(1,1)` and `run_interface_checks` run the design, operating, formation-property and numerical comparisons. For the study grid these drivers read a **local** MAT file specified by `CO2_MODEL_GEOMETRY`, containing scalar `sand_top` and `sand_bottom` in metres below mudline. Fig. 6 uses the simulated inlet temperatures in the workbook.
- `CO2_2(cfg)` runs the wellbore model. The configuration contains `well`, `model`, `thermal`, `injection` and `options`. In older input files, set `options.storage.enabled = false` to run the wellbore alone. `run_wellbore_interface_audit(bhp_MPa,inlet_C,rate_t_h,cfg)` uses the supplied well geometry and environmental inputs. Use a local working copy for drivers that save outputs.
- `run_c1(workbook,parameter_file)` requires the authorized operating schedule, temperature series, layer inputs and site configuration. These are controlled by the operator. The public workbook reports model endpoints and hourly model masses; it cannot drive the field run.

## Data and licence

The MIT [licence](LICENSE) applies to the authors’ software. It does not relicense data, MATLAB, MRST, CoolProp or other third-party materials. The laboratory inputs in `reservoir/li_run1_observations.mat` are from [Li et al., Energy (2024)](https://doi.org/10.1016/j.energy.2024.133224). The rows store sample index, elapsed seconds, pressure (Pa), three axial temperatures (K) and inferred hydrate mass (kg). The property arrays in `reservoir/properties.mat` are computed with CoolProp; `build_properties.py` documents their generation (Python, NumPy, SciPy, CoolProp).

Precise well locations, site locator maps, logs and log-derived profiles, gauge records, and detailed operating inputs are outside the public release. The indexes in the workbook distinguish simulated series from published measurements and identify omitted field records. For enquiries about controlled inputs, contact Liangjie Mao (maoliangjie@qq.com); the operator decides whether access can be granted.
