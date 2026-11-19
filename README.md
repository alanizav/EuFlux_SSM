# Calculating Eulerian Fluxes using the Salish Sea Model

# Description
Code created by Alejandro Alaniz

September 2026

This code is for a capstone project at Western Washington University(WWU) for the Marine and Coastal Science (MACS) program. This code uses hindcasts produced by the Salish Sea Model(SSM), a Finite Volume Community Ocean Model(FVCOM). The SSM and its output is maintained by the Salish Sea Modeling Center (https://ssmc-uw.org/). The purpose of this code is to use the SSM output to create transects that can be used to calculate Eulerian flux. 

Claude AI was used to help create portions of the code and summarize usage of the code. Because of this, some of the code can seem messy and convoluted at times. Future versions of this project will be aimed at simplifying and fixing these issues.

Contact aalaniz20.av@gmail.com for more detail.

# Usage
## Transect flux toolkit: how to use each function

Variable names used below: `hyd` (hydrodynamic run), `wqm` (water-quality run), `grid` (mesh Dataset),
`nodesA` (picked node indices), `A` (a `run_transect` result), `T, Fl, Tr, DS` (outputs of `long_flux`),
`DSf` (low-pass filtered `DS`), `DSs` (subsampled `DSf`), and `DS10`, `DSf10`, `DSs10` (the depth-binned versions
of `DS`, `DSf` and `DSs`), and `DS_stats` (window statistics made by `flux_stats`, section 10).

Run the definition cells first (imports, `META`/`add_metadata`, `unit`, the `hydRun`/`WQMRun` classes, then the
functions). Folders such as `data` and `results` must already exist:

```python
os.makedirs('data', exist_ok=True)
os.makedirs('results', exist_ok=True)
```

### 1. Get the data

```python
# download (or reuse) one day; hyd is a hydRun, wqm is a netCDF4 Dataset
hyd = retrieve_data(resolution='low', hindcast='hyd', year=2023, day=1, save_path='data')
wqm = retrieve_data(resolution='low', hindcast='wqm', year=2023, day=1, save_path='data')

# also get the file path (to delete the file later)
hyd, hyd_path = retrieve_data(hindcast='hyd', year=2023, day=1, save_path='data', return_path=True)

# read an already-downloaded hydrodynamic file directly
hyd = read_hyd('data/ssm.hyd.2023.ssm_0001.nc')

# join consecutive days into one continuous run (list them in date order)
hyd_3days = stitch_files([hyd1, hyd2, hyd3], 'hyd')
wqm_3days = stitch_files([wqm1, wqm2, wqm3], 'wqm')
```

`day` is 1–366. `year` is 2014–2024, and `resolution='high'` only exists for 2023 and 2024.
Tracer names are the keys of `unit` (`list(unit)`). `Heat` and `Salinity` come from `hyd`; every other tracer needs `wqm`.

### 2. Mesh and picking transects

```python
# keep only the mesh from one day, so the big files aren't needed for plots
grid = fvcom_grid(hyd, save_path='results/grid_low.nc')
grid = xr.open_dataset('results/grid_low.nc').load()          # reload it later

# click nodes on the map (use %matplotlib qt); close the window when done
nodesA = pick_transect_nodes(grid)
nodesA = pick_transect_nodes(grid, as_dataarray=True, name='A')   # DataArray with lon/lat

# pick several in a row: one window per name, and each window shows the ones already picked
nodes = pick_transect_nodes(grid, names=['A', 'B', 'C'])           # dict: name -> node indices

# add one more later, with the existing ones drawn and labeled on the map
nodes['D'] = pick_transect_nodes(grid, name='D', existing=nodes)

# map of one or several transects: each is drawn as its control-volume face line, with the elements it uses shaded
# (values can be run_transect Datasets, the saved DS['A'] Datasets, or run_transect dicts)
transects_map(grid, {'A': A['ds'], 'B': B['ds']})
transects_map(grid, {'A': A['ds']}, figsize=(8, 8), pad=0.5)

# one transect with an inset: face line, shaded elements, start and end labels
plot_transect(grid, A['ds'])
plot_transect(grid, A['ds'], save_path='results/A_map.png', margin=1)   # margin: how much map around the transect
```

`pick_transect_nodes` needs the mesh Dataset (`grid`), not the run (`hyd`). Click to add a node, right-click to
undo, scroll to zoom. `existing` takes a dict of arrays or DataArrays. If you pass both `existing` and `names`,
the existing ones are kept and the new names are added to the same dict.

`plot_transect` and `transects_map` read the stored midpoints (`ds['mn']`, see section 3), so they only work on
Datasets made by the current `extract_transect`. Older results raise a `KeyError` and must be rerun.

### 3. Run one transect

```python
# volume flux only
A = run_transect(hyd, nodesA, name='A')

# put the control-volume faces on the right of the picked path instead of the left (default 'left')
A = run_transect(hyd, nodesA, name='A', side='right')

# plus tracers (WQM tracers need wqm=; Heat and Salinity come from hyd)
A = run_transect(hyd, nodesA, name='A', tracers=['NO3', 'Heat'], wqm=wqm)

A['ds']                      # one Dataset with everything
A['Flux']['tf']              # total volume flux, per hour
A['tracers']['NO3']['tf']    # total NO3 flux, per hour (kg N/s)
A['compass']['P']            # direction of positive flux
```

The pieces `run_transect` is built from (you rarely call these yourself):

```python
components = extract_transect(hyd, nodesA)                    # dict: tn, te, dx, dy, side, mn, ds
components = extract_transect(hyd, nodesA, side='right')      # faces on the right of the path
components = extract_transect(hyd, nodesA, fallback=False)    # never switch sides (see below)

# zeta at each segment's element, shape (time, segment)
corners = hyd.grid.triangles[components['te']]
zeta = np.ma.filled(hyd.data.zeta, np.nan)[:, corners].mean(axis=-1)

nx, ny, el, ds = segmented_transect(components)               # per-face normals and lengths
geometry = transect_geometry(hyd, zeta, components)           # dict: el, ct, nx, ny, hc, slc, ds
Area, Flux = transect_flux(hyd, components, geometry)         # dicts: ac, apl, ta / fc, fpl, tf
out = tracer_flux(hyd, wqm, components, Flux, 'NO3')          # dict: fc, fpl, tf, unit, tracer, ds
compass = flux_direction(hyd, components)                     # dict: P, N, pb, ds
```

### 4. Plot

Static plot, by layer:

```python
plot_flux(A['ds'])                                  # flux by sigma layer vs date
plot_flux(A['ds'], tracer='NO3')
plot_flux(A['ds'], save_path='results/A_flux.png')
plot_flux(A['ds'], window=(0, 12))                  # only hours 0 to 12
plot_flux(DS['A'], window=('2023-03-01', '2023-03-15'))     # or give dates
plot_flux(DSf['A'], tracer='NO3', kind='Residual')          # title prefix
plot_flux(DS10['A'], binned=True)                    # depth-binned data (section 8): colorbar in meters
```

Animations (one Dataset per call; `plot_vectors` takes a dict):

```python
plot_crosssection(DSs['A'])                         # animated cross-section
plot_crosssection(DSs['A'], tracer='NO3', save_path='results/A_NO3.gif')

plot_layerflux(DSs['A'], tracer='Heat')             # animated bars by sigma layer
plot_layerflux(DSf10['A'], binned=True)              # animated bars by depth bin

# animated arrows on the map: needs a dict of Datasets (A['ds'], not A); each transect is drawn as its face line,
# with the arrow at its center
plot_vectors(grid, {'A': A['ds'], 'B': B['ds']})
plot_vectors(grid, DSs, tracer='NO3', max_frac=0.2, pad=0.8, window=('2023-03-01', '2023-03-15'))
```

### 5. Run a long record (a year or any date range)

`long_flux` downloads one day at a time, runs every transect, and deletes each day's files, so only the small
results stay in memory.

```python
# whole year (the default date_range is 2023-01-01 to 2023-12-31)
T, Fl, Tr, DS = long_flux({'A': nodesA, 'B': nodesB}, save_path='data',
                          tracers=['NO3', 'DOXG', 'Heat', 'Salinity'])

# any date range (strings, datetime, date or Timestamp); it may cross a year boundary
T, Fl, Tr, DS = long_flux({'A': nodesA, 'B': nodesB}, save_path='data',
                          date_range=('2023-04-10', '2023-05-10'), tracers=['NO3', 'Heat'])
T, Fl, Tr, DS = long_flux({'A': nodesA, 'B': nodesB}, save_path='data',
                          date_range=('2022-12-01', '2024-01-31'))

# faces on the right of the picked paths (default 'left')
T, Fl, Tr, DS = long_flux({'A': nodesA, 'B': nodesB}, save_path='data', side='right')

T['A']                       # name, zeta, components, geometry, compass
Fl['A']['tf']                # volume flux, transect A
Tr['A']['NO3']['tf']         # NO3 flux, transect A
DS['A']                      # combined Dataset for transect A (same layout as A['ds'])
```
Save the results right away:

```python
for n, ds in DS.items():
    ds.to_netcdf(f'results/{n}_2023.nc')

DS = {n: xr.open_dataset(f'results/{n}_2023.nc').load() for n in ['A', 'B']}
```

### 6. Smooth the record (remove the tides)

`lowpass_dataset` applies a Butterworth low-pass filter to every flux in a Dataset: `fc`, `fpl` and `tf`, and
the tracer versions (`NO3_fc`, `NO3_tf`, and so on). Everything else (`zeta`, geometry, direction) is left
alone, so all the plot functions work on the result. It works on one Dataset, so use a dictionary comprehension
for a dict. The output has the same hourly time axis as the input.

```python
DSf = {n: lowpass_dataset(d, cutoff_hours=720, edge_hours=720) for n, d in DS.items()}

plot_flux(DSf['A'], tracer='NO3', kind='Residual')
DSf['A']['tf'].to_series()            # smoothed total flux as a pandas Series
DSf['A'].to_netcdf('results/A_2023_smoothed.nc')
```

### 7. Shorten a long record for the animations

A year of hourly frames is too large to animate. Each frame is about 0.07 MB, and the notebook sets
`plt.rcParams['animation.embed_limit'] = 350` (MB) in the dynamic plots cell. That allows about 5,000 frames,
but long animations still make the notebook slow. `subsample(ds, hours)` keeps the time steps that are
multiples of `hours` and drops everything between them. It takes one Dataset or a dict, and returns the same
kind of object with the same layout, so all the plot functions work on the result.

```python
DSs = subsample(DSf, 18)                                      # whole dict at once
DSs_A = subsample(DSf['A'], 18)                               # one Dataset

plot_crosssection(DSs['A'], tracer='NO3')
plot_vectors(grid, DSs)                                       # same time axis for every transect
```

### 8. Group flux by depth instead of sigma layer

Sigma layers follow the bottom, so one layer is at different depths in different parts of a transect.
`depth_bin_flux` regroups the flux into bins of fixed depth below sea level.

```python
DS10 = depth_bin_flux(DS, n_bins=7)                            # every transect
DS10 = depth_bin_flux(DS, n_bins=7, surface_bin=True)          # plus one bin for water above sea level
B = depth_bin_flux(DS['A'], n_bins=10, max_depth=100)         # one transect, fixed deepest edge

DSf10 = {n: lowpass_dataset(d, cutoff_hours=720, edge_hours=720) for n, d in DS10.items()}
DSs10 = subsample(DSf10, 18)

plot_flux(DS10['A'], binned=True)
plot_layerflux(DSf10['A'], binned=True)
```

### 9. Net flux into an area enclosed by transects

`net_flux` combines the total flux of two or more transects that enclose an area. Pass the transect Datasets
(or one dict of them) and the mesh. Each transect's flux counts as influx when the positive direction points
toward the middle of the area and as outflux when it points away. The middle is the average of the mean
longitude and latitude of each transect's nodes.

```python
net = net_flux(grid, DS['A'], DS['B'])                                  # DataArray, m^3/s, positive = influx
net = net_flux(grid, DS['A'], DS['B'], tracer='NO3')                    # tracer flux
net = net_flux(grid, DSf['A'], DSf['B'], window=('2023-01-01', '2023-12-31'))   # smoothed (residual) net flux
out = net_flux(grid, {'A': DS['A'], 'B': DS['B']}, window=('2023-01-01', '2023-12-31'), detail=True)   # a dict works too
out['influx'], out['outflux'], out['netflux'], out['toward_center']    # detail=True returns a Dataset
```

### 10. Statistics over a window, and static plots

`flux_stats` summarizes flux over a time window and returns a Dataset (`DS_stats`) instead of a table. For every
flux variable it makes seven statistics, named `<variable>_<stat>`: `mean`, `std`, `max`, `time_of_max`, `min`,
`time_of_min` and `n_hours`. The three flux outputs stay separate: `tf` (total), `fpl` (per layer) and `fc`
(per chunk, meaning layer and segment). Each statistic keeps the dimensions that remain after the time
dimension is collapsed, so `tf_*` has none, `fpl_*` has `layer` and `fc_*` has `layer, segment`.

```python
w = ('2023-03-01', '2023-03-31')

DS_stats = flux_stats(DS, window=w)                       # dict, same keys as DS
DS_stats['A']['tf_mean']                                  # total volume flux, mean (one number)
DS_stats['A']['fpl_std']                                  # per-layer standard deviation, dims (layer)
DS_stats['A']['fc_time_of_max']                           # per-chunk time of maximum, dims (layer, segment)

LD = flux_stats(DS['A'], tracer='LDOC', window=w)         # LDOC_tf_mean, LDOC_fpl_std, ...
both = flux_stats(DS['A'], tracer=[None, 'LDOC'], window=w)   # volume and tracer in one Dataset
tot = flux_stats(DS, window=w, include=('tf',))           # totals only (include picks tf, fpl, fc)
sep = flux_stats_all(DS['A'], window=w)                   # {'tf': ..., 'fpl': ..., 'fc': ...}, separate Datasets
```

```python
# cross-section of the mean or the standard deviation of chunk flux (static version of plot_crosssection)
static_crosssection(DS['A'], DS_stats['A'], stat='mean', window=w)
static_crosssection(DS['A'], DS_stats['A'], stat='std', window=w)
static_crosssection(DS['A'], DS_stats['A'], stat='mean', tracer='LDOC', window=w, save_path='A_LDOC_mean.png')

# bars of mean flux per layer with +/- 1 standard deviation error bars (static version of plot_layerflux)
static_layer_flux(DS['A'], DS_stats['A'])
static_layer_flux(DS['A'], DS_stats['A'], tracer='LDOC', save_path='A_LDOC_layers.png')
static_layer_flux(DS10['A'], DS10_stats['A'], binned=True)    # depth-binned data
static_layer_flux(DS['A'], DS_stats['A'], show_std=False)     # means only
```

### 11. Tables

`flux_dataframes` takes one Dataset or a dict of Datasets (`A['ds']`, `DS`, `DSf`) and returns a dict of DataFrames:
`total`, `layers`, `direction`, and `chunks` only when `chunks=True`.

```python
tables = flux_dataframes(DS)                                  # all tracers in the Datasets
tables = flux_dataframes(DS, tracers=['NO3'], start='2023-01-01')
tables = flux_dataframes(DS, start='2023-01-01', end='2023-01-31 23:00')   # end needs start
tables = flux_dataframes(A['ds'], chunks=True)                # single run, includes chunk table

tables['total'].loc['A']['NO3 (kg N/s)'].plot()
tables['direction']
```

