# Gornergletscher sediment provenance and flow-distance analysis
Code and input data for the article "The influence of overdeepenings on sediment flux from an Alpine glacier"

This repository contains a Python workflow for analysing the relationship between **glacier flow distance, altitude, sediment provenance, and sediment characteristics** for Gornergletscher in the Swiss Alps.

The workflow combines gridded topographic datasets, glacier geometry, glacier-part boundaries, geological information, suspended-load provenance estimates, and bedload clast observations to investigate how sediment provenance and sediment characteristics vary across the glacier system.

The analysis produces figures describing:

* glacier flow distance from the snout;
* the relationship between altitude and flow distance;
* suspended-load provenance;
* bedload provenance;
* the distribution of sediment lithologies;
* sediment angularity, grain size, and mass;
* the relationship between clast angularity and flow distance;
* differences between lithological contributions based on clast counts and clast mass.

> **Important:** The suspended-load provenance product used by this repository is generated in a separate companion repository, [Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion). This repository uses the resulting provenance/erosion-rate raster as an input for the glacier-scale analysis described here.

---

# 1. Relationship to the XRD inversion workflow

The suspended-load provenance information used here originates from the **non-linear XRD inversion workflow** developed in the companion repository:

**[Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion)**

That repository contains the code and data used to go from X-ray diffraction (XRD) measurements of source and suspended-sediment samples to spatially distributed erosion-rate maps.

The inversion treats the suspended-sediment signal as a weighted combination of the tracer-mineral concentrations of the source areas, with the weights corresponding to spatially varying erosion rates. Erosion rates are parameterised in log-space, which imposes a positivity constraint and makes the inverse problem non-linear. The companion repository tests both steepest-descent and quasi-Newton approaches.

For the full XRD processing, fingerprint generation, inversion methodology, and uncertainty/resolution analysis, please refer to:

**[fdedonck/Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion)**

---

# 2. Requirements

The workflow was developed in Python 3 and uses the following packages:

```text
rasterio
richdem
geopandas
pandas
numpy
matplotlib
scipy
scikit-learn
```

They can be installed using, for example:

```bash
pip install rasterio geopandas pandas numpy matplotlib scipy scikit-learn richdem
```

Depending on the operating system, installing `rasterio`, `geopandas`, and `richdem` may be easier using `conda` or `mamba`.

---

# 3. Repository structure

The expected directory structure is:

```text
.
├── INPUT/
│   ├── DEM_lowres.tif
│   ├── glacier_mask.tif
│   ├── gornergletscher_2018.shp
│   ├── gorner_parts.shp
│   ├── hillshade_clipped.tif
│   ├── USurf_lowres.tif
│   ├── weighted_median.tif
│   ├── bed_DEM_lowres.tif
│   ├── gorner_parts.tif
│   ├── geology.tif
│   └── BaseData_GS1.csv
│
├── OUTPUT/
│   └── RESULTS/
│
└── <analysis_script>.py
```

The `RESULTS` directory is created automatically if it does not already exist.

---

# 4. Input data

All datasets required by the script are stored in the `INPUT/` directory.

## `DEM_lowres.tif`

Low-resolution digital elevation model used as the main surface elevation dataset.

It is used to:

* determine glacier surface elevation;
* calculate altitude distributions;
* derive the flow-distance relationship;
* provide the surface topography used in the analysis.

The script assumes a raster resolution of **300 × 300 m**, consistent with:

```python
pixel_area_km2 = 0.3 * 0.3
```

---

## `glacier_mask.tif`

Raster mask defining the glacier area.

Pixels where:

```python
mask == 1
```

are considered part of the glacier.

The raster mask is combined with the glacier polygon from `gornergletscher_2018.shp` to define the final analysis area.

---

## `gornergletscher_2018.shp`

Polygon shapefile defining the Gornergletscher outline.

The glacier polygon is used to:

* mask the DEM;
* define the glacier extent;
* plot the glacier boundary on the provenance maps.

The shapefile should use a coordinate reference system compatible with the raster datasets.

---

## `gorner_parts.shp`

Shapefile containing polygons representing different parts of Gornergletscher.

These glacier parts are used to distinguish different areas of the glacier in the altitude–flow-distance plots.

Each polygon is expected to contain an `id` field corresponding to the glacier-part identifiers in `gorner_parts.tif`.

---

## `hillshade_clipped.tif`

Hillshade raster used for visualisation.

It is not used directly in the calculations. Instead, it is overlaid on the provenance maps to provide topographic context.

---

## `USurf_lowres.tif`

Low-resolution surface-velocity raster.

The dataset is loaded by the script as `Usurf`. It is retained as part of the input dataset for the analysis, but **the current version of the script does not subsequently use `Usurf` in the calculations or plotted quantities**.

It should therefore not be interpreted as directly controlling the circle sizes or provenance calculations in the current workflow.

---

## `weighted_median.tif`

Raster containing the **suspended-load provenance / inferred erosion-rate signal** generated by the XRD inversion workflow described in the companion repository:

**[Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion)**

This raster is therefore an **output of the preceding XRD inversion workflow and an input to the present repository**.

The XRD inversion uses mineralogical fingerprints derived from source and suspended-sediment XRD measurements together with a geological source map to infer spatially variable erosion rates.

In the present repository, the resulting raster is used to:

* map the suspended-load provenance signal;
* investigate how the provenance signal varies with altitude and flow distance;
* compare the suspended-load signal with the independently constructed bedload provenance estimate; and
* scale points in the altitude–flow-distance plots.

The raster should therefore be interpreted as an **inferred spatial provenance/erosion-rate signal**, rather than as a direct measurement of erosion at every individual pixel.

For details on how this raster was obtained, see the companion repository:

[Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion)

---

## `bed_DEM_lowres.tif`

Low-resolution bed elevation DEM.

This dataset is used to derive the flow-distance field.

Pixels with negative values are replaced by the corresponding surface DEM values:

```python
bed_filled = np.where(bed_DEM < 0, DEM, bed_DEM)
```

The resulting DEM is then used in the flow-distance calculation.

---

## `gorner_parts.tif`

Raster representation of the glacier-part polygons.

Each glacier part is represented by a numerical identifier.

These identifiers are used to:

* assign colours to glacier parts;
* separate glacier sections in the scatterplots;
* calculate separate linear regressions for each glacier part.

The identifiers should correspond to the polygons in `gorner_parts.shp`.

---

## `geology.tif`

Raster geological map.

Each geological unit is represented by an integer code from 1–7.

The codes are translated into lithological classes using:

| ID | Lithology                               |
| -: | --------------------------------------- |
|  1 | ZSF ophiolites (serpentinites)          |
|  2 | ZSF sediments                           |
|  3 | Stockhorn, Tuftgrat, Gornergrat         |
|  4 | Monte Rosa (granite)                    |
|  5 | ZSF ophiolites (metabasites, eclogites) |
|  6 | Monte Rosa (gneiss, micaschist)         |
|  7 | Furgg series                            |

These geological units are used to construct the spatially distributed bedload provenance estimate.

---

## `BaseData_GS1.csv`

CSV file containing the measured bedload clast dataset.

The analysis expects at least the following columns:

```text
Class
Unresolved
Mass
Angularity
Size (cm)
```

Rows where:

```text
Unresolved == 1
```

are excluded.

The remaining clasts are used to calculate:

* lithological contribution by clast count;
* lithological contribution by clast mass;
* angularity distributions;
* grain-size distributions;
* mass distributions.

---

# 5. Flow-distance calculation

One of the main purposes of this workflow is to estimate the distance travelled by material through the glacier from the upper glacier towards the snout.

The snout location is defined manually in the script:

```python
x_snout = 624173
y_snout = 91430
```

These coordinates are converted into raster row/column coordinates.

The bed DEM is first filled using `RichDEM`:

```python
rd.FillDepressions(...)
```

A D8 flow-routing calculation is then performed.

The script subsequently propagates distance values upstream from the snout using neighbouring raster cells. The resulting raster represents an estimated **flow distance to the glacier snout**.

The flow distance is finally restricted to the combined glacier mask.

> **Note:** The current implementation uses a relatively simple upstream-neighbour propagation based on elevation. It should therefore be considered a derived/modelled flow-path distance rather than a full physical glacier-flow-routing model.

---

# 6. Bedload provenance calculation

The measured bedload dataset provides the relative contribution of each lithological class.

Two measures are calculated.

## Clast-count contribution

The proportion of resolved clasts belonging to each lithological class is calculated as:

```text
number of clasts of lithology
─────────────────────────────
total number of resolved clasts
```

## Mass contribution

The proportion of the total resolved clast mass represented by each lithology is calculated as:

```text
mass of lithology
─────────────────
total mass of resolved clasts
```

The count-based lithological proportions are subsequently combined with the mapped geological areas to construct a spatially distributed bedload provenance signal.

In simplified form, the calculation assumes that the observed lithological proportions in the bedload can be related to the relative contributions of the mapped geological units.

---

# 7. Suspended-load versus bedload provenance

An important aspect of this repository is the comparison between two provenance approaches.

### Suspended load

The suspended-load signal originates from the **non-linear XRD inversion**.

The inversion uses mineralogical tracer concentrations and spatial geological information to infer a spatially varying erosion-rate/provenance signal.

### Bedload

The bedload signal is independently constructed from the observed lithological composition of bedload clasts and the mapped geological units.

Thus, the two signals are derived from different types of observations:

```text
Suspended load
      │
      ▼
XRD mineralogical fingerprints
      │
      ▼
Non-linear inversion
      │
      ▼
Spatial erosion/provenance signal


Bedload
      │
      ▼
Clast lithology + mass/count
      │
      ▼
Geological map
      │
      ▼
Spatial bedload provenance signal
```

The present repository brings these two approaches together to investigate how inferred sediment provenance varies with glacier flow distance and altitude.

---

# 8. Generated output

All figures are saved in:

```text
OUTPUT/RESULTS/
```

The figures are saved as **SVG files**, preserving vector graphics and making them suitable for subsequent editing in applications such as Inkscape or Adobe Illustrator.

## `AngularityKDE.svg`

Probability distributions of clast angularity for the different lithological classes.

The distributions are calculated using kernel density estimation (KDE).

---

## `logMassKDE.svg`

Kernel density distributions of:

```text
log10(Mass)
```

for the different lithological classes.

The logarithmic transformation is used because clast mass can span several orders of magnitude.

---

## `SizeKDE.svg`

Kernel density distributions of clast size for the different lithological classes.

---

## `Angularity_FlowDistance.svg`

Scatterplot showing the relationship between:

* mean flow distance of each lithological class; and
* the probability of observing angular clasts.

Each lithology is shown using its assigned geological colour.

This figure is intended to explore whether clast angularity varies systematically with inferred transport distance through the glacier.

---

## `Suspended_Bedload_test.svg`

Combined figure comparing suspended-load and bedload provenance signals as a function of:

* flow distance; and
* altitude.

Different colours represent different glacier parts.

The point size represents the magnitude of the provenance signal used in each panel.

The figure also contains altitude and flow-distance distributions along the margins of the scatterplots.

Linear regressions are calculated separately for each glacier part when sufficient variation is present.

---

## `Bedload_Suspendedload_maps.svg`

Two-panel map showing:

**Left:** suspended-load provenance signal derived from the XRD inversion.

**Right:** bedload provenance signal derived from the bedload observations and geological map.

Both maps include:

* hillshade;
* Gornergletscher outline;
* glacier-part boundaries.

---

## `Count_Mass.svg`

Comparison of lithological contributions calculated from:

1. the **number of clasts**; and
2. the **total mass of clasts**.

This provides a direct comparison between count-based and mass-based estimates of lithological contribution.

---

# 9. Lithological colours

The analysis uses a fixed colour scheme for the geological classes:

```python
lithology_colors = {
    "ZSF ophiolites (serpentinites)": "#e3dd19",
    "ZSF sediments": "#809847",
    "Stockhorn, Tuftgrat, Gornergrat": "#ff5001",
    "Monte Rosa (granite)": "#fa9b9a",
    "ZSF ophiolites (metabasites, eclogites)": "#32a02d",
    "Monte Rosa (gneiss, micaschist)": "#ff7f41",
    "Furgg series": "#badd68"
}
```

The glacier parts use a separate colour scheme defined by the `colors` list at the beginning of the script.

---

# 10. Running the analysis

Place the input datasets in the `INPUT/` directory and make sure the filenames match those specified in the script.

Then run:

```bash
python <analysis_script>.py
```

The script will:

1. load the raster and vector datasets;
2. construct the glacier analysis mask;
3. calculate flow distance to the snout;
4. load and process the bedload observations;
5. calculate lithological count and mass contributions;
6. construct the bedload provenance signal;
7. load the suspended-load provenance signal generated by the XRD inversion workflow;
8. calculate sediment-property distributions;
9. generate altitude–flow-distance plots;
10. generate provenance maps;
11. save all figures to `OUTPUT/RESULTS/`.

The script also displays the figures interactively using `matplotlib`.

---

# 11. Important assumptions and parameters

Several parameters are defined manually at the beginning of the script.

## Snout coordinates

```python
x_snout = 624173
y_snout = 91430
```

These coordinates determine the endpoint of the flow-distance calculation.

## Pixel area

```python
pixel_area_km2 = 0.3 * 0.3
```

The analysis assumes a raster resolution of:

```text
300 m × 300 m
```

corresponding to:

```text
0.09 km² per pixel
```

If the input rasters have a different resolution, this parameter must be changed.

## Geological class definitions

The mapping between numerical geological codes and lithological names is defined in:

```python
geo_names
```

The geological raster must therefore use the same numerical coding.

---

# 12. Limitations of the suspended-load provenance inversion

The suspended-load provenance raster inherits the assumptions and limitations of the XRD inversion workflow.

The companion repository identifies several important requirements:

* source areas need sufficiently distinct mineralogical signatures;
* accurate spatial delineation of the source areas is important;
* multiple suitable tracer minerals are required;
* sediment storage should be minimal for the provenance signal to represent contemporary source contributions.

The inferred erosion/provenance map should therefore be interpreted in the context of these assumptions.

The XRD inversion additionally provides posterior covariance and resolution information that can be used to assess the uncertainty and spatial resolving power of the inferred erosion-rate maps.

For a detailed description of these aspects, see the companion repository:

[Non-Linear-XRD-Inversion](https://github.com/fdedonck/Non-Linear-XRD-Inversion)

---

# 13. Reproducibility

For reproducibility, the following should be kept consistent:

* raster resolution and alignment;
* raster coordinate reference system;
* glacier mask;
* glacier outline;
* glacier-part identifiers;
* geological class identifiers;
* snout coordinates;
* lithological classification in the bedload dataset;
* input CSV column names;
* suspended-load provenance raster generated by the XRD inversion workflow.

The workflow assumes that the raster datasets are spatially aligned and have compatible dimensions and coordinate reference systems.

Because the suspended-load provenance raster is generated in a separate repository, reproducing the complete analysis from raw XRD measurements requires running the **Non-Linear-XRD-Inversion** workflow first.

---

# 14. Overall workflow

The complete analysis can therefore be viewed as a two-stage workflow.

```text
                 STAGE 1
        XRD PROVENANCE INVERSION
                 │
                 │
   ┌─────────────┴─────────────┐
   │                           │
Source XRD                Sediment XRD
   │                           │
   └─────────────┬─────────────┘
                 │
                 ▼
       Mineralogical fingerprints
                 │
                 +
           Geological map
                 │
                 ▼
        Non-linear inversion
                 │
                 ▼
     Spatial erosion/provenance
              raster
                 │
                 │
                 ▼
              STAGE 2
   GORNERGLETSCHER ANALYSIS
                 │
   ┌─────────────┼─────────────┐
   │             │             │
   ▼             ▼             ▼
Topography   Bedload data   Glacier geometry
   │             │             │
   │             ▼             │
   │       Bedload provenance  │
   │             │             │
   └──────┬──────┴──────┬──────┘
          │             │
          ▼             ▼
    Flow distance   Provenance comparison
          │             │
          └──────┬──────┘
                 ▼
       Sediment characteristics
                 │
                 ▼
        FINAL FIGURES / MAPS
```

In short, the **XRD inversion repository provides the suspended-load provenance information**, while this repository combines that information with glacier topography, flow distance, geological information, and bedload observations to investigate sediment transport and provenance patterns across Gornergletscher.

---

# 15. Related repository

### Non-Linear XRD Inversion

**From XRD data to erosion rate maps**

https://github.com/fdedonck/Non-Linear-XRD-Inversion

This repository contains the XRD processing, tracer-fingerprint generation, forward modelling, and non-linear inversion used to derive the spatially distributed suspended-load erosion/provenance signal used here.
