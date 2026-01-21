#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 18:45:52 2025

@author: fiendedoncker
"""

import rasterio
import rasterio.mask
import richdem as rd

import geopandas as gpd
import pandas as pd

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm


import os

from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde


import numpy as np
import numpy as np
from sklearn.metrics import r2_score
from scipy.cluster.vq import kmeans, vq
import matplotlib.pyplot as plt


plt.rcParams['svg.fonttype'] = 'none'

# -------------------------------
# USER PARAMETERS
# -------------------------------
# Set colours etc
colors = ['#332288',"#117733", "#88CCEE", "#DDCC77", "#CC6677",'#882255','#464646']


x_snout = 624173
y_snout = 91430

pixel_area_km2 = 0.3 * 0.3  # 300m x 300m -> km²

lithology_colors = {
    "ZSF ophiolites (serpentinites)": "#e3dd19",
    "ZSF sediments": "#809847",
    "Stockhorn, Tuftgrat, Gornergrat": "#ff5001",
    "Monte Rosa (granite)": "#fa9b9a",
    "ZSF ophiolites (metabasites, eclogites)": "#32a02d",
    "Monte Rosa (gneiss, micaschist)": "#ff7f41",
    "Furgg series": "#badd68"
}

geo_names = {1: 'ZSF ophiolites (serpentinites)',
 2: 'ZSF sediments',
 3: 'Stockhorn, Tuftgrat, Gornergrat',
 4: 'Monte Rosa (granite)',
 5: 'ZSF ophiolites (metabasites, eclogites)',
 6: 'Monte Rosa (gneiss, micaschist)',
 7: 'Furgg series'}


# -------------------------------
# INPUT FILES
# -------------------------------
INPUT_PATH = './INPUT'
dem_path = f"{INPUT_PATH}/DEM_lowres.tif"
mask_path = f"{INPUT_PATH}/glacier_mask.tif"
glacier_shp = f"{INPUT_PATH}/gornergletscher_2018.shp"
parts_shp = f"{INPUT_PATH}/gorner_parts.shp"
hillshade_path = f"{INPUT_PATH}/hillshade_clipped.tif"
vel_path = f"{INPUT_PATH}/USurf_lowres.tif"
prov_path = f"{INPUT_PATH}/weighted_median.tif"
bed_dem_path = f"{INPUT_PATH}/bed_DEM_lowres.tif"
parts_path = f"{INPUT_PATH}/gorner_parts.tif"
geology_path = f"{INPUT_PATH}/geology.tif"
bedload_path = f"{INPUT_PATH}/BaseData_GS1.csv"

OUTPUT_PATH = "./OUTPUT"
output_folder = os.path.join(OUTPUT_PATH, 'RESULTS')
os.makedirs(output_folder, exist_ok=True)

# -------------------------------
# HELPER FUNCTIONS
# -------------------------------
def read_raster(path):
    with rasterio.open(path) as src:
        data = src.read(1)
        profile = src.profile
    return data, profile

def logistic(x, x0, k):
    return 1 / (1 + np.exp(-k*(x-x0)))

def storage_simple(z_rel, c_low, A1, mu1, sigma1, A2, mu2, sigma2):
    peak1 = A1 * np.exp(-0.5*((z_rel - mu1)/sigma1)**2)
    peak2 = A2 * np.exp(-0.5*((z_rel - mu2)/sigma2)**2)
    return c_low + peak1 + peak2

def plot_pdf_numeric_kde(df, column, xlabel, bw_factor, savename):
    fig, ax = plt.subplots(figsize=(5, 5))  # square figure

    for g_cls, d in df.groupby("Class"):
        values = d[column].dropna().values

        if len(values) < 5:
            continue

        kde = gaussian_kde(values)
        kde.set_bandwidth(kde.factor * bw_factor)

        x = np.linspace(values.min(), values.max(), 500)
        y = kde(x)

        ax.plot(x, y, linewidth=2, label=g_cls)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability density")
    ax.legend()
    ax.set_box_aspect(1)

    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, savename + ".svg"))
    plt.show()

def plot_pdf_categorical(df, column, savename):
    pdf = (
        df.groupby(["Class", column])
          .size()
          .groupby(level=0)
          .apply(lambda x: x / x.sum())
          .unstack(fill_value=0)
    )

    pdf.plot(kind="bar", stacked=False)
    plt.ylabel("Probability")
    plt.title(f"PDF of {column}")
    plt.legend(title=column)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, savename +".svg"))
    plt.show()
    



# -------------------------------
# LOAD RASTERS
# -------------------------------
DEM, dem_prof = read_raster(dem_path)
mask, _ = read_raster(mask_path)
Usurf, _ = read_raster(vel_path)
prov, _ = read_raster(prov_path)
bed_DEM, _ = read_raster(bed_dem_path)
glacier_parts,_ = read_raster(parts_path)
geology,_ = read_raster(geology_path)

# -------------------------------
# MASK TO GLACIER AREA
# -------------------------------
parts_shp = gpd.read_file(parts_shp)
glacier = gpd.read_file(glacier_shp)
with rasterio.open(dem_path) as src:
    glacier_mask_geom = [feature["geometry"] for feature in glacier.__geo_interface__["features"]]
    glacier_masked, _ = rasterio.mask.mask(src, glacier_mask_geom, crop=False)
glacier_masked = glacier_masked[0]

# combine masks (both raster and shapefile)
combined_mask = (mask == 1) & np.isfinite(glacier_masked)

# replace bed dem -9999 pixels by DEM pixels
bed_filled = np.where(bed_DEM<0, DEM, bed_DEM)

parts = glacier_parts[combined_mask]
unique_parts = np.unique(parts); unique_parts = unique_parts[unique_parts!=0]
color_map = dict(zip(unique_parts, colors))

# -------------------------------
# HILLSHADE FOR PLOTTING
# -------------------------------
with rasterio.open(hillshade_path) as src:
    hillshade = src.read(1)
    hillshade[hillshade < 0] = np.nan
    hs_extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
 
print('Data loaded')
# -------------------------------
# FLOW PATH DISTANCE TO SINK
# -------------------------------
print('Calculcating flow distance')
# Convert sink coordinates to raster row/col
with rasterio.open(dem_path) as src:
    sink_row, sink_col = src.index(x_snout, y_snout)

# Build DEM and fill sinks
rdem = rd.rdarray(bed_filled, no_data=np.nan)
rdem.projection = dem_prof['crs']
rd.FillDepressions(rdem, epsilon=True, in_place=True)

# Compute D8 flow directions
fd = rd.FlowAccumulation(rdem, method='D8')  # this returns upstream count

# Compute distance to sink manually using a mask
# Initialize distances with NaN
flow_distance = np.full(DEM.shape, np.nan)
flow_distance[sink_row, sink_col] = 0

# Simple propagation: every downstream neighbor gets parent's distance + cellsize
queue = [(sink_row, sink_col)]

# neighbor offsets for D8
neighbors = [(-1,0), (-1,1), (0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1)]

while queue:
    r, c = queue.pop(0)
    for dr, dc in neighbors:
        rr, cc = r+dr, c+dc
        if 0 <= rr < DEM.shape[0] and 0 <= cc < DEM.shape[1]:
            # If neighbor flows *toward* current cell
            # Simple heuristic: neighbor must be higher than current cell
            if np.isnan(flow_distance[rr, cc]) and DEM[rr, cc] >= DEM[r, c]:
                flow_distance[rr, cc] = flow_distance[r, c] + np.hypot(dr, dc)*dem_prof['transform'][0]
                queue.append((rr, cc))

# Mask to glacier area
dist = flow_distance[combined_mask]

print('Flow distance computed')

# -------------------------------
# COMPUTE BEDLOAD CONTRIBUTION
# -------------------------------
# 1. Load data and remove unresolved
df = pd.read_csv(bedload_path, sep=",", encoding="latin1")
df = df[df["Unresolved"] == 0]

# 2. Get data, either by count or by weight
# Total n
n_total = len(df)
count_per_class = df["Class"].value_counts()

# Build output dictionary (force zeros if missing)
count_pct = {
    key: 100 * count_per_class.get(name, 0) / n_total
    for key, name in geo_names.items()
}
d = np.array(list(count_pct.values()))/100

# Total mass of resolved clasts
mass_total = df["Mass"].sum()
mass_per_class = df.groupby("Class")["Mass"].sum()

# Build output dictionary (force zeros if missing)
mass_pct = {
    key: 100 * mass_per_class.get(name, 0.0) / mass_total
    for key, name in geo_names.items()
}

# 3. Get provenance
geology = geology.astype(int)
unique_geol = list(np.unique(geology)); unique_geol.remove(0)


A = np.zeros((np.shape(d))) # this array will hold the area (km2) for every litho

for i in unique_geol:
    A_i = np.size(np.where(geology == i)) * pixel_area_km2
    A[i-1] = A_i # correct for the fact that number on geol map starts from 1
    
prov_bdld_discr = d * np.sum(A) / A

prov_bedload_map = geology.copy()
prov_bedload_map = prov_bedload_map.astype(float)
prov_bedload_map[geology == 0] = np.nan

for i in unique_geol:
    prov_bedload_map[geology == i] = prov_bdld_discr[i-1]
    
# 4. Plot stuff!
plot_pdf_categorical(df, "Angularity", "AngularityKDE")
plot_pdf_numeric_kde(
    df.assign(logMass=np.log10(df["Mass"])),
    "logMass",
    "log10(Mass)", 
    1,
    "logMassKDE")
plot_pdf_numeric_kde(df, "Size (cm)", "Size (cm)", 1, "SizeKDE")

# 5. Angularity and flow distance
pdf = (df.groupby(["Class", 'Angularity'])
      .size()
      .groupby(level=0)
      .apply(lambda x: x / x.sum())
      .unstack(fill_value=0))

rounded_df = (
    pdf[["Anguleux"]]
      .rename(columns={"Anguleux": "p_anguleux"})
      .assign(Class=pdf.index)
      .reset_index(drop=True)
)

rounded_df["Class"] = rounded_df["Class"].apply(
    lambda x: x[0] if isinstance(x, tuple) else x
)

mean_dist = {}

for litho_id, litho_name in geo_names.items():
    mask = geology == litho_id

    if np.any(mask):
        mean_dist[litho_name] = np.nanmean(flow_distance[mask])
    else:
        mean_dist[litho_name] = np.nan

dist_df = (
    pd.DataFrame.from_dict(mean_dist, orient="index", columns=["mean_flow_distance"])
      .reset_index()
      .rename(columns={"index": "Class"})
)

plot_df = rounded_df.merge(dist_df, on="Class", how="inner")

plt.figure(figsize=(15, 5))

for _, row in plot_df.iterrows():
    g_cls = row["Class"]
    plt.scatter(
        row["mean_flow_distance"],
        row["p_anguleux"],
        color=lithology_colors.get(g_cls, "gray"),
        s=120,
        edgecolor="k",
        linewidth=0.8,
        label=g_cls
    )

plt.ylabel("Probability of angular clasts")
plt.xlabel("Mean flow distance (m)")
plt.grid(alpha=0.3)

# avoid duplicate legend entries
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), frameon=False)
plt.savefig(os.path.join(output_folder, 'Angularity_FlowDistance' +".svg"))
plt.tight_layout()
plt.show()


# -------------------------------
# PLOT SCATTERPLOT OF PROV-FLOW-ALT
# -------------------------------
x = flow_distance[combined_mask]
y = bed_filled[combined_mask]
s = prov[combined_mask]
s = np.nan_to_num(s, nan=0.0)
s[s < 0] = 0
parts = glacier_parts[combined_mask]

# ------------------------------------------
# GridSpec with layout: 
# row1 tall, row2 short; col1 narrow, col2 wide
# ------------------------------------------
x_min, x_max = np.nanmin(x), np.nanmax(x)
y_min, y_max = np.nanmin(y), np.nanmax(y)

x_range = x_max - x_min
y_range = y_max - y_min

ratio = x_range / y_range   # ~5.33


fig = plt.figure(figsize=(20, (20/ratio)*2))

gs = fig.add_gridspec(
    3, 2,
    height_ratios=[ratio, 1, ratio],    # top row tall, bottom row thin
    width_ratios=[1/(ratio), ratio],     # left column thin, right column wide
    wspace=0.15,
    hspace=0.15)


# Axes
ax_hist_y = fig.add_subplot(gs[0, 0])   # altitude histogram (vertical, flipped)
ax_scatter = fig.add_subplot(gs[0, 1])  # main scatter
ax_empty = fig.add_subplot(gs[1, 0])    # empty
ax_hist_x = fig.add_subplot(gs[1, 1])   # flow-distance histogram
ax_hist_2 = fig.add_subplot(gs[2,0])
ax_dist = fig.add_subplot(gs[2,1])


# ------------------------------------------
# 1. Suspended load
# ------------------------------------------
slopes_parts = {}
for p in unique_parts:
    mask = (parts == p)

    # --- Scatter ---
    ax_scatter.scatter(
        x[mask], 
        y[mask], 
        s=s[mask]*30, 
        color=color_map[p], 
        label=str(p),
        alpha=0.7
    )

    # --- Linear regression ---
    if np.sum(mask) > 2:
        try:
            # Only attempt if x and y both have some spread
            if np.ptp(x[mask]) > 1e-6 and np.ptp(y[mask]) > 1e-6:
                a, b = np.polyfit(x[mask], y[mask], 1)

                xx = np.linspace(x[mask].min(), x[mask].max(), 200)
                yy = a * xx + b

                ax_scatter.plot(
                    xx, yy,
                    color=color_map[p],
                    linewidth=2,
                    alpha=0.9
                )

                slopes_parts[str(p)] = a
        except np.linalg.LinAlgError:
            # Skip this part if polyfit fails
            print(f"Skipping regression for part {p} due to ill-conditioned data.")



ax_scatter.set_xlabel("Flow distance")
ax_scatter.set_ylabel("Altitude")

# Set the bounds so histograms align automatically
ax_scatter.set_xlim(np.nanmin(x), np.nanmax(x))
ax_scatter.set_ylim(np.nanmin(y), np.nanmax(y))
ax_scatter.grid(True, alpha=0.3)

# --- Legend for parts
leg1 = ax_scatter.legend(title="Glacier parts", loc="upper right")
ax_scatter.add_artist(leg1)

# --- Legend for circle sizes
positive_vals = prov[combined_mask]; positive_vals = positive_vals[(positive_vals > 0) & (~np.isnan(positive_vals))]
vmin, vmid, vmax = np.percentile(positive_vals, [5, 50, 99])
sizes = [max(vmin*30,0), max(vmid*30,0), max(vmax*30,0)] 
size_labels = [f"{round(vmin,2)} mm/y", f"{round(vmid,2)} mm/y", f"{round(vmax,2)} mm/y"]

size_handles = [
    ax_scatter.scatter([-1], [-1], s=s_i, color="gray", alpha=0.6, label=l)
    for s_i, l in zip(sizes, size_labels)
]

leg2 = ax_scatter.legend(handles=size_handles, title="Circle size = |Prov − ė|", loc="lower right")


# ------------------------------------------
# 2. Altitude histogram (flipped 90°)
# ------------------------------------------
alt_hist, alt_bins = np.histogram(DEM[combined_mask], bins=np.arange(DEM[combined_mask].min(),
                                                                     DEM[combined_mask].max()+40, 40))
alt_hist_area = alt_hist * pixel_area_km2

# Horizontal bars
ax_hist_y.barh(
    alt_bins[:-1],  # y
    alt_hist_area,  # width along x
    height=40,      # bar thickness along y
    color="gray",
    alpha=0.6
)

# Flip so bars grow leftwards
ax_hist_y.invert_xaxis()

# Align y-limits with scatter
ax_hist_y.set_ylim(ax_scatter.get_ylim())
ax_hist_y.set_yticks([])
ax_hist_y.set_xlabel("Area (km²)")
ax_hist_y.grid(True, alpha=0.3)

# ------------------------------------------
# 3. Flow-distance histogram (horizontal)
# ------------------------------------------
dist_hist, dist_bins = np.histogram(flow_distance[combined_mask], bins=np.arange(0, np.nanmax(flow_distance[combined_mask])+200, 200))
dist_hist_area = dist_hist * pixel_area_km2

ax_hist_x.bar(dist_bins[:-1], dist_hist_area, width=200, color="gray", alpha=0.6)

# Align x-limits with scatter
ax_hist_x.set_xlim(ax_scatter.get_xlim())
ax_hist_x.set_ylabel("Area (km²)")
#ax_hist_x.set_yticks([])
ax_hist_x.grid(True, alpha=0.3)

# ------------------------------------------
# 4. Bedload plot
# ------------------------------------------
s = prov_bedload_map[combined_mask]
s = np.nan_to_num(s, nan=0.0)
s[s < 0] = 0

for p in unique_parts:
    mask = (parts == p)

    # --- Scatter ---
    ax_dist.scatter(
        x[mask], 
        y[mask], 
        s=s[mask]*30, 
        color=color_map[p], 
        label=str(p),
        alpha=0.7
    )

    # --- Linear regression ---
    if np.sum(mask) > 2:
        try:
            # Only attempt if x and y both have some spread
            if np.ptp(x[mask]) > 1e-6 and np.ptp(y[mask]) > 1e-6:
                a, b = np.polyfit(x[mask], y[mask], 1)

                xx = np.linspace(x[mask].min(), x[mask].max(), 200)
                yy = a * xx + b

                ax_dist.plot(
                    xx, yy,
                    color=color_map[p],
                    linewidth=2,
                    alpha=0.9
                )

                slopes_parts[str(p)] = a
        except np.linalg.LinAlgError:
            # Skip this part if polyfit fails
            print(f"Skipping regression for part {p} due to ill-conditioned data.")



# --- Legend for parts
leg1 = ax_dist.legend(title="Glacier parts", loc="upper right")
ax_dist.add_artist(leg1)

# --- Legend for circle sizes
positive_vals = prov_bedload_map[combined_mask]; positive_vals = positive_vals[(positive_vals > 0) & (~np.isnan(positive_vals))]
vmin, vmid, vmax = np.percentile(positive_vals, [5, 50, 99])
sizes = [max(vmin*30,0), max(vmid*30,0), max(vmax*30,0)] 
size_labels = [f"{round(vmin,2)} mm/y", f"{round(vmid,2)} mm/y", f"{round(vmax,2)} mm/y"]

size_handles = [
    ax_scatter.scatter([-1], [-1], s=s_i, color="gray", alpha=0.6, label=l)
    for s_i, l in zip(sizes, size_labels)
]

leg2 = ax_dist.legend(handles=size_handles, title="Circle size = |Prov − ė|", loc="lower right")

ax_dist.grid(True, alpha=0.3)
ax_dist.set_xlabel("Flow distance")
ax_dist.set_ylabel("Altitude")

ax_dist.set_xlim(ax_scatter.get_xlim())


ax_empty.axis("off")


# ------------------------------------------
# 5. Altitude plot lower left
# ------------------------------------------
# Horizontal bars
ax_hist_2.barh(
    alt_bins[:-1],  # y
    alt_hist_area,  # width along x
    height=40,      # bar thickness along y
    color="gray",
    alpha=0.6
)

# Flip so bars grow leftwards
ax_hist_2.invert_xaxis()

# Align y-limits with scatter
ax_hist_2.set_ylim(ax_dist.get_ylim())
ax_hist_2.set_yticks([])
ax_hist_2.set_xlabel("Area (km²)")
ax_hist_2.grid(True, alpha=0.3)

plt.savefig(os.path.join(output_folder, "Suspended_Bedload_test.svg"))
plt.show()



# ------------------------------------------
# PLOT PROVENANCE MAPS
# ------------------------------------------
# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)

# --- Left: suspended load provenance
prov_to_plot,_ = read_raster(prov_path)
vmin, vmax = np.nanpercentile(prov, [1, 99])
im1 = axes[0].imshow(prov_to_plot, extent=hs_extent, cmap="Purples", vmin = 0, vmax = 2.5)
axes[0].imshow(hillshade, extent=hs_extent, cmap="gray", alpha=0.4)
glacier.boundary.plot(ax=axes[0], edgecolor="black", linewidth=0.8)
for _, row in parts_shp.iterrows():
    pid = float(row["id"])
    color = color_map.get(pid, "black") #fallback to black in case of problem
    geom = row.geometry.boundary
    x, y = geom.xy
    axes[0].plot(x, y, color=color, linewidth=3)
axes[0].set_title("Suspended Load Provenance signal", fontsize=11)
axes[0].set_xticks([]); axes[0].set_yticks([])
cbar1 = plt.colorbar(im1, ax=axes[0], orientation="horizontal", fraction=0.046, shrink=0.5, pad=0.07)
cbar1.set_label("Provenance (mm/y)")

# Right: bedload provenance
im2 = axes[1].imshow(prov_bedload_map, extent=hs_extent, cmap="Purples", vmin = 0, vmax = 1.5)
axes[1].imshow(hillshade, extent=hs_extent, cmap="gray", alpha=0.4)
glacier.boundary.plot(ax=axes[1], edgecolor="black", linewidth=0.8)
for _, row in parts_shp.iterrows():
    pid = float(row["id"])
    color = color_map.get(pid, "black") #fallback to black in case of problem
    geom = row.geometry.boundary
    x, y = geom.xy
    axes[1].plot(x, y, color=color, linewidth=3)

axes[1].set_title("Bedload Provenance signal", fontsize=11)
axes[1].set_xticks([]); axes[1].set_yticks([])
cbar2 = plt.colorbar(im2, ax=axes[1], orientation="horizontal", fraction=0.046, shrink=0.5, pad=0.07)
cbar2.set_label("Provenance (relative)")

# Save and show
plt.savefig(os.path.join(output_folder, "Bedload_Suspendedload_maps.svg"))
plt.show()



# ------------------------------------------
# PLOT RAW COUNT DATA
# ------------------------------------------
# Order lithologies
keys = list(geo_names.keys())
labels = [geo_names[k] for k in keys]

count_vals = np.array([count_pct.get(k, 0) for k in keys])
mass_vals  = np.array([mass_pct.get(k, 0)  for k in keys])

x = np.arange(len(keys))
width = 0.35

fig, ax1 = plt.subplots(figsize=(8, 4))

# Count %
bars1 = ax1.bar(
    x - width/2,
    count_vals,
    width,
    label="Count %",
    edgecolor="black"
)
ax1.set_ylabel("Count (%)")
ax1.set_ylim(0, 100)

# Mass %
ax2 = ax1.twinx()
bars2 = ax2.bar(
    x + width/2,
    mass_vals,
    width,
    label="Mass %",
    edgecolor="black",
    facecolor = 'red'
)
ax2.set_ylabel("Mass (%)")
ax2.set_ylim(0, 100)

# X axis
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=30, ha="right")

# Legend (combined)
handles = [bars1, bars2]
labels_legend = ["Count %", "Mass %"]
ax1.legend(handles, labels_legend, frameon=False)

plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'Count_Mass' +".svg"))
plt.show()




