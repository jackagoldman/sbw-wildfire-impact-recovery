import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Load the NetCDF file
ds = xr.open_dataset('./presentation/KEN36_2012_1917.nc')

# Check and print the bands
print("Available bands:", list(ds.data_vars))

# Get the extent
x_min = float(ds.x.min().values)
x_max = float(ds.x.max().values)
y_min = float(ds.y.min().values)
y_max = float(ds.y.max().values)

# Create a figure with subplots for Band1, Band2, and Band13
fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={'projection': ccrs.Mercator()})

bands = ['Band1', 'Band2', 'Band13']
titles = ['Year Before Fire', 'Fire Year', '10 years following fire']

band2_im = None
for i, band_name in enumerate(bands):
    ax = axes[i]
    ax.set_extent([x_min, x_max, y_min, y_max], crs=ccrs.PlateCarree())
    im = ds[band_name].plot.imshow(ax=ax, transform=ccrs.PlateCarree(), cmap='plasma', add_colorbar=False, alpha=0.8)
    ax.set_title(titles[i])
    # Remove the box border
    ax.spines['geo'].set_visible(False)
    if band_name == 'Band2':
        band2_im = im

# Add a single colorbar beside Band13
fig.colorbar(band2_im, ax=axes[2], orientation='vertical', shrink=0.8, pad=0.05)

plt.tight_layout()
plt.savefig('./presentation/rbr_map_ken36_no_basemap.png')
print("Map saved as rbr_map_ken36_no_basemap.png")