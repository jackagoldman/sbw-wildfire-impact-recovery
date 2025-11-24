import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import contextily as ctx

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

for i, band_name in enumerate(bands):
    ax = axes[i]
    ax.set_extent([x_min, x_max, y_min, y_max], crs=ccrs.PlateCarree())
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, zoom=10)
    ds[band_name].plot.imshow(ax=ax, transform=ccrs.PlateCarree(), cmap='plasma', add_colorbar=True, alpha=0.8)
    ax.set_title(f'{band_name} Map for KEN36_2012_1917')

plt.tight_layout()
plt.savefig('./presentation/rbr_map_ken36_bands_side_by_side.png')
print("Map saved as rbr_map_ken36_bands_side_by_side.png")