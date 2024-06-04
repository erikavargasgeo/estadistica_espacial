
import rioxarray as rxr

import geopandas as gpd
import pickle
import matplotlib.pyplot as plt
from shapely.geometry import Point
import esda
from pointpats import centrography, PointPattern
import numpy as np
from matplotlib.patches import Ellipse

def geodataframe_from_csv(csv_file):
    gdf = gpd.read_file(csv_file)
    return gdf

def geodataframe_from_pickle(pickle_file):
    picklefile = open(pickle_file, "rb")
    gdf = pickle.load(picklefile)
    picklefile.close()
    return gdf

def save_geodataframe_as_pickle(gdf, pickle_file):
    picklefile = open(pickle_file, "wb")
    pickle.dump(gdf, picklefile)
    picklefile.close()

def convert_dataframe_csv_to_pickle(csv_file, pickle_file):
    gdf = geodataframe_from_csv(csv_file)
    save_geodataframe_as_pickle(gdf, pickle_file)


# Function to reproject the raster to another CRS and save it as a new raster
def reproject_raster(original_raster_path, new_raster_path, desired_crs):
    original_raster = rxr.open_rasterio(original_raster_path, masked=True).squeeze()
    new_raster = original_raster.rio.reproject(desired_crs)
    new_raster.rio.to_raster(new_raster_path)

