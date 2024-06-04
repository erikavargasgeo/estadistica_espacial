import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from shapely.geometry import Point
import rasterio

from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW

from pykrige.uk import UniversalKriging
from sklearn.metrics import r2_score

from copy import deepcopy

from utils import reproject_raster



FIGURE_FORMAT = 'png'
REPROJECT_RASTERS = False
NUMBER_OF_POINTS = 20


# Ruta al gpkg
muestreo_gpkg = "data/muestreo.gpkg"

# Cargar los datos del archivo GeoPackage
no_calibrados = gpd.read_file(muestreo_gpkg, layer='Puntos de muestreo')
calibrados = gpd.read_file(muestreo_gpkg, layer='sensores calibrados')
poligono_area_estudio = gpd.read_file(muestreo_gpkg, layer='area de estudio')

# Reproyectar poligono_area_estudio y no_calibrados a la misma CRS que calibrados
poligono_area_estudio = poligono_area_estudio.to_crs(calibrados.crs)
no_calibrados = no_calibrados.to_crs(calibrados.crs)

# Extraer los datos de temperatura
temperatura_no_calibrados = no_calibrados[['ta_media']]
temperatura_calibrados= calibrados[['ta_2m1']]

#########################################
# Metodo de interpolacion GWR
#########################################

# Variables explicativas
variables_independientes = no_calibrados[['Altitud', 'Forest_P_2010', 'Forest_P_2000']]

# Extraer las coordenadas
coords = np.column_stack((no_calibrados.geometry.x, no_calibrados.geometry.y))

# Interpolar los valores de temperatura usando GWR
gwr_selector = Sel_BW(coords, temperatura_no_calibrados.values, variables_independientes.values)
gwr_bw = gwr_selector.search(bw_min=2)
gwr_model = GWR(coords, temperatura_no_calibrados.values, variables_independientes.values, bw=gwr_bw)

gwr_model.fit().summary()

gwr_model_copy = deepcopy(gwr_model)
####################################################################################################
# Aumentar el GeoDataFrame 'calibrados' con los valores de Altitud, Forest_P_2000 y Forest_P_2010
####################################################################################################
# Convertir las series de pandas calibrados.geometry.x y calibrados.geometry.y a listas de coordenadas x e y
xs = calibrados.geometry.x.tolist()
ys = calibrados.geometry.y.tolist()

#########################################################################
# Reproyectar los rasters a la misma CRS que el GeoDataFrame 'calibrados'
#########################################################################
if REPROJECT_RASTERS:
    reproject_raster(
        original_raster_path = "data/rasters/dem200.tif", 
        new_raster_path = "data/rasters/dem200_reprojected.tif", 
        desired_crs = calibrados.crs
    )
    reproject_raster(
        original_raster_path = "data/rasters/percent_forest_2000.tif", 
        new_raster_path = "data/rasters/percent_forest_2000_reprojected.tif", 
        desired_crs = calibrados.crs
    )
    reproject_raster(
        original_raster_path = "data/rasters/percent_forest_2010.tif", 
        new_raster_path = "data/rasters/percent_forest_2010_reprojected.tif", 
        desired_crs = calibrados.crs
    )
    reproject_raster(
        original_raster_path = "data/rasters/mascara_area_estudio.tif", 
        new_raster_path = "data/rasters/mascara_area_estudio_reprojected.tif", 
        desired_crs = calibrados.crs
    )

####################################################################################################
# Abrir los rasters reproyectados y extraer los valores en las coordenadas de los puntos
####################################################################################################
with rasterio.open("data/rasters/dem200_reprojected.tif", mode='r') as raster:
    samples_dem = raster.sample(zip(xs, ys))
    calibrados['Altitud'] = [elem[0] for elem in samples_dem]

with rasterio.open("data/rasters/percent_forest_2000_reprojected.tif", mode='r') as raster:
    samples_f2000 = raster.sample(zip(xs, ys))
    calibrados['Forest_P_2000'] = [elem[0] for elem in samples_f2000]

with rasterio.open("data/rasters/percent_forest_2010_reprojected.tif", mode='r') as raster:
    samples_f2010 = raster.sample(zip(xs, ys))
    calibrados['Forest_P_2010'] = [elem[0] for elem in samples_f2010]

# Load raster into variable and extract coordinates
mask = rasterio.open("data/rasters/mascara_area_estudio_reprojected.tif")
min_x, min_y, max_x, max_y = mask.bounds



grid_x_area_estudio = np.linspace(min_x, max_x, NUMBER_OF_POINTS)
grid_y_area_estudio = np.linspace(min_y, max_y, NUMBER_OF_POINTS)

# Crear una malla de puntos para interpolar
X_grid_area_estudio, Y_grid_area_estudio = np.meshgrid(grid_x_area_estudio, grid_y_area_estudio)

# Crear un GeoDataFrame con la malla de puntos
grid_gdf_area_estudio = gpd.GeoDataFrame(
    geometry=[Point(x, y) for x, y in zip(X_grid_area_estudio.flatten(), Y_grid_area_estudio.flatten())],
    crs=calibrados.crs
)


####################################################
# Predicción de los valores de temperatura con GWR
####################################################
# Variables explicativas
variables_independientes = calibrados[['Altitud', 'Forest_P_2010', 'Forest_P_2000']]

# Extraer las coordenadas
coords = np.column_stack((calibrados.geometry.x, calibrados.geometry.y))

# Interpolar los valores de temperatura usando GWR
prediction_output = gwr_model.predict(coords, variables_independientes.values)

# Extraer los valores de temperatura predichos del objeto mgwr.gwr.GWRResults
predicted_temperatures = prediction_output.predictions.flatten().tolist()

# Crear una nueva columna en el GeoDataFrame para las temperaturas predichas
calibrados['predicted_temperature_GWR'] = predicted_temperatures

# Calculo del R-squared entre los valores predichos con GWR y los valores reales
R_squared_GWR = r2_score(calibrados[['ta_2m1']], calibrados['predicted_temperature_GWR'])

# Mostrar la temperatura predicha con GWR en los puntos calibrados
fig, ax = plt.subplots(1, 1)
calibrados.plot(column='predicted_temperature_GWR', ax=ax, legend=True, cmap='coolwarm', vmin=7.5, vmax=26, zorder=1)
poligono_area_estudio.plot(ax=ax, color='none', edgecolor='black', zorder=2)
plt.title('Temperaturas predichas con GWR')

# Guardar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/temperaturas_predichas_GWR_calibrados.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/temperaturas_predichas_GWR_calibrados.pdf", dpi=300)

plt.show()


####################################################
# Generar raster de temperaturas predichas con GWR
####################################################
# Obtener valores de altitud y porcentaje de bosque en la malla de puntos
with rasterio.open("data/rasters/dem200_reprojected.tif", mode='r') as raster:
    samples_dem = raster.sample(zip(X_grid_area_estudio.flatten(), Y_grid_area_estudio.flatten()))
    grid_gdf_area_estudio['Altitud'] = [elem[0] for elem in samples_dem]

with rasterio.open("data/rasters/percent_forest_2000_reprojected.tif", mode='r') as raster:
    samples_f2000 = raster.sample(zip(X_grid_area_estudio.flatten(), Y_grid_area_estudio.flatten()))
    grid_gdf_area_estudio['Forest_P_2000'] = [elem[0] for elem in samples_f2000]

with rasterio.open("data/rasters/percent_forest_2010_reprojected.tif", mode='r') as raster:
    samples_f2010 = raster.sample(zip(X_grid_area_estudio.flatten(), Y_grid_area_estudio.flatten()))
    grid_gdf_area_estudio['Forest_P_2010'] = [elem[0] for elem in samples_f2010]

# Crear una matriz 2D de coordenadas a partir de X_grid e Y_grid
coords_grid = np.column_stack((X_grid_area_estudio.flatten(), Y_grid_area_estudio.flatten()))

# Definir las variables independientes en la malla de puntos
variables_independientes_grid = grid_gdf_area_estudio[['Altitud', 'Forest_P_2010', 'Forest_P_2000']].values

# Interpolar los valores de temperatura usando GWR
prediction_output_grid = gwr_model_copy.predict(coords_grid, variables_independientes_grid)

# Extraer los valores de temperatura predichos del objeto mgwr.gwr.GWRResults
predicted_values_grid = prediction_output_grid.predictions.flatten().tolist()

# Añadir los valores predichos al GeoDataFrame
grid_gdf_area_estudio['predicted_temperature_GWR'] = predicted_values_grid

# Graficar el raster de temperaturas predichas con GWR
fig, ax = plt.subplots(1, 1)
# Graficar area de estudio
poligono_area_estudio.plot(ax=ax, color='none', edgecolor='black', zorder=2)
# grid_gdf_area_estudio.plot(column='predicted_temperature_GWR', ax=ax, legend=True, cmap='coolwarm', vmin=7.5, vmax=26)

# Construir matriz de valores de temperatura predichos
pixel_values = grid_gdf_area_estudio['predicted_temperature_GWR'].values.reshape(NUMBER_OF_POINTS, NUMBER_OF_POINTS)

# Invertir los valores de los píxeles verticalmente para que coincidan con la orientación del gráfico
pixel_values = np.flipud(pixel_values)

# Graficar los valores de los píxeles como una imagen
plt.imshow(pixel_values, cmap='coolwarm', extent=[min_x, max_x, min_y, max_y], vmin=7.5, vmax=26, zorder=1)

# Mosrar la barra de color
plt.colorbar()
plt.title('Temperaturas predichas con GWR')

# Guardar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/raster_temperaturas_predichas_gwrim.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/raster_temperaturas_predichas_gwrim.pdf", dpi=300)

plt.show()




###################################################################################
# Metodo de interpolación Kriging Universal
###################################################################################
x_no_cal = no_calibrados.geometry.x.tolist()
y_no_cal = no_calibrados.geometry.y.tolist()

# Crear un objeto UniversalKriging
uk = UniversalKriging(
    x_no_cal,
    y_no_cal, 
    temperatura_no_calibrados, 
    variogram_model='linear'
)

####################################################################################################
# Predicción de los valores de temperatura con Kriging Universal para los puntos calibrados
####################################################################################################
# Interpolar los valores de temperatura usando Kriging Universal
predicted_temperatures_kriging, _ = uk.execute('points', calibrados.geometry.x, calibrados.geometry.y)

# Crear una nueva columna en el GeoDataFrame para las temperaturas predichas
calibrados['predicted_temperature_Kriging'] = predicted_temperatures_kriging.flatten().tolist()

# Calculo del R-squared entre los valores predichos con Kriging Universal y los valores reales
R_squared_Kriging = r2_score(calibrados[['ta_2m1']], calibrados['predicted_temperature_Kriging'])

# Mostrar los valores de R-squared para ambos métodos
print(f"R-squared GWR: {R_squared_GWR}, R-squared Kriging: {R_squared_Kriging}")

# Mostrar la temperatura predicha con Kriging Universal en los puntos calibrados
fig, ax = plt.subplots(1, 1)
poligono_area_estudio.plot(ax=ax, color='none', edgecolor='black', zorder=2)
calibrados.plot(column='predicted_temperature_Kriging', ax=ax, legend=True, cmap='coolwarm', vmin=7.5, vmax=26, zorder=1)

plt.title('Temperaturas predichas con Kriging Universal')

# Guardar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/temperaturas_predichas_kriging_calibrados.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/temperaturas_predichas_kriging_calibrados.pdf", dpi=300)

plt.show()

####################################################################
# Generar raster de temperaturas predichas con Kriging Universal
####################################################################
# Evaluar el metodo de Kriging en la grilla
predicted_temperatures_kriging, _ = uk.execute('grid', grid_x_area_estudio, grid_y_area_estudio)

# Añadir las temperaturas predichas al GeoDataFrame
grid_gdf_area_estudio['predicted_temperature_Kriging'] = predicted_temperatures_kriging.flatten()

# Graficar el raster de temperaturas predichas con Kriging Universal
fig, ax = plt.subplots(1, 1)
# grid_gdf.plot(column='predicted_temperature_Kriging', ax=ax, legend=True, cmap='coolwarm', vmin=7.5, vmax=26)
poligono_area_estudio.plot(ax=ax, color='none', edgecolor='black', zorder=2)

# Construir matriz de valores de temperatura predichos
pixel_values = grid_gdf_area_estudio['predicted_temperature_Kriging'].values.reshape(NUMBER_OF_POINTS, NUMBER_OF_POINTS)

# Invertir los valores de los píxeles verticalmente para que coincidan con la orientación del gráfico
pixel_values = np.flipud(pixel_values)

# Graficar los valores de los píxeles como una imagen
plt.imshow(pixel_values, cmap='coolwarm', extent=[min_x, max_x, min_y, max_y], vmin=7.5, vmax=26, zorder=1)


# Mosrar la barra de color
plt.colorbar()
plt.title('Temperaturas predichas con Kriging Universal')

# Guardar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/raster_temperaturas_predichas_krigingim.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/raster_temperaturas_predichas_krigingim.pdf", dpi=300)

plt.show()







