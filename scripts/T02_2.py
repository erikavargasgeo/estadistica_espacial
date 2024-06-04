import matplotlib.pyplot as plt
from shapely.geometry import Point                                                                                                                                          
import pandas as pd
import geopandas as gpd

FIGURE_FORMAT = 'png'
CSV_FILE = 'data/d_aves.csv'

# Leer el archivo csv
d_aves = pd.read_csv(CSV_FILE, encoding='ISO-8859-1')

####################################################################
# Exploracion de los datos
####################################################################

# Numero de registros en el archivo csv
num_registros = d_aves.shape[0]
print(f"El numero de registros del archivo csv es: {num_registros}")

# Numero de variables en el archivo csv
num_variables = d_aves.shape[1]
print(f"El numero de variables en el archivo csv es: {num_variables}")

# Tipo de datos de las variables
print(d_aves.dtypes)

# Numero de variables de tipo int64
num_variables_int64 = d_aves.select_dtypes(include='int64').shape[1]
print(f"El numero de variables de tipo int64 en el archivo csv es: {num_variables_int64}")

# Numero de variables de tipo object
num_variables_object = d_aves.select_dtypes(include='object').shape[1]
print(f"El numero de variables de tipo object en el archivo csv es: {num_variables_object}")

# Numero de variables de tipo float64
num_variables_float64 = d_aves.select_dtypes(include='float64').shape[1]
print(f"El numero de variables de tipo float64 en el archivo csv es: {num_variables_float64}")

# Conteo de ocurrencias agrupando por año y provincia
conteo_ocurrencias = d_aves.groupby(['year', 'stateProvince']).size().reset_index(name='ocurrencias')
print(conteo_ocurrencias)


####################################################################
# Mapa con los puntos de ocurrencias
####################################################################

fig, axs = plt.subplots(1,1)

poligonos = gpd.read_file("data/limites_provinciales.gpkg")

# Creación de la geometría de los puntos de ocurrencias para el GeoDataFrame
geometry_data = [Point(xy) for xy in zip(d_aves['decimalLongitude'], d_aves['decimalLatitude'])]

# Crear un GeoDataFrame con los puntos de ocurrencias
ocurrencias_gdf = gpd.GeoDataFrame(d_aves, geometry=geometry_data)

ocurrencias_gdf.plot(
    ax = axs, 
    markersize=1.5, 
    color='red', 
    zorder=2 
)

axs.set_xlabel('Longitud')
axs.set_ylabel('Latitud')
axs.set_title('Puntos de ocurrencias - Avistamiento de aves')
axs.set_aspect('equal', adjustable = 'datalim')
poligonos.plot(ax = axs, color='gainsboro', edgecolor = 'black', zorder=1)

if FIGURE_FORMAT == 'png':
    plt.savefig("figures/puntos_ocurrencias_aves.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/puntos_ocurrencias_aves.pdf", dpi=300)

plt.show()



