import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry import box
import matplotlib.colors as colors
from utils import geodataframe_from_pickle

FIGURE_FORMAT = 'png'

CORRECTED_PICKLE_FILE = "data/d_aves_data_corrected.pkl"

GENERO_INTERES = 'Thraupis'

###################################################
# Carga el DataFrame del archivo pickle
###################################################
d_aves_data = geodataframe_from_pickle(CORRECTED_PICKLE_FILE)

poligonos = gpd.read_file("data/limites_provinciales.gpkg")

# Lista de años de interés
años = ['2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020']	

# Filtrar los registros para el genero 'Grallaria'
d_aves_data_filtrada = d_aves_data.query("genus == @GENERO_INTERES and year in @años")

# Asignar la geometría al GeoDataFrame
geometry_data = [Point(xy) for xy in zip(d_aves_data_filtrada['decimalLon'], d_aves_data_filtrada['decimalLat'])]
d_aves_data_filtrada = gpd.GeoDataFrame(d_aves_data_filtrada, geometry=geometry_data)

poligono = {}
# Obtener todos los polígonos en 'poligonos', iterando sobre cada fila
for index, row in poligonos.iterrows():
    poligono[row['ADM1_ES'].replace(' ', '-')] = poligonos.iloc[index:index+1]

# Definir las celdas de la cuadrícula para cada polígono
cells_dict = {}
for key in poligono:
    bounds = poligono[key].bounds

    # Dividir el área del polígono en una cuadrícula de celdas de 10x10
    minx, miny, maxx, maxy = bounds.minx.iloc[0], bounds.miny.iloc[0], bounds.maxx.iloc[0], bounds.maxy.iloc[0]
    cells = []
    for i in range(10):
        for j in range(10):
            minx_cell = minx + (maxx - minx) * i / 10
            miny_cell = miny + (maxy - miny) * j / 10
            maxx_cell = minx + (maxx - minx) * (i + 1) / 10
            maxy_cell = miny + (maxy - miny) * (j + 1) / 10
            cells.append(box(minx_cell, miny_cell, maxx_cell, maxy_cell))

    cells = gpd.GeoDataFrame(geometry=cells)

    cells_dict[key] = cells

# Asignar el número de registros en cada celda
for key in poligono:
    cells = cells_dict[key]
    cells['count'] = 0
    d_aves_filtrada_por_provincia = d_aves_data_filtrada.query("stateProvi == @key")
    for index, row in d_aves_filtrada_por_provincia.iterrows():
        for index2, row2 in cells.iterrows():
            if row['geometry'].intersects(row2['geometry']):
                cells.at[index2, 'count'] = cells.at[index2, 'count'] + 1 if 'count' in cells else 1
                # Detener el bucle si el registro ha sido asignado a una celda
                break

# Imprimir el número de registros en cada celda, imprimirlo como una matriz de 10x10
for key in poligono:
    cells = cells_dict[key]
    print(f"Cells for {key}")
    matrix = np.empty((10, 10))
    for j in range(9, -1, -1):
        for i in range(10):
            print(cells['count'][i*10+j], end='\t')
            matrix[9-j][i] = cells['count'][i*10+j]
        print(f'\n')
    # print(matrix)

fig, axs = plt.subplots(1,2, figsize=(12, 7))

# Graficar los polígonos y la cuadrícula
for i, (key, value) in enumerate(poligono.items()):
    cells = cells_dict[key]
    minx, miny, maxx, maxy = poligono[key].bounds.minx.iloc[0], poligono[key].bounds.miny.iloc[0], poligono[key].bounds.maxx.iloc[0], poligono[key].bounds.maxy.iloc[0]

    poligono[key].plot(ax=axs[i], color='none', edgecolor='blue')
    norm = colors.LogNorm(vmin=1, vmax=2000)
    cells.plot(column='count', ax=axs[i], legend=True, cmap='coolwarm', alpha=0.5, norm=norm)

    # Mostrar el número de registros en cada celda
    for index, row in cells.iterrows():
        x, y = row['geometry'].centroid.x, row['geometry'].centroid.y
        axs[i].text(x, y, int(row['count']), fontsize=6, ha='center', va='center')

    axs[i].set_ylim(miny, maxy)
    axs[i].set_xlim(minx, maxx)
    axs[i].set_xlabel('Longitud')
    axs[i].set_ylabel('Latitud')
    axs[i].set_title(f'Provincia: {key}', fontsize=10)

    # Establecer los ticks para que coincidan con la cuadrícula
    x_ticks = np.linspace(minx, maxx, 11)
    y_ticks = np.linspace(miny, maxy, 11)
    axs[i].set_xticks(x_ticks)
    axs[i].set_yticks(y_ticks)

    # Rotar las etiquetas del eje x 45 grados y establecer las etiquetas de los ticks para incluir solo un decimal
    axs[i].set_xticklabels([f"{x:.2f}" for x in axs[i].get_xticks()], rotation=45)
    axs[i].set_yticklabels([f"{y:.2f}" for y in axs[i].get_yticks()])

fig.suptitle(f'Número de registros de {GENERO_INTERES} por provincia en cuadrantes', fontsize=12)

# Mostrar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/test_cuadrantes.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/test_cuadrantes.pdf", dpi=300)

plt.show()
