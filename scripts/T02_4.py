import os
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
from pointpats import centrography, PointPattern
import numpy as np
from matplotlib.patches import Ellipse
from utils import geodataframe_from_pickle, convert_dataframe_csv_to_pickle, save_geodataframe_as_pickle

FIGURE_FORMAT = 'png'

CSV_FILE = "data/d_aves.csv"
PICKLE_FILE = "data/d_aves_data.pkl"
CORRECTED_PICKLE_FILE = "data/d_aves_data_corrected.pkl"

CORREGIR_ETIQUETADO = False

############################################
# Guardar el GeoDataFrame como un archivo pickle
############################################
if not os.path.exists(PICKLE_FILE):
    convert_dataframe_csv_to_pickle(CSV_FILE, PICKLE_FILE)
else:
    print(f"El archivo {PICKLE_FILE} ya existe")

###################################################
# Carga el DataFrame del archivo pickle
###################################################
d_aves_data = geodataframe_from_pickle(PICKLE_FILE)


###################################################
# Generos de aves más frecuentes
###################################################
# Mostrar los 10 géneros de aves más frecuentes
print(d_aves_data['genus'].value_counts().head(10))


###################################################
# Procesamiento de los datos
###################################################
GENERO_INTERES = 'Thraupis'

# Lista de años de interés
años = ['2010', '2015', '2020', '2022']

####################################################################
# Corrección de etiquetado de provincias de acuerdo a 
# la geolocalización del registro
####################################################################
# Cargar los poligonos de las provincias
poligonos = gpd.read_file("data/limites_provinciales.gpkg")

if CORREGIR_ETIQUETADO:
    # Crear una columna 'geometry' en 'd_aves_data' con los puntos de latitud y longitud
    geometry_data = [Point(xy) for xy in zip(d_aves_data['decimalLon'], d_aves_data['decimalLat'])]
    d_aves_data = gpd.GeoDataFrame(d_aves_data, geometry=geometry_data)

    poligono = {}
    # Obtener todos los poligonos en 'poligonos', iterando sobre cada fila
    for index, row in poligonos.iterrows():
        poligono[row['ADM1_ES'].replace(' ', '-')] = poligonos.iloc[index:index+1]

    # Definir un contador para los registros mal clasificados
    wrongly_classified_counter = 0

    print(f"Iniciando corrección de etiquetado de provincias para {d_aves_data.shape[0]} registros...")
    # Para cada registro en 'd_aves_data' reemplace el valor 'stateProvi' con la clave del poligono para el cual el registro tiene una intersección
    for index, row in d_aves_data.iterrows():
        # Para cada poligono en 'poligono'
        print(f"Procesando registro {index+1} de {d_aves_data.shape[0]}")
        for key in poligono:
            # Si la geometría del registro interseca con el poligono
            if row['geometry'].intersects(poligono[key].iloc[0]['geometry']):
                # Si el valor de 'stateProvi' es diferente de la clave
                if row['stateProvi'] != key:
                    # Incrementar el contador de registros mal clasificados
                    wrongly_classified_counter += 1
                    # Reemplazar el valor de 'stateProvi' con la clave (provincia) actual
                    d_aves_data.at[index, 'stateProvi'] = key
                break

    # Guardar el GeoDataFrame corregido como un archivo pickle
    save_geodataframe_as_pickle(d_aves_data, CORRECTED_PICKLE_FILE)

    print(f"Numero de registros corregidos: {wrongly_classified_counter}")
else:
    # Cargar el GeoDataFrame corregido desde el archivo pickle
    d_aves_data = geodataframe_from_pickle(CORRECTED_PICKLE_FILE)

####################################################################
fig, axs = plt.subplots(1, 4, figsize=(15, 3.5))
plt.subplots_adjust(hspace=0.1, wspace=0.1)

# Color de los puntos
POINT_COLOR = 'purple'
# Color de las provincias
LOJA_COLOR = 'beige'
ZAMORA_COLOR = 'lightgray'

# Color del centro medio
MEAN_CENTER_COLOR = 'blue'
MEAN_CENTER_COLOR_LOJA = 'royalblue'
MEAN_CENTER_COLOR_ZAMORA = 'royalblue'
# Color de la mediana espacial
MEDIAN_CENTER_COLOR = 'red'
MEDIAN_CENTER_COLOR_LOJA = 'red'
MEDIAN_CENTER_COLOR_ZAMORA = 'red'
# Color de las elipses
ELLIPSE_COLOR = 'red'
ELIPSE_COLOR_LOJA = 'green'
ELIPSE_COLOR_ZAMORA = 'green'

# Iterar sobre cada año y calcular la centralidad espacial
for año in años:   

    # Filtrar los registros para el año actual
    d_aves_data_filtrada = d_aves_data.query("genus == @GENERO_INTERES and year == @año")

    # Crear la geometría de los puntos del DataFrame filtrado
    geometry_data = [Point(xy) for xy in zip(d_aves_data_filtrada['decimalLon'], d_aves_data_filtrada['decimalLat'])]

    # Crear el GeoDataFrame con los puntos del DataFrame filtrado
    gdf = gpd.GeoDataFrame(d_aves_data_filtrada, geometry=geometry_data)

    ############################################################################
    # Calcular la centralidad espacial
    ############################################################################

    # Calcular la centralidad espacial con el centro medio y la mediana espacial en la provincia de Loja
    gdf_loja = gdf.query("stateProvi == 'Loja'")
    mean_center_loja = Point(gdf_loja.geometry.x.mean(), gdf_loja.geometry.y.mean())
    median_center_loja = Point(gdf_loja.geometry.x.median(), gdf_loja.geometry.y.median())

    print(f"Centro medio para el año {año} en Loja: {mean_center_loja}")
    print(f"Mediana espacial para el año {año} en Loja: {median_center_loja}")

    # Calcular la centralidad espacial con el centro medio y la mediana espacial en la provincia de Zamora Chinchipe
    gdf_zamora = gdf.query("stateProvi == 'Zamora-Chinchipe'")
    mean_center_zamora = Point(gdf_zamora.geometry.x.mean(), gdf_zamora.geometry.y.mean())
    median_center_zamora = Point(gdf_zamora.geometry.x.median(), gdf_zamora.geometry.y.median())

    print(f"Centro medio para el año {año} en Zamora-Chinchipe: {mean_center_zamora}")
    print(f"Mediana espacial para el año {año} en Zamora-Chinchipe: {median_center_zamora}")




    ############################################################################
    # Calcular la dispersión espacial con el método de las elipses
    ############################################################################
    
    # Extraer coordenadas como una lista de tuplas para la provincia de Loja
    coord_list_loc_loja = np.array(
        [(x, y) for x, y in zip(gdf_loja.geometry.x, gdf_loja.geometry.y)]
    )
    # Extraer las coordenadas de los puntos como un objeto PointPattern para la provincia de Loja
    pp_loc_loja = PointPattern(coord_list_loc_loja)
    # Calcular la elipse de dispersión espacial para la provincia de Loja
    axis_major_loja, axis_minor_loja, rotation_loja = centrography.ellipse(pp_loc_loja.points)
    # Crear la elipse de dispersión espacial para la provincia de Loja
    ellipse_loja = Ellipse(
        xy=(mean_center_loja.x, mean_center_loja.y),
        width=2 * axis_minor_loja,
        height=2 * axis_major_loja,
        angle=np.rad2deg(rotation_loja),
        facecolor='none',
        edgecolor=ELIPSE_COLOR_LOJA,
        linestyle="--",
        label="standard ellipse",
    )


    # Extraer coordenadas como una lista de tuplas para la provincia de Zamora Chinchipe
    coord_list_loc_zamora = np.array(
        [(x, y) for x, y in zip(gdf_zamora.geometry.x, gdf_zamora.geometry.y)]
    )
    # Extraer las coordenadas de los puntos como un objeto PointPattern para la provincia de Zamora Chinchipe
    pp_loc_zamora = PointPattern(coord_list_loc_zamora)
    # Calcular la elipse de dispersión espacial para la provincia de Zamora Chinchipe
    axis_major_zamora, axis_minor_zamora, rotation_zamora = centrography.ellipse(pp_loc_zamora.points)
    # Crear la elipse de dispersión espacial para la provincia de Zamora Chinchipe
    ellipse_zamora = Ellipse(
        xy=(mean_center_zamora.x, mean_center_zamora.y),
        width=2 * axis_minor_zamora,
        height=2 * axis_major_zamora,
        angle=np.rad2deg(rotation_zamora),
        facecolor='none',
        edgecolor=ELIPSE_COLOR_ZAMORA,
        linestyle="--",
        label="standard ellipse",
    )


    #######################################################################
    # Creación de los subplots segun el año
    #######################################################################
    if año == '2010':
        gdf.plot(ax = axs[0], markersize=1.5, color=POINT_COLOR, zorder=2)
        axs[0].set_title(año)
        # axs[0].set_aspect('equal', adjustable = 'datalim')
        # Graficar las provincias con un color distinto para cada una
        for index, row in poligonos.iterrows():
            nombre_provincia = poligonos.iloc[index:index+1]['ADM1_ES'].values[0]
            poligonos.iloc[index:index+1].plot(ax = axs[0], color=LOJA_COLOR if nombre_provincia=='Loja' else ZAMORA_COLOR, edgecolor = 'black', zorder=1, label=nombre_provincia)

        axs[0].plot(mean_center_loja.x, mean_center_loja.y, marker='o', color=MEAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[0].plot(median_center_loja.x, median_center_loja.y, marker='o', color=MEDIAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[0].add_patch(ellipse_loja)
        axs[0].plot(mean_center_zamora.x, mean_center_zamora.y, marker='o', color=MEAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[0].plot(median_center_zamora.x, median_center_zamora.y, marker='o', color=MEDIAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[0].add_patch(ellipse_zamora)

        axs[0].set_ylim([-5.203, -3.159])

  
        
    elif año == '2015':
        gdf.plot(ax = axs[1], markersize = 1.5, color = POINT_COLOR, zorder=2)
        axs[1].set_title(año)
        # axs[1].set_aspect('equal', adjustable='datalim')
        # Graficar las provincias con un color distinto para cada una
        for index, row in poligonos.iterrows():
            nombre_provincia = poligonos.iloc[index:index+1]['ADM1_ES'].values[0]
            poligonos.iloc[index:index+1].plot(ax = axs[1], color=LOJA_COLOR if nombre_provincia=='Loja' else ZAMORA_COLOR, edgecolor = 'black', zorder=1, label=nombre_provincia)

        axs[1].plot(mean_center_loja.x, mean_center_loja.y, marker='o', color=MEAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[1].plot(median_center_loja.x, median_center_loja.y, marker='o', color=MEDIAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[1].add_patch(ellipse_loja)
        axs[1].plot(mean_center_zamora.x, mean_center_zamora.y, marker='o', color=MEAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[1].plot(median_center_zamora.x, median_center_zamora.y, marker='o', color=MEDIAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[1].add_patch(ellipse_zamora)

        axs[1].set_ylim([-5.203, -3.159])
        axs[1].set_yticklabels([])
        axs[1].set_ylabel('')


    elif año == '2020':
        gdf.plot(ax=axs[2], markersize=1.5, color=POINT_COLOR, zorder=2)
        axs[2].set_title(año)
        # axs[2].set_aspect('equal', adjustable='datalim')
        # Graficar las provincias con un color distinto para cada una
        for index, row in poligonos.iterrows():
            nombre_provincia = poligonos.iloc[index:index+1]['ADM1_ES'].values[0]
            poligonos.iloc[index:index+1].plot(ax = axs[2], color=LOJA_COLOR if nombre_provincia=='Loja' else ZAMORA_COLOR, edgecolor = 'black', zorder=1, label=nombre_provincia)

        axs[2].plot(mean_center_loja.x, mean_center_loja.y, marker='o', color=MEAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[2].plot(median_center_loja.x, median_center_loja.y, marker='o', color=MEDIAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[2].add_patch(ellipse_loja)
        axs[2].plot(mean_center_zamora.x, mean_center_zamora.y, marker='o', color=MEAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[2].plot(median_center_zamora.x, median_center_zamora.y, marker='o', color=MEDIAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[2].add_patch(ellipse_zamora)

        axs[2].set_ylim([-5.203, -3.159])
        axs[2].set_yticklabels([])
        axs[2].set_ylabel('')


    elif año == '2022':
        gdf.plot(ax=axs[3], markersize=1.5, color=POINT_COLOR, zorder=2)
        axs[3].set_title(año)
        # axs[3].set_aspect('equal', adjustable='datalim')
        # Graficar las provincias con un color distinto para cada una
        for index, row in poligonos.iterrows():
            nombre_provincia = poligonos.iloc[index:index+1]['ADM1_ES'].values[0]
            poligonos.iloc[index:index+1].plot(ax = axs[3], color=LOJA_COLOR if nombre_provincia=='Loja' else ZAMORA_COLOR, edgecolor = 'black', zorder=1, label=nombre_provincia)

        axs[3].plot(mean_center_loja.x, mean_center_loja.y, marker='o', color=MEAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[3].plot(median_center_loja.x, median_center_loja.y, marker='o', color=MEDIAN_CENTER_COLOR_LOJA, markersize=4, zorder=3)
        axs[3].add_patch(ellipse_loja)
        axs[3].plot(mean_center_zamora.x, mean_center_zamora.y, marker='o', color=MEAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[3].plot(median_center_zamora.x, median_center_zamora.y, marker='o', color=MEDIAN_CENTER_COLOR_ZAMORA, markersize=4, zorder=3)
        axs[3].add_patch(ellipse_zamora)

        axs[3].set_ylim([-5.203, -3.159])
        axs[3].set_yticklabels([])
        axs[3].set_ylabel('')


# Mostrar la leyenda
fig.legend(['Ocurrencias', 'Centro medio', 'Mediana espacial', 'Elipse de desviación estándar'], loc='lower center', ncol=4, fontsize='small')

fig.suptitle(f'Centralidad y dispersión espacial - {GENERO_INTERES}') 

# Mostrar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/centralidad_dispersion_anos.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/centralidad_dispersion_anos.pdf", dpi=300)

plt.show()





