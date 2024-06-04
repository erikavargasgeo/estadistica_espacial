import geopandas as gpd
import matplotlib.pyplot as plt
from utils import geodataframe_from_pickle

FIGURE_FORMAT = 'png'

CORRECTED_PICKLE_FILE = "data/d_aves_data_corrected.pkl"

###################################################
# Carga el DataFrame del archivo pickle
###################################################
d_aves_data = geodataframe_from_pickle(CORRECTED_PICKLE_FILE)

poligonos = gpd.read_file("data/limites_provinciales.gpkg")

###################################################
# Procesamiento de los datos
###################################################
GENERO_INTERES = 'Thraupis'

# Lista de años de interés
años = ['2010', '2015', '2020', '2022']

d_aves_data_filtrada = d_aves_data.query("genus == @GENERO_INTERES and year in @años")

# Obtener el número de ocurrencias por año y provincia
conteo_ocurrencias = d_aves_data_filtrada.groupby(['year', 'stateProvi']).size().reset_index(name='ocurrencias')

# Obtener el número de ocurrencias por año de interés en una lista, para cada provincia
conteo_ocurrencias_prov_year = conteo_ocurrencias.groupby(['stateProvi', 'year'])['ocurrencias'].sum().reset_index()

# Imprimir registros de conteo_ocurrencias_prov_year
print(conteo_ocurrencias_prov_year)

# Graficar el número de ocurrencias por provincia y año, con x siendo el año y y el número de ocurrencias
ax = conteo_ocurrencias_prov_year.pivot(index='year', columns='stateProvi', values='ocurrencias').plot(kind='line', marker='o', figsize=(6, 5))
ax.grid(True)
ax.legend(title='Provincia')
ax.set_xlabel('Año')
ax.set_ylabel('Número de ocurrencias')
ax.set_title(f'Avistamientos de {GENERO_INTERES} por provincia y año')

# Mostrar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/incremento_avistamientos.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/incremento_avistamientos.pdf", dpi=300)

plt.show()