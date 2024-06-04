import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


FIGURE_FORMAT = 'png'
CSV_FILE = 'data/d_aves.csv'

# Leer el archivo csv
d_aves = pd.read_csv(CSV_FILE, encoding='ISO-8859-1')

############################################################################
# Conteo de ocurrencias agrupando por año y provincia
############################################################################

# Conteo de ocurrencias agrupando por año y provincia
conteo_ocurrencias = d_aves.groupby(['year', 'stateProvince']).size().reset_index(name='ocurrencias')

# Para cada registro en conteo_ocurrencias, print el año, la provincia y el numero de ocurrencias
for index, row in conteo_ocurrencias.iterrows():
    print(f"Año: {row['year']},\tProvincia: {row['stateProvince']},\t\tNúmero de ocurrencias: {row['ocurrencias']}")

############################################################################
# Grafico de barras con el número de ocurrencias por provincia
############################################################################

# Conteo de ocurrencias agrupando por provincia
conteo_ocurrencias_prov = conteo_ocurrencias.groupby('stateProvince')['ocurrencias'].sum().reset_index()

# Generar una lista de colores
colors = sns.color_palette('Paired', n_colors=len(conteo_ocurrencias_prov))

# Crear un gráfico de barras con un color distinto para cada barra
ax = conteo_ocurrencias_prov.plot(
    kind='bar', 
    x='stateProvince', 
    y='ocurrencias', 
    color=colors, 
    logy=True,
    ylim=(1, 1e6),
    xlabel='Provincia',
    ylabel='Número de ocurrencias',
    title='Ocurrencias registradas en cada provincia',
    legend=False,
    figsize=(8, 4),
)

# Rotar las etiquetas del eje x para mejor visualización
plt.xticks(rotation=0)

# Mostrar el gráfico
if FIGURE_FORMAT == 'png':
    plt.savefig("figures/barras_ocurrencias_provincias.png", dpi=300)
elif FIGURE_FORMAT == 'pdf':
    plt.savefig("figures/barras_ocurrencias_provincias.pdf", dpi=300)

plt.show()

