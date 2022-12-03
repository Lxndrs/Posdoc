import numpy as np
from astropy.table import Table
import glob
import argparse

# objetivo del programa: Obtener la curva TEC vs tiempo para 27 días previos al evento
# para una misma estación y recopilar los datos de todas las curvas TEC en un mismo archivo

# Obtención de parámetros de línea de comandos.
# - Station: Estación de GPS de la cual recopilaremos las curvas TEC
# - date: Fecha en la cual se detectó el meteoro
# - sets: set de estaciones donde se encuentra la estación que deseamos

parser = argparse.ArgumentParser(
	description=""" Choose a file to work""")

parser.add_argument("--station", type=str,
			default="kvtx",
			help=" Choose station")

parser.add_argument('--date', type=str, default='2000-01-01',
			help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--sets", type=str, default="1", help="Choose a set of stations")


cmd_args = parser.parse_args()
date = cmd_args.date
station = cmd_args.station
set_folder = cmd_args.sets

directory = "./data/{}/set{}/previous/".format(date, set_folder) # Directorio del cual extraemos las curvas TEC de
                                                                 # los días previos
b_data = glob.glob(directory+"{}*.Std".format(station))
outtab = Table(names=("FileID", "time", "TEC", "sigma_TEC")) # Iniciamos tabla de salida, donde recopilaremos todas 
                                                             #las curvas TEC
outdir = "./data/{}/set{}/".format(date, set_folder)         # Directorio de la tabla de salida
# Recopilación de datos
for b in b_data:
    tab = Table.read(b, format="ascii")
    time = tab["col1"]
    TEC = tab["col2"]
    s_TEC = tab["col3"]
    for i, tec in enumerate(TEC):
        if tec=="-": # Reemplazamos guiones en archivo original por nan
            guantec =np.nan
        else:
            guantec = float(tec) # Debido a los guiones, python cree que la variable TEC está compuesta por strings (texto)
                                 # por lo que hay que cambiar los elementos de estas variables por números reales (float)
        if s_TEC[i] == "-":
            s_guantec = np.nan # guantec y s_guantec son variables temporales donde almacenaremos la información que queremos
        else:                  # guantec para el TEC y s_guantec para la desviación estándar
            s_guantec = float(s_TEC[i])
        outtab.add_row([int(i+1), time[i], guantec, s_guantec]) # Añadimos una nueva línea al archivo de salida en cada nueva iteración
outtab.write(outdir+"{}-Avg-prev-days-summary.tab".format(station), format="ascii") # Al finalizar, guardamos
#Calcularemos el dTEC en otro programa.
