import numpy as np
from astropy.table import Table
import glob
import statistics as stats
import argparse

# objetivo del programa: Obtener la curva TEC vs tiempo en diferentes fechas previas al evento
# para una misma estación y obtener ~~un promedio~~ la mediana del comportamiento del TEC y restarlo a los datos
# Observar si hay alguna fluctución inusual que sería nuestra detección.


parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")

parser.add_argument("--station", type=str,
                    default="kvtx",
                    help=" Choose station")

parser.add_argument('--date', type=str, default='2000-01-01',
                    help='Choose date. Format: yyyy-mm-dd')




cmd_args = parser.parse_args()
date = cmd_args.date
station = cmd_args.station

directory = "./data/"+date
b_data = glob.glob(directory+"/background/"+station+"*.Std")



Tec_Matrix = []

for b in b_data:
    Tab = Table.read(b, format="ascii")
    Tec_Matrix.append(Tab["col2"])




m = len(Tec_Matrix) # Number of RINEX files read
n = len(Tec_Matrix[0]) # The number of rows in one Std file
Tec_median = np.zeros(n)

for i in range(0, n):
    guante = []
    for j in range(0, m):
        if Tec_Matrix[j][i]=="-":
            guante.append(np.nan)
        else:
            guante.append(float(Tec_Matrix[j][i]))
    Tec_median[i] = stats.median(guante)

# Create table with the median curve

out_table = Table([Tab["col1"], Tec_median], names=("Time (UT)", "Mean vTEC"))
out_table.write(directory+"/"+station+"-mean-TEC.tab", format="ascii", overwrite=True)
