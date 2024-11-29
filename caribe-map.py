# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This program is meant to get in one piece all the data from the *.shp files and create a single 
# caribbean map and save the info in a single file called caribbean-map.tab

# By Alejandro Tarango-Yong
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load python libraries

from astropy.table import Table
from plotfullmap import plot_map
import shapefile as shp
import numpy as np
import glob


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Go to work !!!

caribe = glob.glob("./Map/*.shp")
X, Y = [], []
for car in caribe:
    sf = shp.Reader(car)
    for shape in sf.shapeRecords():
        for i in shape.shape.points[:]:
#        x = [i[0] for i in shape.shape.points[:]]
#        y = [i[1] for i in shape.shape.points[:]]
            X.append(i[0])     
            Y.append(i[1])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Save info into table

output_table = Table([X, Y], names=("Longitude", "Latitude"))
outfile = "caribbean-map.tab"
output_table.write(outfile, format="ascii", overwrite=True)
print("Map created suucesfully")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
