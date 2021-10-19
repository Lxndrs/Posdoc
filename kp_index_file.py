# The goal of this program is to create a file containing the following data

#     Date   | Kp index median | Max Kp index |     Variation       |
#  ----------+-----------------+--------------+---------------------|
#  2000-01-01|  Number         | Above 4      |Median - min Kp index|
#  2000-01-02|  Number         | 4            |  Number             |
#  2000-01-03|  Number         |  Below 4     |  Number             |

# For use in another program. For this we must recycle some linecodes
# From our other Kp index program

import glob
import seaborn as sns
from astropy.table import Table
import numpy as np
import statistics as stats

def yesterday(date):
    """
    Enters any date and returns the previos day
    (is not as easy as you think). Limited to 1990
    and forward.
    """
    year, month, day = date.split("-") # Get year, month and day from date
    bisiesto = (year =="1992") | (year=="1996") | (year=="2000") | (year=="2004") | (year=="2008") | (year=="2012") | (year=="2016") | (year=="2020") 
    thirty_day_month = (month == "05") | (month == "07") | (month == "10")| (month == "12")
    if day == "01":
        if month == "01": # January 1st case
            pyear = str(int(year)-1)
            pmonth = "12"
            pday = "31"     
        else:
            pyear = year
            pmonth = str(int(month)-1)
            if int(pmonth) <10:
                pmonth="0"+pmonth
            if month == "03": # March 31th case
                if bisiesto== True:
                    pday = "29"
                else:
                    pday = "28"
            elif thirty_day_month ==True: # Previous month has 30 days
                pday = "30"
            else: # Previous month has 31 days
                pday = "31"
    else:
        pyear = year
        pmonth = month
        pday = str(int(day)-1)
        if int(pday) <10:
            pday = "0"+pday
    return(pyear+"-"+pmonth+"-"+pday)

# Create a list with the desired sample (dates)

#sample = ["1995-08-05", "1996-07-12", temporaly remove 1995 and 1996 meteors since their data format is sligthly different
sample = ["1997-10-09", "2000-01-18", "2000-08-25", "2005-11-15", "2015-07-19", "2019-02-01", "2019-05-23", "2019-07-18", "2019-08-10", "2020-04-28", "2019-10-03", "2019-10-09", "2019-11-16", "2019-11-17", "2019-11-19", "2019-11-26", "2019-12-04", "2019-12-15", "2019-12-29", "2020-01-03", "2020-01-06", "2019-06-22", "2021-03-31", "2020-12-29", "2020-01-15", "2020-02-12", "2020-03-03", "2020-03-31", "2020-04-08", "2020-04-18", "2020-04-20", "2020-04-25", "2020-05-08", "2020-07-15", "2020-08-07", "2020-09-13", "2020-09-30", "2020-11-16", "2020-11-17", "2020-12-19", "2020-12-23"] 

# open needed file(s)

#ftpfile = {"1995-08-05":"1995_DGD.txt", "1996-07-12":"1996_DGD.txt", 
ftpfile = {"1997-10-09":"1997_DGD.txt", "2000-01-18":"2000_DGD.txt", "2000-08-25":"2000_DGD.txt", "2005-11-15":"2005_DGD.txt", "2015-07-19":"2015_DGD.txt", "2019-02-01":"2019Q1_DGD.txt", "2019-05-23":"2019Q2_DGD.txt", "2019-07-18":"2019Q3_DGD.txt", "2019-08-10":"2019Q3_DGD.txt", "2020-04-28":"2020Q2_DGD.txt", "2019-10-03":"2019Q4_DGD.txt", "2019-10-09":"2019Q4_DGD.txt", "2019-11-16":"2019Q4_DGD.txt", "2019-11-17":"2019Q4_DGD.txt", "2019-11-19":"2019Q4_DGD.txt", "2019-11-26":"2019Q4_DGD.txt", "2019-12-04":"2019Q4_DGD.txt", "2019-12-15":"2019Q4_DGD.txt", "2019-12-29":"2019Q4_DGD.txt", "2020-01-03":"2020Q1_DGD.txt", "2020-01-06":"2020Q1_DGD.txt", "2019-06-22":"2019Q2_DGD.txt", "2021-03-31":"2021Q1_DGD.txt", "2020-12-29":"2020Q4_DGD.txt", "2020-01-15":"2020Q1_DGD.txt", "2020-02-12":"2020Q1_DGD.txt", "2020-03-03":"2020Q1_DGD.txt", "2020-03-31":"2020Q1_DGD.txt", "2020-04-08":"2020Q2_DGD.txt", "2020-04-18":"2020Q2_DGD.txt", "2020-04-20":"2020Q2_DGD.txt", "2020-04-25":"2020Q2_DGD.txt", "2020-05-08":"2020Q2_DGD.txt", "2020-07-15":"2020Q3_DGD.txt", "2020-08-07":"2020Q3_DGD.txt", "2020-09-13":"2020Q3_DGD.txt", "2020-09-30":"2020Q3_DGD.txt", "2020-11-16":"2020Q4_DGD.txt", "2020-11-17":"2020Q4_DGD.txt", "2020-12-19":"2020Q4_DGD.txt", "2020-12-23":"2020Q4_DGD.txt"}

# Declaring output variables

Kp_median = []
Kp_max = []
Kp_deviation = []


for ftp in ftpfile:
    f = open(ftpfile[ftp], "r")
    
    ## Skip first 12 rows
    
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    
    ## Load data
    
    raw_data = f.readlines()
    
    # Select desired dates from the whole data
    
    kp = []
    Date = ftp
    Datep = yesterday(Date)
    Datepp = yesterday(Datep)
    for d in raw_data:
        k_date = d.split()[0:3]
        kdate = k_date[0]+"-"+ k_date[1]+"-"+k_date[2]
        if((kdate==Date)|(kdate==Datep)|(kdate==Datepp)):
            kp.append(d.split()[-8:])
    
    
    # Reshape array to be unidimensional
    
    Kp = np.array(kp).reshape(24,)
    
    # Convert array elements from strings to integers
    
    Kp = [int(k) for k in Kp]
    kp_median = stats.median(Kp)
    kp_max = max(Kp)
    kp_deviation = kp_median - min(Kp)
        
    Kp_median.append(kp_median)
    Kp_max.append(kp_max)
    Kp_deviation.append(kp_deviation)

#print(len(sample), len(Kp_median), len(Kp_deviation))
t = Table([sample, Kp_median, Kp_max, Kp_deviation], names=("Date", "Median", "Max value", "Variation"))
t.write("Kp_table.csv", format="csv", overwrite=True)
