import numpy as np
import matplotlib.pyplot as plt
import argparse

# Get and plot planetary K index for a determined set of dates

parser = argparse.ArgumentParser(
      description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
		      help='Choose date. Format: yyyy-mm-dd')

parser.add_argument('--datep', type=str, default='2000-01-01',
		      help='Choose date. Format: yyyy-mm-dd')

parser.add_argument('--datepp', type=str, default='2000-01-01',
		      help='Choose date. Format: yyyy-mm-dd')


parser.add_argument("--ftpfile", type=str, default="Q4", 
		    help="choose the file with the corresponding Kp index data")


#Capture data from command line

cmd_args = parser.parse_args()
date = cmd_args.date
datep = cmd_args.datep # Previous day to impact date
datepp = cmd_args.datepp # 2 days before impact date
year = date.split("-")[0]
ftpfile = year+cmd_args.ftpfile+"_DGD.txt"

# Read and load data from Kp index text file

f = open(ftpfile, "r")

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

for d in raw_data:
    k_date = d.split()[0:3]
    kdate = k_date[0]+"-"+ k_date[1]+"-"+k_date[2]
    if((kdate==date)|(kdate==datep)|(kdate==datepp)):
       kp.append(d.split()[-8:])


# Reshape array to be unidimensional

Kp = np.array(kp).reshape(24,)

# Convert array elements from strings to integers

Kp = [int(k) for k in Kp]

# Start plotting. The output will be a bar graph

## Set x coords

x = np.arange(len(Kp))

## Plot bar graph
bar = plt.bar(x, Kp, width=0.5)

## Set graph limits
plt.xlim(0, 24)
plt.ylim(0, 9)

## Set ticks in both axis

plt.xticks([0, 7.5, 15.5], [datepp, datep, date])
plt.yticks(np.arange(10))

## Set vertical lines at the beginning of each day (00:00 UTC)

plt.axvline(x=7.5, ls="--", c="k")
plt.axvline(x=15.5, ls="--", c="k")


## Set different color to bars according to Kp index value

for i in range(24):
    if Kp[i]==4:
       bar[i].set_color("y")
    elif Kp[i] > 4:
       bar[i].set_color("r")
    else:
       bar[i].set_color("g")

## Set label to axis and graph title

plt.ylabel("Kp Index")
plt.title("Estimated Planetary K Index (3 hours data). Begin {} UTC".format(datepp))

# Save graph

plt.savefig(date+" Kp index.pdf")
