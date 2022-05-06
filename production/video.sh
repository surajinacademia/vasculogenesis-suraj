#!/bin/bash

cd codes

#Reads the simulation parameter range file and checks if it exists in the folder or not
INPUT=parameters.csv
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

#Reads the inputs from the parameters file 



while read IF IR spring AA density angle
do

d='IF_'$IF'_IR_'$IR'_spring_'$spring'_AA_'$AA'_density_'$density'_angle_'$angle

if [ -d "${d}" ]; then
echo "Data Exits. Creating animation"
cd $d
cd 1
cd video 


cat>"vid.py"<<END

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.cm as cm
import math
import os
import pandas as pd


files=[]
for (dirpath, dirnames, filenames) in os.walk("frames/"):
	files.extend(["frames/"+filenamest for filenamest in filenames if "DS_Store" not in filenamest])
	break
files = sorted(files)

print('Number of frames:', len(files))
fig = plt.figure(figsize=(20,20), dpi = 100)
ax = plt.gca()
plt.clf()

N = len(files)
def framelist():
    if int(N/100) < 1:
        return files
    else:           
        movie = []
        for i in range(0, len(files), int(N/100)):
            movie.append(files[i])
        return movie



def readData(i):
    plt.clf() ## CLEARS THE PROJECTIONS
    df = pd.read_csv(framelist()[i], sep="\s+",skiprows = 1, names=list(range(300)))
    x = df[2].to_numpy()
    y = df[3].to_numpy()
    u = np.cos(df[4].to_numpy())
    v = np.sin(df[4].to_numpy())
    u_in = -1*np.cos(df[4].to_numpy())
    v_in = -1*np.sin(df[4].to_numpy())
    ## PLOT PARTICLES
    plt.scatter(x, y, c='green', alpha = 0.5, s=500, edgecolors='black')
    plt.quiver(x, y, u, v, pivot = 'mid', minshaft = 4, headwidth = 7)
    plt.quiver(x, y, u_in, v_in, pivot = 'mid', minshaft = 4, headwidth = 7)
    plt.xticks(np.linspace(-28,28,15), fontsize = 20)
    plt.yticks(np.linspace(-28,28,15), fontsize = 20)
    ## DRAW PBC BOX -> RED
    plt.plot([-30,30],[-30,-30],"r-", linewidth = 2.5) ## xmin to xmax at ymin
    plt.plot([-30,30],[30,30],"r-", linewidth = 2.5) ## xmin to xmax at ymax
    plt.plot([-30,-30],[-30,30],"r-", linewidth = 2.5) ## ymin to ymax at xmin
    plt.plot([30,30],[-30,30],"r-", linewidth = 2.5) ## ymin to ymax at xmax

    ## INTERIOR BOX GRIC -> BLACK
    plt.xlim([-31,31])
    plt.ylim([-31,31])

    plt.axis("on")

anim = ani.FuncAnimation(fig,readData,frames=len(framelist()), blit=False)
anim.save("animation.gif", fps=10)#, extra_args=["-vcodec", "libx264"])

END

python vid.py

mv animation.gif 'IF_'${IF}'_IR_'${IR}'_spring_'${spring}'_AA_'${AA}'_density_'${density}'_angle_'${angle}.gif
cp 'IF_'${IF}'_IR_'${IR}'_spring_'${spring}'_AA_'${AA}'_density_'${density}'_angle_'${angle}.gif /Users/surajsahu/Documents/Vasculogensis/Simulation/Surajcodes/22.04.13/animation

cd ..
cd ..
cd ..

fi

done < $INPUT
IFS=$OLDIFS


