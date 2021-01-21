import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.signal import argrelextrema

filenames = sys.argv[1:]
for filename in filenames:
    f = open(filename, 'r', encoding='utf-16')
    lines = f.readlines()

    for ctr, line in enumerate(lines):
        if 'StartOfData' in line:
            startOfDataLine = ctr
        elif 'Sample' in line:
            words = line.split()
            if words[0] == 'Sample':
                title = words[1]
     


    data = np.loadtxt(filename,
                      skiprows = startOfDataLine + 1, 
                      delimiter = '\t', 
                      encoding = 'utf-16',
                      )
    heatFlow = data[:,2]
    temp = data[:,1]
    plt.plot(temp,heatFlow, label=filename)


print("Relative Max:", temp[argrelextrema(heatFlow, np.greater)[0]])
print("Relative Min:", temp[argrelextrema(heatFlow, np.less)[0]])
plt.xlabel(r'Temperature ($\degree$C)', fontname='Arial', fontsize=12)
plt.ylabel('Heat Flow (mW)', fontname='Arial', fontsize=12)
plt.title(title, fontname='Arial', fontsize=14)
plt.tick_params(right=True, top=True, which='both', direction='in')
plt.minorticks_on()
#plt.legend()
#plt.grid()
plt.show()
