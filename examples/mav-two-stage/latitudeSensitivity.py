import pandas as pd
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

N1 = pd.read_excel('LatitudeSensitivity.xlsx', sheetname='Sheet1')

LAT  = N1['Latitude'].values
PROP = N1['PropMass'].values

plt.figure()
plt.style.use('dark_background')
#plt.rcParams.update(plt.rcParamsDefault)
plt.scatter(TOF_y, VINF_kms, c=LC3, cmap='jet',vmin=LC3.min(), vmax=LC3.max())
#plt.scatter(TOF_y, VINF_kms, c=LC3, cmap='jet',vmin=0, vmax=225.0)
cbar= plt.colorbar()
cbar.ax.tick_params(labelsize=16) 
cbar.set_label('Launch C3', labelpad=-30, y=1.05, rotation=0, fontsize=16)
plt.yticks(np.arange(5, 26, step=5),fontsize=16)
plt.xticks(np.arange(4, 16, step=2),fontsize=16)
plt.grid(True,linestyle='-', linewidth=0.2)
plt.ylabel("Arrival "+r'$V_{\infty}$'+' (km/s)',fontsize=20)
plt.xlabel("TOF "+ '(years)' ,fontsize=20)

plt.xticks(fontsize=16)
plt.xlim([4.0,16.0])
#plt.ylim([5.0,25.0])
plt.title('Gravity Assist Trajectories to Neptune 2025-2034',fontsize=16)
plt.show()








