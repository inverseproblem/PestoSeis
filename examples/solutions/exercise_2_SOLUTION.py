# # Exercise 2 - Solution

# %matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
import pestoseis.ttimerays as tr

# # Problem Setup
#
# The objective of this exrcise is to trace a series of rays across a few different input velocity models.
#
# First, we start off by importing our data:

i = 1 #1 2 3 4

# +
filename = f'../inputdata/exe2_input_velmod_{i}.npy'
print(filename)

inpdat=np.load(filename,allow_pickle=True).item()
gridpar=inpdat['gridpar']
sources=inpdat['srcs']
receivers=inpdat['recs']
velmod=inpdat['velmod']
# -

# # Compute the Traveltimes
#
# Similar to in the previous exercise, we want to compute the traveltimes for the different source-receiver positions.

ttpick,ttime = tr.traveltime(velmod, gridpar, sources, receivers)

# # Trace the Rays
#
# Now that we have the traveltimes we can trace the rays to see the ray paths throughout the domain.

rays = tr.traceallrays(gridpar, sources, receivers, ttime)

# # Plot the Results
#
# We can not plot the following results:
# - The traveltimes for each receiver receiver position
# - The ray paths travelled through the domain

# +
import matplotlib.gridspec as gridspec

plt.figure(figsize=(9,6))
gs = gridspec.GridSpec(2, 2)
plt.subplot(gs[0, 0]) 
plt.title('Traveltimes')
for i in range(ttpick.size):
    plt.plot(ttpick[i][:],'.-',label='src {}'.format(i))
plt.legend()

plt.subplot(gs[0,1]) 
tr.plotrays(sources,receivers,rays)
tr.plotgrid(gridpar)
#plotttimemod(gridpar,ttime[0])
#extent_ttime = [gridpar['xttmin'],gridpar['xttmax'],
#                gridpar['yttmin'],gridpar['yttmax'] ]
#plt.contour(ttime[0].T,50,extent=extent_ttime,colors='black',linewidth=0.5)


plt.subplot(gs[1,:]) 
tr.plotvelmod(gridpar,velmod)
tr.plotrays(sources,receivers,rays)

name = filename.split("/")[-1].split(".")[0]
plt.savefig("figs/{}.pdf".format(name))
plt.show()
# -


