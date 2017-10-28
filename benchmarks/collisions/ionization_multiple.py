import happi

import numpy as np
import matplotlib.pyplot as plt
ln = np.log
plt.ion()

D = []
colors = ["k", "r", "g", "b", "m"]
for elm in ["C", "Al", "Zn", "Sn", "Au"]:
	S1=happi.Open("ionization_multiple"+elm+"1")
	S2=happi.Open("ionization_multiple"+elm+"2")
	
	if S1.valid and S2.valid:
		
		color = colors.pop()
		
		timestep1 = np.round(np.double(S1.namelist.Main.timestep), decimals=1)
		D += [
			S1.ParticleBinning(0,sum={"ekin":[0,1]}, linestyle="-", color=color, label=elm)
		]
		
		timestep2 = int(np.double(S2.namelist.Main.timestep))
		D += [
			S2.ParticleBinning(0,sum={"ekin":[0,1]}, linestyle="", marker=".", color=color )
		]


# Plot simulation result
plt.figure(1, figsize=(6,3.5))
happi.multiPlot(*D)
fig =plt.gcf()
ax = plt.gca()
# Make nicer plot
legend = plt.legend( loc=2, prop={'size':10})
ax.add_artist(legend)
l1, = plt.plot([0],'-k' , label="$\Delta t="+str(timestep1)+"\;\omega_0^{-1}$")
l2, = plt.plot([0],linestyle="", marker=".",color="k", label="$\Delta t="+str(timestep2)+"\;\omega_0^{-1}$")
legend = plt.legend(handles=[l1, l2], loc=9, prop={'size':10})
for t in legend.texts: t.set_verticalalignment('center')
ax.xaxis.labelpad = 3
ax.yaxis.labelpad = 0
ax.set_xlabel("Time /$\omega_0^{-1}$")
ax.set_ylabel("Density of new electrons /$n_c$")
ax.set_title("")


