import ROOT
import numpy as np
from ROOT import TChain
from math import isnan,pi
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

print "Reading tree"

# Opening root file
#fname = "ElectronSelectionAna.root"
fname = "../FitSimPhotons.root"

# Creating TChain
chain    = TChain("myana/FitSimPhotons")
chain.Add(fname)

# Printing the number of entries
entries    = chain.GetEntries()
print "Number of entries in the ANA tree: ", entries

MaxEntries = 0
if MaxEntries > 0:
  entries = MaxEntries

nrPMT = 32

hits=[]
sps =[]
simp=[]
time=[]

print "Running over ", entries, " entries."
for jentry in xrange( entries ):                
    chain.GetEntry( jentry )
    hits.append(chain.q_z_hit)
    print "Charge from spacepoints: ",chain.q_z_sps
    sps.append(chain.q_z_sps)
    simp.append(len(chain.simphot_channel))
    print "Number of sim photons: ",simp[-1]
    for t in range(len(chain.simphot_time)):
    	time.append(t)


xy = np.vstack([simp,sps])
z = gaussian_kde(xy)(xy)

idx = z.argsort()
simp, sps, z =  np.asarray(simp)[idx], np.asarray(sps)[idx], z[idx]

plt.scatter(simp, sps, c=z, s=50, edgecolor='')
plt.title("Charge in function of Simulated PE")
plt.xlabel("Simulated Photons [PE]")
plt.ylabel("Sum Hit Integrals obtained through spacepoints [ADC]")
plt.xlim(-100,5000)
plt.ylim(-2000,100000)
plt.tight_layout()

plt.show()


