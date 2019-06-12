import sys
sys.path.append('..')

import numpy as np
from metatb.lattice import Lattice
from metatb.tbmodel import TBmodel
import matplotlib.pyplot as plt

lats = [[[1,0],[np.cos(np.pi*2/3),np.sin(np.pi*2/3)]]]
positions = [[1/3,2/3]]
groupnums = [17]
lattice = Lattice(lats,positions,groupnums)

tbmodel = TBmodel(lattice)
graph = tbmodel.build(method='guassian',tolerance=0.5)
tbmodel.plotgraph(span=3,fileformat='png')

ks =[[-0.5,0],[0.0, 0.0],[1./3, 1./3],[0.5,0]]
k_vec,k_dist,k_node,evals = tbmodel.bandstructure(ks)

fig, ax = plt.subplots()
for i in range(np.size(evals,0)):
    ax.plot(k_dist,evals[i,:])
ax.set_xticks(k_node)
ax.set_xticklabels(["M","$\Gamma$","K","M"])
ax.set_xlim(k_node[0],k_node[-1])
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")
#plt.savefig('band.png')
plt.show()

bandgaps,diracstatus,diracenergy = tbmodel.calcbandgap()
print('Band gaps info: ',bandgaps)
print('Hold Dirac point: ', diracstatus)
print('Dirac point Enery: ', diracenergy)