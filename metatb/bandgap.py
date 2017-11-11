# This file is part of MetaTB.  MetaTB is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# MetaTB is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# A copy of the GNU General Public License should be available
# alongside this source in a file named gpl-3.0.txt.  If not,
# see <http://www.gnu.org/licenses/>.

import numpy as np 

def calcbandgap(model,nk=200,tolerance=1e-6):
    diracstatus = False
    diracenergy = 0

    ks = [[0,0],[1,0],[0,1],[0,0]]
    (k_vec,k_dist,k_node)=model.k_path(ks,nk,report=False)
    bands=model.solve_all(k_vec)
    N = np.size(bands,axis=0)
    maxvalues = np.zeros(N)
    minvalues = np.zeros(N)
    for i in range(N):
        maxvalues[i] = np.max(bands[i,:])
        minvalues[i] = np.min(bands[i,:])

    for i in range(N):
        for j in range(i+1,N):
            if (abs(maxvalues[i] - minvalues[j])<tolerance):
                diracstatus = True
                diracenergy = maxvalues[i]
            if (abs(minvalues[i] - maxvalues[j])<tolerance):
                diracstatus = True
                diracenergy = minvalues[i]
                
    maxvalue = np.max(maxvalues)
    minvalue = np.min(minvalues)

    status = np.ones(int((maxvalue-minvalue)/tolerance))
    for i in range(N):
        lowbound = int((minvalues[i]-minvalue)/tolerance)
        upbound = int((maxvalues[i]-minvalue)/tolerance)
        status[lowbound:upbound] = 0
    
    bandgaps = list()
    possets, = np.where(status==1)

    end = 0
    start = 0
    while end < len(possets)-2:
        if possets[end]+1 == possets[end+1]:
            end += 1
        else:
            startenergy = possets[start]*tolerance + minvalue
            endenergy = possets[end]*tolerance + minvalue
            center = (endenergy + startenergy)/2
            width = endenergy-startenergy
            bandgaps.append([center,width])

            start = end+1

    if len(possets) > 0:
        startenergy = possets[start]*tolerance + minvalue
        endenergy = possets[end]*tolerance + minvalue
        center = (endenergy + startenergy)/2
        width = endenergy-startenergy
        bandgaps.append([center,width])

    bandgaps = np.array(bandgaps)
    return bandgaps,diracstatus,diracenergy