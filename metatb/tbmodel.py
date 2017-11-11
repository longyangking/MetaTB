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
from .pythtb import *  
from .lattice import Lattice
import matplotlib.pyplot as plt

class TBmodel:
    def __init__(self,lattice):
        self.lattice = lattice
        self.model = None
        self.subspacenum = lattice.getsubspacenum()
        self.graph = None

    def getpositions(self):
        return self.lattice.getpositions()

    def build(self,onsite=0,hopping=1.0,method='exp',width=1.0,tolerance=0.0,eps=1e-3):
        lats = self.lattice.getlattices()
        positions = self.lattice.getpositions()
        sites = self.lattice.getsites()
        graph = list()

        if len(lats)==1:
            lat = lats[0]
            orb = sites
            self.model=tb_model(2,2,lat,orb)

            lat = np.array(lat)
            N = len(positions)
            for i in range(N):
                for j in range(i,N):
                    positioni = positions[i,:]
                    positionj0 = positions[j,:]
                    values = list()
                    directions = list()
                    for m in range(-1,2):
                        for n in range(-1,2):
                            if (i==j) and (m==0) and (n==0):
                                continue
                            positionj = positionj0 + m*lat[0] + n*lat[1]
                            if method == 'guassian':
                                t = hopping*np.exp(-(np.sqrt(np.sum(np.square(positioni-positionj)))/width)**2)
                            else:
                                t = hopping*np.exp(-(np.sqrt(np.sum(np.square(positioni-positionj)))/width))

                            if (len(values) == 0) or (t - values[0] > hopping*eps) or (abs(t-values[0])<hopping*eps):
                                if (len(values) > 0) and (t - values[0] > hopping*eps):
                                    values = list()
                                    directions = list()
                                if [-m,-n] not in directions:
                                    values.append(t)
                                    directions.append([m,n])
                    
                    for m in range(len(values)):
                        if values[m] < tolerance:
                            continue
                        self.model.set_hop(values[m], i, j, directions[m])
                        mi,ni = directions[m]
                        graph.append([positioni,positionj0 + mi*lat[0] + ni*lat[1],values[m]])
                    #graph.append([positioni,positionj0,value])

        else:
            raise Exception("Not realize it yet!")

        self.graph = graph
        return graph

    def bandstructure(self,ks,pointsnum=100):
        (k_vec,k_dist,k_node)=self.model.k_path(ks, pointsnum, report=False)
        evals=self.model.solve_all(k_vec)
        return k_vec,k_dist,k_node,evals

    def plotgraph(self,span=5,graph=None,filename=None,fileformat='pdf'):
        if graph is None:
            graph = self.graph

        fig = plt.figure()
        positions = self.lattice.supercell(span,span)

        lats = self.lattice.getlattices()
        if len(lats)==1:
            lat = np.array(lats[0])
            style = 'y--'
            xs = [0,lat[0,0]]
            ys = [0,lat[0,1]]
            plt.plot(xs,ys,style)
            xs = [0,lat[1,0]]
            ys = [0,lat[1,1]]
            plt.plot(xs,ys,style)
            xs = [lat[0,0],lat[0,0]+lat[1,0]]
            ys = [lat[0,1],lat[0,1]+lat[1,1]]
            plt.plot(xs,ys,style)
            xs = [lat[1,0],lat[0,0]+lat[1,0]]
            ys = [lat[1,1],lat[0,1]+lat[1,1]]
            plt.plot(xs,ys,style)

        plt.scatter(positions[:,0],positions[:,1],s=100)
        plt.xlim([np.min(positions[:,0]),np.max(positions[:,0])])
        plt.ylim([np.min(positions[:,1]),np.max(positions[:,1])])
        ls = list()
        for i in range(len(graph)):
            start,end,value = np.array(graph[i])
            xs = [start[0],end[0]]
            ys = [start[1],end[1]] 
            plt.text(np.mean(xs), np.mean(ys), np.around(value,decimals=2), fontsize=8,color='m')
            plt.plot(xs,ys,c='r',alpha=value)

        cellpositions = self.lattice.getpositions()
        plt.scatter(cellpositions[:,0],cellpositions[:,1],c='black',s=100)
        if filename is None:
            filename = 'model'
        plt.savefig(filename+'.'+fileformat)

    def calcbandgap(self,nk=350,tolerance=1e-3):
        diracstatus = False
        diracenergies = list()
        degeneratepoints = list()

        ks = [[0,0],[0.5,0],[1.0,0.0],[0.5,0.5],[0,1],[0,0.5],[0,0],[0.5,0.5]]
        (k_vec,k_dist,k_node) = self.model.k_path(ks,nk,report=False)
        bands = self.model.solve_all(k_vec)
        N = np.size(bands,axis=0)
        maxvalues = np.zeros(N)
        minvalues = np.zeros(N)
        for i in range(N):
            maxvalues[i] = np.max(bands[i,:])
            minvalues[i] = np.min(bands[i,:])

        # Find Dirac point
        for i in range(N):
            for j in range(i+1,N):
                if (abs(maxvalues[i] - minvalues[j])<tolerance):
                    diracstatus = True
                    diracenergies.append(maxvalues[i])
                if (abs(minvalues[i] - maxvalues[j])<tolerance):
                    diracstatus = True
                    diracenergies.append(minvalues[i])

        #for k in range(np.size(bands,axis=1)):
        #    for i in range(N):
        #        for j in range(i+1,N):
        #            if (abs(bands[i,k]-bands[j,k])<tolerance) and (bands[i,k] not in degeneratepoints):
        #                degeneratepoints.append(bands[i,k]) 

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
                end += 1

        if len(possets) > 0:
            startenergy = possets[start]*tolerance + minvalue
            endenergy = possets[end]*tolerance + minvalue
            center = (endenergy + startenergy)/2
            width = endenergy-startenergy
            bandgaps.append([center,width])

        bandgaps = np.array(bandgaps)
        return bandgaps,diracstatus,diracenergies