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

class Group:
    def __init__(self,lattices,groupnum,position):
        self.sites = None
        self.groupnum = groupnum
        self.lattices = lattices

        if position is None:
            x,y = np.random.random(2)
        else:
            x,y = position
        self.setsites(x,y)

    def getpositions(self,sites=None,tolerance=0.05):
        if sites is None:
            sites = self.sites
        positions = sites.dot(self.lattices)
        return positions

    def supercell(self,xspan,yspan):
        sites = None
        for i in range(-xspan,xspan+1):
            for j in range(-yspan,yspan+1):
                if sites is not None:
                    trans = np.ones(self.sites.shape)
                    trans[:,0] = i*trans[:,0]
                    trans[:,1] = j*trans[:,1]
                    sites = np.concatenate((sites, self.sites + trans),axis=0)
                else:
                    trans = np.ones(self.sites.shape)
                    trans[:,0] = i*trans[:,0]
                    trans[:,1] = j*trans[:,1]
                    sites = self.sites + trans
        return self.getpositions(sites)

    def getsites(self):
        return self.sites

    def getlattices(self):
        return self.lattices

    def setlattices(self,lattices):
        self.lattices = lattices

    def setsites(self,x,y,tolerance=0.05):
        if self.groupnum == 1: # p1
            self.p1(x,y)
        if self.groupnum == 2: # p2
            self.p2(x,y)
        if self.groupnum == 3: # p1m1
            self.p1m1(x,y)
        if self.groupnum == 4: # p1g1
            self.p1g1(x,y)
        if self.groupnum == 5: # c1m1
            self.c1m1(x,y)
        if self.groupnum == 6: # p2mm
            self.p2mm(x,y)
        if self.groupnum == 7: # p2mg
            self.p2mg(x,y)
        if self.groupnum == 8: # p2gg
            self.p2gg(x,y)
        if self.groupnum == 9:
            self.c2mm(x,y)
        if self.groupnum == 10:
            self.p4(x,y)
        if self.groupnum == 11:
            self.p4mm(x,y)
        if self.groupnum == 12:
            self.p4gm(x,y)
        if self.groupnum == 13:
            self.p3(x,y)
        if self.groupnum == 14:
            self.p3m1(x,y)
        if self.groupnum == 15:
            self.p31m(x,y)
        if self.groupnum == 16:
            self.p6(x,y)
        if self.groupnum == 17:
            self.p6mm(x,y)

        # Merge the point that can not be distinguished
        r0 = tolerance
        positions = self.getpositions()
        N = np.size(positions,axis=0)
        status = np.ones(N)
        for i in range(N):
            for j in range(i+1,N):
                if np.sqrt(np.sum(np.square(positions[i,:]-positions[j,:]))) < r0:
                    status[i] = 0
                    break  # Find next point
        pos = np.where(status==1)
        self.sites = self.sites[pos,:][0]

    def p1(self,x,y):
        self.sites = np.array([[x,y]])%1

    def p2(self,x,y):
        self.sites = np.array([[x,y],[-x,-y]])%1

    def p1m1(self,x,y):
        self.sites = np.array([[x,y],[-x,y]])%1

    def p1g1(self,x,y):
        self.sites = np.array([[x,y],[-x,y+1/2]])%1

    def c1m1(self,x,y):
        self.sites = np.array([[x,y],[-x,y],
            [x+1/2,y+1/2],[-x+1/2,y+1/2]])%1

    def p2mm(self,x,y):
        self.sites = np.array([[x,y],[-x,-y],
            [-x,y],[x,-y]])%1

    def p2mg(self,x,y):
        self.sites = np.array([[x,y],[-x,-y],
            [-x+1/2,y],[x+1/2,-y]])%1

    def p2gg(self,x,y):
        self.sites = np.array([[x,y],[-x,-y],
            [-x+1/2,y+1/2],[x+1/2,-y+1/2]])%1

    def c2mm(self,x,y):
        self.sites = np.array([
            [x,y],[-x,-y],
            [-x,y],[x,-y],
            [x+1/2,y+1/2],[-x+1/2,-y+1/2],
            [-x+1/2,y+1/2],[x+1/2,-y+1/2]])%1

    def p4(self,x,y):
        self.sites = np.array([
            [x,y],[-x,-y],
            [-y,x],[y,-x]
        ])%1

    def p4mm(self,x,y):
        self.sites = np.array([
            [x,y],[-x,-y],
            [-y,x],[y,-x],
            [-x,y],[x,-y],
            [y,x],[-y,-x]
        ])%1

    def p4gm(self,x,y):
        self.sites = np.array([
            [x,y],[-x,-y],
            [-y,x],[y,-x],
            [-x+1/2,y+1/2],[x+1/2,-y+1/2],
            [y+1/2,x+1/2],[-y+1/2,-x+1/2]
        ])%1

    def p3(self,x,y):
        self.sites = np.array([
            [x,y],[-y,x-y],
            [-x+y,-x]
        ])%1

    def p3m1(self,x,y):
        self.sites = np.array([
            [x,y],[-y,x-y],
            [-x+y,-x],[-y,-x],
            [-x+y,y],[x,x-y]
        ])%1

    def p31m(self,x,y):
        self.sites = np.array([
            [x,y],[-y,x-y],
            [-x+y,-x],[y,x],
            [x-y,-y],[-x,-x+y]
        ])%1

    def p6(self,x,y):
        self.sites = np.array([
            [x,y],[-y,x-y],
            [-x+y,-x],[-x,-y],
            [y,-x+y],[x-y,x]
        ])%1

    def p6mm(self,x,y):
        self.sites = np.array([
            [x,y],[-y,x-y],
            [-x+y,-x],[-x,-y],
            [y,-x+y],[x-y,x],
            [-y,-x],[-x+y,y],
            [x,x-y],[y,x],
            [x-y,-y],[-x,-x+y]
        ])%1