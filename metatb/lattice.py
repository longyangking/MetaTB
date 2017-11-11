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
from .symmetry import Group

class Lattice:
    def __init__(self,lattices,positions,groupnums):
        self.lattices = lattices
        self.positions = positions
        self.groupnums = groupnums

        self.groups = list()
        for i in range(len(groupnums)):
            lattice = lattices[i]
            groupnum = groupnums[i]
            position = positions[i]
            self.groups.append(Group(lattice,groupnum,position))

        self.sites = None
        self.mainsites = None
        for group in self.groups:
            sites = group.getsites()
            if self.sites is not None:
                self.mainsites = np.concatenate((self.mainsites,sites[0]),axis=0)
                self.sites = np.concatenate((self.sites,sites),axis=0)
            else:
                self.mainsites = sites[0]
                self.sites = sites

    def getlattices(self):
        return self.lattices

    def getpositions(self):
        positions = None
        for group in self.groups:
            if positions is not None:
                position = group.getpositions()
                positions = np.concatenate((positions,position),axis=0)
            else:
                positions = group.getpositions()
        return positions
        
    def getsubspacenum(self):
        return len(self.groupnums)
        
    def supercell(self,xspan,yspan):
        positions = None
        for group in self.groups:
            if positions is not None:
                position = group.supercell(xspan,yspan)
                positions = np.concatenate((positions,position),axis=0)
            else:
                positions = group.supercell(xspan,yspan)
        return positions

    def getsites(self):
        return self.sites

    def deform(self,lattices,positions):
        for i in range(len(positions)):
            x,y = positions[i]
            self.groups[i].setlattices(lattices[i])
            self.groups[i].setsites(x,y)

