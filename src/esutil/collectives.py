#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#  
#  This file is part of ESPResSo++.
#  
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. 


"""
*******************************
**espressopp.esutil.collectives**
*******************************

"""
import _espressopp
from espressopp import pmi

ResultNone = _espressopp.esutil_Collectives_ResultNone

def locateItem(here):
    """locate the node with here=True (e.g. indicating that data of a
    distributed storage is on the local node). This is a collective
    SPMD function.

    here is a boolean value, which should be True on at most one
    node. Returns on the controller the number of the node with
    here=True, or an KeyError exception if no node had the item,
    i.e. had here=True.
    """
    res = _espressopp.esutil_Collectives_locateItem(here, pmi.CONTROLLER)
    if pmi.isController:
        if res == ResultNone:
            raise IndexError("collectives.locateItem could not find anything")
        return res