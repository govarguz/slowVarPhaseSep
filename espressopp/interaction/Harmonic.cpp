/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "python.hpp"
#include "Harmonic.hpp"
#include "FixedPairListInteractionTemplate.hpp"

namespace espressopp {
  namespace interaction {
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void 
    Harmonic::registerPython() {
      using namespace espressopp::python;

      class_< Harmonic, bases< Potential > >
    	("interaction_Harmonic", init< real, real, real >())
	.def(init< real, real, real, real >())
	.add_property("K", &Harmonic::getK, &Harmonic::setK)
	.add_property("r0", &Harmonic::getR0, &Harmonic::setR0)
    	;

      typedef class FixedPairListInteractionTemplate< Harmonic >
        FixedPairListHarmonic;
      class_< FixedPairListHarmonic, bases< Interaction > >
        ("interaction_FixedPairListHarmonic",
           init< shared_ptr<System>, shared_ptr<FixedPairList>, shared_ptr<Harmonic> >())
        .def(init< shared_ptr<System>, shared_ptr<FixedPairListAdress>, shared_ptr<Harmonic> >())
        .def("setPotential", &FixedPairListHarmonic::setPotential)
        .def("getPotential", &FixedPairListHarmonic::getPotential)
        .def("setFixedPairList", &FixedPairListHarmonic::setFixedPairList)
        .def("getFixedPairList", &FixedPairListHarmonic::getFixedPairList);
     ;
    }

  }
}