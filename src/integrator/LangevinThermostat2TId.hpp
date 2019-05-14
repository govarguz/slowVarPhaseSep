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

// ESPP_CLASS
#ifndef _INTEGRATOR_LangevinThermostat2TId_HPP
#define _INTEGRATOR_LangevinThermostat2TId_HPP

#include "types.hpp"
#include "logging.hpp"
#include "Particle.hpp"
#include "SystemAccess.hpp"

#include "Extension.hpp"
#include "VelocityVerlet.hpp"


#include "boost/signals2.hpp"


namespace espressopp {
  namespace integrator {

    /** Langevin thermostat */

    class LangevinThermostat2TId : public Extension {

      public:

        LangevinThermostat2TId(shared_ptr<System> system);
        virtual ~LangevinThermostat2TId();

        void setType1(int type1);
        int getType1();
        void setType2(int type2);
        int getType2();
        void setMpc(int mpc);
        int getMpc();


        void setGamma1(real gamma1);
        real getGamma1();
        void setGamma2(real gamma2);
        real getGamma2();
        void setTemperature1(real temperature1);
        real getTemperature1();
        void setTemperature2(real temperature2);
        real getTemperature2();

        void setAdress(bool _adress);
        bool getAdress();

        void initialize();

        /** update of forces to thermalize the system */
        void thermalize();
        void thermalizeAdr(); // same as above, for AdResS

        /** very nasty: if we recalculate force when leaving/reentering the integrator,
            a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
            numbers are drawn twice, resulting in a different variance of the random force.
            This is corrected by additional heat when restarting the integrator here.
            Currently only works for the Langevin thermostat, although probably also others
            are affected.
        */
        void heatUp();

        /** Opposite to heatUp */
        void coolDown();

        /** Register this class so it can be used from Python. */
        static void registerPython();

      private:

        boost::signals2::connection _initialize, _heatUp, _coolDown,
                                       _thermalize, _thermalizeAdr;

        void frictionThermo(class Particle&);

        // this connects thermalizeAdr
        void enableAdress();
        bool adress;

        void connect();
        void disconnect();

        int type1;
        int type2;
        real gamma1;        //!< friction coefficient of type1 particles
        real gamma2;        //!< friction coefficient of type2 particles
        real temperature1;  //!< desired user temperature of type1 particles
        real temperature2;  //!< desired user temperature of type2 particles
        int mpc;

        real pref1;  //!< prefactor, reduces complexity of thermalize
        real pref2;  //!< prefactor, reduces complexity of thermalize
        real pref12;  //!< prefactor, reduces complexity of thermalize
        real pref22;  //!< prefactor, reduces complexity of thermalize

        real pref2buffer; //!< temporary to save value between heatUp/coolDown
        real pref22buffer; //!< temporary to save value between heatUp/coolDown

        shared_ptr< esutil::RNG > rng;  //!< random number generator used for friction term

    };
  }
}

#endif
