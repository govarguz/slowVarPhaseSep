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
#include "LangevinThermostat2TId.hpp"

#include "types.hpp"
#include "System.hpp"
#include "storage/Storage.hpp"
#include "iterator/CellListIterator.hpp"
#include "esutil/RNG.hpp"

namespace espressopp {

  namespace integrator {

    using namespace espressopp::iterator;


    LangevinThermostat2TId::LangevinThermostat2TId(shared_ptr<System> system)
    :Extension(system) {

      type = Extension::Thermostat;

      gamma1  = 0.0;
      gamma2  = 0.0;
      temperature1 = 0.0;
      temperature2 = 0.0;
      type1 = 0;
      type2 = 1;
      mpc = 1;

      adress = false;

      if (!system->rng) {
        throw std::runtime_error("system has no RNG");
      }

      rng = system->rng;

      LOG4ESPP_INFO(theLogger, "Langevin2T constructed");


    }
    void LangevinThermostat2TId::setType1(int _type1)
    {
      type1 = _type1;
    }

    int LangevinThermostat2TId::getType1()
    {
      return type1;
    }
    void LangevinThermostat2TId::setType2(int _type2)
    {
      type2 = _type2;
    }

    int LangevinThermostat2TId::getType2()
    {
      return type2;
    }

    void LangevinThermostat2TId::setGamma1(real _gamma1)
    {
      gamma1 = _gamma1;
    }

    real LangevinThermostat2TId::getGamma1()
    {
      return gamma1;
    }
    void LangevinThermostat2TId::setGamma2(real _gamma2)
    {
      gamma2 = _gamma2;
    }

    real LangevinThermostat2TId::getGamma2()
    {
      return gamma2;
    }

    void LangevinThermostat2TId::setAdress(bool _adress){
        adress = _adress;
    }

    bool LangevinThermostat2TId::getAdress(){
        return adress;
    }

    void LangevinThermostat2TId::setTemperature1(real _temperature1)
    {
      temperature1 = _temperature1;
    }

    real LangevinThermostat2TId::getTemperature1()
    {
      return temperature1;
    }

    void LangevinThermostat2TId::setTemperature2(real _temperature2)
    {
      temperature2 = _temperature2;
    }

    real LangevinThermostat2TId::getTemperature2()
    {
      return temperature2;
    }

    LangevinThermostat2TId::~LangevinThermostat2TId() {
        disconnect();
    }

    void LangevinThermostat2TId::setMpc(int _mpc)
    {
      mpc = _mpc;
    }

    int LangevinThermostat2TId::getMpc()
    {
      return mpc;
    }



    void LangevinThermostat2TId::disconnect() {

        _initialize.disconnect();
        _heatUp.disconnect();
        _coolDown.disconnect();
        _thermalize.disconnect();
        _thermalizeAdr.disconnect();

    }

    void LangevinThermostat2TId::connect() {

        // connect to initialization inside run()
        _initialize = integrator->runInit.connect(
                boost::bind(&LangevinThermostat2TId::initialize, this));

        _heatUp = integrator->recalc1.connect(
                boost::bind(&LangevinThermostat2TId::heatUp, this));

        _coolDown = integrator->recalc2.connect(
                boost::bind(&LangevinThermostat2TId::coolDown, this));

        if (adress) {
            _thermalizeAdr = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat2TId::thermalizeAdr, this));
        }
        else {
            _thermalize = integrator->aftCalcF.connect(
                boost::bind(&LangevinThermostat2TId::thermalize, this));
        }
    }


    void LangevinThermostat2TId::thermalize()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();
      
      CellList cells = system.storage->getRealCells();

      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }
    }

    // for AdResS
    void LangevinThermostat2TId::thermalizeAdr()
    {
      LOG4ESPP_DEBUG(theLogger, "thermalize");

      System& system = getSystemRef();

      // thermalize CG particles
      /*CellList cells = system.storage->getRealCells();
      for(CellListIterator cit(cells); !cit.isDone(); ++cit) {
        frictionThermo(*cit);
      }*/

      // TODO: It doesn't make that much sense to thermalize both CG and AT particles, since CG particles get velocities of AT particles anyway.
      
      // thermalize AT particles
      ParticleList& adrATparticles = system.storage->getAdrATParticles();
      for (std::vector<Particle>::iterator it = adrATparticles.begin();
              it != adrATparticles.end(); it++) {
            frictionThermo(*it);
            
        // Only in hybrid region!          
        /*Particle &at = *it;
        real w = at.lambda();  
        if(w!=1.0 && w!=0.0) {
            //std::cout << "w: " << w << std::endl;
            //std::cout << "pos_x: " << at.position()[0] << std::endl;
            
            frictionThermo(*it);
        }*/           
            
      }
    }
      
    void LangevinThermostat2TId::frictionThermo(Particle& p)
    {
      real massf = sqrt(p.mass());

      // get a random value for each vector component

      Real3D ranval((*rng)() - 0.5, (*rng)() - 0.5, (*rng)() - 0.5);
      if((p.id()/ mpc) % 2 ==type1){
          p.force() += pref1 * p.velocity() * p.mass() +
                       pref2 * ranval * massf;
      }
      else{
          if((p.id()/ mpc) % 2==type2){
              p.force() += pref12 * p.velocity() * p.mass() +
                           pref22 * ranval * massf;
          }
      }
      LOG4ESPP_TRACE(theLogger, "new force of p = " << p.force());
    }

    void LangevinThermostat2TId::initialize()
    { // calculate the prefactors

        real timestep = integrator->getTimeStep();

      LOG4ESPP_INFO(theLogger, "init, timestep = " << timestep <<
            ", gamma1 = " << gamma1 <<
            ", temperature1 = " << temperature1 <<
            ", gamma2 = " << gamma2 <<
            ", temperature2 = " << temperature2);

      pref1 = -gamma1;
      pref2 = sqrt(24.0 * temperature1 * gamma1 / timestep);
      pref12 = -gamma2;
      pref22 = sqrt(24.0 * temperature2 * gamma2 / timestep);


    }

    /** very nasty: if we recalculate force when leaving/reentering the integrator,
	a(t) and a((t-dt)+dt) are NOT equal in the vv algorithm. The random
	numbers are drawn twice, resulting in a different variance of the random force.
	This is corrected by additional heat when restarting the integrator here.
	Currently only works for the Langevin thermostat, although probably also others
	are affected.
    */

    void LangevinThermostat2TId::heatUp()
    {
      LOG4ESPP_INFO(theLogger, "heatUp");

      pref2buffer = pref2;
      pref2       *= sqrt(3.0);
      pref22buffer = pref22;
      pref22       *= sqrt(3.0);
    }

    /** Opposite to heatUp */

    void LangevinThermostat2TId::coolDown()
    {
      LOG4ESPP_INFO(theLogger, "coolDown");

      pref2 = pref2buffer;
      pref22 = pref22buffer;
    }

    /****************************************************
    ** REGISTRATION WITH PYTHON
    ****************************************************/

    void LangevinThermostat2TId::registerPython() {


      using namespace espressopp::python;


      class_<LangevinThermostat2TId, shared_ptr<LangevinThermostat2TId>, bases<Extension> >
        ("integrator_LangevinThermostat2TId", init<shared_ptr<System> >())
        .def("connect", &LangevinThermostat2TId::connect)
        .def("disconnect", &LangevinThermostat2TId::disconnect)
        .add_property("adress", &LangevinThermostat2TId::getAdress, &LangevinThermostat2TId::setAdress)
        .add_property("type1", &LangevinThermostat2TId::getType1, &LangevinThermostat2TId::setType1)
        .add_property("type2", &LangevinThermostat2TId::getType2, &LangevinThermostat2TId::setType2)
        .add_property("gamma1", &LangevinThermostat2TId::getGamma1, &LangevinThermostat2TId::setGamma1)
        .add_property("gamma2", &LangevinThermostat2TId::getGamma2, &LangevinThermostat2TId::setGamma2)
        .add_property("temperature1", &LangevinThermostat2TId::getTemperature1, &LangevinThermostat2TId::setTemperature1)
        .add_property("temperature2", &LangevinThermostat2TId::getTemperature2, &LangevinThermostat2TId::setTemperature2)
        .add_property("mpc", &LangevinThermostat2TId::getMpc, &LangevinThermostat2TId::setMpc)

        ;


    }

  }
}

