#ifndef _NEURON_H
#define _NEURON_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <list>
#include <vector>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */


#include "synapse.h"
#include "electrode.h"

#define efunc(X,Y,Z) (1.0/(1.0+exp(((X)-(Y))/(Z))))
#define forall(l, it) for (it= (l).begin(); it != (l).end(); it++)
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define ATOL RCONST(1.0e-14)
#define INF 1.0e14

class Synapse;
class Electrode;

using namespace std;

/* let this be the base class from which all other neuron subclasses will derive methods
 * define parameters: conductances, equilibrium potentials
 * let neuron class define own function derivative y' = f(t,y)
 */
class Neuron {
  protected:
    static int idx_generator;
    int enabled, idx;
    realtype rm, cm;
    int k;
    
  public:
    static int varCount;
    list<Synapse *> den;       /* list of all incoming synapses and electrodes */
    list<Synapse *>::const_iterator den_it;
    vector <Electrode *> electrodes;
    vector<Electrode *>::const_iterator iapp_it;
    vector<double> spiketimes;
    int spiking;
    double spiketime, oldSpiketime;
    double previousTime;
    int increasing;
    double voltage, oldVoltage;
    int iVarNo;
    string name;
    
    int mode;
        //iniVarNo is number of state variables
    Neuron(const string name, int iniVarNo, int mode);
    virtual ~Neuron();
        /* will be redefined by derived neruon classes, to create a
         * pure virtual function we assign it the value of 0. It acts
         * as an interface to describe the the available methods */
    virtual int derivative(realtype t, N_Vector *, N_Vector *, void *user_data) = 0;
        /*use currents function to save all the currents and sum to a file using os stream */
    virtual void currents(realtype t, N_Vector, N_Vector, ostream&) = 0;
    inline int getIdx() { return idx;}
    virtual void init(N_Vector y, int, double *, int);
    virtual void setTol(N_Vector, int, int);
    int getIVarNo();
    
    double getVoltage() {
      return this->voltage;
    }
    void setVoltage(double v) {
      this->voltage = v;
    }
    void detectSpikes(double, N_Vector);
};
#endif
