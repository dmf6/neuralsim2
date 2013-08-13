#ifndef _NEURONMODEL_H
#define _NEURONMODEL_H

#include <iostream>
#include <list>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */

using namespace std;

class Neuron;
class Synapse;

class NeuronModel
{
  private:
    list<Neuron *> *neurs;
    list<Neuron *>::iterator niter;
    list<Synapse *> *syns;
    list<Synapse *>::iterator siter;
    
  public:
    int iVarCnt;
    NeuronModel(list<Neuron *> *, list<Synapse *> *,
                int &, ostream &); 
    ~NeuronModel(){}
        /* loops over neuron and synapse objects in the list so that
         * RK4 can solve each one at time t */
    
    virtual void derivative(realtype t, N_Vector, N_Vector, void *);
    void currents(realtype t, N_Vector y);
    
};

#endif
