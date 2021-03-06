#ifndef _SYNAPSE_H_INCLUDED_
#define _SYNAPSE_H_INCLUDED_

#include "neuron.h"
class Neuron;

class Synapse {
  protected:
    double ISyn;
    double gmax;
    double erev;
    int idx;
    Neuron *source;
    Neuron *target;
    
  private:
    static int idx_generator;
    int enabled;
    
  public:
    static int varCount;
    int type;
        /* have access to pre- and post-synaptic voltages through Neuron objects */

    int iVarNo;

    Synapse(Neuron *source, Neuron *target, int type, int iniVarNo);
        /* synaptic (or electrode) current */
    virtual void setISyn() = 0;
    virtual ~Synapse();
    
    int getIdx();
    void setIdx(int);
    virtual void setGmax(double gmax)=0;
        // for a depressing/facilitating synapse defined a derivative
        // vector to be passed to the integrator setting a virtual
        // function = 0 means that we MUST derive a class (and implement
        // said function) in order to use it. It's like a base interface class
    virtual int derivative(realtype, N_Vector *, N_Vector *)=0;
    virtual double getISyn()=0;
    void init(N_Vector y, double *iniVars);
    void setTol(N_Vector tol);
};

#endif/* _SYNAPSE_H_INCLUDED_*/
