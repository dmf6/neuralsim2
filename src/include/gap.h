#include "synapse.h"

class Gap : public Synapse {
    
  public:
    Gap(Neuron *source, Neuron *target, int type, int iniVarNo);
    ~Gap();
    
    void setISyn();
    virtual int derivative(realtype, N_Vector *, N_Vector *);
    void setGmax(double gmax);
    double getISyn();
};

