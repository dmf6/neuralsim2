#include "synapse.h"

class ChemSyn : public Synapse {

 public:
  ChemSyn(Neuron *source, Neuron *target, int type, int iniVarNo);
  ~ChemSyn();
  void setISyn();
  double getISyn();
  virtual int derivative(realtype, N_Vector *, N_Vector *);
  void setGmax(double gmax);
  void setErev(double erev);
};
  
