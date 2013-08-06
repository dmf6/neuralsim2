#ifndef _CHEMSYN_H_INCLUDED
#define _CHEMSYN_H_INCLUDED

#include "synapse.h"

class ChemSyn : public Synapse {
 public:
  int k;
  double tauFacil;
  double tauRec;
  double tauDecay;

  ChemSyn(Neuron *source, Neuron *target, int type, int iniVarNo, double, double, double);
  ~ChemSyn();
  void setISyn();
  virtual int derivative(realtype, N_Vector *, N_Vector *);
  void setGmax(double gmax);
  double getISyn();
  void setErev(double);
};

#endif 
