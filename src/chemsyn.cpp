#include "chemsyn.h"

ChemSyn::ChemSyn(Neuron *source, Neuron *target, int type, int iniVarNo, 
		 double tauFacil, double tauRec, double tauDecay)
  : Synapse(source, target, type, iniVarNo) {
    this->source = source;
    this->target = target;  
    this->tauFacil=tauFacil;
    this->tauRec=tauRec;
    this->tauDecay=tauDecay;
    ISyn = 0;
    erev = -70;
    k = varCount - iVarNo;
}



ChemSyn::~ChemSyn() {
    cout << "Destroying Chemical Synapse\n";
};

int ChemSyn::derivative(realtype t, N_Vector *y, N_Vector *ydot) {
   
  realtype *yi, *yd;
  N_Vector yn = *y;
  N_Vector ydn = *ydot;
  double ut, xt, yt, zt;
  int flag;

  yi = NV_DATA_S(yn);
  yd = NV_DATA_S(ydn);
  
  ut = yi[k]; xt=yi[k+1]; yt=yi[k+2];
  zt = yi[k+3];
  
  /* del(t-t0) */
  if(source->spiking) {
    flag = 1;
  } else {
    flag = 0;
  }
  yd[k] = -(ut/tauFacil) + 2*(1-ut)*flag;
  yd[k+1] = zt/tauRec - ut*xt*flag;
  yd[k+2] = -(yt/tauDecay) + ut*xt*flag;
  yd[k+3] = (yt/tauDecay) - zt/tauRec;
  return(0);
}

/* compue gap junction current */
void ChemSyn::setISyn() {
        
}

void ChemSyn::setGmax(double g) {
    gmax = g;
}

double ChemSyn::getISyn(){
  return 0.0;
}

void ChemSyn::setErev(double e_rev){
  erev=e_rev;
}
