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
  //cout << "idx of chemsyn " << idx << "\n";
  ut = yi[idx]; 
  
  /* detect the maximum outside of the derivative method 
     Once we detect a maximum, we have a discontinuity

     do the following:
     (1) Integrate up to the point of discontinuity 
     (2) Incorporate the discontinuity
     (3) Reinitialize the solver at the point of discontinuity 
     (4) Integrate from the point of discontinuity 

     But we need to know when the maxima or spikes occur ahead of time 
     so we can tell the integrator when it should stop integrating...

     Use 4/5 adaptive Runge-Kutta for synapse models that involve
     discontinuities
  */
      //cout  << t << "\t" << source->max << "\n";
  /* del(t-t0) */

  if(source->max) {
    flag = 1;
        //cout << "Found Max: Incorporate discontinuity\n";
    yd[idx] = -5*ut;
  } else {
    flag = 0;
    yd[idx] = -ut;
  }
  
      // cout << "flag = " << flag << "\n";
      // yd[idx] = 0;//-(ut/tauFacil) + 1;//0.01*(1-ut)*flag;
  // yd[idx+1] = zt/tauRec - ut*xt*flag;
  // yd[idx+2] = -(yt/tauDecay) + ut*xt*flag;
  // yd[idx+3] = (yt/tauDecay) - zt/tauRec;

  //cout << "AT time: " << t << "\t" << ut << "\n";
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
