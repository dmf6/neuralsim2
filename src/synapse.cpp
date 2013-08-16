#include <cassert>
#include "synapse.h"

int Synapse::idx_generator = 1;

Synapse::Synapse(Neuron *source, Neuron *target, int type, int iniVarNo) {
        //assert(target != NULL);
    source= source;
    target= target;
    if (target != NULL) {
        target->den.push_back(this);
    }
    iVarNo= iniVarNo; /*number of state variables */
    gmax = 0;
    ISyn = 0;
}

Synapse::~Synapse() {
    cout << "Destroying Synapse\n";
}

void Synapse::setIdx(int inidx) {
    assert(!enabled);
    idx= inidx;
    enabled= 1;
}

int Synapse::getIdx() {
    return idx;
}

void Synapse::init(N_Vector y, double *iniVars) {
  // copy the contents of from iniVars to y
    cout << "start is " << idx << " and number of variables is " << iVarNo << endl;
  
    realtype *yi;
    yi = NV_DATA_S(y);
        //cout << idx << "\n";
    for(int i = 0; i < iVarNo; i++) {
            //cout << idx << "\t" << i << "\t" <<  iniVars[i] << "\n";
        yi[idx+i] = iniVars[i];
    }
}

void Synapse::setTol(N_Vector tol) {
  realtype *ytol;
  ytol = NV_DATA_S(tol);
  
  for(int i = 0; i < iVarNo; i++) {
    ytol[idx+i] = ATOL;
  }
}

