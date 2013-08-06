#include <cassert>
#include "neuron.h"

int Neuron::idx_generator = 1;
int Neuron::varCount=0;

Neuron::Neuron(const string name, int iniVarNo, int mode) {
    this->name = name;
    idx= idx_generator++;    //index number of neuron
    enabled= 0;
    iVarNo= iniVarNo; /*number of state variables */
    varCount+=iVarNo;
    this->mode = mode;
        //cout << "The index of this neuron is " << idx << "\n";
    voltage = -60;
    oldVoltage = voltage;
    spiketime=-1;
    oldSpiketime=-1;
}

Neuron::~Neuron() {
        //cout << "Deallocated neuron memory" << endl;
    
    for (den_it = den.begin(); den_it != den.end(); ++den_it) {
        (*den_it)->target= NULL;
    }
}

void Neuron::init(N_Vector y, int start, double *iniVars, int n) {
  // copy the contents of from iniVars to y
  //cout << "start is " << start << " and number of variables is " << n << endl;
  
  realtype *yi;
  yi = NV_DATA_S(y);
  
  for(int i = 0; i < n; i++) {
    //cout << idx << "\t" << i << "\t" <<  iniVars[i] << "\n";
    
    yi[start+i] = iniVars[i];
  }
}

void Neuron::setTol(N_Vector tol, int start, int n) {
  realtype *ytol;
  ytol = NV_DATA_S(tol);
  
  for(int i = 0; i < n; i++) {
    ytol[start+i] = ATOL;
  }
}

int Neuron::getIVarNo() {
    return iVarNo;
}

void Neuron::detectSpikes(double t, N_Vector v){
 
  realtype *vi;
  vi = NV_DATA_S(v);
  double newVoltage = vi[0];
  if(newVoltage > oldVoltage) {
    spiking=0;
    increasing = 1;
  } 
  else if ((newVoltage < oldVoltage) && increasing) {
    increasing = 0;
    /* make sure we only get real spikes */
    if ((t-oldSpiketime) > 20) { 
      spiketime = t;
      oldSpiketime = spiketime;
      spiketimes.push_back(spiketime);
      //cout << t << '\t' << newVoltage << '\n';
      spiking = 1;
    }
  }
  oldVoltage = newVoltage;
}
