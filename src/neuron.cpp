#include <cassert>
#include "neuron.h"

int Neuron::idx_generator = 1;

Neuron::Neuron(const string name, int iniVarNo, int mode) {
    this->name = name;
    idx= idx_generator++;    //index number of neuron
    enabled= 0;
    iVarNo= iniVarNo; /*number of state variables */
    this->mode = mode;
        //cout << "The index of this neuron is " << idx << "\n";
    voltage = -60;
    oldVoltage = voltage;
    spiketime=-1;
    oldSpiketime=-1;
    previousTime=-1;
    v1 = -60;
    v2 = -70;
    v3 =-80;
    spiking = 0;
    vi = new realtype[iVarNo];
    
}

Neuron::~Neuron() {
        //cout << "Deallocated neuron memory" << endl;  
    // for (den_it = den.begin(); den_it != den.end(); ++den_it) {
    //     (*den_it)->target= NULL;
    // }
    delete [] vi;
    
}

void Neuron::init(N_Vector y, double *iniVars) {
    realtype *yi;
    yi = NV_DATA_S(y);
    
    for(int i = 0; i < iVarNo; i++) {
        yi[idx+i] = iniVars[i];
    }
}

void Neuron::setTol(N_Vector tol) {
  realtype *ytol;
  ytol = NV_DATA_S(tol);
  
  for(int i = 0; i < iVarNo; i++) {
    ytol[idx+i] = ATOL;
  }
}

int Neuron::getIVarNo() {
    return iVarNo;
}

int Neuron::detectSpikes(double t, N_Vector v){ 

    vi = NV_DATA_S(v);
    
    if(voltage > oldVoltage) {
        spiking=0;
        increasing = 1;
    } 
    else if ((voltage < oldVoltage) && increasing) {
        increasing = 0;
        
            /* make sure we only get real spikes */
        if ((t-oldSpiketime) > 20) { 
            spiketime = t;
            cout << "SPIKE\n";
            oldSpiketime = spiketime;
            spiketimes.push_back(oldSpiketime);
            spiking = 1;
        }
    }
    else {
        spiking = 0;
    }
    
    oldVoltage = voltage;
    previousTime = t;
    return 0;
    
}

int Neuron::detectMaximum(double t){
    if (v2 > v1 && v2 > v3) {
        cout << "Time of max: " << t << "\n";
        max = 1;
    }
    else {
        max = 0;
    }
    
    v3 = v2;
    v2 = v1;
    v1 = voltage;
    return max;
    
        //cout << v3 << "\t" << v2 << "\t" << v1 << "\n";
    
    // double newVoltage = voltage;
    // if(newVoltage > oldVoltage) {
    //     max=0;
    //     increasing = 1;
    // }
    // else if ((newVoltage < oldVoltage) && increasing) {
    //     increasing = 0;
    //         max = 1;
    // }
    //oldVoltage = newVoltage;
    //previousTime = t;
    //return max;
}

void Neuron::clampVoltage(realtype  t, N_Vector y, N_Vector ydot, ostream& out) {
    realtype Iext;
    Iext = 0;    
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        Iext += (*iapp_it)->getIapp();
    }
    
    voltage = Iext;
}

void Neuron::writeStateVector(N_Vector y) {
    realtype *yi;
    yi = NV_DATA_S(y);
    cout << name << " State Vector\n";
    
    for(int i = idx; i < idx + iVarNo; i++) {
        cout << i << "\t" << yi[i] << "\n";
    }
}

