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
    spiking = 0;
    this->mode = mode;
    spikeCount = 0; /* number of spikes */
    period = INF;
        //cout << "The index of this neuron is " << idx << "\n";
    voltage = -60;
    
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

/* neuron class can keep track of its own spike times via a vector
 * member variable
 *
 * If we just want to measure period then we can stop the numerical
 * integration at spikeCount > 2 and output the period
 *
 * burst duration is time spent above the voltage threshold
 */
void Neuron::detectSpike(double t, double v) {
    if (v > SPIKE_VTHR) {
        if (spiking == 0) {
                /* save the time of the crossing in vector */
            spikeTimes.push_back(t);
            spikeCount +=1;
            spiking = 1;
                /* save time of start of active phase */
            t_on = t;
        }
        else {
                /* active phase has ended so save the time as t_off */
            t_off = t;
            spiking = 0;
        }
    }
}

double Neuron::getPeriod() {
    if (spikeCount > 1) {
        period = spikeTimes[1] - spikeTimes[0];
    }
        /* Neuron is constructed with an infinite period so we don't
         * need to set it */
    return period;
}

double Neuron::getBurstDur() {
    return (t_off - t_on);
}

int Neuron::getIVarNo() {
    return iVarNo;
}
