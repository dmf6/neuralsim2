#include "gap.h"

Gap::Gap(Neuron *source, Neuron *target, int type, int iniVarNo)
    : Synapse(source, target, type, iniVarNo) {
    this->source = source;
    this->target = target;  
    ISyn = 0;
}

Gap::~Gap() {
    cout << "Destroying Gap\n";
};

/* compue gap junction current */
void Gap::setISyn() {
        
}

void Gap::setGmax(double g) {
    gmax = g;
}

/* a gap junction has no dynamics */
int Gap::derivative(realtype a, N_Vector *b, N_Vector *c) {
    return 0;
}

double Gap::getISyn(){
    return gmax*(source->getVoltage()-target->getVoltage());
}
