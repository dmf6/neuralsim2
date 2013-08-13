#include "gap.h"

Gap::Gap(Neuron *source, Neuron *target, int type, int iniVarNo)
    : Synapse(source, target, type, iniVarNo) {
    ISyn = 0;
        cout << target->voltage << "\n";
        this->target = target;
        this->source = source;
        
}

Gap::~Gap() {
    cout << "Destroying Gap\n";
};

/* compue gap junction current */
void Gap::setISyn() { 
}

void Gap::setGmax(double g) {
    gmax = g;
    cout << target->voltage << "\n";
}

/* a gap junction has no dynamics */
int Gap::derivative(realtype a, N_Vector *b, N_Vector *c) {
    return 0;
}

double Gap::getISyn(){

        //cout <<  source->name << ": " << source->voltage << "\t" << target->name << ": " << target->voltage << "\n";
        /* has to be target - source voltag because this synapse is
         * associated with the target neuron */
    return gmax*(target->voltage-source->voltage);
}
