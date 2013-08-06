#include <cassert>
#include "synapse.h"

int Synapse::idx_generator = 1;
int Synapse::varCount=0;

Synapse::Synapse(Neuron *source, Neuron *target, int type, int iniVarNo) {
    assert(target != NULL);
    source= source;
    target= target;
    target->den.push_back(this);
    gmax = 0;
    ISyn = 0;
    iVarNo= iniVarNo; /*number of state variables */
    varCount+=iVarNo;
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
