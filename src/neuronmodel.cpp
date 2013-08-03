#include "neuronmodel.h"
#include "neuron.h"
#include "synapse.h"

NeuronModel::NeuronModel(list<Neuron *> *neurs, list<Synapse *> *syns,
			 int &N, ostream &msgos) {
    iVarCnt= 0;
    Neuron *n;
    Synapse *s;

    this->neurs = neurs;
    this->syns = syns;

    forall(*neurs, niter) {
        n= *niter;
        iVarCnt += n->iVarNo;
    }
    forall(*syns, siter) {
        s= *siter;
        iVarCnt+= s->iVarNo;
    }
    N= iVarCnt;
    msgos << "We have " << N << " state variables ..." << endl;
}

void NeuronModel::derivative(realtype t, N_Vector y, N_Vector ydot, void *user_data) {  
    forall(*neurs, niter) {
        (*niter)->derivative(t , &y, &ydot, user_data);
    }
    forall((*syns), siter) {
        (*siter)->derivative(t, &y, &ydot);
    }
}

