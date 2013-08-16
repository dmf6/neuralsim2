#ifndef __ABNEURON_H
#define __ABNEURON_H

#include <cstring>
#include "neuron.h"

#define IAPP 0// (1 nA)
#define CM 9 //(uF)

/* maximal conductances  (uS) */
#define GL 0.045
#define GH 0.054
#define GA 21.6
#define  GKD 1890
#define GCAT 55.2
#define GCAS 9
#define GKCA 6000
#define GNAP 2.7
#define GPROC 570

/* reversal potentials */
#define EL  -50
#define EH -20
#define EK -80
#define ENA 50
#define EPROC 0

#define GAXIAL 0.3//nS
#define CM_AXON 1.5//uF
#define EL_AXON -60
#define AXON_GL 0.0018
#define AXON_GK 52.5
#define AXON_GNA 300

#define F 0.418
#define CO 0.5
#define TAUCA 303
#define CAO 13000 //(uM)

class ABNeuron: public Neuron {
    int type;
    int pno;
    double *p;
    ofstream out;
    char outfile[80];
    int k;
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    ABNeuron(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~ABNeuron();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector);
};

#endif
