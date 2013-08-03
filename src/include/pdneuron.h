#ifndef __PDNEURON_H
#define __PDNEURON_H

#include <cstring>
#include "neuron.h"

#define IAPP 0// (1 nA)
#define CM 12 //(uF)

/* maximal conductances  (uS) */
#define GL 0.105
#define GH 0.219
#define GA 390.42
#define  GKD 1576.8
#define GCAT 10
#define GCAS 54
#define GKCA 251.85
#define GNAP 4.38
#define GPROC 0

/* reversal potentials */
#define EL  -55
#define EH -20
#define EK -80
#define ENA 50
#define EPROC 570

#define GAXIAL 1.05//nS
#define CM_AXON 6.0//uF
#define EL_AXON -55
#define AXON_GL 0.00018
#define AXON_GK 150
#define AXON_GNA 1110

#define F 0.515
#define CO 0.5
#define TAUCA 300.0
#define CAO 13000 //(uM)

class PDNeuron: public Neuron {
    int type;
    int m_dValue;
    int pno;
    double *p;
    ofstream out;
    char outfile[80];
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    PDNeuron(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~PDNeuron();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, N_Vector, ostream&);
};

#endif
