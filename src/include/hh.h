#ifndef __HH_H
#define __HH_H

#include <cstring>
#include "neuron.h"

#define VCLAMP 1

#define CM 150.0

#define GL 0.0003
#define GNA 0.12
#define GK 0.036
#define EL  -50.6
#define ENA 55.0
#define EK -75.0

using namespace std;

class HH: public Neuron {
    int type;
    int pno;
    double *p;
    double iapp;
    double olddvdt;
    

  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    HH(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~HH();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, N_Vector, ostream&);
};

#endif
