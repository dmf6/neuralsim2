#ifndef __SIMPLECELL_H
#define __SIMPLECELL_H

#include <cstring>
#include "neuron.h"

#define VCLAMP 1

#define CM 1

#define GL p[0]
#define GH p[1]
#define GCA p[2]
#define EL  -60
#define EH -20
#define ECA 120
#define ENA 50
#define EK -80
/* h-current activation pars */
#define VH_H p[3]
#define K_H p[4]

#define TAU1H p[5]
#define TAU2H p[6]

/* Ca activation pars */
#define VH_CA_M p[7]
#define K_CA_M p[8]

#define TAU1_CA_M p[9]
#define TAU2_CA_M p[10]

#define VH_CA_H p[11]
#define K_CA_H p[12]

#define TAU1_CA_H p[13]
#define TAU2_CA_H p[14]

using namespace std;

class SimpleCell: public Neuron {
    int type;
    int pno;
    double *p;
    realtype iapp;
    int flag;
    ofstream out;
    char outfile[80];
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    SimpleCell(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~SimpleCell();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, ostream&);
};

#endif
