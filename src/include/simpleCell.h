#ifndef __SIMPLECELL_H
#define __SIMPLECELL_H

#include <cstring>
#include "neuron.h"

#define CM 1000

#define GL p[0]
#define GH p[1]
#define GCA p[2]
#define EL  p[3]
#define EH p[4]
#define ECA p[5]

/* h-current activation pars */
#define VH_H p[6]
#define K_H p[7]

#define TAU1H p[8]
#define TAU2H p[9]
#define TAU_VH_H p[10]
#define TAU_K_H p[11]

/* Ca activation pars */
#define VH_CA_M p[12]
#define K_CA_M p[13]

#define TAU1_CA_M p[14]
#define TAU2_CA_M p[15]
#define TAU_VH_CA_M p[16]
#define TAU_K_CA_M p[17]
 
#define TAU1_CA_H p[18]
#define TAU2_CA_H p[19]
#define TAU_VH_CA_H p[20]
#define TAU_K_CA_H p[21]

#define VH_CA_H p[22]
#define K_CA_H p[23]


class SimpleCell: public Neuron {
    int type;
    int m_dValue;
    int pno;
    double *p;
    ofstream out;
    char outfile[80];
    double iapp;
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    SimpleCell(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~SimpleCell();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, N_Vector);
};

#endif
