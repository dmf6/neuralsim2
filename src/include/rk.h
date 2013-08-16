#ifndef _RK_H_
#define _RK_H_

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include "neuronmodel.h"

#define MIN(X,Y) (X < Y ? X : Y)
#define MAX(X,Y) (X > Y ? X : Y)
#define RK_SUCCESS 1
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */

class RK  {
  public:
    double **k;
    N_Vector yt; /*temporary y vector for each RK function evaluation */
     realtype hh;
    int n; /* problem dimension */
    int o;  /* integrator order */
    int i, j; /*loop variables */
    
        /*constructor should take , yerr_vec, h, time, dimension N, NeuronModel */
    RK(int dim, int order);
    ~RK();
    void fadvance(realtype t, N_Vector, N_Vector, N_Vector, NeuronModel *model, double h);

    inline int dim() {
        return n;
    }
};

#endif
