#ifndef _CVODESOLVER_H_INCLUDED_
#define _CVODESOLVER_H_INCLUDED_

/* cvode libraries */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "neuronmodel.h"

#define T0 0.095367432            /* start time */
#define RTOL  RCONST(1.0e-8)   /* scalar relative tolerance            */

/* definition of the function to be passed to CVodeSolver */
typedef int (*RhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data);

class NeuronModel;

class CVodeSolver {
  private:
    N_Vector yi, abstol;
    realtype reltol;
    void *cvode_mem;
    RhsFn rhsFn;
    
  public:
    
/*pass in y, ydot, yout, and abstol vectors to constructor */
        /* Also pass in the NeuronModel, which is our user_data */
        /* in the constructor we should do the following:
         *  (1) Set vector of intial values with problem size
         *  (2) Create CVode object
         *  (3) Initialize CVode solver
         *  (4) Specify Integration tolerances
         *  (5) Set user data
         *  (6) Attach Linear solver module (CVDense)
         */
    CVodeSolver(int, N_Vector, N_Vector, NeuronModel *, RhsFn);
    
        /* Deallocate memory for solution vector */
    ~CVodeSolver();
    
        /* function to be solved in time */
    int f(realtype t, N_Vector y, N_Vector dy, void *user_data);
    
        /* advance solution in time */
    int fadvance(realtype tout, realtype t);

    int setStopTime(realtype t);
};
#endif /* _CVODESOLVER_H_INCLUDED */
