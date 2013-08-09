#include "CVodeSolver.h"

/* User has to pass in a function pointer defining f(t, y) to the
 * constructor of the CVodeSolver class */
CVodeSolver::CVodeSolver(int n, N_Vector yi, N_Vector abstol, NeuronModel *model, RhsFn rhsfn ) {
    this->yi = yi;
    this->abstol = abstol;
    rhsFn = rhsfn;
    reltol = RTOL;
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeSetUserData(cvode_mem, static_cast<void *>(model));
    CVodeInit(cvode_mem, rhsFn, T0, yi);
    CVodeSVtolerances(cvode_mem, reltol, abstol);
    CVDense(cvode_mem, n);
}

CVodeSolver::~CVodeSolver() {
        // cout << "Deallocated CVode memory" << endl;
    
   N_VDestroy_Serial(yi);
   N_VDestroy_Serial(abstol);
   CVodeFree(&cvode_mem);
}

/* make the call to perform the integration of the IVP */
/* 
   flag = CVode(cvode mem, tout, yout, &tret, itask);
   itask can be either CV_NORMAL or CV_ONE_STEP 
*/
int CVodeSolver::fadvance(realtype tout, realtype t) {
    return CVode(cvode_mem, tout, yi, &t, CV_NORMAL);
}

int CVodeSolver::setStopTime(realtype t) {
  return CVodeSetStopTime(cvode_mem, t);
}

int CVodeSolver::rootInit(int n, gFn g) {
  return CVodeRootInit(cvode_mem, n, g);
}

int CVodeSolver::getRootInfo(int *rootsfound) {
  return CVodeGetRootInfo(cvode_mem, rootsfound);
}
