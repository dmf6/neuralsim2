#include "CVode.h"


CVode::CVode(int n, N_Vector yi, NeuronModel *model) {
    this->yi = yi;
    this->abstol = abstol;
    
    reltol = RTOL;
    
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeSetUserData(cvode_mem, static_cast<void *>(model));

    CVodeInit(cvode_mem, this->f, T0, yi);

    CVodeSVtolerances(cvode_mem, reltol, abstol);  
    CVDense(cvode_mem, n);
}

int CVode::f(realtype t, N_Vector y, N_Vector dy, void *user_data) {
        /* for multiple cells cast user_data to NeuronModel */
    NeuronModel *model = static_cast<NeuronModel *>(user_data);
    model->derivative(t, y, dy, NULL);
    return 0;   
} 

/* make the call to perform the integration of the IVP */
void CVode::fadvance(realtype tout, realtype t) {
    int flag;
    flag = CVode(cvode_mem, tout, yi, &t, CV_NORMAL);
}
