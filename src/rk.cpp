#include "rk.h"

RK::RK(int dim, int order){
    n= dim;
    o = order;
    
    yt = N_VNew_Serial(dim);
    k = new double*[o];

     for(i = 0; i < o; i++) {
         k[i] =  new double[n];
     }
}

RK::~RK() {
    for(int i = 0; i < o; i++) {
        delete [] k[i];
    }
    delete [] k;
    N_VDestroy_Serial(yt);
}

void RK::fadvance(realtype t, N_Vector y, N_Vector yout, N_Vector ydot, NeuronModel *model, realtype h) {
    realtype *yi, *yd, *yo, *yt_data;
    yi = NV_DATA_S(y);
    yo = NV_DATA_S(yout);
    yd = NV_DATA_S(ydot);

    model->derivative(t, y, ydot, NULL); 
    for (i=0; i < n; i++) {
        k[0][i] = yd[i]; //k1= f(yn, tn)
    }
    
    for (i=0; i < n; i++) {
        yi[i] = yi[i] + h*k[0][i];
    }
    
    model->derivative(t+h, y, ydot, NULL);  
    for (i=0; i < n; i++) {
        k[1][i] = yd[i]; //k2 = h* f(tn+h, yn + h*f(tn, yn))
    
        
    }
    /* yn+1 = yn + (h/2)*[f(tn, yn) + f(tn+h, yn + h*f(tn, yn))] */
    for (i=0; i < n; i++) {
        yo[i] = yi[i] + (h/2)*(k[0][i] + k[1][i]);
    }
}