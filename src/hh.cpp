#include "hh.h"

HH::HH(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0) 
    : Neuron(name, iniVarNo, mode), pno(inpno) {

    if (pno > 0) {
        p = new double[pno];
        set_p(inp);
    }
    iapp = 0;
    voltage = -70.0;
}


HH::~HH() {
    delete [] p;
}

void HH::set_p(double *inp)
{
    for (int i= 0; i < pno; i++) {
        p[i]= inp[i];
    }
}

int HH::derivative(realtype t, N_Vector *y, N_Vector *ydot, void *user_data) {
    realtype am, bm, ah, bh, an, bn;
    realtype minf, hinf, ninf, taum, tauh, taun;
    realtype v, m, h, n;
    realtype Il, INa, Ik, iapp;
    realtype *yi, *yd;
        
    N_Vector yn = *y;
    N_Vector ydn = *ydot;
    
    yi = NV_DATA_S(yn);
    yd = NV_DATA_S(ydn);

    v = yi[0]; m = yi[1]; h = yi[2]; n = yi[3];

    // am =  -0.1 * (25.0-v) / (exp(-0.1*(25.0+v)) - 1.0);
    // bm =  4.0 * exp(-v/18.0);

    // ah =  0.07 * exp(-v/20.0);
    // bh = 1.0 / (exp(-(30.0-v)/10.0) + 1.0);
    
    // an = 0.01 * (10.0+v) / (exp(0.1*(10.0-v)) - 1.0);
    // bn = 0.125 * exp(-v/80.0);
    minf = 1/(1+exp(-(v+40)/9));
    ninf =  1/(1+exp(-(v+53)/16));
    hinf =  1/(1+exp((v+62)/10));
    taum = 0.3;
    tauh= 1 + (11/(1+exp((v+62)/10)));
    taun= 1 + (6/(1+exp((v+53)/16)));
    
    iapp = 0.0;
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        iapp += (*iapp_it)->getIapp();
    }
    
    Il = GL*(v-EL);
    Ik = GK*pow(n, 4)*(v -EK);
    INa = GNA*pow(m, 3)*h*(v-ENA);

        //cout << v << "\t" <<  m << "\t" << h << "\t" << n << "\n";
    
        /* just one compartment in voltage clamp ...for now */
    if (mode == VCLAMP) {
        yd[0] = (-iapp-(Il+Ik+INa))/CM;
    }
    else {
         yd[0] = (iapp - (Il + Ik + INa))/CM;
    }
               
    // yd[1] = am*(1-m)-bm*m;
    // yd[2] = ah*(1-h)-bh*h;
    // yd[3] = an*(1-n)-bn*n;
    yd[1] = (minf-m)/taum;
    yd[2] = (hinf-h)/tauh;
    yd[3] = (ninf-n)/taun;
       
    return (0);
}

void HH::currents(realtype  t, N_Vector y, N_Vector ydot, ostream& out) {
    realtype *yi, *yd;
    realtype v, m, h, n;
    realtype Il, INa, Ik, iapp, Icap, dvdt;
    realtype I_ion;
    
    yi = NV_DATA_S(y);
    yd = NV_DATA_S(ydot);
    v = yi[0]; m = yi[1]; h = yi[2]; n = yi[3];
    dvdt = yd[0];

    Il = GL*(v-EL);
    Ik = GK*pow(n, 4)*(v - EK);
        //INa = GNA*pow(m, 3)*h*(v-ENA);
    Il = GL*(v-EL);
    Icap = CM*dvdt;
 
    // Iapp = sum of I_ion + Icap 
    I_ion = (Il + INa + Ik);
    iapp = (I_ion + 0.001*Icap);
    out << t  << '\t' << iapp << '\t' << v << '\n';
}
