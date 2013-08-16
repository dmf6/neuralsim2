#include "simplecell.h"

SimpleCell::SimpleCell(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0) 
    : Neuron(name, iniVarNo, mode), pno(inpno) {

    if (pno > 0) {
        p = new double[pno];
        set_p(inp);
    }
    iapp = 0;
    voltage = -60.0;
    strcpy (outfile, name.c_str());
    strcat (outfile,".asc");
    out.open(outfile);
}


SimpleCell::~SimpleCell() {
    delete [] p;
    out.close();
}

void SimpleCell::set_p(double *inp)
{
    for (int i= 0; i < pno; i++) {
        p[i]= inp[i];
    }
}

int SimpleCell::derivative(realtype t, N_Vector *y, N_Vector *ydot, void *user_data) {
    realtype h_minf, h_mtau, ca_minf, ca_mtau, ca_hinf, ca_htau;
    realtype v, h_m, ca_m, ca_h;
    realtype Ih, Il, ICa, Iext;
    realtype *yi, *yd;

    N_Vector yn = *y;
    N_Vector ydn = *ydot;
    
    yi = NV_DATA_S(yn);
    yd = NV_DATA_S(ydn);

    v = yi[idx]; h_m = yi[idx+1]; ca_m = yi[idx+2]; ca_h = yi[idx+3];
     
    Iext = 0;    
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        Iext += (*iapp_it)->getIapp();
    }
    if (mode == VCLAMP) {
        v = Iext;
        voltage = v;
        
    }
    
        /* CHECK SIGN OF SIGMOID FUNCTION */
    h_minf    = 1/(1+exp((v-VH_H)/K_H));
    h_mtau   = TAU2H+(TAU1H - TAU2H)*h_minf;
    
    ca_minf = 1/(1+exp(-(v-VH_CA_M)/K_CA_M));  
    ca_mtau   = TAU2_CA_M+(TAU1_CA_M - TAU2_CA_M)*(1-ca_minf);
    
    ca_hinf = 1/(1+exp((v-VH_CA_H)/K_CA_H));
    ca_htau   = TAU2_CA_H+(TAU1_CA_H - TAU2_CA_H)*(1-ca_hinf);
    
    Il = GL*(v-EL);
    Ih=GH*h_m*(v-EH);
    ICa = GCA*pow(ca_m, 3)*ca_h*(v-ECA);

    if (mode == VCLAMP) {
        yd[0] = 0;
        
    }
    else {
        yd[idx] = (Iext-(Il + Ih + ICa))/CM;
    }
    
    yd[idx+1] = (h_minf - h_m)/h_mtau;
    yd[idx+2] = (ca_minf - ca_m)/ca_mtau;
    yd[idx+3] = (ca_hinf - ca_h)/ca_htau;
    
    return (0);
}

void SimpleCell::currents(realtype  t, N_Vector y, ostream &out) {
    realtype *yi;
    realtype dvdt, v, h_m, ca_m, ca_h, Il, Ih, ICa, Icap;
    realtype I_ion;
    
    yi = NV_DATA_S(y);
    v = voltage; h_m = yi[1]; ca_m = yi[2]; ca_h = yi[3];
    
    Il = GL*(v-EL);
    Ih=GH*h_m*(v-EH);
    ICa = GCA*pow(ca_m, 3)*ca_h*(v-ECA);

    Icap = CM*((v-oldVoltage)/0.05);
    oldVoltage = v;

    // Iapp = sum of I_ion + Icap 
    I_ion = (Il + Ih + ICa);
    
    iapp = (I_ion+Icap);
    out << t  << '\t' << iapp << '\t' << voltage <<'\n';
}
