#include "simplecell.h"

SimpleCell::SimpleCell(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0) 
    : Neuron(name, iniVarNo, mode), pno(inpno) {

    if (pno > 0) {
        p = new double[pno];
        set_p(inp);
    }
    iapp = 0;
    voltage = -60.0;
    flag = 0;
    oldV=voltage;
    
}


SimpleCell::~SimpleCell() {
    delete [] p;
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

    realtype na_minf, na_hinf, na_mtau, na_htau, k_minf, k_mtau;
    realtype v_axon, na_m, na_h, k_m;
    realtype Ina_axon, Ik_axon, Il_axon;
    realtype iaxial_s, iaxial_a;
        
    N_Vector yn = *y;
    N_Vector ydn = *ydot;
    
    yi = NV_DATA_S(yn);
    yd = NV_DATA_S(ydn);

    v = yi[0]; h_m = yi[1]; ca_m = yi[2]; ca_h = yi[3];
   
    // if (v > V_THR && flag == 0) {
    //     flag = 1;
    //     v -= 3200.0;
    // }
    // else if (v > V_THR && flag == 1) {
    //     flag = 0;
    //     v+=3600.0;
    // }
     
    Iext = 0;    
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        Iext += (*iapp_it)->getIapp();
    }

    voltage = v;
 
    
        /* CHECK SIGN OF SIGMOID FUNCTION */
    h_minf    = 1/(1+exp((v-VH_H)/K_H));
    h_mtau   = TAU2H+(TAU1H - TAU2H)*h_minf;
    
    ca_minf = 1/(1+exp(-(v-VH_CA_M)/K_CA_M));  
    ca_mtau   = TAU2_CA_M+(TAU1_CA_M - TAU2_CA_M)*(1-ca_minf);
    
    ca_hinf = 1/(1+exp((v-VH_CA_H)/K_CA_H));
    ca_htau   = TAU2_CA_H+(TAU1_CA_H - TAU2_CA_H)*(ca_hinf);
    
    Il = GL*(v-EL);
    Ih=GH*h_m*(v-EH);
    ICa = GCA*pow(ca_m, 3)*ca_h*(v-ECA);

        /* see Soto-Trevino et al., 2005 */
        /* axon compartment Na+/K+ currents */
    v_axon = yi[4]; na_m = yi[5]; na_h = yi[6]; k_m=yi[7];
    
    na_minf = 1/(1+exp(-(v_axon + 24.7)/5.29));
    na_mtau = 1.32 - 1.26/(1+exp(-(v_axon+120)/25));

    na_hinf = 1/(1+exp((v_axon+48.9)/5.18));
    na_htau = (0.67/(1+exp(-(v_axon+62.9)/10)))*(1.5+1/(1+exp((v_axon + 34.9)/3.6)));

    k_minf = 1/(1+exp(-(v_axon+14.2)/11.8));
    k_mtau = 7.2-6.4/(1+exp(-(v_axon+28.3)/19.2));
    
    Ina_axon = GNA_AXON*pow(na_m, 3)*na_h*(v_axon-ENA);
    Ik_axon = GK_AXON*pow(k_m, 4)*(v_axon-EK);
    Il_axon = GL_AXON*(v_axon-EL_AXON);
    
        /* axial current in the soma and axon compartment - axial
         * currents are symmetrical */
    iaxial_s = GAXIAL1*(v-v_axon);
    iaxial_a = GAXIAL2*(v_axon-v);
 
        /* just one compartment in voltage clamp ...for now */
    if (mode == VCLAMP) {
        yd[0] = -(Il + Ih + ICa + iaxial_s)/CM;
        yd[4] = -(Ina_axon + Ik_axon + Il_axon + iaxial_a)/CM_AXON;
    }
    else {
        yd[0] = (Iext-(Il + Ih + ICa + iaxial_s))/CM;
        yd[4] = (-(Ina_axon + Ik_axon + Il_axon + iaxial_a))/CM_AXON;
    }
    
    
    yd[1] = (h_minf - h_m)/h_mtau;
    yd[2] = (ca_minf - ca_m)/ca_mtau;
    yd[3] = (ca_hinf - ca_h)/ca_htau;
        /* axon state variables */
    yd[5] = (na_minf-na_m)/na_mtau;
    yd[6] = (na_hinf-na_h)/na_htau;
    yd[7] = (k_minf - k_m)/k_mtau;
    
    return (0);
}

void SimpleCell::currents(realtype  t, N_Vector y, N_Vector ydot, ostream& out) {
    realtype *yi, *yd;
    realtype dvdt, v, h_m, ca_m, ca_h, Il, Ih, ICa, Icap;
    realtype I_ion;
    realtype v_axon, iaxial_s, iaxial_a;
    
    yi = NV_DATA_S(y);
    yd = NV_DATA_S(ydot);
    v = yi[0]; h_m = yi[1]; ca_m = yi[2]; ca_h = yi[3];
    v_axon = yi[4];
    dvdt = yd[0];
    
    Il = GL*(v-EL);
    Ih=GH*h_m*(v-EH);
    ICa = GCA*pow(ca_m, 3)*ca_h*(v-ECA);
    Icap = CM*dvdt;
        //out << t << '\t' << (olddvdt -dvdt)/0.05 << '\n';
    Icap = CM*(-(oldV-v)/0.05);
    oldV = v;
    
    iaxial_s = GAXIAL1*(v - v_axon);
    iaxial_a = GAXIAL2*(v_axon-v);
    // Iapp = sum of I_ion + Icap 
    I_ion = (Il + Ih +Il+ ICa);
    
    iapp = (I_ion+Icap);
    out << t  << '\t' << iapp << '\t' << voltage <<'\n';
}
