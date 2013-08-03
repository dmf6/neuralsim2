#include "simplecell3.h"

SimpleCell3::SimpleCell3(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0) 
    : Neuron(name, iniVarNo, mode), pno(inpno) {

    if (pno > 0) {
        p = new double[pno];
        set_p(inp);
    }
    iapp = 0;
}


SimpleCell3::~SimpleCell3() {
    delete [] p;
}

void SimpleCell3::set_p(double *inp)
{
    for (int i= 0; i < pno; i++) {
        p[i]= inp[i];
    }
}
/* 3 ODEs */
int SimpleCell3::derivative(realtype t, N_Vector *y, N_Vector *ydot, void *user_data) {
    realtype m_minf, m_mtau, mi_minf, mi_mtau;
    realtype v, m_m, mi_m;
    realtype Im, Il, Imi, iapp;
    realtype *yi, *yd;
    
    N_Vector yn = *y;
    N_Vector ydn = *ydot;
    
    yi = NV_DATA_S(yn);
    yd = NV_DATA_S(ydn);

    v = yi[0]; m_m = yi[1]; mi_m = yi[2];
    
    iapp = 0.0;
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        iapp += (*iapp_it)->getIapp();
    }
    
    m_minf    = 1/(1+exp(-(v-VH_M)/K_M));
    m_mtau   = TAU1M+TAU2M/(1+exp(-(v-TAU_VH_M)/TAU_K_M));

    mi_minf = 1/(1+exp(-(v-VH_MI_M)/K_MI_M));
    mi_mtau   = TAU1_MI+TAU2_MI/(1+exp(-(v-VH_MI_M)/K_MI_M));
        
    Il = GL*(v-EL);
    Im=GM*m_m*(v-EM);
    Imi = GMI*mi_m*(v-EMI);
        
    yd[0] = (-iapp-(Il+Im+Imi))/CM;
    yd[1] = (m_minf - m_m)/m_mtau;
    yd[2] = (mi_minf - mi_m)/mi_mtau;
    
    return (0);
}

void SimpleCell3::currents(realtype  t, N_Vector y, N_Vector ydot, ostream& out) {
    realtype *yi, *yd;
    realtype dvdt, v, m_m, mi_m, Il, Im, Imi, Icap;
    realtype I_ion;

    yi = NV_DATA_S(y);
    yd = NV_DATA_S(ydot);
    v = yi[0]; m_m = yi[1]; mi_m = yi[2];
    dvdt = yd[0];
    
    Il = GL*(v-EL);
    Im=GM*m_m*(v-EM);
    Imi = GMI*mi_m*(v-EMI);
    Icap = CM*dvdt;
    // Iapp = sum of I_ion + Icap 
    I_ion = (Il + Im + Imi);
    iapp = (I_ion + 0.001*Icap);
    out << t  << '\t' << iapp << '\t' << v <<  '\n';
}
