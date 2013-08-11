#include "abneuron.h"

ABNeuron::ABNeuron(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0)
    : Neuron(name, iniVarNo, mode), pno(inpno) {
    if (pno > 0) {
         p = new double[pno];
         set_p(inp);       
    }
  
    strcpy (outfile, name.c_str());
    strcat (outfile,".asc");
    out.open(outfile);
        //cout << varCount << "\t" << idx << "\t" << iVarNo << "\t" << k<< endl;
}

ABNeuron::~ABNeuron() {
    delete [] p;
    out.close();
}

void ABNeuron::set_p(double *inp)
{
    for (int i= 0; i < pno; i++) {
        p[i]= inp[i];
    }
}
int ABNeuron::derivative(realtype t, N_Vector *y, N_Vector *ydot, void *user_data) {
    realtype h_minf, h_mtau, a_minf, a_mtau, a_hinf, a_htau;
    realtype kd_minf, kd_mtau, cat_minf, cat_hinf, cat_mtau, cat_htau, cas_minf, cas_mtau;
    realtype kca_minf, kca_mtau, nap_minf, nap_mtau, nap_hinf, nap_htau, k_minf, k_mtau;
    realtype na_minf, na_mtau, na_hinf, na_htau, proc_minf;
    realtype v, h_m, a_m, a_h, kd_m, cat_m, cat_h, cas_m, cai, kca_m, nap_m, nap_h, v_a, k_m, na_m, na_h, proc_m;
    realtype Ih, Il, Ia, Ikd, ICaT, ICaS, IKCa, INaP, IProc, Isyn, Ina_axon, Ik_axon, Il_axon, icoup0, icoup1, Iapp;
    realtype eca;
    realtype *yi, *yd;
    
        /* sum current from all incoming synapses */
    // Isyn= 0.0;
    // for (den_it = den.begin(); den_it != den.end(); ++den_it) {
    //     Isyn += (*den_it)->getISyn();
    // }
    //cout << "Synaptic current is " <<Isyn << "\n";
    
    Iapp = 0.0;
    for ( iapp_it = electrodes.begin();  iapp_it != electrodes.end(); ++iapp_it) {
        Iapp += (*iapp_it)->getIapp();
    }

    N_Vector yn = *y;
    N_Vector ydn = *ydot;

    yi = NV_DATA_S(yn);
    yd = NV_DATA_S(ydn);

    //cout << "voltage index is " << v << "\n";
    
    v = yi[idx]; h_m=yi[idx+1]; a_m = yi[idx+2]; a_h = yi[idx+3]; kd_m = yi[idx+4]; cat_m = yi[idx+5]; cat_h = yi[idx+6];
    cas_m = yi[idx+7]; cai = yi[idx+8]; kca_m = yi[idx+9]; nap_m = yi[idx+10]; nap_h = yi[idx+11]; v_a = yi[idx+12]; k_m = yi[idx+13];
    na_m=yi[idx+14]; na_h=yi[idx+15]; proc_m=yi[idx+16];
    
        /*compute calcium reversal potential using Nernst equation  (RT/zF)*(ln(cao/cai))*/
        //2.303*(RT/zF)*(cao/cai)*1000
        // T = 283, R = 8.135, F = 96840
    eca = 2.303*((8.314*310)/(2*(96400)))*log(CAO/cai)*1000;

         /* calculate voltage and time dependent rates */
    h_minf    = 1/(1+exp((v+70.0)/6.0));
    h_mtau   = 272+1499/(1+exp((v+42.2)/8.73));

    a_minf = 1/(1+exp(-(v+27.0)/8.7));
    a_mtau = 11.6-10.4/(1+exp(-(v+32.9)/15.2));
    a_hinf = 1/(1+exp((v+56.9)/4.9));
    a_htau = 38.6-29.2/(1+exp(-(v+38.9)/26.5));

    kd_minf = 1/(1+exp(-(v+14.2)/11.8));
    kd_mtau = 7.2-6.4/(1+exp(-(v+28.3)/19.2));

    cat_minf = 1/(1+exp(-(v+25)/7.2));
    cat_mtau =  55-49.5/(1+exp((v+58)/17));
    cat_hinf = 1/(1+exp((v+36)/7));
    cat_htau = 87.5-75/(1+exp(-(v+50)/16.9));

    cas_minf = 1/(1+exp(-(v+22)/8.5));
    cas_mtau =  16-13.1/(1+exp(-(v+25.1)/26.4));

    kca_minf = (cai/(cai+30))*(1/(1+exp(-(v+51)/4)));
    kca_mtau = 90.3-75.09/(1+exp(-(v+46)/22.7));
    
    nap_minf = 1/(1+exp(-(v+26.8)/8.2));
    nap_mtau =  19.8-10.7/(1+exp(-(v+26.5)/8.6));
    nap_hinf = 1/(1+exp((v+48.5)/4.8));
    nap_htau = 666-379/(1+exp(-(v+33.6)/11.7));

    proc_minf = 1/(1+exp(-(v+12)/3.05));
    
        /* axon rate equations */
    k_minf = 1/(1+exp(-(v_a+14.2)/11.8));
    k_mtau = 7.2-6.4/(1+exp(-(v_a+28.3)/19.2));

        /* for VCLAMP mode set GNA to 0 since we don't want spikes */
    na_minf = 1/(1+exp(-(v_a+24.7)/5.29));
    na_mtau = 1.32-1.26/(1+exp(-(v_a+120)/25));
    na_hinf = 1/(1+exp((v_a+48.9)/5.18));
    na_htau = (0.67/1+exp(-(v_a+62.9)/10))*(1.5+1/(1+exp((v_a+34.9)/3.6)));
        
        /* somatic currents */
    Il = GL*(v-EL);
    Ih=GH*h_m*(v-EH);
    Ia =GA*pow(a_m, 3)*a_h*(v-EK);
    Ikd = GKD*pow(kd_m, 4)*(v-EK);
    ICaT = GCAT*pow(cat_m, 3)*cat_h*(v-eca);
    ICaS = GCAS*pow(cas_m, 3)*(v-eca);
    IKCa =GKCA*pow(kca_m, 4)*(v-EK);
    INaP = GNAP*pow(nap_m, 3)*nap_h*(v-ENA);
    IProc = GPROC*proc_m*(v-EPROC);
    

        /*axon currents */
    Ik_axon=AXON_GK*pow(k_m, 4)*(v_a-EK);
    Il_axon=AXON_GL*(v_a-EL_AXON);
    Ina_axon=AXON_GNA*pow(na_m, 3)*na_h*(v_a-ENA);
                                    
    icoup0 = GAXIAL*(v-v_a);
    icoup1 = GAXIAL*(v_a-v);
    
    if(mode == 0) {
        if (electrodes[0]->isEnabled()) {
            yd[idx] =  (Iapp-(Ih + Il + Ikd + INaP+Ia+ICaT+ICaS+IKCa+icoup0))/CM;
        }
        else {
             yd[idx] =  (-(Ih + Il + Ikd + INaP+Ia+ICaT+ICaS+IKCa+icoup0))/CM;
        }
    }
    else {
        yd[idx+0] = 0;
    }

    voltage = v;
    //cout << "AB Voltage = " << voltage << '\n';
        //cout << this->getVoltage() << "\n";
    yd[idx+1] = (h_minf - h_m)/h_mtau;
    yd[idx+2] = (a_minf - a_m)/a_mtau;
    yd[idx+3] = (a_hinf - a_h)/a_htau;
    yd[idx+4] = (kd_minf - kd_m)/kd_mtau;
    yd[idx+5] = (cat_minf-cat_m)/cat_mtau;
    yd[idx+6] = (cat_hinf-cat_h)/cat_htau;
    yd[idx+7] = (cas_minf-cas_m)/cas_mtau;
    yd[idx+8] = (-F*(ICaT+ICaS)-cai+CO)/TAUCA;
    yd[idx+9] = (kca_minf-kca_m)/kca_mtau;
    yd[idx+10] = (nap_minf-nap_m)/nap_mtau;
    yd[idx+11] = (nap_hinf-nap_h)/nap_htau;
    yd[idx+12] = (-(Ina_axon + Ik_axon + Il_axon+icoup1))/CM_AXON;
    yd[idx+13] = (k_minf - k_m)/k_mtau;
    yd[idx+14] = (na_minf - na_m)/na_mtau;
    yd[idx+15] = (na_hinf - na_h)/na_htau;
    yd[idx+16] = (proc_minf - proc_m)/0.5;
   
    return (0);
}

void ABNeuron::currents(realtype  t, N_Vector y, N_Vector ydot, ostream& out) {
    // realtype *yi;
    // yi = NV_DATA_S(y);
     
    // double v, h_m, a_m, a_h, kd_m, cat_m, cat_h, cas_m, cai, kca_m, nap_m, nap_h, v_a, k_m, na_m, na_h;
    // double Ih, Il, Ia, Ikd, ICaT, ICaS, IKCa, INaP, Ik_axon, Il_axon, Ina_axon, icoup0, icoup1, isum_soma;
    // double eca;

    // v = yi[k+0]; h_m=yi[k+1]; a_m = yi[k+2]; a_h = yi[k+3]; kd_m = yi[k+4]; cat_m = yi[k+5]; cat_h = yi[k+6];
    // cas_m = yi[k+7]; cai = yi[k+8]; kca_m = yi[k+9]; nap_m = yi[k+10]; nap_h = yi[k+11]; v_a = yi[k+12]; k_m = yi[k+13];
    // na_m=yi[k+14]; na_h=yi[k+15];
    
    // eca = 2.303*((8.314*310)/(2*(96400)))*log(CAO/cai)*1000;

    //      /* compute and send current value to file stream os passed to function */
    //      /* somatic currents */
    // Il = GL*(v-EL);
    // Ih=GH*h_m*(v-EH);
    // Ia = GA*pow(a_m, 4)*a_h*(v-EK);
    // Ikd = GKD*pow(kd_m, 4)*(v-EK);
    // ICaT = GCAT*pow(cat_m, 3)*cat_h*(v-eca);
    // ICaS = GCAS*pow(cas_m, 3)*(v-eca);
    // IKCa =GKCA*pow(kca_m, 4)*(v-EK);
    // INaP =GNAP*pow(nap_m, 3)*nap_h*(v-ENA);

    //      /* axon currents */
    // Ik_axon=AXON_GK*pow(k_m, 4)*(v_a-EK);
    // Il_axon=AXON_GL*(v_a-EL);
    // Ina_axon=AXON_GNA*pow(na_m, 3)*na_h*(v_a-ENA);
    
    //    icoup0 = GAXIAL*(v-v_a);
    // icoup1 = GAXIAL*(v_a-v);

    //      /* total somatic current */
    // isum_soma =Il + Ih + Ia + Ikd + INaP + ICaT + ICaS + IKCa + icoup0;

    // out << t  << "\t" << isum_soma <<  "\t" << v << endl;
}
