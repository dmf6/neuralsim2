#include "simplecell2.h"

SimpleCell2::SimpleCell2(const string name=NULL, double *inp=NULL, int inpno=0, int iniVarNo=0, int mode=0) 
    : Neuron(name, iniVarNo, mode), pno(inpno) {

    if (pno > 0) {
        p = new double[pno];
        set_p(inp);
    }
    iapp = 0;
    voltage=-60;
}


SimpleCell2::~SimpleCell2() {
    delete [] p;
}

void SimpleCell2::set_p(double *inp)
{
    for (int i= 0; i < pno; i++) {
        p[i]= inp[i];
    }
}

int SimpleCell2::derivative(realtype t, N_Vector *y, N_Vector *ydot, void *user_data) {
  // No ODES
  return (0);
}

void SimpleCell2::currents(realtype  t, N_Vector y) {}
