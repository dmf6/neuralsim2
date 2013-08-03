#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
//#include <ostream>
//#include <boost/iostreams/device/file.hpp>
//#include <boost/iostreams/stream.hpp>
#include "mm.h"
//#include <boost/foreach.hpp>
#include <algorithm>
#include <iomanip>
#include <sstream>

//#define foreach_         BOOST_FOREACH

typedef struct {
    realtype p[3];
} *UserData;

using namespace std;

char *program_name;
static void show_usage(string name) {
    cerr << "Usage: " << name << " option(s) -i -p\n"
         << "Options:\n"
         << "\t-h, --help\t\t Show help message\n"
         << "\t-p, --params\t\t Specify parameter string"
         << endl;
}


int main(int argc, char *argv[]) {
        /* parse command line arguments to extract parameters */
    program_name = argv[0];
        /* pass this array to the PD constructor once filled with command-line arguments */
  
    string id;
    double upper_bound, lower_bound;
    realtype t, tout, reltol;
    void *cvode_mem;
    N_Vector yi, yo, yd, abstol;
    int flag, iout;
    int numVars;
    ofstream os;
    list<Neuron *> neurs;
    list<Synapse *> syns;
    
    int mode = ICLAMP;
    double *par_ptr = new double[NUM_PARS];
         
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if  ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        }
        else if ((arg == "-i")  || (arg == "identifier")) {
                /* save id */
            id = argv[++i];
                //cout << "The thread id is " << id << endl;
        }
        else if ((arg == "-vub") || (arg == "upper_bound")) {
            upper_bound = atof(argv[++i]);
        }
        else if ((arg == "-vlb") || (arg == "lower_bound")) {
            lower_bound = atof(argv[++i]);
        }
        else if ((arg == "-p") || (arg == "params")) {
            // cout << "The parameters provided are:\n";
            for (int  j = i+1, k = 0; k < NUM_PARS; j++, k++) {
                par_ptr[k] = atof(argv[j]);
            }
        }
    }
        /* create neurons */
    Neuron *ab = new ABNeuron("AB", par_ptr, NUM_PARS, NEQ, mode);
    Neuron *pd = new PDNeuron("PD", par_ptr, NUM_PARS, NEQ, mode);
    
        /*initialize the y vector using the number of state variables
         * for each neuron */
    int abVarNo, pdVarNo;
    int nVars =0;
    abVarNo = ab->getIVarNo();
    pdVarNo = pd->getIVarNo();
    numVars = abVarNo + pdVarNo;
    cout << "Number of Variables " << numVars << endl;

    yi = N_VNew_Serial(numVars);
             //yo = N_VNew_Serial(numVars);
             //yd = N_VNew_Serial(numVars);
    abstol = N_VNew_Serial(numVars);
    reltol = RTOL;
    
    ab->init(yi, nVars, PD_INIVARS, abVarNo);
    ab->setTol(abstol, nVars, abVarNo);
    nVars+=abVarNo;
    pd->init(yi, nVars, PD_INIVARS, pdVarNo);
    pd->setTol(abstol, nVars, pdVarNo);
    
           /* add neurons to the list */
    neurs.push_back(ab);
    neurs.push_back(pd);
    
        /* create synapses and add to list */
    Synapse *gapAB = new Gap(ab, pd, 1, 1);
    gapAB->setGmax(0);
    Synapse *gapPD = new Gap(pd, ab, 1, 1);
    gapPD->setGmax(0);

    syns.push_back(gapAB);
    syns.push_back(gapPD);
    
        /* create electrode */
    Electrode *e = new Electrode(ab, 0.0);

        /* switch on numerical method given in command-line arguments */
    NeuronModel *model = new NeuronModel(&neurs, &syns, numVars, os);

        /*make integration classes for RK and CVode */
    RhsFn myfunc =&f;
        /*myfunc is the function to be solved by CVode */
    CVodeSolver *cvode = new CVodeSolver(numVars, yi, abstol, model, myfunc);

        /* class for RK fadvance */
        /*int order = 2;
         * RK *rk = new RK(numVars, order);
         */

        /* append thread id to output file name */
    string filename = "SC";
    filename += id;
    filename += ".asc";
    ofstream out;
    out.open(filename.c_str());
     
    realtype isyn = 0.0; tout = 0.0;
    int count = 0;

     out.unsetf(ios::floatfield);            // floatfield not set
     out.precision(15);
     setprecision(15);
         /* write the data in 10k blocks */
     char buf[10241];
     out.rdbuf()->pubsetbuf(buf, sizeof(buf));

     
         /*****************************************************************
         *                     start CVode loop (only use for current clamp)
         *
         ******************************************************************/  
    iout = 0;  tout = T1;
    while(1) {
            //rk->fadvance();
        
        flag = cvode->fadvance(tout, t);

        if (check_flag(&flag, "CVode", 1)) {
            break;
        }
        if (flag == CV_SUCCESS) {
            iout++;
            tout += TSTEP;
        }
        out << tout << '\t' << e->getIapp()<< '\t' << ab->getVoltage() << '\t' <<  pd->getVoltage() <<   '\n';
        if (tout > TSTOP) {
            break;
        }
    }
     /* free up memory */
    delete e;
    delete ab; delete pd;
    delete gapAB; delete gapPD;
    delete model;
    delete [] par_ptr;
    delete cvode;
    out.close();
}

void swap(N_Vector y, N_Vector yout) {
    realtype *yi, *yo;
   
    yi = NV_DATA_S(y);
    yo = NV_DATA_S(yout);
    
    for (int i = 0; i < NEQ; i++) {
        yi[i] = yo[i];
    }
}

int f(realtype t, N_Vector y, N_Vector dy, void *user_data) {
        /* for multiple cells cast user_data to NeuronModel */
        //cout << "My function has been advaned in time\n";
    NeuronModel *model = static_cast<NeuronModel *>(user_data);
    model->derivative(t, y, dy, NULL);
    return 0;   
} 


static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


     // while(1) {
     //     if (mode == VCLAMP) {
     //             //go into voltage clamp mode
     //             //e->setIapp(tout, 16.5, VCLAMP);
     //         e->setIapp(tout, (upper_bound+fabs(lower_bound))/2, VCLAMP);
     //         Ith(yi, 1) = e->getIapp();
     //         sc->currents(tout, yi, yd, out);
     //     }
     //     else {
     //             //go into current clamp mode
     //             //e->setIapp(tout, 1.8, ICLAMP);
     //             //e->setIapp(tout, 1.35, ICLAMP);
     //             //sc->currents(tout, yi, yd, out);
     //         out << tout << '\t' << e->getIapp()<< '\t' << Ith(yi, 1) << '\n';
     //     }
     //     rk->integrate(tout, yi, yo, yd, model, TSTEP);
     //         //cout << tout << '\t' << e->getIapp() << '\t' << TSTOP << '\n';
         
     //     if (tout > TSTOP) {
     //         break;
     //     }
     //     tout+=TSTEP;
     //     swap(yi, yo);
     // }
