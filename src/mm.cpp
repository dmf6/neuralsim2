#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include "mm.h"

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
    string id;
    double upper_bound, lower_bound;
    realtype t, tout, reltol;
    N_Vector yi, yo, yd, abstol;
    int solver;
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
    
    /* create neuron and synapse objects */
    Neuron *ab = new ABNeuron("AB", par_ptr, NUM_PARS, NEQ, ICLAMP);
    ab->setIdx(0);
    
    /* set up y vector and tolerances */
    numVars = NEQ;
    cout << numVars <<"\n";
    yi = N_VNew_Serial(numVars);
    yo = N_VNew_Serial(numVars);
    yd = N_VNew_Serial(numVars);
    abstol = N_VNew_Serial(numVars);
    reltol = RTOL;

    /* initialize neuron and synapse state variables */
    ab->init(yi, PD_INIVARS);
    ab->setTol(abstol);

    /* add neurons to the list */
    neurs.push_back(ab);

        /* create electrode to drive sc neuron membrane potential with ZAP*/
    Electrode *e = new Electrode(ab, 0.0);
    //e->setWaveform(ZAP);
    //e->setDuration(TSTOP);
    //e->setBias(-60);
    
    /* switch on numerical method given in command-line arguments */
    NeuronModel *model = new NeuronModel(&neurs, &syns, numVars, os);

    /*make integration classes for RK and CVode */
    RhsFn myfunc =&f;
    /*myfunc is the function to be solved by CVode */
    CVodeSolver *cvode = new CVodeSolver(numVars, yi, abstol, model, myfunc);
    cvode->setStopTime(TSTOP);
 
        /* class for RK fadvance */
    RK *rk = new RK(numVars, 2);

    double rk65_MINDT= 0.05;
    double rk65_eps= 1e-12;
    double rk65_relEps= 1e-9;
    double rk65_absEps= 1e-16;
    //double mindt=1e-6;
    
    RK65n machine(numVars, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
    
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
    tout = T1;
    double dt= 0.0001;
    double dtx= 0.005;

    /* set solver to either CVODE or RK65 If we need to handle
       spike-mediated events in synaptic mechanisms then RK65 is the
       choice 
    */
    
    solver = RK65;
    while(1) {
      if (solver == RK65) {
	dtx= machine.integrate(tout, yi, yo, model, dt);
	dtx= min(dtx, 2.0*dt);
	swap(yi, yo);
	dt= dtx;
	tout+=dt; 
	if (tout > TSTOP) {
	  break;
	}
	out << tout << '\t' << e->getIapp()<< '\t' << ab->getVoltage() <<   '\n';
	//sc->clampVoltage(tout, yi, NULL, out);
	//cout << "MAX: " << sc->detectMaximum(tout) << "\n";
      }
      else if (solver == CVODE) {
	flag = cvode->fadvance(tout, t);
	if (check_flag(&flag, "CVode", 1)) {
	  break;
	}
	if (flag == CV_SUCCESS) {
	  tout += TSTEP;
	}
	out << tout << '\t' << e->getIapp()<< '\t' << ab->getVoltage() <<   '\n';
	if (tout > TSTOP) {
	  break;
	}   
      }
    }
    
    /* free up memory */
    delete e;
    delete ab;
    delete model;
    delete [] par_ptr;
    delete cvode;
    delete rk, machine;
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
      return(1); }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); 
  }

  return(0);
}


  
