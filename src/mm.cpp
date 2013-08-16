#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <Python.h>
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
        /* pass this array to the PD constructor once filled with command-line arguments */
  
    string id;
    double upper_bound, lower_bound;
    realtype reltol;
    N_Vector yi, yo, yd, abstol;
    int solver;
    int flag, flagm, flagr, iout;
    int numVars;
    ofstream os;
    list<Neuron *> neurs;
    list<Synapse *> syns;
    int rootsfound[2];
    int tInj[2] = {2000, 3500};
 
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
    // Neuron *ab = new PDNeuron("AB", par_ptr, NUM_PARS, 17, ICLAMP);
    // ab->setIdx(0);
    Neuron *sc = new SimpleCell("SC", par_ptr, 15, 4, VCLAMP);
    sc->setIdx(0);
    
    // /* set up y vector and tolerances */
    numVars = 17;
    
    yi = N_VNew_Serial(numVars);
    yo = N_VNew_Serial(numVars);
    yd = N_VNew_Serial(numVars);
    abstol = N_VNew_Serial(numVars);
    reltol = RTOL;

    /* initialize neuron and synapse state variables */
    // ab->init(yi, AB_INIVARS);
    // ab->setTol(abstol);
    sc->init(yi, SC_INIVARS);
    sc->setTol(abstol);
    // /* add neurons to the list */
    // neurs.push_back(ab);
    neurs.push_back(sc);
    
    //     /* create electrode to drive sc neuron membrane potential with ZAP*/
    Electrode *e = new Electrode(sc, 0.0);
    e->setWaveform(ZAP);
    e->setStart(0);
    e->setDuration(TSTOP);
    e->setBias(-60);
    
    // /* switch on numerical method given in command-line arguments */
    NeuronModel *model = new NeuronModel(&neurs, &syns, numVars, os);

    // /*make integration classes for RK and CVode */
    RhsFn myfunc =&f;
    /*myfunc is the function to be solved by CVode */
    CVodeSolver *cvode = new CVodeSolver(numVars, yi, abstol, model, myfunc);
    gFn myRootFunc = &g;
    cvode->rootInit(2, myRootFunc); 

    //     /* class for RK fadvance */
    RK *rk = new RK(numVars, 2);

    double rk65_MINDT= 0.05;
    double rk65_eps= 1e-9;
    double rk65_relEps= 1e-9;
    double rk65_absEps= 1e-9;
    
    RK65n machine(numVars, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
    
        /* append thread id to output file name */
    string filename = "SC";
    filename += id;
    filename += ".asc";
    ofstream out;
    out.open(filename.c_str());
    
/* write the data in 10k blocks */
    char buf[10241];
    out.rdbuf()->pubsetbuf(buf, sizeof(buf));
    out.unsetf(ios::floatfield);     
    out.precision(15);
    setprecision(15);
    
    double dt= 0.001;
    double dtx= 0.05;
    realtype t, tout;
    tout = 0.0; iout = 0;  tout = 0;
    
    tout = T1;
    solver =  RK65;
    while(1) {
       if (solver == RK65) {
           e->setIapp(tout, 15, VCLAMP);
           
           // dtx= machine.integrate(tout, yi, yo, model, dt);
           // dtx= min(dtx, 2.0*dt);
           // swap(yi, yo);
           // dt= dtx;
           // tout+=dt;
           rk->fadvance(tout, yi, yo, yd, model, 0.05);
           swap(yi, yo);
           sc->currents(tout, yi, out);
           tout +=0.05;
           if (tout > TSTOP) {
               break;
           }
           
       }

      else if (solver == CVODE) {
              //e->setIapp(tout, 15, VCLAMP);
          flag = cvode->fadvance(tout, t);
          
    //       if (flag == CV_ROOT_RETURN) {
    //               /* CVODE succeeded and found 1 or more roots
    //                  This would be the time of a spike because the function
    //                  del(t - t_spk) would evaluate to 0
    //                  first test this by checking when voltage are passing through
    //                  a point on the rising phase
    //               */
    //           flagr = cvode->getRootInfo(rootsfound);
    //           if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
    //               //cout << rootsfound[0] << "\t" << rootsfound[1] << "\n";
    //           if (rootsfound[0] == 1) {
    //                   //cout << "The soma voltage is rising through -40\n";
    //           }
    //           if (rootsfound[1] == 1) {
    //                   // cout << "The axon voltage is rising through 0\n";
    //           }
    //       }
              //out << tout << '\t' << ab->getVoltage() << "\n";
          sc->currents(tout, yi, out);
          if (check_flag(&flag, "CVode", 1)) {
              break;
          }
 
          if (flag == CV_SUCCESS) {
              iout++;
              tout += 0.05;
          }
               
          if (tout > TSTOP) {
              break;
          }   
      }
    }

    // ab->writeStateVector(yi);
    // pd->writeStateVector(yi);

    
    /* free up memory */
    N_VDestroy_Serial(yi);
    N_VDestroy_Serial(yo);
    N_VDestroy_Serial(yd);
    N_VDestroy_Serial(abstol);
    delete  sc;
    
        //delete ab; // delete pd;
    // delete gapAB; delete gapPD;
    delete model;
    delete [] par_ptr;
    delete cvode;
    delete rk;
    delete e;
    
    out.close();
    return (0);
}

void swap(N_Vector y, N_Vector yout) {
    realtype *yi, *yn, tmp;
    yi = NV_DATA_S(y);
    yn = NV_DATA_S(yout);
    for(int i=0; i<NEQ; i++){
        tmp = yi[i];
        yi[i] = yn[i];
        yn[i] = tmp;
    }
}

int f(realtype t, N_Vector y, N_Vector dy, void *user_data) {
        /* for multiple cells cast user_data to NeuronModel */
        //cout << "My function has been advaned in time\n";
    NeuronModel *model = static_cast<NeuronModel *>(user_data);
    model->derivative(t, y, dy, NULL);
    return 0;   
} 

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data) {
  realtype vSoma, vAxon;
  /* to get the spike time, do exactly what we did in the function f
     user data is a pointer to user data, the same as the user data
     parameter passed to CVodeSetUserData.
     NeuronModel *model = static_cast<NeuronModel *>(user_data);
     get the required neuron from the list using its index
     and check whether its spiking or max variables are 1
     if so record the time
     the func will be
     gout[0] = t - tspk
  */
  vSoma = Ith(y,1); vAxon = Ith(y,13);
  gout[0] = vSoma - RCONST(-40);
  gout[1] = vAxon - RCONST(0.0);
  return(0);
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
