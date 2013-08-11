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
    Neuron *ab = new ABNeuron("AB", par_ptr, NUM_PARS, NEQ, ICLAMP);
    ab->setIdx(0);
    Neuron *sc = new SimpleCell2("SC", par_ptr, NUM_PARS, 0, VCLAMP);
    Synapse *resonantSyn = new ChemSyn(sc, NULL, 1, NEQSYN, 125, 3, 100);
    resonantSyn->setIdx(NEQ);
    
    /* set up y vector and tolerances */
    numVars = NEQ + NEQSYN;
    cout << numVars <<"\n";
    yi = N_VNew_Serial(numVars);
    yo = N_VNew_Serial(numVars);
    yd = N_VNew_Serial(numVars);
    abstol = N_VNew_Serial(numVars);
    reltol = RTOL;

    /* initialize neuron and synapse state variables */
    ab->init(yi, PD_INIVARS);
    ab->setTol(abstol);
    resonantSyn->init(yi, SYN_INIVARS);
    resonantSyn->setTol(abstol);

    /* add neurons to the list */
    neurs.push_back(ab);
    neurs.push_back(sc);
    /* add synapses to the list */
    syns.push_back(resonantSyn);

        /* create synapses and add to list */
    /* Synapse *gapAB = new Gap(ab, pd, 1, 1);
       gapAB->setGmax(0);
       Synapse *gapPD = new Gap(pd, ab, 1, 1);
       gapPD->setGmax(0);
       syns.push_back(gapAB);
       syns.push_back(gapPD); 
    */
   
        /* create electrode to drive sc neuron membrane potential with ZAP*/
    Electrode *e = new Electrode(ab, 0.0);
    e->setWaveform(PULSE);
    e->setStart(2000);
    e->setDuration(1500);
    e->setBias(0);
    
    /* switch on numerical method given in command-line arguments */
    NeuronModel *model = new NeuronModel(&neurs, &syns, numVars, os);

    /*make integration classes for RK and CVode */
    RhsFn myfunc =&f;
    /*myfunc is the function to be solved by CVode */
    CVodeSolver *cvode = new CVodeSolver(numVars, yi, abstol, model, myfunc);
    gFn myRootFunc = &g;
    cvode->rootInit(2, myRootFunc); 

        /* class for RK fadvance */
    RK *rk = new RK(numVars, 2);

    double rk65_MINDT= 0.05;
    double rk65_eps= 1e-12;
    double rk65_relEps= 1e-9;
    double rk65_absEps= 1e-16;
    
    RK65n machine(numVars, rk65_MINDT, rk65_eps, rk65_absEps, rk65_relEps);
    
        /* append thread id to output file name */
    string filename = "SC";
    filename += id;
    filename += ".asc";
    ofstream out;
    out.open(filename.c_str());

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

    double dt= 0.0001;
    double dtx= 0.005;
    realtype t, tout;
    tout = 0.0; iout = 0;  tout = T1;
    /* set solver to either CVODE or RK65 If we need to handle
       spike-mediated events in synaptic mechanisms then RK65 is the
       choice 
    */
    /* for CVode specifiy a rootfinding function that when evaluated to 0 
       can signal that the integrator should be reinitialized
       provide a function g whose output is 0 when (t - tspk) = 0
       typedef int (*CVRootFn)(realtype t, N Vector y, realtype *gout,
       void *user data);

       To stop when the location of the discontinuity is determined by
       the solution (implicitly-defined discontinuity), use the
       rootfinding feature (available both in CVODE/CVODES and
       IDA/IDAS). In other words, the location of the discontinuity is
       the zero of some function g of the solution (i.e.,g(tdisc,
       y(tdisc))=0). Note that, in this situation, you cannot use
       CVodeSetStopTime to prevent the discontinuity from being "seen"
       during this integration phase. Indeed, the rootfinding
       algorithm relies on detecting sign changes in g and therefore
       needs to be allowed to evaluate it (and implicitly the RHS) on
       both sides of tdisc. Therefore, during this phase only, you
       need to have a smooth extension over the discontinuity so that
       the step across it and the subsequent rootfinding can be done
       efficiently.

       Call CVodeRootInit to specify the root function g with 2 components
       flag = CVodeRootInit(cvode_mem, 2, g);
       if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

    */
    int e0, e1;
    e0= 0; e1 = 0;
    
    
    solver = CVODE;
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
	e->setIapp(tout, 15, VCLAMP);
	sc->clampVoltage(tout, yi, NULL, out);
	out << tout << '\t' << e->getIapp()<< '\t' << ab->getVoltage() << "\t" << Ith(yi, 12) << "\t" << sc->getVoltage()<< '\t' << Ith(yi, 18) << "\n";
	flag = sc->detectMaximum(tout);
      }
      else if (solver == CVODE) {
	flag = cvode->fadvance(tout, t);

	if (flag == CV_ROOT_RETURN) {
	  /* CVODE succeeded and found 1 or more roots
	     This would be the time of a spike because the function
	     del(t - t_spk) would evaluate to 0
	     first test this by checking when voltage are passing through
	     a point on the rising phase
	  */
            flagr = cvode->getRootInfo(rootsfound);
            if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
            cout << rootsfound[0] << "\t" << rootsfound[1] << "\n";
            if (rootsfound[0] == 1) {
                    //cout << "The soma voltage is rising through -40\n";
            }
            if (rootsfound[1] == 1) {
                    //cout << "The axon voltage is rising through -40\n";
            }
        }
	if (check_flag(&flag, "CVode", 1)) {
            break;
	}
	if (flag == CV_SUCCESS) {
	  iout++;
	  tout += TSTEP;
	}

        e->setIapp(tout, -2, ICLAMP);
	//sc->clampVoltage(tout, yi, NULL, out);
        out << tout << '\t' << e->getIapp()<< '\t' << ab->getVoltage() << "\t" << Ith(yi, 12) << "\t" << sc->getVoltage()<< '\t' << Ith(yi, 18) << "\n";
	flagm = sc->detectMaximum(tout);
            /* if electrode was turned on or off, signal to the
               integrator to reinitialize
               only reinit when current changes
            */
        if (e1=(e->isEnabled()) != e0 && (count < 2)) {
            e0 = e1;
            cout << "Electrode pulse\n";
            cvode->reinit(tInj[count]);
            count++;
            cout << count << "sahsdfhshss\n";
            
        }
        
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


