#ifndef _ELECTRODE_H
#define _ELECTRODE_H

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <math.h>
#include "neuron.h"
#include "random.h"

#define PI RCONST(3.1415926535898) 
#define ZAP 1
#define PULSE 2
#define SINE 3

#define FMIN 0.1
#define FMAX 4.0
#define DUR 100000

#define VCLAMP 1

class Neuron;

class Electrode {
    
  public:
    realtype amp;
    realtype current;
    realtype freq;
    int duration;
    int start;
    Neuron *target;
    int num_precycles;
    int waveform;
    realtype bias;
    Random rand;
    int enabled;
    
        /* target neuron and amplitude (amplitude is common to both zap and pulse) */
    Electrode(Neuron *, realtype inI);
    
    realtype getIapp();
    void setIapp(realtype t, realtype amp, int mode);
        /* electrode derives abstract base class Synapse so we have to provide implementations of pure virtual functions */
    realtype getFreq();
    
    void set_gsyn(double) { }
    int derivative(realtype, N_Vector, N_Vector) { return 1;}
    void setDuration(int dur);
    void setStart(int start);
    void setWaveform(int w);
    void setBias(realtype b);
    void setAmplitude(realtype);
    int isEnabled();
};

#endif
