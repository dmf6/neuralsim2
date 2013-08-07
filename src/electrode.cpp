#include "electrode.h"
#include "random.h"


Electrode::Electrode(Neuron *intarget, realtype inI)  {
    current = inI;
    num_precycles = 3;
    freq = 0;
    target= intarget;
    target->electrodes.push_back(this);
    waveform = 0;
    bias = 0.0;
}

void Electrode::setWaveform(int w) {
    waveform = w;
}

realtype Electrode::getIapp() {  
     return current;  
}

void Electrode::setIapp(realtype t, realtype amp, int mode)  {
    double L = log(FMAX/FMIN)/duration;
        //double precycle_dur = 0.1*num_precycles*duration;
    switch(waveform) {
        case SINE:
            current = bias + amp*sin(2*PI*0.75*t/1000);
            break;
            
        case ZAP:
                // freq = (t < 30000) ? 0.1*t : (FMIN*((exp(L*(t-30000))-1)/(L*(t-30000))))*(t-30000);
            freq = (FMIN*((exp(L*(t))-1)/(L*(t))));
                /* for voltage clamp mode */
            if (mode==VCLAMP) {
                current = bias + amp+amp*sin(2*PI*freq/1000*t + 3*(PI/2));
            }
                /* for current clamp mode */
            else {  
                current =  bias +amp*sin(2*PI*freq*t/1000 + 3*(PI/2));
            }
            break;

        case PULSE:
                /* command pulse of amplitude amp between START and
                 * START + DURATION */
            if (((t - start) > 0) && (t < (start + duration))) {
                current = bias + amp;
                    // cout << current << "\t" << bias << "\t" << amp << "\n";
                
            }
            else {
                current = bias;
            }
            
            break;

        default:
            cout << "Command waveform type not recognized!" << endl;
    }
    
}

realtype Electrode::getFreq() {
    return freq;
}

void Electrode::setDuration(int dur) {
    duration = dur;
}

void Electrode::setStart(int s) {
    start = s;
}

void Electrode::setBias(realtype b) {
    bias = b;
}
void Electrode::setAmplitude(realtype a) {
    amp = a;
}
