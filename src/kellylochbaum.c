#include <math.h>
#include <stdlib.h>
#include "m_pd.h"


static t_class *kellylochbaum_class;

typedef struct _kellylochbaum
{
    t_object  x_obj;
    t_inlet  *x_grc; 
    t_inlet  *x_dmp; 
    t_outlet  *x_outlet;
    int     x_sr;
    t_float x_in; //dummy
    //pointers to the delay bufs
    double  * x_f0;
    double * x_f0old;
    double * x_f1;
    double * x_f1old;
    double * x_b0;
    double * x_b0old;
    double * x_b1;
    double * x_b1old;
    double * x_cylrad;
    double * x_kcoeff;
    unsigned int x_ntubesec;
    unsigned int x_njunct;
} t_kellylochbaum;

static void kellylochbaum_clear(t_kellylochbaum *x){
    unsigned int i;
    for(i=0; i<x->x_njunct; i++){
        x->x_f0[i] = 0.;
        x->x_f0old[i] = 0.;
        x->x_f1[i] = 0.;
        x->x_f1old[i] = 0.;
        x->x_b0[i] = 0.;
        x->x_b0old[i] = 0.;
        x->x_b1[i] = 0.;
        x->x_b1old[i] = 0.;
    };
}

static t_int *kellylochbaum_perform(t_int *w)
{
    t_kellylochbaum *x = (t_kellylochbaum *)(w[1]);
    int n = (int)(w[2]);
    t_float *_xin = (t_float *)(w[3]);
    t_float *_grc = (t_float *)(w[4]);
    t_float *_dmp = (t_float *)(w[5]);
    t_float *out = (t_float *)(w[6]);
    unsigned int nj = x->x_njunct;
    unsigned int nts = x->x_ntubesec;
    int i;
    for(i=0; i<n;i++){
        unsigned int j;
        t_float xin = _xin[i];
        t_float grc = _grc[i];
        t_float dmp = _dmp[i];
        for(j=0; j<nj; j++) {
            // deep copy arrays
            x->x_f0old[j] = x->x_f0[j];
            x->x_f1old[j] = x->x_f1[j];
            x->x_b0old[j] = x->x_b0[j];
            x->x_b1old[j] = x->x_b1[j];
        };

        x->x_f0[0] = x->x_b0[0] * grc + xin;
        x->x_f1[0] = x->x_f0[0];
        x->x_b1[0] = (x->x_b1[1] * (1.0 + x->x_kcoeff[0]) + x->x_f1[0] * x->x_kcoeff[0]) * dmp;

        for(j=1;j<(nj-1);j++){
            //percolate input
            x->x_b0[j] = x->x_b1old[j];
            x->x_f0[j] = (x->x_f1[j-1] * (1.0 - x->x_kcoeff[j-1]) - x->x_b0[j] * x->x_kcoeff[j-1])*dmp;
            x->x_f1[j] = x->x_f0old[j-1];
            x->x_b1[j] = (x->x_b1[j+1] * (1.0 + x->x_kcoeff[j]) + x->x_f1[j] * x->x_kcoeff[j])*dmp;
        };
        x->x_b0[nj-1] = 0.0;
        x->x_f0[nj-1] = x->x_f1[nj-2]*(1.0 - x->x_kcoeff[nj-1])*dmp;
        out[i] = x->x_f0[nj-1];
    };
    
    return (w + 7);
}

static void kellylochbaum_dsp(t_kellylochbaum *x, t_signal **sp)
{
    int sr = sp[0]->s_sr;
    if(sr != x->x_sr){
        //if new sample rate isn't old sample rate, need to realloc
        x->x_sr = sr;
    };
    dsp_add(kellylochbaum_perform, 6, x, sp[0]->s_n,
	    sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
	    sp[3]->s_vec);
}

static void kellylochbaum_compute_kcoeff(t_kellylochbaum *x) {
    unsigned int i;
    unsigned int nj = x->x_njunct;
    for(i=0; i < nj; i++){
        x->x_kcoeff[i] = (x->x_cylrad[i]*x->x_cylrad[i]-x->x_cylrad[i+1]*x->x_cylrad[i+1])/(x->x_cylrad[i]*x->x_cylrad[i]+x->x_cylrad[i+1]*x->x_cylrad[i+1]);
    };
}

static void kellylochbaum_setcylrad(t_kellylochbaum *x, t_symbol *s, int argc, t_atom * argv){
  
    unsigned int idx = 0;
    t_float _idx = 0;
    t_float rad = 0.;
    unsigned int ntubesec = x->x_ntubesec;
    int argnum = 0; //current argument
    while(argc){
        if(argv -> a_type == A_FLOAT){
            t_float curf = atom_getfloatarg(0, argc, argv);
            switch(argnum){
                case 0:
                    _idx = curf > 0 ? curf : 0; 
                    break;
                case 1:
                    rad = curf > 0 ? curf : 0; 
                    break;
                default:
                    break;
            };
            argnum++;
        };
        argc--;
        argv++;
    };
    idx = (unsigned int)idx;
    idx = idx < ntubesec ? idx : ntubesec;
    x->x_cylrad[idx] = rad;
    kellylochbaum_compute_kcoeff(x);
}

static void kellylochbaum_setcylrads(t_kellylochbaum *x, t_symbol *s, int argc, t_atom * argv){
  
    int ntubesec = (int) x->x_ntubesec;
    int argnum = 0; //current argument
    while(argc){
        if(argv -> a_type == A_FLOAT){
            t_float curf = atom_getfloatarg(0, argc, argv);
            if(argnum < ntubesec) {
                x->x_cylrad[argnum] = curf > 0 ? curf : 0;
            };
            argnum++;
        };
        argc--;
        argv++;
    };
    kellylochbaum_compute_kcoeff(x);
}


static void *kellylochbaum_new(t_symbol *s, int argc, t_atom * argv){
    t_kellylochbaum *x = (t_kellylochbaum *)pd_new(kellylochbaum_class);
   
    //defaults
    t_float grc = 0.969;
    t_float dmp = 1.0;
    unsigned int i;
    unsigned int njunct = 0;
    unsigned int ntubesec = 1;
    t_float _ntubesec = 1;
    x->x_sr = sys_getsr();
    
    int argnum = 0; //current argument
    while(argc){
        if(argv -> a_type == A_FLOAT){
            t_float curf = atom_getfloatarg(0, argc, argv);
            switch(argnum){
                case 0:
                    _ntubesec = curf;
                    break;
                default:
                    break;
            };
            argnum++;
        };
        argc--;
        argv++;
    };
    
    
    _ntubesec = _ntubesec > 0 ? _ntubesec : 0;
    if(_ntubesec < 2) {
        goto errstate;
    }
    else {
        ntubesec = (unsigned int)_ntubesec;
        njunct = (unsigned int)ntubesec - 1;
        x->x_ntubesec = ntubesec;
        x->x_njunct = njunct;
    };
    //boundschecking
    //this is 1/44.1 (1/(sr*0.001) rounded up, good enough?

    x->x_f0 = (double *)malloc(sizeof(double)*njunct);
    x->x_f0old = (double *)malloc(sizeof(double)*njunct);
    x->x_f1 = (double *)malloc(sizeof(double)*njunct);
    x->x_f1old = (double *)malloc(sizeof(double)*njunct);
    x->x_b0 = (double *)malloc(sizeof(double)*njunct);
    x->x_b0old = (double *)malloc(sizeof(double)*njunct);
    x->x_b1 = (double *)malloc(sizeof(double)*njunct);
    x->x_b1old = (double *)malloc(sizeof(double)*njunct);
    x->x_kcoeff = (double *)malloc(sizeof(double)*njunct);
    x->x_cylrad = (double *)malloc(sizeof(double)*ntubesec);

    for(i=0; i < ntubesec; i++){
        x->x_cylrad[i] = 1.0;
    };

    kellylochbaum_compute_kcoeff(x);

    //inlets outlets
    x->x_grc = inlet_new((t_object *)x, (t_pd *)x, &s_signal, &s_signal);
    pd_float((t_pd *)x->x_grc, grc);
    x->x_dmp = inlet_new((t_object *)x, (t_pd *)x, &s_signal, &s_signal);
    pd_float((t_pd *)x->x_dmp, dmp);
    x->x_outlet = outlet_new((t_object *)x, &s_signal);
    return (x);

    errstate:
        pd_error(x, "kellylochbaum~: improper args");
        return NULL;
}

static void * kellylochbaum_free(t_kellylochbaum *x){
    free(x->x_f0);
    free(x->x_f0old);
    free(x->x_f1);
    free(x->x_f1old);
    free(x->x_b0);
    free(x->x_b0old);
    free(x->x_b1);
    free(x->x_kcoeff);
    free(x->x_cylrad);
    inlet_free(x->x_grc);
    inlet_free(x->x_dmp);
    outlet_free(x->x_outlet);
    return (void *)x;
}

void kellylochbaum_tilde_setup(void)
{
    kellylochbaum_class = class_new(gensym("kellylochbaum~"),
			   (t_newmethod)kellylochbaum_new,
			   (t_method)kellylochbaum_free,
			   sizeof(t_kellylochbaum), 0,
			   A_GIMME, 0);
    CLASS_MAINSIGNALIN(kellylochbaum_class, t_kellylochbaum, x_in);
    //class_addmethod(kellylochbaum_class, nullfn, gensym("signal"), 0);
    class_addmethod(kellylochbaum_class, (t_method)kellylochbaum_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(kellylochbaum_class, (t_method)kellylochbaum_clear, gensym("clear"), 0);
    class_addmethod(kellylochbaum_class, (t_method)kellylochbaum_setcylrad, gensym("cylrad"), A_GIMME, 0);
    class_addmethod(kellylochbaum_class, (t_method)kellylochbaum_setcylrads, gensym("cylrads"), A_GIMME, 0);

}
