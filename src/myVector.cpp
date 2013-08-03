#include <iostream>
#include <cstdlib>
#include <math.h>
#include "myVector.h"

using namespace std;

void M_VScale(double c, My_Vector x) {
    int i, n;
    double *xv;
    xv = NULL;

    n  = V_LENGTH(x);
    xv = (double *) V_DATA(x);
    
    for (i = 0; i < n; i++) {
        xv[i] *= c;
    }
}

void product(My_Vector x, My_Vector y, My_Vector z) 
{
    int i, n;
    double *xv, *yv, *zv;
    xv = yv = zv = NULL;

    n = V_LENGTH(x);
    xv = V_DATA(x);
    yv = V_DATA(y);
    zv = V_DATA(z);

    for (i=0; i < n; i++) {
        zv[i] = xv[i] * yv[i];
    }
}

 double *M_VGetArrayPointer(My_Vector v) {
  return (V_DATA(v));
}

void VPrint(My_Vector x)
{
    int i, n;
    double *xv;

    xv = NULL;

    n  = V_LENGTH(x);
    xv = V_DATA(x);

    for (i = 0; i < n; i++) {
        cout << xv[i] << endl;
    }
    return;
}

/* Function to create a new empty vector using C-style structs*/
My_Vector M_VNewEmpty(int length) {
    My_Vector v;
    My_VectorOps ops;
    My_VectorContent content;

    v = NULL;
    v = (My_Vector) malloc(sizeof*v);

    ops = (My_VectorOps) malloc(sizeof(My_VectorOps));
    ops->scale           = M_VScale;
    ops->getarraypointer = M_VGetArrayPointer;
    
    
    content = (My_VectorContent) malloc(sizeof(My_VectorContent));

    content ->length = length;
    content->data = NULL;

    v->content = content;
    v->ops         = ops;

    return v;
}

My_Vector M_VCreate(int length)
{
    My_Vector v;
    double *data;

    v = NULL;
    v = M_VNewEmpty(length);
    if (v == NULL) {
        return(NULL);
    }
    
        /* Create data */
    if (length > 0) {

            /* Allocate memory */
        data = NULL;
        data = (double *) malloc(length * sizeof(double));
            //if data is null then free memory allocated and return using cerr
        
            /* Attach data */
        V_DATA(v)     = data;
    }
    
    return(v);
}

void VDestroy(My_Vector v)
{
    free(V_DATA(v));
    V_DATA(v) = NULL;
  
    free(v->content); v->content = NULL;
    free(v->ops); v->ops = NULL;
    free(v); v = NULL;

    return;
}
