#ifndef _MYVECTOR_H
#define _MYVECTOR_H

typedef struct _generic_Vector *My_Vector; //define new type My_Vector that points to the _generic_Vector structure 

struct _generic_Vector {
    void *content; //implementation-dependent content field
    struct _My_VectorOps *ops; //operations field pointing to a structure with all possible generic vector operations
};

struct _My_VectorContent  //structure to hold the length and content of a _generic_Vector
{
    int length;
    double *data; //this could be of any type, not necessarily double (int, float or any typedef, class or another struct */
};

typedef struct _My_VectorContent *My_VectorContent; //define new type My_VectorContent that points to the _My_VectorContent structure

/* define a structure that holds all the operations that can be performed using My_Vector */
struct _My_VectorOps { 
    double*        (*getarraypointer)(My_Vector);
    void               (*product)(My_Vector, My_Vector, My_Vector);
    void               (*scale)(double, My_Vector);
    void               (*norm)(My_Vector, My_Vector);
};

typedef struct _My_VectorOps *My_VectorOps; //define a new type My_VectorOps that points to the _My_VectorOps structure

/* These macros give access to the contents of a _generic_Vector structure */
#define V_CONTENT(v)            ((My_VectorContent)(v->content)) /* retrieve content part of _generic_Vector casting to type _MyVector_Content */
#define V_DATA(v)                     ( V_CONTENT(v)->data ) /* returns a pointer to the first component of v */
#define V_Ith(v,i)                         (V_DATA(v)[i]) /* retrieve the ith element of v */
#define V_LENGTH(v)               (V_CONTENT(v)->length) /* retrieve the length of My_Vector */

void M_VScale(double c, My_Vector x) ;
double *M_VGetArrayPointer(My_Vector v);
My_Vector M_VNewEmpty(int length);
My_Vector M_VCreate( int length);
void VPrint(My_Vector x);
void VDestroy(My_Vector);

#endif
