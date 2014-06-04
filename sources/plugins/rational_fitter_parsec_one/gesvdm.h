#ifndef _gesvdm_h_
#define _gesvdm_h_

BEGIN_C_DECLS

typedef struct gesvdm_args_s {
    dplasma_qrtree_t *qrtree;
    int ib;
    int M;
    int np;
    int nq;
    int ny;
    int rc;

    /* For now the matrix with CI coefficients */
    tiled_matrix_desc_t *A;
    /* The matric of ci vectors */
    tiled_matrix_desc_t *ci;

    const void *data;
    void *rptr;

    void (*fillp)(int i0, int M, int N,
                  double *CI, int ldci,
                  const void *data,
                  void* rptr);

    void (*fillq)(int i0, int M, int N,
                  double *CI, int ldci,
                  const void *data,
                  void* rptr, int ny);

    int (*solve)(int np, int nq, int M,
                 double *ptr_CI, double *ptr_ci,
                 void *rptr);

} gesvdm_args_t;

END_C_DECLS

#endif /* _gesvdm_h_ */
