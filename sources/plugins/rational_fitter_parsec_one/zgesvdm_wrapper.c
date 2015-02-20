/*
 * Copyright (c) 2011-2012 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "dague_internal.h"
#include <plasma.h>
#include "dplasma.h"
#include "dplasma/lib/dplasmatypes.h"
#include "dplasma/lib/dplasmaaux.h"
#include "dplasma/lib/memory_pool.h"

#include "gesvdm.h"
#include "zgesvd_multiple.h"

dague_handle_t* dplasma_zgesvd_multiple_New( gesvdm_args_t *args )
{
    dague_zgesvd_multiple_handle_t* handle;
    tiled_matrix_desc_t *A = args->A;

    /*
     * TODO: We consider ib is T->mb but can be incorrect for some tricks with GPU,
     * it should be passed as a parameter as in getrf
     */
    handle = dague_zgesvd_multiple_new( args, (dague_ddesc_t*)(A),
                                        (dague_ddesc_t*)(args->ci),
                                        *(args->qrtree), args->ib,
                                        NULL, NULL/*, NULL*/ );

    handle->p_work = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( handle->p_work, args->ib * A->nb * sizeof(dague_complex64_t) );

    handle->p_tau = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( handle->p_tau, A->nb * sizeof(dague_complex64_t) );

    /* handle->p_superb = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t)); */
    /* dague_private_memory_init( handle->p_superb, 5 * dague_imin(A->nb, A->m) * sizeof(dague_complex64_t) ); */

    /* Default type */
    dplasma_add2arena_tile( handle->arenas[DAGUE_zgesvd_multiple_DEFAULT_ARENA],
                            A->mb*A->nb*sizeof(dague_complex64_t),
                            DAGUE_ARENA_ALIGNMENT_SSE,
                            MPI_DOUBLE_COMPLEX, A->mb );

    /* Upper triangular part of tile with diagonal */
    dplasma_add2arena_upper( handle->arenas[DAGUE_zgesvd_multiple_UPPER_TILE_ARENA],
                             A->mb * A->nb * sizeof(dague_complex64_t),
                             DAGUE_ARENA_ALIGNMENT_SSE,
                             MPI_DOUBLE_COMPLEX, A->mb, 1 );

    /* /\* Lower triangular part of tile without diagonal *\/ */
    /* dplasma_add2arena_lower( handle->arenas[DAGUE_zgesvd_multiple_LOWER_TILE_ARENA], */
    /*                          A->mb*A->nb*sizeof(dague_complex64_t), */
    /*                          DAGUE_ARENA_ALIGNMENT_SSE, */
    /*                          MPI_DOUBLE_COMPLEX, A->mb, 0 ); */

    /* Vector */
    dplasma_add2arena_rectangle( handle->arenas[DAGUE_zgesvd_multiple_VECTOR_ARENA],
                                 args->ci->mb * args->ci->nb * sizeof(double),
                                 DAGUE_ARENA_ALIGNMENT_SSE,
                                 MPI_DOUBLE, 1, A->nb, -1);

    return (dague_handle_t*)handle;
}

int dplasma_zgesvd_multiple( dague_context_t *dague,
                             gesvdm_args_t *args )
{
    dague_handle_t *dague_zgesvd_multiple = NULL;

    dague_zgesvd_multiple = dplasma_zgesvd_multiple_New( args );
    dague_enqueue(dague, (dague_handle_t*)dague_zgesvd_multiple);
    dplasma_progress(dague)
    dplasma_zgesvd_multiple_Destruct( dague_zgesvd_multiple );
    return 0;
}

void
dplasma_zgesvd_multiple_Destruct( dague_handle_t *o )
{
    dague_zgesvd_multiple_handle_t *dague_zgesvd_multiple = (dague_zgesvd_multiple_handle_t *)o;

    dplasma_datatype_undefine_type( &(dague_zgesvd_multiple->arenas[DAGUE_zgesvd_multiple_DEFAULT_ARENA   ]->opaque_dtt) );
    /* dplasma_datatype_undefine_type( &(dague_zgesvd_multiple->arenas[DAGUE_zgesvd_multiple_LOWER_TILE_ARENA]->opaque_dtt) ); */
    dplasma_datatype_undefine_type( &(dague_zgesvd_multiple->arenas[DAGUE_zgesvd_multiple_UPPER_TILE_ARENA]->opaque_dtt) );
    dplasma_datatype_undefine_type( &(dague_zgesvd_multiple->arenas[DAGUE_zgesvd_multiple_VECTOR_ARENA    ]->opaque_dtt) );

    /* dague_private_memory_fini( dague_zgesvd_multiple->p_superb ); */
    dague_private_memory_fini( dague_zgesvd_multiple->p_work );
    dague_private_memory_fini( dague_zgesvd_multiple->p_tau  );
    /* free( dague_zgesvd_multiple->p_superb ); */
    free( dague_zgesvd_multiple->p_work );
    free( dague_zgesvd_multiple->p_tau  );

    DAGUE_INTERNAL_HANDLE_DESTRUCT(dague_zgesvd_multiple);
}
