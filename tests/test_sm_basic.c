#define SMType float
#define SMProd(a, b) (a * b)
#define SMSum(a, b) (a + b)
#define SM_PRINT_FMT "%f"

#ifndef SM_DEBUG
#include <assert.h>
#include <stdio.h>
#endif // SM_DEBUG

#include "../sparse-matrix.h"

#include <math.h>

int main() {
    int nrows = 3, ncols = 3;
    SparseMatrix *psm = sparse_matrix_init(nrows, ncols);
    assert(psm != NULL);
    assert(psm->nrows == nrows && psm->ncols == ncols);

    /**
     *Matrix: 1.10   0   0
     *           0 1.2 1.7
     *           0   0   0
     */
    assert(sparse_matrix_add_entry(psm, 1.1, 0, 0) == 0);
    assert(sparse_matrix_add_entry(psm, 1.2, 1, 1) == 0);
    assert(sparse_matrix_add_entry(psm, 1.7, 1, 2) == 0);

    float epsilon = 1e-6f;
    float e1, e2, e3;
    assert(sparse_matrix_get_entry(psm, 0, 0, &e1) == 0);
    assert(sparse_matrix_get_entry(psm, 1, 1, &e2) == 0);
    assert(sparse_matrix_get_entry(psm, 1, 2, &e3) == 0);

    assert(fabs(e1 - 1.1) < epsilon);
    assert(fabs(e2 - 1.2) < epsilon);
    assert(fabs(e3 - 1.7) < epsilon);

#if defined(SM_DEBUG)
    sparse_matrix_debug_print(psm);
#endif

    int dim = 3;
    Vec *pv = vec_init(dim);
    assert(pv != NULL);
    assert(pv->dim == dim);

    pv->entries[0] = 1.2f;
    pv->entries[1] = 1.7f;
    pv->entries[2] = 1.8f;

#if defined(SM_DEBUG)
    sparse_matrix_debug_print_vec(pv);
#endif

    Vec *pr = vec_init(dim);
    assert(pr != NULL);
    assert(pr->dim == dim);

    // perform left multiplication with the vector;
    assert(sparse_matrix_mult_vec_left(pv, psm, pr) == 0);

    /*
    ** Result: [1.32 2.04 2.89]
    */
    assert(fabs(pr->entries[0] - 1.32) < epsilon);
    assert(fabs(pr->entries[1] - 2.04) < epsilon);
    assert(fabs(pr->entries[2] - 2.89) < epsilon);

#if defined(SM_DEBUG)
    sparse_matrix_debug_print_vec(pr);
#endif

    sparse_matrix_deinit(psm);
    vec_deinit(pv);
    vec_deinit(pr);

    printf("\n[PASS] test_sm_basic\n\n");
    return 0;
}
