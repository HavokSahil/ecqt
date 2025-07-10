typedef struct Complex {
    float re;
    float im;
} Complex;

#define SMType Complex

float cabs(Complex a);

Complex cadd(Complex a, Complex b);
#define SMSum(a, b) cadd(a, b)

Complex cprod(Complex a, Complex b);
#define SMProd(a, b) cprod(a, b)

void print_complex(Complex a);
#define SM_DBG_PRINT_FN(a) print_complex(a)

#define SM_DEBUG

#include "../sparse-matrix.h"

#include <math.h>

#define EPSILON 1e-6f
#define CEQ(a, b) fabs(cabs(a) - cabs(b)) < EPSILON

int main() {
    int nrows = 3, ncols = 3;
    SparseMatrix *psm = sparse_matrix_init(nrows, ncols);
    assert(psm != NULL);
    assert(psm->nrows == nrows);
    assert(psm->ncols == ncols);
    assert(psm->nnz == 0);

    Complex e1 = {0.0f, 1.0f};
    Complex e2 = {1.0f, 0.0f};
    Complex e3 = {1.0f, 2.0f};
    Complex e4 = {1.0f, -0.3f};

    assert(sparse_matrix_add_entry(psm, e1, 0, 0) == 0);
    assert(sparse_matrix_add_entry(psm, e2, 0, 1) == 0);
    assert(sparse_matrix_add_entry(psm, e3, 1, 1) == 0);
    assert(sparse_matrix_add_entry(psm, e4, 2, 2) == 0);

    sparse_matrix_debug_print(psm);

    int dim = 3;
    Vec *pv = vec_init(dim);
    assert(pv != NULL);
    assert(pv->dim == dim);

    Complex v1 = {0.0f, 1.0f};
    Complex v2 = {0.2f, 0.3f};
    Complex v3 = {0.1f, 0.0f};

    pv->entries[0] = v1;
    pv->entries[1] = v2;
    pv->entries[2] = v3;

    sparse_matrix_debug_print_vec(pv);

    Vec *pr = vec_init(dim);
    assert(pr != NULL);
    assert(pr->dim == dim);

    Complex rl1 = {-1.0f, 0.0f};
    Complex rl2 = {-0.4f, 1.7f};
    Complex rl3 = {0.1f, -0.03f};
    // perform the left product on the sparse matrix
    assert(sparse_matrix_mult_vec_left(pv, psm, pr) == 0);
    sparse_matrix_debug_print_vec(pr);
    assert(CEQ(pr->entries[0], rl1));
    assert(CEQ(pr->entries[1], rl2));
    assert(CEQ(pr->entries[2], rl3));

    Complex rr1 = {-0.8f, 0.3f};
    Complex rr2 = {-0.4f, 0.7f};
    Complex rr3 = {0.1f, -0.03f};
    // perform the right product on the sparse matrix
    assert(sparse_matrix_mult_vec_right(psm, pv, pr) == 0);
    sparse_matrix_debug_print_vec(pr);
    assert(CEQ(pr->entries[0], rr1));
    assert(CEQ(pr->entries[1], rr2));
    assert(CEQ(pr->entries[2], rr3));

    printf("\n[PASS] test_sm_complex\n\n");
    return 0;
}

Complex cadd(Complex a, Complex b) {
    Complex res;
    res.re = a.re + b.re;
    res.im = a.im + b.im;
    return res;
}

Complex cprod(Complex a, Complex b) {
    Complex res;
    res.re = a.re * b.re - a.im * b.im;
    res.im = a.re * b.im + a.im * b.re;
    return res;
}

float cabs(Complex a) { return a.re * a.re + a.im * a.im; }

void print_complex(Complex a) { printf("%.2f+%.2fi ", a.re, a.im); }
