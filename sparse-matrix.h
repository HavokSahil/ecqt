/**
 * MIT License
 *
 * Copyright (c) 2025 Sahil
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file sparsematrix.h
 * @brief Sparse matrix module with configurable precision and vector
 * multiplication support.
 *  - This module provides a basic sparse matrix representation using coordinate
 * format (COO), along with vector multiplication (left and right), transpose,
 * and debugging support.
 */

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

/**
 * @def Type
 * @brief Define the precision type used in the sparse matrix and vector
 * entries. Default is double. You can override it before including this header.
 */
#ifndef SMType
#define SMType double
#define SMProd(a, b) (a * b)
#define SMSum(a, b) (a + b)
#endif

/**
 * @def SMProd
 * @brief There should be a user defined product function for the type.
 * @usage Example:
 * Type product(Type a, Type b) { return a * b; }
 * #define SMProd(a, b) product(a, b);
 */
#ifndef SMProd
#error "product not defined."
#endif

/**
 * @def SMSum
 * @brief There should be a user defiend sum function for the type.
 * @usage Example:
 * Type sum(Type a, Type b) { return a + b; }
 * #define SMSum(a, b) sum(a, b);
 */
#ifndef SMSum
#error "sum not defined."
#endif

/**
 * @def SM_PRINT_FMT
 * @brief Print format for matrix entries in debug mode.
 * Adjusts automatically depending on the `Type` macro.
 */
#ifndef SM_PRINT_FMT
#define SM_PRINT_FMT "%.3f"
#endif // SM_PRINT_FMT

#define MX_15_BIT 0x7FFF

// Limit checks for configuration macros
#if defined(SM_MAX_NZERO) && SM_MAX_NZERO > MX_15_BIT
#error "Number of non zero entries exceed the limit"
#endif

#if defined(SM_MAX_ROWS) && SM_MAX_ROWS > MX_15_BIT
#error "Number of rows exceed the max limit"
#endif

#if defined(SM_MAX_COLS) && SM_MAX_COLS > MX_15_BIT
#error "Number of cols exceed the max limit"
#endif

// Default configuration if not already defined
#if !defined(SM_MAX_NZERO)
#define SM_MAX_NZERO MX_15_BIT
#endif

#if !defined(SM_MAX_ROWS)
#define SM_MAX_ROWS MX_15_BIT
#endif

#if !defined(SM_MAX_COLS)
#define SM_MAX_COLS MX_15_BIT
#endif

// Utility macros
#if !defined(MIN)
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#if !defined(MAX)
#define MAX(a, b) (((a) < (b)) ? (b) : (a))
#endif

#include <stdlib.h>
#include <string.h>

#ifdef SM_DEBUG
#include <assert.h>
#include <stdio.h>
#endif // SM_DEBUG

/**
 * @struct Vec
 * @brief Represents a vector with configurable precision.
 */
typedef struct Vec {
    SMType *entries; /**< Pointer to vector elements */
    int dim;         /**< Dimension of the vector */
} Vec;

/**
 * @struct SparseMatrix
 * @brief Represents a sparse matrix in coordinate format (COO).
 */
typedef struct SparseMatrix {
    SMType *nzentries; /**< Non-zero entries */
    int *prows;        /**< Row indices corresponding to each entry */
    int *pcols;        /**< Column indices corresponding to each entry */
    int nnz;           /**< Current number of non-zero entries */
    int nrows;         /**< Number of matrix rows */
    int ncols;         /**< Number of matrix columns */
} SparseMatrix;

/**
 * @brief Allocate and initialize a vector of given dimension.
 * @param dim Dimension of the vector
 * @return Pointer to the initialized Vec or NULL on failure
 */
Vec *vec_init(int dim);

/**
 * @brief Allocate and initialize a vector of given dimension and allocated
 * memory safely (turns the initial pointer to NULL, taking the ownership)
 * @param dim Dimension of the vector
 * @param pmem Pointer to the allocated memory
 * @return Pointer to the initialized Vec or NULL on failure
 */
Vec *vec_init_mem(int dim, SMType **pmem);

/**
 * @brief Deallocate a vector and its memory.
 * @param pv Pointer to the vector to be deallocated
 */
void vec_deinit(Vec *pv);

/**
 * @brief Initialize a sparse matrix with given dimensions.
 * @param nrows Number of rows
 * @param ncols Number of columns
 * @return Pointer to the initialized SparseMatrix or NULL on failure
 */
SparseMatrix *sparse_matrix_init(int nrows, int ncols);

/**
 * @brief Deallocate all memory associated with the sparse matrix.
 * @param psm Pointer to the sparse matrix to be deallocated
 */
void sparse_matrix_deinit(SparseMatrix *psm);

/**
 * @brief Transpose the given sparse matrix in-place.
 * Swaps row and column indices for each non-zero entry.
 * @param psm Pointer to the sparse matrix
 */
void sparse_matrix_transpose(SparseMatrix *psm);

/**
 * @brief Add a non-zero entry to the sparse matrix.
 * @param psm Pointer to the sparse matrix
 * @param entry Non-zero value
 * @param irow Row index
 * @param icol Column index
 * @return 0 on success, -1 on failure (invalid input or matrix full)
 */
int sparse_matrix_add_entry(SparseMatrix *psm, SMType entry, int irow,
                            int icol);

/**
 * @brief Multiply a vector from the left: result = vecᵗ × matrix.
 * @param pv Input vector (dimension = matrix rows)
 * @param psm Sparse matrix
 * @param res Result vector (dimension = matrix cols)
 * @return 0 on success, -1 on error (e.g., dimension mismatch)
 */
int sparse_matrix_mult_vec_left(Vec *pv, SparseMatrix *psm, Vec *res);

/**
 * @brief Multiply a array from the left: result = memᵗ × matrix.
 * @param pmem Input memory array (dimension = matrix rows)
 * @param size Dimension of the input array
 * @param psm Sparse matrix
 * @param res Result vector (dimension = matrix cols)
 * @return 0 on success, -1 on error (e.g., dimension mismatch)
 */
int sparse_matrix_mult_mem_left(SMType *pmem, int size, SparseMatrix *psm,
                                Vec *res);

/**
 * @brief Multiply a vector from the right: result = matrix × vec.
 * @param psm Sparse matrix
 * @param pv Input vector (dimension = matrix cols)
 * @param res Result vector (dimension = matrix rows)
 * @return 0 on success, -1 on error (e.g., dimension mismatch)
 */
int sparse_matrix_mult_vec_right(SparseMatrix *psm, Vec *pv, Vec *res);

/**
 * @brief Multiply a array from the right: result = matrix × vec.
 * @param psm Sparse matrix
 * @param pmem Input memory array (dimension = matrix cols)
 * @param size Dimension of the input array
 * @param res Result vector (dimension = matrix rows)
 * @return 0 on success, -1 on error (e.g., dimension mismatch)
 */
int sparse_matrix_mult_mem_right(SparseMatrix *psm, SMType *pmem, int size,
                                 Vec *res);

/*
** Implementations
*/
Vec *vec_init(int dim) {
    Vec *pv = (Vec *)malloc(sizeof(Vec));
    if (pv == NULL)
        return NULL;
    pv->dim = dim;
    pv->entries = (SMType *)malloc(sizeof(SMType) * pv->dim);
    if (pv->entries == NULL) {
        free(pv);
        return NULL;
    }
    // initialize the zero vector
    memset(pv->entries, 0, sizeof(SMType) * dim);
    return pv;
}

Vec *vec_init_mem(int dim, SMType **pmem) {
    if (pmem == NULL || *pmem == NULL)
        return NULL;
    Vec *pv = (Vec *)malloc(sizeof(Vec));
    if (pv == NULL)
        return NULL;
    pv->dim = dim;
    pv->entries = *pmem;
    // zero out the memory region
    memset(*pmem, 0, sizeof(SMType) * dim);
    *pmem = NULL; // take the ownership
    return pv;
}

void vec_deinit(Vec *pv) {
    if (pv) {
        if (pv->entries)
            free(pv->entries);
        free(pv);
    }
}

SparseMatrix *sparse_matrix_init(int nrows, int ncols) {
    SparseMatrix *psm = (SparseMatrix *)malloc(sizeof(SparseMatrix));
    if (psm == NULL)
        return NULL;
    psm->nrows = MIN(nrows, SM_MAX_ROWS);
    psm->ncols = MIN(ncols, SM_MAX_COLS);

    // initialize the pointers with null value
    psm->prows = NULL;
    psm->pcols = NULL;
    psm->nzentries = NULL;

    psm->nnz = 0; // no entries present initially;
    int mxnnz = MIN(nrows * ncols, SM_MAX_NZERO);
    psm->nzentries = (SMType *)malloc(sizeof(SMType) * mxnnz);
    if (psm->nzentries == NULL) {
        sparse_matrix_deinit(psm);
        return NULL;
    }
    psm->prows = (int *)malloc(sizeof(int) * mxnnz);
    if (psm->prows == NULL) {
        sparse_matrix_deinit(psm);
        return NULL;
    }
    psm->pcols = (int *)malloc(sizeof(int) * mxnnz);
    if (psm->pcols == NULL) {
        sparse_matrix_deinit(psm);
        return NULL;
    }

    return psm;
}

void sparse_matrix_deinit(SparseMatrix *psm) {
    if (psm != NULL) {
        if (psm->nzentries != NULL)
            free(psm->nzentries);
        if (psm->prows != NULL)
            free(psm->prows);
        if (psm->pcols != NULL)
            free(psm->pcols);
        free(psm);
        psm = NULL;
    }
}

void sparse_matrix_transpose(SparseMatrix *psm) {
    if (psm && psm->nzentries && psm->prows && psm->pcols) {
        for (int i = 0; i < psm->nnz; i++) {
            int *row = psm->prows + i;
            int *col = psm->pcols + i;
            // swap the rows and columns
            int tmp = *row;
            *row = *col;
            *col = tmp;
        }
    }
}

int sparse_matrix_add_entry(SparseMatrix *psm, SMType entry, int irow,
                            int icol) {
    if (psm && psm->nzentries && psm->prows && psm->pcols) {
        if (irow >= 0 && irow < psm->nrows && icol >= 0 && icol < psm->ncols) {
            // For valid row and column index
            int mxnnz = MIN(psm->nrows * psm->ncols, SM_MAX_NZERO);
            if (psm->nnz < mxnnz) {
                // if space is available for the entry
                psm->nzentries[psm->nnz] = entry;
                psm->prows[psm->nnz] = irow;
                psm->pcols[psm->nnz] = icol;
                psm->nnz++;
                return 0;
            }
        }
    }
    return -1;
}

int sparse_matrix_mult_vec_left(Vec *pv, SparseMatrix *psm, Vec *res) {
    if (pv && pv->entries && res && res->entries && psm && psm->prows &&
        psm->pcols && psm->nzentries) {
        if (pv->dim != psm->nrows) {
            return -1; // dimension mismatch for input
        }
        if (res->dim != psm->ncols) {
            return -1; // dimension mismatch for output
        }
        memset(res->entries, 0, res->dim * sizeof(SMType));
        for (int i = 0; i < psm->nnz; i++) {
            int irow = psm->prows[i];
            int icol = psm->pcols[i];
            res->entries[icol] =
                SMSum(SMProd(psm->nzentries[i], pv->entries[irow]),
                      res->entries[icol]);
        }
        return 0; // success
    }
    return -1; // invalid pointer
}

int sparse_matrix_mult_mem_left(SMType *pmem, int size, SparseMatrix *psm,
                                Vec *res) {
    if (pmem && psm && psm->prows && psm->pcols && psm->nzentries && res &&
        res->entries) {
        if (size != psm->nrows) {
            return -1; // dimension mismatch for input
        }
        if (res->dim != psm->ncols) {
            return -1; // dimemsion mismatch for output
        }
        memset(res->entries, 0, res->dim * sizeof(SMType));
        for (int i = 0; i < psm->nnz; i++) {
            int irow = psm->prows[i];
            int icol = psm->pcols[i];
            res->entries[icol] = SMSum(SMProd(psm->nzentries[i], pmem[irow]),
                                       res->entries[icol]);
        }
        return 0; // success
    }
    return -1; // invalid pointer
}

int sparse_matrix_mult_vec_right(SparseMatrix *psm, Vec *pv, Vec *res) {
    if (pv && pv->entries && res && res->entries && psm && psm->prows &&
        psm->pcols && psm->nzentries) {
        if (pv->dim != psm->ncols) {
            return -1; // dimension mismatch for input
        }
        if (res->dim != psm->nrows) {
            return -1; // dimension mismatch for output
        }
        memset(res->entries, 0, res->dim * sizeof(SMType));
        for (int i = 0; i < psm->nnz; i++) {
            int irow = psm->prows[i];
            int icol = psm->pcols[i];
            res->entries[irow] =
                SMSum(SMProd(psm->nzentries[i], pv->entries[icol]),
                      res->entries[irow]);
        }
        return 0; // success
    }
    return -1; // invalid pointer
}

int sparse_matrix_mult_mem_right(SparseMatrix *psm, SMType *pmem, int size,
                                 Vec *res) {
    if (pmem && res && res->entries && psm && psm->prows && psm->pcols &&
        psm->nzentries) {
        if (size != psm->ncols) {
            return -1; // dimension mismatch for input
        }
        if (res->dim != psm->nrows) {
            return -1; // dimension mismatch for output
        }
        memset(res->entries, 0, res->dim * sizeof(SMType));
        for (int i = 0; i < psm->nnz; i++) {
            int irow = psm->prows[i];
            int icol = psm->pcols[i];
            res->entries[irow] = SMSum(SMProd(psm->nzentries[i], pmem[icol]),
                                       res->entries[irow]);
        }
        return 0; // success
    }
    return -1; // invalid pointer
}

#if defined(SM_DEBUG)
/**
 * @brief Print the vector for debugging.
 * Requires `SM_DEBUG` to be defined.
 * @param pv Pointer to the vector
 */
void vec_debug_print(Vec *pv) {
    if (pv == NULL) {
        printf("Null pointer provided.\n");
        return;
    }
    if (pv->entries == NULL) {
        printf("pv->entries is null.\n");
        return;
    }

    assert(pv->dim >= 0);
    printf("[ ");
    for (int i = 0; i < pv->dim; i++) {
#ifdef SM_DBG_PRINT_FN
        SM_DBG_PRINT_FN(pv->entries[i]);
#else
        printf(SM_PRINT_FMT " ", pv->entries[i]);
#endif
    }
    printf("]\n");
}

/**
 * @brief Print the full matrix in dense format for debugging.
 * Requires `SM_DEBUG` to be defined.
 * @param psm Pointer to the sparse matrix
 */
void sparse_matrix_debug_print(SparseMatrix *psm) {
    if (psm == NULL) {
        printf("Null pointer provided.\n");
        return;
    }
    if (psm->nzentries == NULL) {
        printf("psm->nzentries is null.\n");
    }
    if (psm->prows == NULL) {
        printf("psm->prows is null.\n");
    }
    if (psm->pcols == NULL) {
        printf("psm->pcols is null.\n");
    }

    assert(psm->nrows >= 0);
    assert(psm->ncols >= 0);

    int nrows = psm->nrows, ncols = psm->ncols;
    SMType matrix[nrows][ncols];

    memset(matrix, 0, sizeof(SMType) * nrows * ncols);

    for (int i = 0; i < psm->nnz; i++) {
        matrix[psm->prows[i]][psm->pcols[i]] = psm->nzentries[i];
    }

    for (int i = 0; i < nrows; i++) {
        for (int ii = 0; ii < ncols; ii++) {
#ifdef SM_DBG_PRINT_FN
            SM_DBG_PRINT_FN(matrix[i][ii]);
#else
            printf(SM_PRINT_FMT " ", matrix[i][ii]);
#endif
        }
        printf("\n");
    }
}

/**
 * @brief Prints detailed information about a SparseMatrix.
 *
 * This function prints the number of non-zero entries (nnz), the matrix
 * dimensions, and for each non-zero entry, its value and corresponding (row,
 * col) position.
 *
 * @param mat Pointer to the SparseMatrix structure.
 */
void sparse_matrix_debug_info(const SparseMatrix *mat) {
    if (mat == NULL) {
        printf("SparseMatrix: NULL pointer\n");
        return;
    }

    printf("=== SparseMatrix Info ===\n");
    printf("Dimensions: %d rows x %d cols\n", mat->nrows, mat->ncols);
    printf("Non-zero entries (nnz): %d\n", mat->nnz);

    for (int i = 0; i < mat->nnz; i++) {
        int row = mat->prows[i];
        int col = mat->pcols[i];
        SMType val = mat->nzentries[i];

#define COMPLEX_FLOAT
#ifdef COMPLEX_FLOAT
        printf("Entry %d: (%d, %d) = %.6f + %.6fi\n", i, row, col, val.re,
               val.im);
#else
        printf("Entry %d: (%d, %d) = %.6f\n", i, row, col, val);
#endif
    }

    printf("=========================\n");
}

#endif // SM_DEBUG
#endif // SPARSEMATRIX_H_
