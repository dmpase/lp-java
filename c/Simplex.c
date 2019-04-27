/*******************************************************************************
 * Copyright (c) 1988, 2019 Douglas M. Pase                                    *
 * All rights reserved.                                                        *
 * Redistribution and use in source and binary forms, with or without          *
 * modification, are permitted provided that the following conditions          *
 * are met:                                                                    *
 * o       Redistributions of source code must retain the above copyright      *
 *         notice, this list of conditions and the following disclaimer.       *
 * o       Redistributions in binary form must reproduce the above copyright   *
 *         notice, this list of conditions and the following disclaimer in     *
 *         the documentation and/or other materials provided with the          *
 *         distribution.                                                       *
 * o       Neither the name of the copyright holder nor the names of its       *
 *         contributors may be used to endorse or promote products derived     *
 *         from this software without specific prior written permission.       *
 *                                                                             *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" *
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   *
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         *
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        *
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    *
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF      *
 * THE POSSIBILITY OF SUCH DAMAGE.                                             *
 *******************************************************************************/

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "LinearProgram.h"

#include "Simplex.h"

Simplex_t* new_simplex(LinearProgram_t* s)
{
    Simplex_t* ths = NULL;

    ths = malloc(sizeof *ths);
    ths->system = s;

    ths->rows = s->rows;
    ths->cols = s->cols;

    ths->A = NULL;
    ths->B = NULL;
    ths->C = NULL;
    ths->x = NULL;
    ths->Z = 0;

    ths->augmented_rows = 0;
    ths->augmented_cols = 0;
    ths->extras         = 0;

    ths->art_var_ct           = 0;
    ths->artificial_variables = NULL;
    ths->basic_variables      = NULL;

    return ths;
}


void delete_simplex(Simplex_t* ths)
{
    if (ths != NULL) {
	if (ths->A != NULL) {
	    if (ths->A[0] != NULL) free(ths->A[0]);
	    free(ths->A);
	}
	if (ths->B != NULL) free(ths->B);
	if (ths->C != NULL) free(ths->C);
	if (ths->x != NULL) free(ths->x);

	if (ths->artificial_variables != NULL) free(ths->artificial_variables);
	if (ths->basic_variables      != NULL) free(ths->basic_variables);

	free(ths);
    }
}


LP_State_t optimize_system(Simplex_t* ths, bool minimize, FILE* out)
{
    setup_system(ths, minimize);

    fprintf(out, "Setup minimize=%s\n",minimize?"true":"false");
    fprintf(out, "extras = %d\n", ths->extras);
    print_augmented_system(ths, stdout);
    fprintf(out, "\n");

    LP_State_t soln_type = simplex(ths);

    int i;
    if (soln_type == LP_SOLUTION || soln_type == LP_UNBOUNDED) {
	for (i=1; i <= ths->system->rows; i++) {
	    if (ths->basic_variables[i-1] < 0 && 0.0 < ths->B[i-1]) {
		if (0 < ths->art_var_ct) {
		    soln_type = LP_NO_SOLUTION;
		} else {
		    soln_type = LP_UNBOUNDED;
		}
		break;
	    }
	}
    }

    ths->Z = 0;
    for (i=1; i <= ths->system->rows; i++) {
	int j = ths->basic_variables[i-1];
	if (1 <= j && j <= ths->system->cols) {
	    if (soln_type == LP_NO_SOLUTION) {
		ths->x[j-1] = 0.0;
	    } else {
		ths->x[j-1] = ths->B[i-1];
		ths->Z += ths->system->c[j-1] * ths->x[j-1];
	    }
	}
    }

    fprintf(out, "Optimized\n");
    print_augmented_system(ths, stdout);
    fprintf(out, "\n");

    return soln_type;
}


// LHS    RHS
// Ax  <= b
void setup_system(Simplex_t* ths, bool minimize)
{
    int i,j;

    // make sure all values of b[i] are zero or positive,
    // not convinced this step is necessary.
    for (i=1; i <= ths->system->rows; i++) {
	// if b[i] is negative, reverse the signs of b[i] and a[i][*],
	// and reverse the direction of the inequality (if <= or >=).
	if (ths->system->b[i-1] < 0) {
	    ths->system->b[i-1] = - ths->system->b[i-1];
	    if (ths->system->e[i-1] == LP_LE) {
		ths->system->e[i-1] = LP_GE;
	    } else if (ths->system->e[i-1] == LP_GE) {
		ths->system->e[i-1] = LP_LE;
	    }
	    for (j=1; j <= ths->system->cols; j++) {
		ths->system->a[i-1][j-1] = - ths->system->a[i-1][j-1];
	    }
	}
    }

    // are all of the inequalities of the form sum(a[i][j]*x[j]) <= b[i] ?
    // if so, leave bigM false, otherwise set it to true.
    ths->extras = 0;
    bool bigM = false;
    for (i=1; i <= ths->system->rows; i++) {
	switch (ths->system->e[i-1]) {
	case LP_LE :
	    ths->extras += 1;	// add a slack variable
	    break;
	case LP_EQ :
	    ths->extras += 1;	// add an artificial variable
	    bigM = true;
	    break;
	case LP_GE :
	    ths->extras += 2;	// add a surplus and an artificial variable
	    bigM = true;
	    break;
	}
    }

    // find a value for M (used only if bigM is true)
    double M = 0;
    if (bigM) {
	double MaxA, MaxB, MaxC;

	// using abs of the value; original formulation used only the value
	MaxA = MaxB = MaxC = 1.0;
	for (i=1; i <= ths->system->rows; i++) {
	    for (j=1; j <= ths->system->cols; j++) {
		if (fabs(ths->system->a[i-1][j-1]) > MaxA) {
		    MaxA = fabs(ths->system->a[i-1][j-1]);
		}
	    }

	    if (fabs(ths->system->b[i-1]) > MaxB) {
		MaxB = fabs(ths->system->b[i-1]);
	    }
	}
	for (j=1; j <= ths->system->cols; j++) {
	    if (fabs(ths->system->c[j-1]) > MaxC) {
		MaxC = fabs(ths->system->c[j-1]);
	    }
	}
	// M must be larger than the largest value in A, B or C
	M = max(MaxA, max(MaxB, MaxC)) * 2.0;
    }

    // size and allocate the augmented matrixes
    ths->augmented_rows = ths->system->rows;
    ths->augmented_cols = ths->system->cols + ths->extras;

    ths->A = malloc(ths->augmented_rows * sizeof(double*));
    double* buf = malloc(ths->augmented_rows * ths->augmented_cols * sizeof(double));
    for (i=0; i < ths->augmented_rows; i++) {
	ths->A[i] = buf + ths->augmented_cols * i;
    }

    ths->B = malloc(sizeof(double) * ths->augmented_rows);
    ths->C = malloc(sizeof(double) * ths->augmented_cols);
    ths->x = malloc(sizeof(double) * ths->system->cols);

    // allocate the tracking matrixes
    ths->basic_variables      = malloc(sizeof(int) * ths->system->rows);
    ths->artificial_variables = malloc(sizeof(int) * (ths->system->cols + ths->extras));

    // initialize C, x and av
    for (j=1; j <= ths->system->cols; j++) {
	ths->C[j-1] = minimize ? ths->system->c[j-1] : - ths->system->c[j-1];
	ths->x[j-1] = 0.0;
	ths->artificial_variables[j-1] = j;            // normal variable
    }

    for (j=ths->system->cols+1; j <= ths->augmented_cols; j++) {
	ths->C[j-1] = 0.0;
    }

    // set up the slack, artificial and surplus variables.
    ths->art_var_ct = 0;
    int k = ths->system->cols + 1;
    for (i=1; i <= ths->system->rows; i++) {

	// first approximation of the row contents
	for (j=1; j <= ths->system->cols; j++) {
	    ths->A[i-1][j-1] = ths->system->a[i-1][j-1];
	}

	for (j=ths->system->cols + 1; j <= ths->system->cols + ths->extras; j++) {
	    ths->A[i-1][j-1] = 0.0;
	}

	// modify the row contents (aux variables) based on the constraint type
	switch (ths->system->e[i-1]) {
	case LP_LE :
	    // Ax <= b 
	    // set up the slack variable
	    ths->A[i-1][k-1] = 1.0;
	    ths->basic_variables[i-1] = k;
	    ths->artificial_variables[k-1] = k;		// slack variable
	    k++;
	    break;
	case LP_EQ :
	    // Ax == b
	    // set up the artificial variable
	    ths->A[i-1][k-1] = 1.0;
	    for (int j=1; j <= ths->system->cols; j++) {
		ths->C[j-1] -= ths->A[i-1][j-1] * M;
	    }
	    ths->basic_variables[i-1] = -k;
	    ths->artificial_variables[k-1] = -k;
	    ths->art_var_ct++;
	    k++;
	    break;
	case LP_GE :
	    // Ax >= b
	    // set up the artificial variable
	    ths->A[i-1][k-1] = -1.0;
	    for (int j=1; j <= ths->system->cols; j++) {
		ths->C[j-1] -= ths->A[i-1][j-1] * M;
	    }
	    ths->C[k-1] = M;	// note: this is NOT duplicated in EQ
	    ths->basic_variables[i-1] = -k;
	    ths->artificial_variables[k-1] = k;
	    ths->art_var_ct++;
	    k++;

	    // set up the surplus variable
	    ths->A[i-1][k-1] = 1.0;
	    k++;
	    break;
	}

	ths->B[i-1] = ths->system->b[i-1];
    }
}

    
LP_State_t simplex(Simplex_t* ths)
{
    LP_State_t result = LP_SOLUTION;

    if (ths->augmented_cols <= ths->augmented_rows) {
	return LP_NO_SOLUTION;
    }
	    
    // while there exists another pivot column
	// select pivot column
	// select pivot row
	// update basic, a and c

    int pc;
    for (pc=pivot_col(ths); pc != 0; pc=pivot_col(ths)) {
	int pr = pivot_row(ths,pc);

	fprintf(stdout, "pivot=(%d,%d)\n", pr, pc);

	if (0 != pr) {
	    if ((result=lp_update(ths, pr, pc)) != LP_SOLUTION) {
		break;
	    }
	} else {
	    // Z is unbounded
	    result = LP_UNBOUNDED;
	    break;
	}

	print_augmented_system(ths, stdout);
	fprintf(stdout, "\n");
    }

    return result;
}


int pivot_col(Simplex_t* ths)
{
    int j = 0;                              // search across the columns
    double e = 0.0;
    for (int i=1; i <= ths->augmented_cols; i++) {
	double f = ths->C[i-1];
	if (f < e) {
	    j = i;
	    e = ths->C[i-1];
	}
    }

    return j;
}


int pivot_row(Simplex_t* ths, int pc)
{
    int j = 0;                              // search down the rows
    double e = 0.0;
    for (int i=1; i <= ths->augmented_rows; i++) {
	double f = ths->B[i-1];
	double g = ths->A[i-1][pc-1];
	if (g > 0.0) {
	    double h = f/g;
	    if (h < e || j == 0) {
		j = i;
		e = h;
	    }
	}
    }

    return j;
}


LP_State_t lp_update(Simplex_t* ths, int pr, int pc)
{
    int i,j;

    if (pr == 0 || pc == 0) {
	// Z is unbounded
	return LP_UNBOUNDED;
    }

    double pivot = ths->A[pr-1][pc-1];

    // update pivot row
    for (j=1; j <= ths->augmented_cols; j++) {
	ths->A[pr-1][j-1] = ths->A[pr-1][j-1] / pivot;
    }
    ths->B[pr-1] = ths->B[pr-1]/pivot;

				// update c vector
    double aux = ths->C[pc-1];
    for (j=1; j <= ths->augmented_cols; j++) {
	ths->C[j-1] = ths->C[j-1]-(aux * ths->A[pr-1][j-1]);
    }

				// update matrix a
    for (i=1; i <= ths->augmented_rows; i++) {
	aux = ths->A[i-1][pc-1];
	if (i == pr || aux == 0.0) {
	    continue;
	}

	for (j=1; j <= ths->augmented_cols; j++) {
	    ths->A[i-1][j-1] = (ths->A[i-1][j-1]-aux * ths->A[pr-1][j-1]);
	}
	ths->B[i-1] = (ths->B[i-1]-aux * ths->B[pr-1]);
    }

    if (0 < ths->art_var_ct && ths->basic_variables[pr-1] < 0 && ths->artificial_variables[pc-1] > 0) {
	ths->art_var_ct--;
    } else if (0 < ths->art_var_ct && ths->basic_variables[pr-1] > 0 && ths->artificial_variables[pc-1] < 0) {
	ths->art_var_ct++;
    }

    ths->basic_variables[pr-1] = ths->artificial_variables[pc-1];

    return LP_SOLUTION;
}


#if defined(UNDEFINED)
#endif


void print_augmented_system(Simplex_t* ths, FILE* out)
{
    for (int j=0; j < ths->augmented_cols; j++) {
	fprintf(out, "%8.2f ", ths->C[j]);
    }
    fprintf(out, " = %8.2f\n", ths->Z);

    for (int i=0; i < ths->augmented_rows; i++) {
	for (int j=0; j < ths->augmented_cols; j++) {
	    fprintf(out, "%8.2f ", ths->A[i][j]);
	}
	fprintf(out, " = %8.2f\n", ths->B[i]);
    }

    for (int j=0; j < ths->system->cols; j++) {
	fprintf(out, "%8.2f ", ths->x[j]);
    }
    fprintf(out, "     = x\n");

    for (int j=0; j < ths->augmented_rows; j++) {
	fprintf(out, "%8d ", ths->basic_variables[j]);
    }
    fprintf(out, " basic vars\n");

    for (int j=0; j < ths->augmented_cols; j++) {
	fprintf(out, "%8d ", ths->artificial_variables[j]);
    }
    fprintf(out, " artificial vars\n");
}


double max(double a, double b)
{
    return (a < b) ? b : a;
}
