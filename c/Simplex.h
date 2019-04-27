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

#if !defined(Simplex_h)
#define Simplex_h

#include "LinearProgram.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>


typedef struct Simplex_s {
    LinearProgram_t* system;

    // original system size
    int rows;
    int cols;

    // augmented matrices
    double** A;		// a[rows][cols]
    double*  B;		// b[rows]
    double*  C;		// c[cols]
    double*  x;		// x[cols]
    double   Z;
    
    // augmented system size
    int augmented_rows;
    int augmented_cols;
    int extras;
	
    // artificial and basic variables
    int  art_var_ct;
    int* artificial_variables;
    int* basic_variables;

} Simplex_t;


extern Simplex_t* new_simplex(LinearProgram_t* s);
extern void delete_simplex(Simplex_t* ths);
extern LP_State_t optimize_system(Simplex_t* ths, bool minimize, FILE* out);
extern void setup_system(Simplex_t* ths, bool minimize);
extern LP_State_t simplex(Simplex_t* ths);
extern int pivot_col(Simplex_t* ths);
extern int pivot_row(Simplex_t* ths, int pc);
extern LP_State_t lp_update(Simplex_t* ths, int pr, int pc);
extern void print_augmented_system(Simplex_t* ths, FILE* out);
extern double max(double a, double b);

#endif
