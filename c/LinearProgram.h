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

#if !defined(LinearProgram_h)
#define LinearProgram_h

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#define DEBUG(s) {fprintf(stdout, "%s: %d\n", s, __LINE__); fflush(stdout);}

typedef enum {LP_LE, LP_EQ, LP_GE} LP_Equality_t;
typedef enum {LP_SOLUTION, LP_NO_SOLUTION, LP_UNBOUNDED} LP_State_t; 
typedef struct LinearProgram_s {
    bool minimize;

    int rows;
    int cols;

    double**       a;		// a[rows][cols]
    double*        b;		// b[rows]
    double*        c;		// c[cols]
    LP_Equality_t* e;		// e[rows]
    double*        x;		// x[cols]
    double         Z;

    char** row_labels;		// row_labels[rows]
    char** col_labels;		// col_labels[cols]
    char*  obj_label;
} LinearProgram_t;

extern LinearProgram_t* new_lp(int rows, int cols);
extern void delete_lp(LinearProgram_t* ths);
extern size_t read_line(char** line, int* n, FILE* file);
extern LinearProgram_t* read_lp(FILE* file);
extern void write_lp(LinearProgram_t* ths, FILE* out);
extern void print_lp(LinearProgram_t* ths, FILE* out);

#endif
