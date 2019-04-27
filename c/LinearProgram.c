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

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

#include "Simplex.h"

#include "LinearProgram.h"

/* 
 * b <= Ax
 * Z  = c*x
 */


LinearProgram_t* new_lp(int rows, int cols)
{
    LinearProgram_t* ths = NULL;

    if (0 < rows && 0 < cols) {
	int i, j;
	ths = malloc(sizeof(LinearProgram_t));
	ths->minimize = false;
	ths->rows = rows;
	ths->cols = cols;
	ths->a = malloc(rows*sizeof(double*));
	double* buf = malloc(rows*cols*sizeof(double));
	for (i=0; i < rows; i++) {
	    ths->a[i] = buf+cols*i;
	}
	ths->b = malloc(rows*sizeof(double));
	ths->c = malloc(cols*sizeof(double));
	ths->e = malloc(rows*sizeof(LP_Equality_t));
	ths->x = malloc(cols*sizeof(double));
	ths->Z = 0;

	ths->row_labels = malloc(rows*sizeof(char*));
	ths->col_labels = malloc(cols*sizeof(char*));

	for (i=0; i < rows; i++) {			// rows
	    for (j=0; j < cols; j++) {			// cols
		ths->a[i][j] = 0;
	    }
	    ths->b[i]          = 0;
	    ths->e[i]          = LP_LE;
	    ths->row_labels[i] = NULL;
	}

	for (int j=0; j < cols; j++) {			// cols
	    ths->c[j]          = 0;
	    ths->x[j]          = 0;
	    ths->col_labels[j] = NULL;
	}
    }

    return ths;
}


void delete_lp(LinearProgram_t* ths)
{
    if (ths != NULL) {
	if (ths->col_labels != NULL) {
	    int i;
	    for (i=0; i < ths->cols; i++) {
		if (ths->col_labels[i] != NULL) {
		    free(ths->col_labels[i]);
		}
	    }
	    free(ths->col_labels);
	}

	if (ths->row_labels != NULL) {
	    int i;
	    for (i=0; i < ths->rows; i++) {
		if (ths->row_labels[i] != NULL) {
		    free(ths->row_labels[i]);
		}
	    }
	    free(ths->row_labels);
	}

	if (ths->x != NULL) free(ths->x);
	if (ths->e != NULL) free(ths->e);
	if (ths->c != NULL) free(ths->c);
	if (ths->b != NULL) free(ths->b);
	if (ths->a != NULL) {
	    if (ths->a[0] != NULL) free(ths->a[0]);
	    free(ths->a);
	}

	free(ths);
    }
}


size_t get_line(char** line, size_t* n, FILE* file)
{
    if (line == NULL || n == NULL || file == NULL) return -1;

    if (*line == NULL || *n <= 0) {
	*line = malloc(4096);
	*n = 4096;
    }

    char*  buf = *line;
    size_t len = 0;
    do {
	int c = fgetc(file);
	buf[len++] = c;
	buf[len] = '\0';
    } while (buf[len-1] != '\n' && len < (*n - 1));
    buf[len] = '\0';

    return len;
}


// return the number of useful bytes read, 
// ignoring comments, blank lines and white space
//
// line:	a pointer to the line buffer
// n:		a pointer to the line size
// file:	a FILE pointer
size_t read_line(char** line, int* n, FILE* file)
{
    if (line == NULL || n == NULL || file == NULL) return -1;

    char*  buf = NULL;		// temporary buffer, discarded later
    size_t bn  = 0;		// size of buf
    size_t len = 0;		// length of data read into buf

    int bidx = 0;
    int lidx = 0;

    // get the next non-comment line
    while (0 < (len=getline(&buf, &bn, file))) {
	// adjust the caller's buffer size, if necessary
	if (*line == NULL || *n < 1) {
	    // the caller provided no buffer, allocate one
	    *line = malloc(bn);
	    *n = bn;
	} else if (*line != NULL && *n < bn) {
	    // buffer was too small, make it bigger
	    *line = realloc(*line, bn);
	    *n = bn;
	}

	// remove any trailing white space
	while (0 < len && buf[len-1] <= ' ') {
	    buf[--len] = '\0';
	}

	// remove any comments
	char* comment = index(buf, '#');
	if (comment != NULL) {
	    *comment = '\0';
	    len = comment - buf;
	}

	// start at the beginning
	bidx = 0;
	lidx = 0;

	// skip over any leading white space
	while (bidx < len && buf[bidx] <= ' ') {
	    bidx += 1;
	}

	// copy the line content
	while (bidx < len) {
	    // copy the byte
	    (*line)[lidx] = buf[bidx];

	    // advance to the next byte
	    lidx += 1;
	    bidx += 1;
	}
	(*line)[lidx] = '\0';

	if (0 < lidx) {
	    // the line has content, return it
	    break;
	}

	// the line has no content, start over
    }

    free(buf);

    return lidx;
}


LinearProgram_t* read_lp(FILE* file)
{
    LinearProgram_t* r = NULL;

    if (file != NULL) {
	char* line = NULL;
	int   size = 0;
	int   i,j;

	int len = read_line(&line, &size, file);
	bool minimize = strncasecmp(line, "minimize", len) == 0;

	// read the problem dimensions (rows,cols)
	len = read_line(&line, &size, file);
	int rows = atoi(line);
	int cols = atoi(index(line, ',')+1);

	r = new_lp(rows, cols);
	r->minimize = minimize;

	// read col labels -- "a","b",...,"z"
	len = read_line(&line, &size, file);
	char* p = line;
	for (j=0; j < cols && p != NULL; j++) {
	    char* b = index(p, '"');
	    if (b == NULL) {
		// ill-formed label, no beginning '"'
		break;
	    }
	    char* e = index(b+1, '"');
	    *e = '\0';
	    r->col_labels[j] = strdup(b+1);
	    p = index(e+1, ',');	// find the next ','
	    p = (p != NULL) ? p+1 : NULL;
	}


	// read "row label", a, e, b
	for (i=0; i < rows && 0 < (len=read_line(&line, &size, file)); i++) {
	    p = line;

	    // read the row label
	    char* b = index(p, '"');
	    if (b == NULL) {
		// ill-formed label, no beginning '"'
		break;
	    }
	    char* e = index(b+1, '"');
	    *e = '\0';
	    r->row_labels[i] = strdup(b+1);
	    p = index(e+1, ',');
	    if (p == NULL) {
		// no next ','
		break;
	    }
	    p += 1;			// move past the ','

	    // read a
	    for (j=0; j < cols; j++) {
		e = index(p, ',');
		if (e == NULL) {
		    break;
		}
		*e = '\0';
		r->a[i][j] = atof(p);
		p = e+1;
	    }

	    // read e
	    e = index(p, ',');
	    if (e != NULL) {
		*e = '\0';
		if ((b=index(p, '<')) != NULL) {
		    r->e[i] = LP_LE;
		} else if ((b=index(p, '>')) != NULL) {
		    r->e[i] = LP_GE;
		} else if ((b=index(p, '=')) != NULL) {
		    r->e[i] = LP_EQ;
		} else {
		}
		p = e+1;
	    }

	    // read b
	    r->b[i] = atof(p);
	}

	// read Z and c
	len=read_line(&line, &size, file);
	p = line;

	// read the objective function label
	char* b = index(p, '"');
	if (b != NULL) {
	    char* e = index(b+1, '"');
	    *e = '\0';
	    r->obj_label = strdup(b+1);
	    p = index(e+1, ',');
	}
	p = (p != NULL) ? p+1 : NULL;

	for (j=0; j < cols && p != NULL; j++) {
	    r->c[j] = atof(p);
	    p = index(p, ',');
	    if (p == NULL) {
		// no next ','
		break;
	    }
	    p += 1;
	}
    }

    return r;
}




void write_lp(LinearProgram_t* ths, FILE* out)
{
    if (ths != NULL) {
	int i,j;

	// print the operation
	fprintf(out, "%s\n", ths->minimize ? "minimize" : "maximize");

	// print the dimensions
	fprintf(out, "%d,%d\n", ths->rows, ths->cols);

	// print the column labels
	fprintf(out, "\"%s\"", ths->col_labels[0]);
	for (j=1; j < ths->cols; j++) {
	    fprintf(out, ",\"%s\"", ths->col_labels[j]);
	}
	fprintf(out, "\n");

	// print the constraints - "label[i]", a[i][*], "<=", b[i]
	for (i=0; i < ths->rows; i++) {
	    fprintf(out, "\"%s\",", ths->row_labels[i]);
	    for (int j=0; j < ths->cols; j++) {
		fprintf(out, "%f,\n", ths->a[i][j]);
	    } 
	    fprintf(out, "%s,", (ths->e[i] == LP_LE ? "<=" : (ths->e[i] == LP_EQ ? "==" : ">=")));
	    fprintf(out, "%f\n", ths->b[i]);
	}

	// print the objective function
	fprintf(out, "\"%s\"", ths->obj_label);
	for (j=0; j < ths->cols; j++) {
	    fprintf(out, ",%f", ths->c[j]);
	}
	fprintf(out, "\n");
    }
}


void print_lp(LinearProgram_t* ths, FILE* out)
{
    if (ths != NULL) {
	int i,j;

	// print the operation
	fprintf(out, "%s\n", ths->minimize ? "minimize" : "maximize");

	// print the column labels
	fprintf(out, "           %10.10s", ths->col_labels[0]);
	for (j=1; j < ths->cols; j++) {
	    fprintf(out, ",%10.10s", ths->col_labels[j]);
	}
	fprintf(out, "\n");

	// print the constraints - "label[i]", a[i][*], "<=", b[i]
	for (i=0; i < ths->rows; i++) {
	    fprintf(out, "%10.10s", ths->row_labels[i]);
	    for (int j=0; j < ths->cols; j++) {
		fprintf(out, " %10.2f", ths->a[i][j]);
	    } 
	    fprintf(out, " %s", (ths->e[i] == LP_LE ? "<=" : (ths->e[i] == LP_EQ ? "==" : ">=")));
	    fprintf(out, " %10.2f\n", ths->b[i]);
	}

	// print the objective function
	fprintf(out, "%10.10s", ths->obj_label);
	for (j=0; j < ths->cols; j++) {
	    fprintf(out, " %10.2f", ths->c[j]);
	}
	fprintf(out, "\n");
    }
}


int main(int argc, char* argv[]) 
{
    LinearProgram_t* lp = NULL;

    bool minimize = false;

    for (int i=1; i < argc; i++) {
	if (0 <= access(argv[i], R_OK)) {
	    FILE* file = fopen(argv[i], "r");
	    lp = read_lp(file);
	    fclose(file);
	    minimize = lp->minimize;
	} else if (strcasecmp(argv[i],"-min") == 0 || strcasecmp(argv[i],"-minimize")) {
	    minimize = true;
	} else if (strcasecmp(argv[i],"-max") == 0 || strcasecmp(argv[i],"-maximize")) {
	    minimize = false;
	}
    }
    lp->minimize = minimize;

    print_lp(lp, stdout);
    fprintf(stdout, "\n");

    Simplex_t* simplex = new_simplex(lp);

    LP_State_t soln = optimize_system(simplex, lp->minimize, stdout);

    if (soln == LP_SOLUTION) {
	fprintf(stdout, "%s = %10.2f\n", lp->obj_label, simplex->Z);
	fprintf(stdout, "x = ");
	for (int i=0; i < simplex->cols; i++) {
	    fprintf(stdout, "%10.2f ", simplex->x[i]);
	}
	fprintf(stdout, "\n");
    } else if (soln == LP_NO_SOLUTION) {
	fprintf(stdout, "No Solution\n");
    } else if (soln == LP_UNBOUNDED) {
	fprintf(stdout, "Unbounded Solution\n");
    }
}
