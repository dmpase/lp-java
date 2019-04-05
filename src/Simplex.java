/*******************************************************************************
 * Copyright (c) 1988,2019 Douglas M. Pase                                     *
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

public class Simplex {
    public enum State {SOLUTION, NO_SOLUTION, UNBOUNDED};


    LinearProgram     system = null;

    public double[][] A = null;
    public double[]   B = null;
    public double[]   C = null;
    public double[]   x = null;
    public double     Z = 0;
    
    public int augmented_rows = 0;
    public int augmented_cols = 0;
    public int extras         = 0;
	
    public int   art_var_ct           = 0;
    public int[] artificial_variables = null;
    public int[] basic_variables      = null;

    public Simplex(LinearProgram s)
    {
	system = s;
    }
	

    public State optimize_system(boolean minimize)
    {
	setup_system(minimize);
	    
	System.out.println("Setup minimize="+minimize);
	System.out.printf("extras = %d\n", extras);
	print_system();
	System.out.println();

	State soln_type = simplex();

	if (soln_type == State.SOLUTION || soln_type == State.UNBOUNDED) {
	    for (int i=1; i <= system.rows; i++) {
		if (basic_variables[i-1] < 0 && 0.0 < B[i-1]) {
		    if (0 < art_var_ct) {
			soln_type = State.NO_SOLUTION;
		    } else {
			soln_type = State.UNBOUNDED;
		    }
		    break;
		}
	    }
	}

	Z = 0;
	for (int i=1; i <= system.rows; i++) {
	    int j = basic_variables[i-1];
	    if (1 <= j && j <= system.cols) {
		if (soln_type == State.NO_SOLUTION) {
		    x[j-1] = 0.0;
		} else {
		    x[j-1] = B[i-1];
		    Z += system.c[j-1]*x[j-1];
		}
	    }
	}
	    
	System.out.println("Optimized");
	print_system();
	System.out.println();
	
	return soln_type;
    }


    // LHS    RHS
    // Ax  <= b
    private void setup_system(boolean minimize)
    {
	// make sure all values of b[i] are zero or positive,
	// not convinced this step is necessary.
	for (int i=1; i <= system.rows; i++) {
	    // if b[i] is negative, reverse the signs of b[i] and a[i][*],
	    // and reverse the direction of the inequality (if <= or >=).
	    if (system.b[i-1] < 0) {
		system.b[i-1] = -system.b[i-1];
		if (system.e[i-1] == LinearProgram.Equality.LE) {
		    system.e[i-1] = LinearProgram.Equality.GE;
		} else if (system.e[i-1] == LinearProgram.Equality.GE) {
		    system.e[i-1] = LinearProgram.Equality.LE;
		}
		for (int j=1; j <= system.cols; j++) {
		    system.a[i-1][j-1] = -system.a[i-1][j-1];
		}
	    }
	}

	// are all of the inequalities of the form sum(a[i][j]*x[j]) <= b[i] ?
	// if so, leave bigM false, otherwise set it to true.
	extras = 0;
	boolean bigM = false;
	for (int i=1; i <= system.rows; i++) {
	    switch (system.e[i-1]) {
	    case LE :
		extras += 1;	// add a slack variable
		break;
	    case EQ :
		extras += 1;	// add an artificial variable
		bigM = true;
		break;
	    case GE :
		extras += 2;	// add a surplus and an artificial variable
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
	    for (int i=1; i <= system.rows; i++) {
		for (int j=1; j <= system.cols; j++) {
		    if (Math.abs(system.a[i-1][j-1]) > MaxA) {
			MaxA = Math.abs(system.a[i-1][j-1]);
		    }
		}

		if (Math.abs(system.b[i-1]) > MaxB) {
		    MaxB = Math.abs(system.b[i-1]);
		}
	    }
	    for (int j=1; j <= system.cols; j++) {
		if (Math.abs(system.c[j-1]) > MaxC) {
		    MaxC = Math.abs(system.c[j-1]);
		}
	    }
	    // M must be larger than the largest value in A, B or C
	    M = max(MaxA, max(MaxB, MaxC)) * 2.0;
	}

	// size and allocate the augmented matrixes
	augmented_rows = system.rows;
	augmented_cols = system.cols + extras;

	A = new double[augmented_rows][augmented_cols];
	B = new double[augmented_rows];
	C = new double[augmented_cols];
	x = new double[system.cols];

	// allocate the tracking matrixes
	basic_variables      = new int[system.rows];
	artificial_variables = new int[system.cols+extras];

	// initialize C, x and av
	for (int j=1; j <= system.cols; j++) {
	    C[j-1] = minimize ? system.c[j-1] : -system.c[j-1];
	    x[j-1] = 0.0;
	    artificial_variables[j-1] = j;            // normal variable
	}

	for (int j=system.cols+1; j <= augmented_cols; j++) {
	    C[j-1] = 0.0;
	}

	// set up the slack, artificial and surplus variables.
	art_var_ct = 0;
	int k = system.cols+1;
	for (int i=1; i <= system.rows; i++) {
	    // first approximation of the row contents
	    for (int j=1; j <= system.cols; j++) {
		A[i-1][j-1] = system.a[i-1][j-1];
	    }
	    for (int j=system.cols+1; j <= system.cols+extras; j++) {
		A[i-1][j-1] = 0.0;
	    }

	    // modify the row contents (aux variables) based on the constraint type
	    switch (system.e[i-1]) {
	    case LE :
		// Ax <= b 
		// set up the slack variable
		A[i-1][k-1] = 1.0;
		basic_variables[i-1] = k;
		artificial_variables[k-1] = k;		// slack variable
		k++;
		break;
	    case EQ :
		// Ax == b
		// set up the artificial variable
		A[i-1][k-1] = 1.0;
		for (int j=1; j <= system.cols; j++) {
		    C[j-1] -= A[i-1][j-1]*M;
		}
		basic_variables[i-1] = -k;
		artificial_variables[k-1] = -k;
		art_var_ct++;
		k++;
		break;
	    case GE :
		// Ax >= b
		// set up the artificial variable
		A[i-1][k-1] = -1.0;
		for (int j=1; j <= system.cols; j++) {
		    C[j-1] -= A[i-1][j-1]*M;
		}
		C[k-1] = M;	// note: this is NOT duplicated in EQ
		basic_variables[i-1] = -k;
		artificial_variables[k-1] = k;
		art_var_ct++;
		k++;

		// set up the surplus variable
		A[i-1][k-1] = 1.0;
		k++;
		break;
	    }

	    B[i-1] = system.b[i-1];
	}
    }

	
    private State simplex()
    {
	if (! system.is_valid() || augmented_cols <= augmented_rows) {
	    System.err.println("No Solution: cols="+augmented_cols+" rows="+augmented_rows);
	    return State.NO_SOLUTION;
	}
		
	// while there exists another pivot column
	    // select pivot column
	    // select pivot row
	    // update basic, a and c

	State result = State.SOLUTION;
	for (int pc=pivot_col(); pc != 0; pc=pivot_col()) {
	    int pr = pivot_row(pc);

	    System.out.println("pivot=("+pr+","+pc+")");

	    if (0 != pr) {
		if ((result=lp_update(pr, pc)) != State.SOLUTION) {
		    break;
		}
	    } else {
		// Z is unbounded
		result = State.UNBOUNDED;
		break;
	    }
		
	    print_system();
	    System.out.println();
	}
		
	return result;
    }
	
	
    private int pivot_col()
    {
	int j = 0;                              // search across the columns
	double e = 0.0;
	for (int i=1; i <= augmented_cols; i++) {
	    double f = C[i-1];
	    if (f < e) {
		j = i;
		e = C[i-1];
	    }
	}

	return j;
    }

	
    private int pivot_row(int pc)
    {
	int j = 0;                              // search down the rows
	double e = 0.0;
	for (int i=1; i <= augmented_rows; i++) {
	    double f = B[i-1];
	    double g = A[i-1][pc-1];
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

	
    private State lp_update(int pr, int pc)
    {
	if (pr == 0 || pc == 0) {
	    // Z is unbounded
	    return State.UNBOUNDED;
	}

	double pivot = A[pr-1][pc-1];

	// update pivot row
	for (int j=1; j <= augmented_cols; j++) {
	    A[pr-1][j-1] = A[pr-1][j-1]/pivot;
	}
	B[pr-1] = B[pr-1]/pivot;

				    // update c vector
	double aux = C[pc-1];
	for (int j=1; j <= augmented_cols; j++) {
	    C[j-1] = C[j-1]-(aux*A[pr-1][j-1]);
	}

				    // update matrix a
	for (int i=1; i <= augmented_rows; i++) {
	    aux = A[i-1][pc-1];
	    if (i == pr || aux == 0.0) {
		continue;
	    }

	    for (int j=1; j <= augmented_cols; j++) {
		A[i-1][j-1] = (A[i-1][j-1]-aux*A[pr-1][j-1]);
	    }
	    B[i-1] = (B[i-1]-aux*B[pr-1]);
	}

	if (0 < art_var_ct && basic_variables[pr-1] < 0 && artificial_variables[pc-1] > 0) {
	    art_var_ct--;
	} else if (0 < art_var_ct && basic_variables[pr-1] > 0 && artificial_variables[pc-1] < 0) {
	    art_var_ct++;
	}

	basic_variables[pr-1] = artificial_variables[pc-1];

	return State.SOLUTION;
    }

	
    private void print_system()
    {
	for (int j=0; j < C.length; j++) {
	    System.out.printf("%8.2f ", C[j]);
	}
	System.out.printf(" = %8.2f", Z);
	System.out.println();

	for (int i=0; i < A.length; i++) {
	    for (int j=0; j < A[i].length; j++) {
		System.out.printf("%8.2f ", A[i][j]);
	    }
	    System.out.printf(" = %8.2f", B[i]);
	    System.out.println();
	}

	for (int j=0; j < x.length; j++) {
	    System.out.printf("%8.2f ", x[j]);
	}
	System.out.printf("     = x");
	System.out.println();

	for (int j=0; j < basic_variables.length; j++) {
	    System.out.printf("%8d ", basic_variables[j]);
	}
	System.out.printf(" basic vars");
	System.out.println();

	for (int j=0; j < artificial_variables.length; j++) {
	    System.out.printf("%8d ", artificial_variables[j]);
	}
	System.out.printf(" artificial vars");
	System.out.println();
    }

    
    public static double max(double a, double b)
    {
	return (a < b) ? b : a;
    }
	
	

    // v is the basic(?) variable
    // n is the number of variables
    public boolean is_basic(int v, int n)
    {
	int j;

	for (int i=0; i < n; i++) {
	    if (basic_variables[i] < 0) {
		j = -basic_variables[i];
	    } else {
		j = basic_variables[i];
	    }

	    if (v == j) {
		return true;
	    }
	}

	return false;
    }
}
