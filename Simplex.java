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

package matrix;

public class Simplex {
	public enum State {SOLUTION, NO_SOLUTION, UNBOUNDED};


	LinearProgram     system = null;

	public double[][] A = null;
	public double[]   B = null;
	public double[]   C = null;
	public double[]   x = null;
	public double     Z = 0;
	
	int aug_rows = 0;
	int aug_cols = 0;
    int extras   = 0;
	
	int   avct  = 0;
	int[] art_vars    = null;
	int[] basic = null;

	public Simplex(LinearProgram s)
	{
		system = s;
	}
	

	public State minimize_system()
	{
	    setup_system(true);
	    
		System.out.println("Setup");
	    print_system();
		System.out.println();

	    State soln_type = simplex();
	    
		System.out.println("Simplex");
	    print_system();
		System.out.println();

	    Z = -Z;

	    if (soln_type == State.SOLUTION || soln_type == State.UNBOUNDED) {
	        for (int i=1; i <= system.rows; i++) {
	            if (basic[i-1] < 0 && 0.0 < B[i-1]) {
	                if (0 < avct) {
	                    soln_type = State.NO_SOLUTION;
	                } else {
	                    soln_type = State.UNBOUNDED;
	                }
	                break;
	            }
	        }
	    }

	    for (int i=1; i <= system.rows; i++) {
	        int j = basic[i-1];
	        if (1 <= j && j <= system.cols) {
	            if (soln_type == State.NO_SOLUTION) {
	                x[j-1] = 0.0;
	            } else {
	                x[j-1] = B[i-1];
	            }
	        }
	    }
	    
		System.out.println("Minimize");
	    print_system();
		System.out.println();
	    
	    return soln_type;
	}


	public State maximize_system()
	{
	    setup_system(false);
	    
		System.out.println("Setup");
	    print_system();
		System.out.println();

	    State soln_type = simplex();
	    
		System.out.println("Simplex");
	    print_system();
		System.out.println();
	
	    if (soln_type == State.SOLUTION || soln_type == State.UNBOUNDED) {
	        for (int i=1; i <= system.rows; i++) {
	            if (basic[i-1] < 0 && 0.0 < B[i-1]) {
	                if (0 < avct) {
	                    soln_type = State.NO_SOLUTION;
	                } else {
	                    soln_type = State.UNBOUNDED;
	                }
	                break;
	            }
	        }
	    }
	
	    for (int i=1; i <= system.rows; i++) {
	        int j = basic[i-1];
	        if (1 <= j && j <= system.cols) {
	            if (soln_type == State.NO_SOLUTION) {
	                x[j-1] = 0.0;
	            } else {
	                x[j-1] = B[i-1];
	            }
	        }
	    }
	    
		System.out.println("Maximize");
	    print_system();
		System.out.println();

	    return soln_type;
	}

	
	private void setup_system(boolean minimize)
	{
	    extras = system.rows;
	    boolean bigM = false;
	    for (int i=1; i <= system.rows; i++) {
	        switch (system.e[i-1]) {
	        case LE :
	            break;
	        case EQ :
	            bigM = true;
	            break;
	        case GE :
	            extras += 1;
	            bigM = true;
	            break;
	        }
	    }

	    double M = 0;

	    if (bigM) {
	    	double MaxA, MaxB, MaxC;

	        MaxA = MaxB = MaxC = 1.0;
	        for (int i=1; i <= system.rows; i++) {
	            for (int j=1; j <= system.cols; j++) {
	                if (system.a[i-1][j-1] > MaxA) {
	                    MaxA = system.a[i-1][j-1];
	                }
	            }

	            if (system.b[i-1] > MaxB) {
	                MaxB = system.b[i-1];
	            }
	        }
	        for (int j=1; j <= system.cols; j++) {
	            if (system.c[j-1] > MaxC) {
	                MaxC = system.c[j-1];
	            }
	        }
	        M = max(MaxA, max(MaxB,MaxC)) * 128.0;
	    }
	    
	    aug_rows = system.rows;
	    aug_cols = system.cols + extras;

	    A = new double[aug_rows][aug_cols];
	    B = new double[aug_rows];
	    C = new double[aug_cols];
	    x = new double[system.cols];
	    Z = 0.0;

	    basic = new int[system.rows];
	    art_vars    = new int[system.cols+extras];

	    for (int j=1; j <= system.cols; j++) {
	        C[j-1] = (minimize?1:-1)*system.c[j-1];
	        x[j-1] = 0.0;
	        art_vars[j-1] = j;            // normal variable
	    }

	    for (int j=system.cols+1; j <= aug_cols; j++) {
	        C[j-1] = 0.0;
	    }

	    avct = 0;
	    int k = system.cols+1;
	    for (int i=1; i <= system.rows; i++) {
	        for (int j=1; j <= system.cols; j++) {
	            A[i-1][j-1] = system.a[i-1][j-1];
	        }
	        for (int j=system.cols+1; j <= system.cols+extras; j++) {
	            A[i-1][j-1] = 0.0;
	        }

	        switch (system.e[i-1]) {
	        case LE :
	            A[i-1][k-1] = 1.0;
	            basic[i-1] = k;
	            art_vars[k-1] = k;                // slack variable
	            k++;
	            break;
	        case GE :
	            C[k-1] = M;
	            A[i-1][k-1] = -1.0;
	            art_vars[k-1] = k;                // surplus variable
	            k++;
	            break;
	        case EQ :
	            Z -= B[i-1] * M;
	            for (int j=1; j <= system.cols; j++) {
	        	    double oldC = C[j-1];
	        	    double newC = oldC - A[i-1][j-1]*M;
	                C[j-1] = newC;
	            }
	            A[i-1][k-1] = 1.0;
	            basic[i-1] = -k;
	            art_vars[k-1] = -k;               // artificial variable
	            avct++;
	            k++;
	            break;
	        }

	        B[i-1] = system.b[i-1];
	    }
	}
	
	private void print_system()
	{
		for (int i=0; i < A.length; i++) {
			System.out.printf("%8.2f = ", B[i]);
			for (int j=0; j < A[i].length; j++) {
				System.out.printf("%8.2f ", A[i][j]);
			}
			System.out.println();
		}

		System.out.printf("       C = ");
		for (int j=0; j < C.length; j++) {
			System.out.printf("%8.2f ", C[j]);
		}
		System.out.println();

		System.out.printf("       x = ");
		for (int j=0; j < x.length; j++) {
			System.out.printf("%8.2f ", x[j]);
		}
		System.out.println();

		System.out.printf("   basic = ");
		for (int j=0; j < basic.length; j++) {
			System.out.printf("%8d ", basic[j]);
		}
		System.out.println();

		System.out.printf("art_vars = ");
		for (int j=0; j < art_vars.length; j++) {
			System.out.printf("%8d ", art_vars[j]);
		}
		System.out.println();

		System.out.printf("  extras = %8d\n", extras);
		System.out.printf("       Z = %8.2f\n", Z);
	}

	
	private State simplex()
	{
		if (! system.is_valid() || aug_cols <= aug_rows) {
			System.err.println("No Solution: cols="+aug_cols+" rows="+aug_rows);
			return State.NO_SOLUTION;
		}
		
	    // while there exists another pivot column
	        // select pivot column
	        // select pivot row
	        // update basic, a and c

		State result = State.SOLUTION;
		for (int pc=pivot_col(); pc != 0; pc=pivot_col()) {
			int pr = pivot_row(pc);

	        if (0 != pr) {
	            if ((result=lp_update(pr, pc)) != State.SOLUTION) {
	                break;
	            }
	        } else {
	            // Z is unbounded
	            result = State.UNBOUNDED;
	            break;
	        }
		}
		
		return result;
	}
	
	
	private int pivot_col()
	{
	    int j = 0;                              // search across the columns
	    double e = 0.0;
	    for (int i=1; i <= aug_cols; i++) {
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
	    for (int i=1; i <= aug_rows; i++) {
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
	    for (int j=1; j <= aug_cols; j++) {
	        A[pr-1][j-1] = A[pr-1][j-1]/pivot;
	    }
	    B[pr-1] = B[pr-1]/pivot;

	                                // update c vector
	    double aux = C[pc-1];
	    for (int j=1; j <= aug_cols; j++) {
	        C[j-1] = C[j-1]-(aux*A[pr-1][j-1]);
	    }
	    Z = Z-(aux*B[pr-1]);

	                                // update matrix a
	    for (int i=1; i <= aug_rows; i++) {
	        aux = A[i-1][pc-1];
	        if (i == pr || aux == 0.0) {
	            continue;
	        }

	        for (int j=1; j <= aug_cols; j++) {
	            A[i-1][j-1] = (A[i-1][j-1]-aux*A[pr-1][j-1]);
	        }
	        B[i-1] = (B[i-1]-aux*B[pr-1]);
	    }

	    if (0 < avct && basic[pr-1] < 0 && art_vars[pc-1] > 0) {
	        avct--;
	    } else if (0 < avct && basic[pr-1] > 0 && art_vars[pc-1] < 0) {
	        avct++;
	    }

	    basic[pr-1] = art_vars[pc-1];

	    return State.SOLUTION;
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
			if (basic[i] < 0) {
				j = -basic[i];
			} else {
				j = basic[i];
			}

			if (v == j) {
				return true;
			}
		}

		return false;
	}
}
