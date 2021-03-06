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

import java.io.*;

public class LinearProgram {
	public enum Equality {LE, EQ, GE};
	public enum State {SOLUTION, NO_SOLUTION, UNBOUNDED};

	public boolean minimize = false;

	public int rows = 0;
	public int cols = 0;
	public double[][] a = null;
	public double[]   b = null;
	public double[]   c = null;
	public Equality[] e = null;
	public double[]   x = null;
	public double     Z = 0;

	public String[]   row_labels = null;
	public String[]   col_labels = null;
	public String     obj_label  = null;
	
	/* 
	 * b <= Ax
	 * Z  = c*x
	 */

	public LinearProgram(int rows, int cols)
	{
		this.rows  = rows;
		this.cols  = cols;
		a          = new double[rows][cols];
		b          = new double[rows];
		c          = new double[cols];
		e          = new Equality[rows];
		x          = new double[cols];
		Z          = 0;
		row_labels = new String[rows];
		col_labels = new String[cols];
		
		for (int i=0; i < rows; i++) {			// rows
			for (int j=0; j < cols; j++) {		// cols
				a[i][j] = 0;
			}
			b[i]          = 0;
			e[i]          = Equality.LE;
			row_labels[i] = null;
		}

		for (int j=0; j < cols; j++) {			// cols
			c[j]          = 0;
			x[j]          = 0;
			col_labels[j] = null;
		}
	}

	
	public LinearProgram(double[][] a0, double[] b0, double[] c0, Equality[] e0, String[] rlabels, String[] clabels)
	{
		if (a0 == null || b0 == null || c0 == null || e0 == null || rlabels == null || clabels == null) {
			System.err.println("LinearProgram: null arguments are not allowed.");
		}

		rows = a0.length;
		if (b0.length != rows || e0.length != rows || rlabels.length != rows) {
			System.err.println("LinearProgram: row lengths are inconsistent.");
		}

		cols = c0.length;
		if (clabels.length != cols) {
			System.err.println("LinearProgram: row lengths are inconsistent.");
		}

		a          = a0;
		b          = b0;
		c          = c0;
		e          = e0;
		x          = new double[cols];
		Z          = 0;
		row_labels = rlabels;
		col_labels = clabels;
		
		for (int j=0; j < cols; j++) {			// cols
			x[j]          = 0;
		}
	}

	
	public static double set(double[][] a, int i, int j, double elt)
	{
		double result = 0;
		if (1 <= i && i <= a.length && 1 <= j && j <= a[j-1].length) {
			result = a[i-1][j-1];
			a[i-1][j-1] = elt;
		}

		return result;
	}
	
	
	public static double get(double[][] a, int i, int j)
	{
		double result = 0;
		if (1 <= i && i <= a.length && 1 <= j && j <= a[j-1].length) {
			result = a[i-1][j-1];
		}

		return result;
	}
	
	
	public static double set(double[] a, int i, double elt)
	{
		double result = 0;
		if (1 <= i && i <= a.length) {
			result = a[i-1];
			a[i-1] = elt;
		}

		return result;
	}

	
	public static double get(double[] a, int i)
	{
		double result = 0;
		if (1 <= i && i <= a.length) {
			result = a[i-1];
		}

		return result;
	}
	
	
	public Equality set(Equality[] e, int i, Equality elt)
	{
		Equality result = Equality.LE;
		if (1 <= i && i <= rows) {
			result = e[i-1];
			e[i-1] = elt;
		}

		return result;
	}
	
	
	public Equality set(Equality[] e, int i, byte elt)
	{
		Equality result = Equality.LE;
		if (1 <= i && i <= e.length) {
			result = e[i-1];
			if (elt == '<') {
				e[i-1] = Equality.LE;
			} else if (elt == '=') {
				e[i-1] = Equality.EQ;
			} else if (elt == '>') {
				e[i-1] = Equality.GE;
			}
		}

		return result;
	}
	
	
	public Equality get(Equality[] e, int i)
	{
		Equality result = Equality.LE;
		if (1 <= i && i <= rows) {
			result = e[i-1];
		}

		return result;
	}
	
	
	public double set_a(int i, int j, double elt)
	{
		double result = 0;
		if (1 <= i && i <= rows && 1 <= j && j <= cols) {
			result = a[i-1][j-1];
			a[i-1][j-1] = elt;
		}

		return result;
	}
	
	
	public double get_a(int i, int j)
	{
		double result = 0;
		if (1 <= i && i <= rows && 1 <= j && j <= cols) {
			result = a[i-1][j-1];
		}

		return result;
	}
	
	
	public double set_b(int i, double elt)
	{
		double result = 0;
		if (1 <= i && i <= rows) {
			result = b[i-1];
			b[i-1] = elt;
		}

		return result;
	}
	
	
	public double get_b(int i)
	{
		double result = 0;
		if (1 <= i && i <= rows) {
			result = b[i-1];
		}

		return result;
	}
	
	
	public double set_c(int j, double elt)
	{
		double result = 0;
		if (1 <= j && j <= cols) {
			result = c[j-1];
			c[j-1] = elt;
		}

		return result;
	}
	
	
	public double get_c(int j)
	{
		double result = 0;
		if (1 <= j && j <= cols) {
			result = c[j-1];
		}

		return result;
	}
	
	
	public Equality set_e(int i, Equality elt)
	{
		Equality result = Equality.LE;
		if (1 <= i && i <= rows) {
			result = e[i-1];
			e[i-1] = elt;
		}

		return result;
	}
	
	
	public Equality get_e(int i)
	{
		Equality result = Equality.LE;
		if (1 <= i && i <= rows) {
			result = e[i-1];
		}

		return result;
	}
	
	
	public void add_row(boolean before, int row, double[] a_r, double b_r, Equality e_r, String label)
	{
		if (row < 1 || rows < row) {
			return;
		}

		double[][] new_a = new double[rows+1][cols];
		double[]   new_b = new double[rows+1];
		Equality[] new_e = new Equality[rows+1];
		String[]   new_l = new String[rows+1];
		
		if (before) {
			// insert the row before row "row"

			// copy everything up to, but not including, row
			for (int i=1; i < row; i++) {
				for (int j=1; j <= cols; j++) {
					new_a[i-1][j-1] = a[i-1][j-1];
				}
				new_b[i-1] = b[i-1];
				new_e[i-1] = e[i-1];
				new_l[i-1] = row_labels[i-1];
			}
			
			// copy the new row
			for (int j=1; j <= cols; j++) {
				new_a[row-1][j-1] = a_r[j-1];
			}
			new_b[row-1] = b_r;
			new_e[row-1] = e_r;
			new_l[row-1] = label;
			
			// copy everything from row on
			for (int i=row; i <= rows; i++) {
				for (int j=1; j <= cols; j++) {
					new_a[i][j-1] = a[i][j-1];
				}
				new_b[i] = b[i-1];
				new_e[i] = e[i-1];
				new_l[i] = row_labels[i-1];
			}
		} else {
			// insert the row after row "row"

			// copy everything up to, and including, row
			for (int i=1; i <= row; i++) {
				for (int j=1; j <= cols; j++) {
					new_a[i-1][j-1] = a[i-1][j-1];
				}
				new_b[i-1] = b[i-1];
				new_e[i-1] = e[i-1];
				new_l[i-1] = row_labels[i-1];
			}
			
			// copy the new row
			for (int j=1; j <= cols; j++) {
				new_a[row][j-1] = a_r[j-1];
			}
			new_b[row] = b_r;
			new_e[row] = e_r;
			new_l[row] = label;
			
			// copy everything from row+1 on
			for (int i=row+1; i <= rows; i++) {
				for (int j=1; j <= cols; j++) {
					new_a[i][j-1] = a[i][j-1];
				}
				new_b[i] = b[i-1];
				new_e[i] = e[i-1];
				new_l[i] = label;
			}
		}
		
		a = new_a;
		b = new_b;
		e = new_e;
		row_labels = new_l;
		rows += 1;
	}
	
	
	public void add_col(boolean before, int col, double[] a_c, double c_c, String label)
	{
		if (col < 1 || cols < col) {
			return;
		}

		double[][] new_a = new double[rows][cols+1];
		double[]   new_c = new double[cols+1];
		String[]   new_l = new String[cols+1];
		
		if (before) {
			// insert the column before column "col"
			
			// copy everything up to, but not including, col
			for (int j=1; j < col; j++) {
				for (int i=1; i <= rows; i++) {
					new_a[i-1][j-1] = a[i-1][j-1];
				}
				new_c[j-1] = c[j-1];
				new_l[j-1] = col_labels[j-1];
			}
			
			// copy the new column
			for (int i=1; i <= rows; i++) {
				new_a[i-1][col-1] = a_c[i-1];
			}
			new_c[col-1] = c_c;
			new_l[col-1] = label;
			
			// copy everything from col on
			for (int j=col; j <= cols; j++) {
				for (int i=1; i <= rows; i++) {
					new_a[i-1][j] = a[i-1][j-1];
				}
				new_c[j] = c[j-1];
				new_l[j] = col_labels[j-1];
			}
		} else {
			// insert the column after column "col"
			
			// copy everything up to, and including, col
			for (int j=1; j <= col; j++) {
				for (int i=1; i <= rows; i++) {
					new_a[i-1][j-1] = a[i-1][j-1];
				}
				new_c[j-1] = c[j-1];
				new_l[j-1] = col_labels[j-1];
			}
			
			// copy the new column
			for (int i=1; i <= rows; i++) {
				new_a[i-1][col] = a_c[i-1];
			}
			new_c[col] = c_c;
			new_l[col] = label;
			
			// copy everything from col on
			for (int j=col+1; j <= cols; j++) {
				for (int i=1; i <= rows; i++) {
					new_a[i-1][j] = a[i-1][j-1];
				}
				new_c[j] = c[j-1];
				new_l[j] = col_labels[j-1];
			}
		}
		
		a = new_a;
		c = new_c;
		col_labels = new_l;
		cols += 1;
	}
	
	
	public LinearProgram clone()
	{
		LinearProgram copy = new LinearProgram(rows, cols);
 
		for (int i=0; i < rows; i++) {			// rows
			for (int j=0; j < cols; j++) {		// cols
				copy.a[i][j] = a[i][j];
			}
			copy.b[i]          = b[i];
			copy.e[i]          = e[i];
			copy.row_labels[i] = row_labels[i];
		}

		for (int j=0; j < cols; j++) {			// cols
			copy.c[j]          = c[j];
			copy.col_labels[j] = col_labels[j];
		}

		return copy;
	}
	
	
	public boolean is_valid()
	{
		return a != null && b != null && c != null && e != null;
	}
	
	public static String read_line(RandomAccessFile raf) throws Exception
	{
		String line = null;
		
    	// get the next non-comment line
    	while ((line=raf.readLine()) != null) {
    		// throw out leading and trailing white space
    		line = line.trim();
    		
    		// if the line is empty, throw it out
    		if (line.length() == 0) continue;
    		
    		// does the line contain a comment?
    		int c = line.indexOf('#');
    		if (c == 0) {
    			// '#' starts the line, so throw the whole line out
    			continue;
    		} else if (0 < c) {
    			// the line contains a comment, so throw out just the comment
        		line = line.substring(0, c);
    		}
    		// the line has some content, so use it
    		break;
    	}
    	
		return line;
	}
	
	
	public static LinearProgram read(File file)
	{
		LinearProgram r = null;
		
		if (file != null && file.exists()) {
		    try {
		    	RandomAccessFile raf = new RandomAccessFile(file, "r");

		    	String line = read_line(raf);
		    	boolean minimize = line.equalsIgnoreCase("minimize");

		    	// read the problem dimensions (rows,cols)
		    	line = read_line(raf);
		    	String[] size = line.split(",");
		    	int rows = Integer.parseInt(size[0]);
		    	int cols = Integer.parseInt(size[1]);
		    	
		    	r = new LinearProgram(rows, cols);
		    	r.minimize = minimize;

		    	// read col labels
		    	line = read_line(raf);
		    	r.col_labels = line.split(",");
		    	for (int j=0; j < cols && j < r.col_labels.length; j++) {
		    		r.col_labels[j] = r.col_labels[j].trim();
		    		r.col_labels[j] = r.col_labels[j].substring(1, r.col_labels[j].length()-1);
		    	}

		    	// read a, e, b
		    	for (int i=0; i < rows && (line=read_line(raf)) != null; i++) {
		    		String[] laeb = line.split(",");
		    		laeb[0] = laeb[0].trim();
		    		r.row_labels[i] = laeb[0].substring(1,laeb[0].length()-1);
		    		for (int j=0; j < cols && j < laeb.length; j++) {
			    		laeb[j+1] = laeb[j+1].trim();
		    			r.a[i][j] = Double.parseDouble(laeb[j+1]);
		    		}
		    		laeb[cols+1] = laeb[cols+1].trim();
		    		char eq = laeb[cols+1].charAt(0);
		    		r.e[i] = (eq == '<') ? Equality.LE : (eq == '=') ? Equality.EQ : Equality.GE;
		    		laeb[cols+2] = laeb[cols+2].trim();
		    		r.b[i] = Double.parseDouble(laeb[cols+2]);
		    	}

		    	// read Z and c
		    	line = read_line(raf);
		    	String[] Zs = line.split(",");
		    	r.obj_label = Zs[0].trim();
	    		r.obj_label = r.obj_label.substring(1, r.obj_label.length()-1);
		    	for (int j=0; j < cols && (j+1) < Zs.length; j++) {
		    		r.c[j] = Double.parseDouble(Zs[j+1]);
		    	}

		    	raf.close();
		    } catch(Exception exception) {
		    	System.out.println(exception);
		    }
		}
		
		return r;
	}


	public void write(PrintStream out)
	{
		System.out.println(minimize?"minimize":"maximize");
		out.println(rows+","+cols);

		out.print("\""+col_labels[0]+"\"");
		for (int j=1; j < cols; j++) {
			out.print(",\""+col_labels[j]+"\"");
		}
		out.println();

		for (int i=0; i < rows; i++) {
			out.print("\""+row_labels[i]+"\",");
			out.print(b[i]+"," + ((e[i] == Equality.LE) ? "<=" : (e[i] == Equality.EQ) ? "==" : ">="));
			for (int j=0; j < cols; j++) {
				out.print(","+a[i][j]);
			}
			out.println();
		}

		out.print(obj_label+","+c[0]);
		for (int j=1; j < cols; j++) {
			out.print(","+c[j]);
		}
		out.println();
	}


	public void write(File file)
	{
		try {
			PrintStream out = new PrintStream(file);
			
			write(out);
			
			out.close();
		} catch (Exception exception) {
			System.err.println(exception);
		}
	}
	
	private void print_system()
	{
		System.out.println(minimize?"minimize":"maximize");
		System.out.printf("           ");
		for (int j=0; j < col_labels.length; j++) {
			System.out.printf("%10.10s ", col_labels[j]);
		}
		System.out.println();

		for (int i=0; i < a.length; i++) {
			System.out.printf("%10.10s ", row_labels[i]);
			for (int j=0; j < a[i].length; j++) {
				System.out.printf("%10.2f ", a[i][j]);
			}
			System.out.print(((e[i] == Equality.LE) ? "<= " : ((e[i] == Equality.EQ) ? "== " : ">= ")));
			System.out.printf("%10.2f", b[i]);
			System.out.println();
		}

		System.out.printf("           ");
		for (int j=0; j < c.length; j++) {
			System.out.printf("%10.2f ", c[j]);
		}
		System.out.printf("== %10s", obj_label);
		System.out.println();
	}

	
	public static void main(String[] args) 
	{
		LinearProgram lp = read(new File("windor.csv"));
		
		boolean minimize = false;
		for (int i=0; i < args.length; i++) {
			if ((new File(args[i]).exists())) {
				lp = read(new File(args[i]));
				minimize = lp.minimize;
			} else if (args[i].equalsIgnoreCase("-min") || args[i].equalsIgnoreCase("-minimize")) {
				minimize = true;
			} else if (args[i].equalsIgnoreCase("-max") || args[i].equalsIgnoreCase("-maximize")) {
				minimize = false;
			}
		}
		lp.minimize = minimize;
		
		System.out.println(System.getProperty("user.dir"));

		lp.print_system();
		System.out.println();
		
		Simplex simplex = new Simplex(lp);
		Simplex.State soln = simplex.optimize_system(lp.minimize);
		
		if (soln == Simplex.State.SOLUTION) {
			System.out.println(lp.obj_label+" = "+simplex.Z);
			System.out.print("x = ");
			for (int i=0; i < simplex.x.length; i++) {
				System.out.print(simplex.x[i]+" ");
			}
			System.out.println();
		} else if (soln == Simplex.State.NO_SOLUTION) {
			System.out.println("No Solution");
		} else if (soln == Simplex.State.UNBOUNDED) {
			System.out.println("Unbounded Solution");
		}
	}
}
