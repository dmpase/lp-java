# lp-java
lp is a linear programming system, originally written in 1988 while the author (Douglas Pase) was attending the 
Oregon Graduate Center. The original code was written in the C programming language using the "curses" CRT (e.g., VT-100) 
terminal control package. The solver uses the "Big M" version of the Simplex method, as described in chapter 2 of 
"Operations Research" by Hillier and Lieberman, 2nd Ed., Holden Day, Copyright 1968 and 1974.

lp-java is a spin-off from that package, converted to Java, which eliminates the session editing capabilities that required 
terminal control. Instead, problems are communicated to the solver using a simple, structured text file. The structure of the
text file is simple enough and will be explained shortly.

Linear programming problems consist of three components, namely, the optimization to be performed (i.e., maximize or minimize), 
a system of linear constraints, and a linear objective function. The system of linear constraints takes the form:

<PRE>
Ax <= b
</PRE>

Here A is an _m_ x _n_ matrix representing _m_ equations or inequalities in _n_ unknowns. In other words, they are the form:

<PRE>
a[1][1]*x[1] + ... + a[1][n]*x[n] <= b[1]
a[2][1]*x[1] + ... + a[2][n]*x[n] <= b[2]
... 
a[i][1]*x[1] + ... + a[i][n]*x[n] <= b[i]
...
a[m][1]*x[1] + ... + a[m][n]*x[n] <= b[m]
</PRE>

In this example, a[i][j] and b[i] can be any value, and the constraint '<=' can be any of the three relational operators 
'<=', '==' or '>='. However, the x[j] values will only be non-negative (i.e., zero or positive).

The objective function takes the form:

<PRE>
Z = c[1]*x[1] + c[2]*x[2] + ... + c[n]*x[n]
</PRE>

The task is to maximize or minimize the value Z subject to the constraint Ax <= b and x[i] >= 0.

File structure is all based on a line, with values separated by commas and text fields enclosed in quotes. The lines must be 
presented in a specific order, but blank lines and lines that begin the '#' character are ignored. The order is:
* **operation**
* **dimensions**
* **column labels**
* **constraints**
* **objective function**

The **operation** is either 'maximize' or 'minimize'. The operations, however, can be abbreviated to 'max' or 'min'. 
These are the only two optimization operations that can be performed.

The **dimension** field is a comma-separated pair of integers -- **m**,**n** -- that represent the number of rows and 
columns in the constraint matrix **A**.

**Column labels** is a comma-separated list of **n** labels, one for each column of the constraint matrix **A**. Labels
are enclosed in quotes and must not, themselves, contain commas. The quotes are not actually checked but they are assumed
to be there and are removed.

**Constraints** are the **m** rows of the Ax <= b constraints. The first field of each row is a row label, followed by 
**n** real values of a[i], followed by the one of the three relational operators listed above, followed by the real
value b[i].

The **objective function** is itself a list of real values starting with a label in quotes.

An example file follows:

<PRE>
  1 # Optimize the cost of feeding fish at a fish farm.
  2 # Fish must have minimum amounts of protein and calcium.
  3 # Three brands - A, B, and C - to choose from.
  4 # Costs per unit are $50, $30 and $75, respectively.
  5
  6 minimize
  7 2,3
  8           "Brand A","Brand B","Brand C"
  9 "Protein",        2,        1,        4, >=, 10
 10 "Calcium",       12,       14,       14, >=, 84
 11    "Cost",       50,       30,       75
</PRE>

