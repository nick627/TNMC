# TNMC
Source codes for discipline number-theoretic methods in cryptography (IBKS, SPbSTU, St. Petersburg, fall 2017)

I use MPIR and MPFR library 

Link for help install
https://www.youtube.com/watch?v=S06mX5dwIJ0
https://www.youtube.com/watch?v=je5ei4rCFiw


Task number 1
For each pair of numbers, find the greatest common divisor and its linear representation using the following algorithms:
1. Advanced Euclid's algorithm.
2. Extended Euclidean binary algorithm.
3. Advanced Euclidean algorithm with "truncated" remnants.
In the report, lead all the iterations of the algorithms (sequence of residues and two sequences of linear representation coefficients). With the number of iterations greater than 20, give the first 5 and the last 5 elements of each sequence. Iterations must be numbered.
When obtaining different linear representations of the same greatest common divisor, justify the correctness of the results obtained.
In the report, compare the speed of implemented algorithms.


Task number 2
Each of these numbers is checked for simplicity with the help of tests: Farma, Solovea-Strassen, Rabin-Miller.
In the report lead:
1. The result of the test for each of the tests for each number is "probably simple" or "composite."
2. For a compound number, a base for which the condition of simplicity is violated; show what kind of condition is violated.
3. For a prime number, there are 5 bases for which the condition of simplicity is satisfied; show the feasibility of this condition.
4. For each of the tests, give an example of two or three Carmichael numbers (at least 30 decimal places each, indicate the source) and signs that this Carmichael number is "simple" and composite.


Task number 3
Each of these numbers is factorized:
1. Pollard ro-method.
2. (p-1) -the Pollard method.
3. By the method of a quadratic sieve.
4. The method of continued fractions.

In the report lead:
1. In the case of a successful completion of the program: the result of decomposition.
2. For the Pollard ro-method:
2.0. Parameters of the algorithm (used display, initial value (several if they had to be changed)).
2.1. When the work is completed successfully, the first 5 and the last 5 values of a, b, GCD (a-b, n), the number of iterations, the program time, conclusions (explanation of the effective completion).
2.2. When the program is running for more than several hours, the first 5 and the last (at the time of the program interruption) 5 values of a, b, GCD (a-b, n), the number of iterations performed, the time taken to execute them, the estimated time remaining until the completion of work (through an estimation of complexity of the algorithm), conclusions (explanation of non-productive completion).
3. For the (p-1) -method Pollard:
3.1 Decomposition base (original, modified (if a change was required, justify the need for a change)).
3.2 Values of indicators l_i.
3.3 With the successful completion of work - the basis of a, under which the decomposition is carried out (several if it was required to change the base).
3.4 If it is impossible to find the decomposition more than a few hours - the base a, at which the decomposition is performed (several if it was required to change the base), the number of iterations performed (one iteration - running the algorithm with one base a), the time taken to execute them, remaining until the completion of the work (through the evaluation of the complexity of the algorithm), conclusions (explanation of the non-productive completion).
4. For the method of continued fractions:
4.1 With the successful completion of work:
4.1.1 The decomposition database (original, modified (if needed, justify the need for a change)), with a large base volume - the number of elements in the database and its last element, the criteria for excluding elements from the database.
4.1.2 The numerators P_i of the appropriate fractions used in the expansion for the D-smooth values (P_i) ^ 2 (mod n), for a large data volume - 5 values of Pi and (P_i) ^ 2 (mod n).
4.1.3 The vectors of the exponents for the D-smooth values (P_i) ^ 2 (mod n), for a large data volume - 5 vectors corresponding to the values (P_i) ^ 2 (mod n) from the previous point.
4.1.4 The method used to search for linearly dependent vectors (if a ready procedure was used, specify from which library).
4.1.5 The values of s and t.
4.2 If it is not possible to find the decomposition more than a few hours:
4.2.1 Data on subparagraphs 4.1.1-4.1.3 at the time of the program interruption.
4.2.2 The calculated required number of smooth numbers (P_i) ^ 2 (mod n) and the time required to search for them.
4.3 Conclusions (explanation of the results of the work).
5. For the quadratic sieve method
5.1 With the successful completion of work:
5.1.1 The decomposition database (original, modified (if needed, justify the need for a change)), with a large base volume - the number of elements in the database and its last element, the criteria for excluding elements from the database.
5.1.2 The used range of values of x, the result of the screening procedure.
5.1.3. The values of x and the corresponding D-smooth values of f (x), for a large data volume, are 5 values of x and f (x).
5.1.4 The metric vectors for the D-smooth values of f (x), with a large data volume, are 5 vectors corresponding to the values of f (x) from the previous point.
5.1.5 The method used to search for linearly dependent vectors (if a ready procedure was used, specify from which library).
5.1.6 The values of s and t.
5.2 If it is not possible to find the decomposition more than a few hours:
5.2.1 4.2.1. Data on pp. 5.1.1-5.1.4 at the time of the program interruption.
5.2.2 4.2.2. The estimated required number of smooth numbers f (x) and the time required to search for them.
5.3 Conclusions (explanation of the results of the work).


Task number 4

The element a is of order q modulo p. Find the discrete logarithm of x - an integer 1 <x <q such that a ^ x = b (mod p):
1. Pollard ro-method.
2. The Gelfond method;
3. The decomposition base method.

In the report lead:
1. In the case of a successful completion of the program: the result of discrete logarithm is the number x.
2. For the Gel'fond method:
2.1. The value of s.
2.2. A database of the form {k, b * a ^ (-k s) (mod p)}, sorted by the first coordinate. With a large base volume - the first 5 and the last 5 of its elements.
2.3. The value of t for which a ^ t = b * a ^ (- k s) (mod p) for some k. Equation for calculating the logarithm of x. The method of solving the equation.
2.4. If it is not possible to find the discrete logarithm for more than a few hours - the value of t, to which the sequence a ^ t (mod p) is constructed, the time spent on constructing the second sequence, the estimated time remaining until the completion of the work (through an algorithm complexity estimation), conclusions (explanation of the non-productive completion).
3. For the decomposition base method:
3.1 With the successful completion of work:
3.1.1 The decomposition base B (original, modified (if needed, justify the need for change)), with a large base volume - the number of elements in the database and its last element.
3.1.2. The B-smooth values of a ^ (u_i) (mod p), their corresponding exponents u_i, for a large data volume - for the five specified values
3.1.3 The vectors of the exponents for the B-smooth values of a ^ (u_i) (mod p), for a large data volume, are 5 vectors corresponding to the values of a ^ (u_i) (mod p) from the previous point.
3.1.4 The exponent v for which the value b ^ v (mod p) is B-smooth, the corresponding index vector.
3.1.5 The method used to exclude variables (if a ready procedure was used, specify from which library).
3.1.6 The method used to solve the linear comparison with respect to the unknown x.
3.2 If it is not possible to find the discrete logarithm for more than a few hours:
3.2.1 Information on subparagraphs 4.1.1-4.1.3 at the time of the program interruption.
3.2.2 Estimated number of smooth numbers a ^ (u_i) (mod p) and the time required to search for them.
3.2.3 The exponent v for which the value b ^ v (mod p) is B-smooth, the corresponding index vector.

4. For the Pollard ro-method:
4.1 Algorithm parameters (used display, initial value (several if they had to be changed)).
4.2 When the result is completed, the first 5 and the last 5 values of c, d, log_a_ (c), log_a_ (d), the number of iterations, the program time, conclusions (explanation of the effective completion).
4.3 When the program is running for more than several hours, the first 5 and the last (at the time of the program interruption) 5 values of c, d, log_a_ (c), log_a_ (d), the number of iterations completed, the time taken to execute them, the estimated time remaining until completion of work (through an estimation of the complexity of the algorithm), conclusions (explanation of non-productive completion).
