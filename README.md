# C-library-for-matrix-operation

Compile the file using following command in Linux using g++ compiler:
g++ neog_sampriti.cpp -o neog_sampriti
Run the file using following command in Linux terminal:
./neog_sampriti


Code explanation
I have performed matrix multiplication for all kinds of matrices(Non square, square, odd or even size) of different sizes using Strassen algorithm. Strassen algorithm does modification over traditional divide and conquer algorithm. In divide and conquer method, we have 8 recursive calls while in Strassen it is only 7 recursive calls to do a matrix multiplication. It is a huge improvement over traditional iterative algorithm to calculate matrix product. i.e. 
 
The time complexity of iterative method is O(N^3), divide and conquer is O(N^3) but for Strassen, it is O(N^2.8074) where N denotes size of a input square matrix.
In Strassen algorithm, the steps are:-
	Divide the matrix into square sub-matrices of size (N/2) using recursion on Strassen function where N is the size of the initial matrix. Continue till the submatrix size is equal to 1.
	For rectangular matrix, we pad the sub matrices with 0’s to make it a square matrix. (Made this modification from traditional Strassen algorithm to use for rectangular matrices too)
	Then calculate 7 products using recursion:
P1 = (A11 + A22) (B11 + B22) 
P2 = (A21 + A22) B11 
P3 = A11(B12 − B22) 
P4 = A22(B21 − B11) 
P5 = (A11 + A12)B22 
P6 = (A21 − A11)(B11 + B12) 
P7 = (A12 − A22)(B21 + B22)
Where A11,A12,A21,A22 are coefficient of 1st submatrix and coefficient of 2nd submatrix.
	Then the above products are used to calculate the coefficients of final matrix:
C11 = P1 + P4 − P5 + P7
 C12 = P3 + P5 
C21 = P2 + P4 
C22 = P1 − P2 + P3 + P6
	The Final matrix is 
(■(C11&C12@C21&C22))

Matrix transpose is calculated by changing columns to rows and rows to columns. Its time complexity is O(N^2) 

The overall time complexity of our algorithm is O(N^2.8074).


