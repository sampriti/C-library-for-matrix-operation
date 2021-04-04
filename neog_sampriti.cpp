
//C++ implementation of matrix multiplication and transpose

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

/* Function Declaration*/

void subtract(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C, int n, int l);
void strassen(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C,
	int n,int l,int m);
void sum(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C, int n,int l);
void ikjalgorithm(vector< vector<int> > A,
	vector< vector<int> > B,
	vector< vector<int> > &C, int n, int l, int m);
void transpose(vector< vector<int> > &A, vector< vector<int> > &B, int n, int m);
void printMatrix(vector< vector<int> > matrix, int n,int m);

/* Calculating subtraction of two matrices */

void subtract(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C, int n, int l) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < l; j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
}

/* Iterative matrix multiplication  */

void ikjalgorithm(vector< vector<int> > A,
	vector< vector<int> > B,
	vector< vector<int> > &C, int n,int l,int m){
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m;j++) {
			
			for (int k = 0; k < l; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}


/* Strassen Matrix multiplication 
   where A(nxl) is the 1st matrix ,
   B(lxm) is the 2nd matrix,
   C(nxm) is the product  */
   

void strassen(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C, int n,int l,int m) 
{
	if (n == 1||l==1||m==1)
	{
		ikjalgorithm(A, B, C, n,l,m);
		return;
	}

	
	else 
	{
		int nval = (n >> 1) + (n & 1);  // nval=(n/2) when n is even or (n/2)+1 when n is odd
		int lval = (l >> 1) + (l & 1);  // lval=(l/2) when l is even or (l/2)+1 when l is odd
		int mval = (m >> 1) + (m & 1);  // mval=(m/2) when m is even or (m/2)+1 when m is odd

		vector<int> inner(lval);   
		vector<int> inner1(mval);
		vector<int> inner2(nval);
		vector< vector<int> >
			a11(nval, inner), a12(nval, inner), a21(nval, inner), a22(nval, inner),
			b11(lval, inner1), b12(lval, inner1), b21(lval, inner1), b22(lval, inner1),
			c11(nval, inner1), c12(nval, inner1), c21(nval, inner1), c22(nval, inner1),
			p1(nval, inner1), p2(nval, inner1), p3(nval, inner1), p4(nval, inner1),
			p5(nval, inner1), p6(nval, inner1), p7(nval, inner1),
			aResult(nval, inner), bResult(lval, inner1);
		
		int i, j;
		
		//Dividing the 1st matrix in 4 sub-matrices:
		
		for (i = 0; i < nval; i++) {
			for (j = 0; j < lval; j++) 
			{
				a11[i][j] = A[i][j];
				if (i + nval < n)
					a21[i][j] = A[i + nval][j];
				else
					a21[i][j] = 0;
				if ((j+lval)<l)
				{
					a12[i][j] = A[i][j + lval];
					
				}
				else
				{
					a12[i][j] = 0;
					
				} 
				if (i + nval < n && j+lval< l)
					a22[i][j] = A[i + nval][j + lval];
				else
					a22[i][j] = 0;

			}
		}
		

		//Dividing the 2nd matrix in 4 sub-matrices:
		for (i = 0; i < lval; i++) {
			for (j = 0; j < mval; j++)
			{

				b11[i][j] = B[i][j];
				if (j + mval < m)
					b12[i][j] = B[i][j + mval];
				else
					b12[i][j] = 0;
				if ((i + lval)<l)
				{

					b21[i][j] = B[i + lval][j];
					
				}
				else
				{
					b21[i][j] = 0;
					
				}
				if (i + nval < n && j + lval< l)
					b22[i][j] = B[i + lval][j + mval];
				else
					b22[i][j] = 0;

			}
		}
		
		// Calculating p1 to p7 :

		sum(a11, a22, aResult, nval,lval);                               // a11 + a22
		sum(b11, b22, bResult, lval,mval);                               // b11 + b22
		strassen(aResult, bResult, p1, nval, lval, mval);                // p1 = (a11+a22) * (b11+b22)
		
		sum(a21, a22, aResult, nval,lval);                               // a21 + a22
		strassen(aResult, b11, p2, nval, lval, mval);                    // p2 = (a21+a22) * (b11)
		
		subtract(b12, b22, bResult, lval,mval);                          // b12 - b22
		strassen(a11, bResult, p3, nval, lval, mval);                    // p3 = (a11) * (b12 - b22)
		
		subtract(b21, b11, bResult, lval,mval);                          // b21 - b11
		strassen(a22, bResult, p4, nval, lval, mval);                    // p4 = (a22) * (b21 - b11)
		
		sum(a11, a12, aResult, nval,lval);                               // a11 + a12
		strassen(aResult, b22, p5, nval, lval, mval);                    // p5 = (a11+a12) * (b22)
		
		subtract(a21, a11, aResult, nval,lval);                          // a21 - a11
		sum(b11, b12, bResult, lval, mval);                              // b11 + b12
		strassen(aResult, bResult, p6, nval, lval, mval);                // p6 = (a21-a11) * (b11+b12)


		subtract(a12, a22, aResult, nval,lval);                          // a12 - a22
		sum(b21, b22, bResult, lval,mval);                               // b21 + b22
		strassen(aResult, bResult, p7, nval, lval, mval);                // p7 = (a12-a22) * (b21+b22)
		
		
		// calculating c21, c21, c11 c22:
		
		sum(p3, p5, c12, nval, mval);                                    // c12 = p3 + p5
		sum(p2, p4, c21, nval, mval); // c21 = p2 + p4
	
		sum(p1, p4, aResult, nval, mval); // p1 + p4
		sum(aResult, p7, bResult, nval, mval); // p1 + p4 + p7
		subtract(bResult, p5, c11, nval, mval); // c11 = p1 + p4 - p5 + p7

		sum(p1, p3, aResult, nval, mval); // p1 + p3
		sum(aResult, p6, bResult, nval, mval); // p1 + p3 + p6
		subtract(bResult, p2, c22, nval, mval); // c22 = p1 + p3 - p2 + p6
		
		
		// Grouping the results obtained in the final result matrix:
		
		for (i = 0; i < nval; i++)
		{
			for (j = 0; j < mval; j++) 
			{
				C[i][j] = c11[i][j];
				if(j+mval<m)
					C[i][j+mval] = c12[i][j];
				if(i+nval<n)
					C[i + nval][j] = c21[i][j];
				if (i + nval<n  && j+mval<m)
					C[i + nval][j + mval] = c22[i][j];
			}
		}
		
		
	}
}


/* Calculating sum of two matrices*/

void sum(vector< vector<int> > &A,
	vector< vector<int> > &B,
	vector< vector<int> > &C, int n,int l) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < l; j++) {
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}


/* PRINT A MATRIX*/

void printMatrix(vector< vector<int> > matrix, int n,int m) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (j != 0) {
				cout << "\t";
			}
			cout << matrix[i][j];
		}
		cout << endl;
	}
}

/* Calculating transpose of matrix*/

void transpose(vector< vector<int> > &A, vector< vector<int> > &B, int n,int m)
{
	for (int i = 0;i < m;i++)
		for (int j = 0;j < n;j++)
			B[i][j] = A[j][i];
}

int main()
 {

	/* Input the rows and columns of matrices to multiply */

	int Arows;
	cout << "Enter the Number of rows in 1st Matrix A: ";
	cin >> Arows;
	int AcolBrow;
	cout << "Number of Rows in 1st matrix equals Number of columns in 2nd Matrix for Matrix multiplication\n";
	cout << "Enter the Number of columns in 1st Matrix A or Number of rows in 2nd Matrix B: ";
	cin >> AcolBrow;
	int Bcol;
	cout << "Enter the Number of columns in 2nd Matrix B: ";
	cin >> Bcol;

	/* Create vectors for 1st and 2nd Matrix */

	vector<int> inner(AcolBrow);
	vector<int> inner1(Bcol);
	vector<int> inner2(Arows);
	vector< vector<int> > A(Arows, inner), B(AcolBrow, inner1), C(Arows, inner1);


	/* Taking user input for the values of both the matrices to multiply */

	cout << "Enter the 1st matrix values: \n";
	for (int i = 0; i < Arows; i++)
	{
		cout << "the contents of " << i << "th row\n";

		for (int j = 0; j < AcolBrow; j++)
		{
			cin >> A[i][j];
		}

		cin.ignore(numeric_limits<streamsize>::max(), '\n');                                    // ignore any limit to the number of characters
	}
	cout << "\n";
	cout << "The 1st Matrix to multiply is : \n";
	printMatrix(A, Arows, AcolBrow);
	cout << "\n";

	cout << "Enter the 2nd Matrix values: \n";
	for (int i = 0; i < AcolBrow; i++)
	{
		cout << "Enter the contents of " << i << "th row\n";

		for (int j = 0; j < Bcol; j++)
		{
			cin >> B[i][j];
		}

		cin.ignore(numeric_limits<streamsize>::max(), '\n');                                    // ignore any limit to the number of characters
	}
	cout << "\n";
	cout << "The 2nd Matrix to multiply is : \n";
	printMatrix(B, AcolBrow, Bcol);
	cout << "\n";


	/* Performing strassen multiplication of the two matrices */

	strassen(A, B, C, Arows, AcolBrow, Bcol);
	cout << "Product of 1st and 2nd matrix is: \n";
	printMatrix(C, Arows, Bcol);
	cout << "\n";




	

	/****** Taking user input for size of matrix to transpose *******/

	int Xrow, Xcol;
	cout << "Enter the number of Rows and Columns for the Matrix to transpose: ";
	cin >> Xrow >> Xcol;


	vector<int> inner3(Xcol);
	vector<int> inner4(Xrow);
	vector< vector<int> > X(Xrow, inner3), D(Xcol, inner4);

	/* User input the values of matrix for transpose*/

	cout << "Enter the values of matrix to transpose:\n";
	for (int i = 0; i < Xrow; i++)
	{
		cout << "Enter the contents of " << i << "th row:";

		for (int j = 0; j < Xcol; j++)
		{
			cin >> X[i][j];
		}

		cin.ignore(numeric_limits<streamsize>::max(), '\n');
	}
	cout << "User Input Matrix for Transpose is: \n";
	printMatrix(X, Xrow, Xcol);
	cout << "\n";

	/* Calculating the transpose of the user provided matrix */

	cout << "Matrix transpose is \n";
	transpose(X, D, Xrow, Xcol);
	printMatrix(D, D.size(), D[0].size());

	cout << "End of Output";

	return 0;
}
