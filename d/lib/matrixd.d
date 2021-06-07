module matrixd;

// D
import std.conv: to; 

// Enums
enum uint MAXROWS = 100;
enum uint MAXCOLUMNS = 100;

// Output is csv
string PrintMatrix(const Matrix matrix_IN) {
	string result;
	for(ulong r = 0; r < matrix_IN._nr; ++r) {
		for(ulong c = 0; c < matrix_IN._nc; ++c) {
			result ~= to!string(matrix_IN._m[r][c]);
			if(c == matrix_IN._nc - 1) result ~= "\n";
			else result ~= ",";
		}
	}
	return result;
}

// Lightweight
class Matrix {
public:
// nxm filled with n
this(
	const ulong rowLength_IN, 
	const ulong columnLength_IN, 
	const double n) 
pure {
	initialize(rowLength_IN, columnLength_IN, n);
}

// A = double[][]
this(const double[][] matrixRHS_IN) pure {
	initialize(matrixRHS_IN.length, matrixRHS_IN[0].length, 0.0);
	for(ulong r = 0; r < _nr; ++r) {
		for(ulong c = 0; c < _nc; ++c) {
			_m[r][c] = matrixRHS_IN[r][c];
		}
	}   	
}

// A = double[]
this(const double[] rowvectorRHS_IN) pure {
	initialize(1, rowvectorRHS_IN.length, 0.0);
	for(ulong r = 0; r < _nr; ++r) {
		for(ulong c = 0; c < _nc; ++c) {
			_m[r][c] = rowvectorRHS_IN[c];
		}
	}   	
}

// A = [[x1, x2 ...], [y1, y2, ...]]
void opAssign(const double[][] matrixRHS_IN) pure {
	initialize(matrixRHS_IN.length, matrixRHS_IN[0].length, 0.0);
	for(ulong r = 0; r < _nr; ++r) {
		for(ulong c = 0; c < _nc; ++c) {
			_m[r][c] = matrixRHS_IN[r][c];
		}
	}    	
}		

// += matrix
// -= matrix
// *= matrix
void opOpAssign(string operation_IN)(const Matrix rhs_IN) pure {
	Matrix result = new Matrix(_nr, _nc, 0.0);
	if(operation_IN == "+"){	
		result = this + rhs_IN;	
		_m = result._m;		
	}
	else if(operation_IN == "-"){	
		result = this - rhs_IN;	
		_m = result._m;		
	}    	
	else if(operation_IN == "*"){	
		result = this * rhs_IN;	
		_m = result._m;
		_nr = result.Size()[0];
		_nc = result.Size()[1];
	}    	
}    

// += scalar
// -= scalar
// *= scalar
// /= scalar
void opOpAssign(string operation_IN)(const double rhs_IN) pure {
	Matrix result = new Matrix(_nr, _nc, 0.0);
	if(operation_IN == "+"){	
		result = this + rhs_IN;		
		_m = result._m;	
	}
	else if(operation_IN == "-"){	
		result = this - rhs_IN;		
		_m = result._m;	
	}
	else if(operation_IN == "*"){	
		result = this * rhs_IN;		
		_m = result._m;	
	}
	else if(operation_IN == "/"){	
		result = this / rhs_IN;		
		_m = result._m;	
	}  		 	
}      

// matrix + matrix
// matrix - matrix
// matrix * matrix
Matrix opBinary(string operation_IN)(const Matrix rhs_IN) pure const {
	Matrix result = new Matrix(_nr, _nc, 0.0);
	bool sum = operation_IN == "+";
	bool subtract = operation_IN == "-";
	bool multiply = operation_IN == "*";
	ulong rhs_nr = rhs_IN._nr;
	ulong rhs_nc = rhs_IN._nc;
	if(sum || subtract){
		if(rhs_nr == _nr && rhs_nc == _nc) {
			for(ulong r = 0; r < _nr; ++r) {
				for(ulong c = 0; c < _nc; ++c) {
					double rhs = rhs_IN._m[r][c];
					if(sum) {
						result._m[r][c] = _m[r][c] + rhs;
					}
					else if (subtract) {
						result._m[r][c] = _m[r][c] - rhs;
					}
				}
			}		
		}
		else {
			string matrix_addition_err = "matrix add/subtract: wrong dimensions.";
			throw new Exception(matrix_addition_err);
		}	  
	}  	
	else if(multiply) {			
		// Verify mXn * nXp condition
		if(_nc == rhs_nr) {
			result = new Matrix(_nr, rhs_nc, 0.0); // reshape to mXp
			for (ulong r = 0; r<_nr; r++) {
				for (ulong c = 0; c<rhs_nc; c++) {
					for (ulong k = 0; k<_nc; k++) {
						result._m[r][c] += _m[r][k] * rhs_IN._m[k][c];
					}
				}
			}				
		}
		else {
			string matrix_multip_err = "matrix multiplication: wrong dimensions.";
			throw new Exception(matrix_multip_err);
		}				
	}
	return result;
}    

// matrix + scalar
// matrix - scalar
// matrix * scalar
// matrix / scalar
Matrix opBinary(string operation_IN)(const double rhs_IN) pure const {
	Matrix result = new Matrix(_nr, _nc, 0.0);
	if(operation_IN == "+"){	
		for(ulong r = 0; r < _nr; r++) {
			for(ulong c = 0; c < _nc; c++) {
				result._m[r][c] = _m[r][c] + rhs_IN;
			}
		}
	}
	else if(operation_IN == "-"){	
		for(ulong r = 0; r < _nr; r++) {
			for(ulong c = 0; c < _nc; c++) {
				result._m[r][c] = _m[r][c] - rhs_IN;
			}
		}
	}   
	else if(operation_IN == "*"){	
		for(ulong r = 0; r < _nr; r++) {
			for(ulong c = 0; c < _nc; c++) {
				result._m[r][c] = _m[r][c] * rhs_IN;
			}
		}
	} 
	else if(operation_IN == "/"){	
		for(ulong r = 0; r < _nr; r++) {
			for(ulong c = 0; c < _nc; c++) {
				result._m[r][c] = _m[r][c] / rhs_IN;
			}
		}
	} 			 	
	return result;
}        

// A == B
override bool opEquals(Object o) pure const {
	auto rhs = cast(const Matrix)o;
	if(rhs.Size()[0] == Size()[0] && rhs.Size[1] == Size[1]){
		for(ulong r = 0; r<rhs._nr; ++r) {
			for(ulong c = 0; c<rhs._nc; ++c) {
				if(rhs[r,c] != _m[r][c]) return false;
			}
		}
	}
	else return false;
	return true;
}

// A[1,2] = x
void opIndexAssign(double val, ulong r, ulong c) pure {
	_m[r][c] = val;
}

// x = A[1,2]
double opIndex(ulong r, ulong c) pure const {
	return _m[r][c];
}

// [x1, x2, ...] = A[1]
Matrix opIndex(ulong r) pure const {
	Matrix row_vector = new Matrix(1, _nc, 0.0);
	for (ulong c = 0; c<_nc; ++c) {
		row_vector[0, c] = _m[r][c];
	}
	return row_vector;
}

// Transpose
Matrix T() pure const {
	Matrix transpose = new Matrix(_nc, _nr, 0.0);
	for (ulong r=0; r<transpose.Size()[0]; r++) {
		for (ulong c=0; c<transpose.Size()[1]; c++) {
			transpose._m[r][c] = _m[c][r];
		}
	}
	return transpose;
}

Matrix Inv() pure const
{
	ulong nr = Size()[0];
	ulong nc = Size()[1];

	if(nr == nc) {	
	    // Find determinant of A[][] 
	    double determinant = Det(nr); 
	    if (determinant == 0) { 
	        throw new Exception("Inverse error: determinant must be non-zero\n"); 
	    } 
	  
	    // Inverse
	    Matrix inverse = new Matrix(nr,nr,0.0);

	    // Find adjoint 
	    Matrix adj = new Matrix(nr, nr, 0.0);
	    adjoint(adj); 

	    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	    for (ulong i=0; i<nr; i++) 
	        for (ulong j=0; j<nr; j++) 
	        	inverse[i,j] = adj[i,j]/determinant;
	  
	    return inverse; 
    }	
	else {
		throw new Exception("Inverse error: not square\n");
	}    
}

// Obtained from geeks-for-geeks: https://www.geeksforgeeks.org/determinant-of-a-matrix/
double Det(ulong n) pure const {
	ulong nr = Size()[0];
	ulong nc = Size()[1];

	if(nr == nc) {
		//  Base case : if matrix contains single element 
		double D = 0; // Initialize result 
		if (n == 1) return _m[0][0]; 

		//std::vector<std::vector<double>> temp(n, std::vector<double>(n)); // To store cofactors 
		Matrix temp = new Matrix(n,n,0.0); // To store cofactors 

		int sign = 1;  // To store sign multiplier 

		// Iterate for each element of first row 
		for (ulong f = 0; f < n; f++) {
			// Getting Cofactor of mat[0][f] 
			cofactor(temp, 0, f, n); 
			D += sign * _m[0][f] * temp.Det(n - 1); 
			// terms are to be added with alternate sign 
			sign = -sign; 
		}

		return D;	
	}
	else {
		throw new Exception("Determinant error: not square\n");
	}
}

// Sums all elements
double Sum() pure const {
	double all_sum = 0.0;
	for(ulong r = 0; r < _nr; r++) {
		for(ulong c = 0; c < _nc; c++) {
			all_sum += _m[r][c];
		}
	}
	return all_sum;		
}    

// Sums all elements in row r.
double Sum(const ulong r) pure const {
	if(r < _nr) {
		double row_sum = 0.0;
		for(ulong c = 0; c < _nc; c++) {
			row_sum += _m[r][c];
		}
		return row_sum;	
	}
	else {
		string sum_row_err = "Sum row: index out of bounds";
		throw new Exception(sum_row_err);
	}
}	

// [rows, cols]
ulong[2] Size() pure const {
	return [_nr, _nc];
}    

private:
void initialize(
	const ulong r_IN, 
	const ulong c_IN, 
	const double n_IN) 
pure {
	string initialize_err_maxrow = "initialize: Increase MAXROWS.";
	string initialize_err_maxcolumn = "initialize: Increase MAXCOLUMNS.";		
	if(r_IN > MAXROWS) throw new Exception(initialize_err_maxrow);
	if(c_IN > MAXCOLUMNS) throw new Exception(initialize_err_maxcolumn);
	if(r_IN != 0 && c_IN != 0) {
		for(ulong r = 0; r<r_IN; ++r) {
			for(ulong c = 0; c<c_IN; ++c) {
				_m[r][c] = n_IN;
			}
		}
		_nr = r_IN;
		_nc = c_IN;
	}
	else {
		throw new Exception("Matrix dimensions can't be zero.");
	}
}

// Obtained from geeks-for-geeks: https://www.geeksforgeeks.org/determinant-of-a-matrix/
void cofactor(ref Matrix temp, ulong p, ulong q, ulong n) pure const {
    ulong i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (ulong row = 0; row < n; row++) { 
        for (ulong col = 0; col < n; col++) { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) { 
            	//temp.Set(i,j++,_m[row][col]);
                temp[i,j++] = _m[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 	
}

// Obtained from geeks-for-geeks: https://www.geeksforgeeks.org/adjoint-inverse-matrix/
void adjoint(ref Matrix adj) pure const
{ 
	ulong N = _nr;
    if (N == 1) { 
        adj[0,0] = 1.0; 
        return; 
    } 
  
    // temp is used to store cofactors of A[][] 
    int sign = 1;
    //std::vector<std::vector<double>> temp(N, std::vector<double>(N)); // To store cofactors 
  	Matrix temp = new Matrix(N, N, 0.0);

    for (ulong i=0; i<N; i++) { 
        for (ulong j=0; j<N; j++) { 
            // Get cofactor of A[i][j] 
            cofactor(temp, i, j, N); 

            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j,i] = (sign)*(temp.Det(N-1));
        } 
    } 
} 

ulong _nr = 0;
ulong _nc = 0;
double[MAXROWS][MAXCOLUMNS] _m;
}