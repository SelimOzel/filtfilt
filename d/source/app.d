import std.algorithm: max;
import std.math: PI;
import std.stdio;

import plt = matplotlibd.pyplot;

import matrixd;

Matrix filter(Matrix B, Matrix A, const Matrix X, const Matrix Zi) {
	Matrix Y; // output

	// input sanity
	if(A.empty) {
		throw new Exception("The feedback filter coefficients are empty.");
	}
	if(A.Size[0] > 1) {
		throw new Exception("The feedback filter coefficients is not a row vector.");
	}
	bool all_zeros = true;
	for(ulong c = 0; c<A.Size()[1]; ++c) {
		if(A[0,c] != 0.0) {
			all_zeros = false;
			continue;
		}
	}
	if(all_zeros) {
		string msg = "At least one of the feedback filter coefficients has to be non-zero.";
		throw new Exception(msg);
	}
	if(A[0,0] == 0.0) {
		throw new Exception("First feedback coefficient has to be non-zero.");
	}

	double a0 = A[0,0];
	if(a0 != 1.0) {
		for(ulong i = 0; i< A.Size()[0]; ++i) {
			A[0,i] = A[0,i] / a0;
		}
		for(ulong i = 0; i< B.Size()[0]; ++i) {
			B[0,i] = B[0,i] / a0;
		}		
	}

	ulong input_size = X.Size()[1];
	ulong filter_order = max(A.Size()[1], B.Size[1]);
	Matrix b = new Matrix(1, filter_order, 0.0);
	Matrix a = new Matrix(1, filter_order, 0.0);
	Matrix z = new Matrix(1, filter_order, 0.0);
	for(ulong i = 0; i<filter_order; ++i) {
		if(i < B.Size()[1]) {
			b[0,i] = B[0,i];
		}
		else b[0,i] = 0.0;		
		if(i < A.Size()[1]) {
			a[0,i] = A[0,i];
		}
		else a[0,i] = 0.0;	
		if(i < Zi.Size()[1]) {
			z[0,i] = Zi[0,i];
		}
		else z[0,i] = 0.0;				
	}	
	Y = new Matrix(1, input_size, 0);

	for(ulong i = 0; i<input_size; ++i) {
		ulong order = filter_order - 1;
        while (order) {
            if (i >= order) {
                z[0, order - 1] = b[0, order] * X[0, i - order] - a[0, order] * Y[0, i - order] + z[0, order];
            }
            --order;
        }
        Y[0, i] = b[0, 0] * X[0, i] + z[0, 0];
	}

	return Y;
}

// filter: https://www.mathworks.com/help/matlab/ref/filter.html
void main() {
	// Output of cpp hard coded for comparison
	Matrix y_filter_out_cpp = new Matrix([-0.02, -0.0274419, 0.0176247, 0.035101, 0.064839, 0.126642, 0.227709, 0.247798, 0.306673, 0.3641, 0.459854, 0.473714, 0.525467, 0.57491, 0.621848, 0.626094, 0.667475, 0.705827, 0.780999, 0.852852, 0.921261, 0.946112, 0.967308, 0.984764, 0.958413, 0.928201, 0.974091, 0.976056, 0.934092, 0.928203, 0.958415, 0.904764, 0.927307, 0.946111, 0.96126, 0.932852, 0.900999, 0.825827, 0.787475, 0.706094, 0.661848, 0.61491, 0.565467, 0.473714, 0.419854, 0.3241, 0.306673, 0.287798, 0.267709, 0.206642, 0.144839, 0.0425427, -0.0200001, -0.082543, -0.144839, -0.206643, -0.267709, -0.367798, -0.426673, -0.4841, -0.499854, -0.553714, -0.605467, -0.65491, -0.701848, -0.746094, -0.747475, -0.745827, -0.820999, -0.812852, -0.88126, -0.906111, -0.927307, -0.944763, -0.958413, -0.928201, -0.93409, -0.976055, -0.93409, -0.928201, -0.958413, -0.944763, -0.887307, -0.866111, -0.84126, -0.772851, -0.780998, -0.785827, -0.787475, -0.786094, -0.741848, -0.65491, -0.605467, -0.553714, -0.499854, -0.4441, -0.386673, -0.327798, -0.267709, -0.166642, -0.104839]);

	// generate input signal
	Matrix X = new Matrix([0:100])*2.0*PI/100.0; // create x-axis ending at pi
	Matrix Y = matrixd.sin(X); // full period sine wave
	Y = Y + matrixd.noise(X, 0.0, 0.1); // add noise vector

	// create filter
	Matrix b_coeff = new Matrix([1./5., 1./5., 1./5., 1./5., 1./5.]); 
	Matrix a_coeff = new Matrix([1.]); 
	Matrix zi = new Matrix([ 0. ]);

	// filter
	Matrix input_signal = Y; 
	Matrix y_filter_out_d; 
	y_filter_out_d = filter(b_coeff, a_coeff, input_signal, zi);

	// Fig: input filter, cpp filter out, dlang filter out
	plt.plot(Y.toDouble_v, "b-");
	plt.plot(y_filter_out_cpp.toDouble_v, "r-");
	plt.plot(y_filter_out_d.toDouble_v, "g-");
	plt.xlabel("Samples");
	plt.ylabel("Magnitude");	
	plt.legend(["X", "Y_cpp", "Y_dlang"]);
	plt.grid();
	plt.savefig("result_tilter.png");
	plt.clear();	
}