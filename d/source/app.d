import std.algorithm: max, reverse, maxElement;
import std.math: PI;
import std.stdio;

import plt = matplotlibd.pyplot;

import matrixd;

void add_index_range(ref ulong[] indices, ulong beg, ulong end, ulong inc = 1) {
    for (ulong i = beg; i <= end; i += inc) {
       indices ~= (i);
    }
}

void add_index_const(ref ulong[] indices, ulong value, ulong numel) {
    while (numel--) {
        indices ~= value;
    }
}

void append_vector(ref double[] vec, ref double[] tail) {
    vec ~= tail;
}

double[] subvector_reverse(const ref Matrix vec, ulong idx_end, ulong idx_start) {
	double[] result;
	for(ulong i = idx_start; i<idx_end+1; ++i) {
		result ~= vec[0, i];
	}
    return result.reverse;
}

Matrix filter(Matrix B, Matrix A, const Matrix X, const Matrix Zi) pure {
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

Matrix filtfilt(Matrix B, Matrix A, const Matrix X) {
	Matrix Y; // output

    ulong len = X.Size()[1]; // length of input
    ulong na = A.Size()[1];
    ulong nb = B.Size()[1];
    ulong nfilt = (nb > na) ? nb : na;
    ulong nfact = 3 * (nfilt - 1); // length of edge transients

    if (len <= nfact) {
    	string msg = "Input data too short! Data must have length more than 3 times filter order.";
        throw new Exception(msg);
    }

    // set up filter's initial conditions to remove DC offset problems at the
    // beginning and end of the sequence
	Matrix b = new Matrix(1, nfilt, 0.0);
	Matrix a = new Matrix(1, nfilt, 0.0);    
	for(ulong i = 0; i<nfilt; ++i) {
		if(i < B.Size()[1]) {
			b[0,i] = B[0,i];
		}
		if(i < A.Size()[1]) {
			a[0,i] = A[0,i];
		}
	}   

    ulong[] rows, cols;
    //rows = [1:nfilt-1           2:nfilt-1             1:nfilt-2];
    add_index_range(rows, 0, nfilt - 2);
    if (nfilt > 2) {
        add_index_range(rows, 1, nfilt - 2);
        add_index_range(rows, 0, nfilt - 3);
    }
    //cols = [ones(1,nfilt-1)         2:nfilt-1          2:nfilt-1];
    add_index_const(cols, 0, nfilt - 1);
    if (nfilt > 2) {       
        add_index_range(cols, 1, nfilt - 2);
        add_index_range(cols, 1, nfilt - 2);
    }
    // data = [1+a(2)         a(3:nfilt)        ones(1,nfilt-2)    -ones(1,nfilt-2)];    

    ulong klen = rows.length;
    double[] data;
    for(int i = 0; i<klen; ++i) {
    	data ~= 0.0;
    }

    data[0] = 1.0 + a[0,1];  ulong j = 1;
    if (nfilt > 2) {
        for (ulong i = 2; i < nfilt; i++)
            data[j++] = a[0, i];
        for (ulong i = 0; i < nfilt - 2; i++)
            data[j++] = 1.0;
        for (ulong i = 0; i < nfilt - 2; i++)
            data[j++] = -1.0;
    }

    double[] leftpad = subvector_reverse(X, nfact, 1);
    double _2x0 = 2 * X[0, 0];
	for(ulong i = 0; i< leftpad.length; ++i) {
		leftpad[i] = _2x0 - leftpad[i];
	}

	double[] rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
    double _2xl = 2 * X[0, len-1];
	for(ulong i = 0; i< rightpad.length; ++i) {
		rightpad[i] = _2xl - rightpad[i];
	}

    double y0;
    double[] signal1, signal2, zi;

    double[] Xd = toDouble_v(X);
    append_vector(signal1, leftpad);
    append_vector(signal1, Xd);
    append_vector(signal1, rightpad);    

    // Calculate initial conditions
    Matrix sp = new Matrix(rows.maxElement+1, cols.maxElement+1, 0.0);
    for (ulong k = 0; k < klen; ++k) {
        sp[rows[k], cols[k]] = data[k];
    }    
	Matrix bb = new Matrix(1, nfilt, 0.0);
	for(ulong i = 0; i<nfilt; ++i) {
		bb[0,i] = b[0,i]; 
	}
	Matrix aa = new Matrix(1, nfilt, 0.0);  
	for(ulong i = 0; i<nfilt; ++i) {
		aa[0,i] = a[0,i]; 
	}	
	Matrix bb_segment = new Matrix(1, nfilt-1, 0.0);
	for(ulong i = 0; i<nfilt-1; ++i) {
		bb_segment[0,i] = bb[0,i+1]; 
	}	
	Matrix aa_segment = new Matrix(1, nfilt-1, 0.0);
	for(ulong i = 0; i<nfilt-1; ++i) {
		aa_segment[0,i] = aa[0,i+1]; 
	}		
	Matrix zzi = sp.Inv()*(bb_segment - (aa_segment * bb[0,0])).T();
    for(int i = 0; i<zzi.Size()[1]; ++i) {
    	zi ~= 0.0;
    }

    // Do the forward and backward filtering
    y0 = signal1[0];
	for(ulong i = 0; i< zzi.Size()[1]; ++i) {
		zi[i] = zzi[0, i] * y0;
	}    
	Matrix signal1_m = new Matrix(signal1);
	Matrix zi_m = new Matrix(zi);
    signal2 = toDouble_v(filter(B, A, signal1_m, zi_m));
    signal2 = signal2.reverse;
    y0 = signal2[0];
	for(ulong i = 0; i< zzi.Size()[1]; ++i) {
		zi[i] = zzi[0, i] * y0;
	} 
	Matrix signal2_m = new Matrix(signal2); 
	signal1_m = filter(B, A, signal2_m, zi_m);
	Y = new Matrix(subvector_reverse(signal1_m, signal1_m.Size()[1] - nfact - 1, nfact));
	return Y;
}

// filter: https://www.mathworks.com/help/matlab/ref/filter.html
void main() {
	// Output of cpp hard coded for comparison
	Matrix y_filter_out_cpp = new Matrix([-0.01, -0.013721, 0.00881235, 0.0175505, 0.0324195, 0.0533212, 0.100134, 0.132712, 0.170887, 0.21447, 0.293248, 0.350712, 0.386633, 0.440792, 0.492974, 0.542974, 0.570594, 0.615647, 0.677955, 0.73735, 0.773678, 0.806793, 0.836567, 0.882882, 0.905633, 0.924731, 0.960101, 0.971682, 0.959428, 0.943308, 0.943308, 0.939428, 0.951682, 0.940102, 0.944732, 0.945633, 0.902882, 0.876567, 0.866793, 0.833677, 0.79735, 0.757955, 0.695647, 0.630594, 0.562974, 0.492974, 0.460792, 0.426633, 0.370711, 0.313248, 0.23447, 0.174608, 0.133899, 0.0925831, 0.0309016, -0.0309018, -0.112583, -0.193899, -0.254608, -0.31447, -0.353248, -0.410712, -0.486633, -0.540792, -0.592974, -0.622974, -0.650594, -0.675647, -0.737955, -0.75735, -0.813677, -0.826793, -0.836567, -0.882881, -0.885633, -0.904731, -0.920101, -0.951681, -0.939427, -0.943307, -0.943307, -0.939427, -0.931681, -0.9001, -0.88473, -0.865632, -0.862881, -0.836567, -0.826793, -0.813677, -0.757349, -0.717954, -0.695647, -0.670594, -0.642974, -0.592974, -0.520792, -0.466633, -0.410712, -0.333248, -0.27447]);
	Matrix y_filtfilt_out_cpp = new Matrix([-0.1, -0.0292282, 0.0413037, 0.107357, 0.172694, 0.235082, 0.296289, 0.35209, 0.406265, 0.4586, 0.510888, 0.558931, 0.604539, 0.649532, 0.693741, 0.735007, 0.773183, 0.812134, 0.847737, 0.875885, 0.89648, 0.913443, 0.926707, 0.938218, 0.94394, 0.94785, 0.94994, 0.944218, 0.934707, 0.925443, 0.91448, 0.899884, 0.881737, 0.856134, 0.825183, 0.787007, 0.741741, 0.697532, 0.652539, 0.602931, 0.550888, 0.4946, 0.436265, 0.38009, 0.326289, 0.273082, 0.220694, 0.163357, 0.101304, 0.0387717, -0.0240001, -0.0827719, -0.141304, -0.203357, -0.266694, -0.329082, -0.388289, -0.44209, -0.490265, -0.5386, -0.582888, -0.628931, -0.670539, -0.705532, -0.739741, -0.769007, -0.797183, -0.824133, -0.851737, -0.871884, -0.89048, -0.903443, -0.914706, -0.924218, -0.925939, -0.925849, -0.921939, -0.916217, -0.904706, -0.893443, -0.880479, -0.861884, -0.839736, -0.816133, -0.793182, -0.769007, -0.741741, -0.707532, -0.670539, -0.628931, -0.580888, -0.5326, -0.482265, -0.42409, -0.358289, -0.287082, -0.214694, -0.141357, -0.0633037, 0.0192282, 0.1]);

	// generate input signal
	Matrix X = new Matrix([0:100])*2.0*PI/100.0; // create x-axis ending at pi
	Matrix Y = matrixd.sin(X); // full period sine wave
	Y = Y + matrixd.noise(X, 0.0, 0.1); // add noise vector

	// create filter
	Matrix b_coeff = new Matrix([1./10., 1./10., 1./10., 1./10., 1./10.,1./10., 1./10., 1./10., 1./10., 1./10.]); 
	Matrix a_coeff = new Matrix([1.]); 
	Matrix zi = new Matrix([ 0. ]);

	// filter, filtfilt
	Matrix input_signal = Y; 
	Matrix y_filter_out_d; 
	Matrix y_filtfilt_d;
	y_filter_out_d = filter(b_coeff, a_coeff, input_signal, zi);
	y_filtfilt_d = filtfilt(b_coeff, a_coeff, input_signal);

	// Fig - filter: input, cpp out, dlang out
	plt.plot(Y.toDouble_v, "b-");
	plt.plot(y_filter_out_cpp.toDouble_v, "r-");
	plt.plot(y_filter_out_d.toDouble_v, "g-");
	plt.xlabel("Samples");
	plt.ylabel("Magnitude");	
	plt.legend(["X", "Y_cpp", "Y_dlang"]);
	plt.grid();
	plt.savefig("result_filter.png");
	plt.clear();	

	// Fig - filtfilt: input, cpp out, dlang out
	plt.plot(Y.toDouble_v, "b-");
	plt.plot(y_filtfilt_out_cpp.toDouble_v, "r-");
	plt.plot(y_filtfilt_d.toDouble_v, "g-");
	plt.xlabel("Samples");
	plt.ylabel("Magnitude");	
	plt.legend(["X", "Y_cpp", "Y_dlang"]);
	plt.grid();
	plt.savefig("result_filtfilt.png");
	plt.clear();	
}