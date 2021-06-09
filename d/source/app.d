import std.math: PI;
import std.stdio;

import plt = matplotlibd.pyplot;

import matrixd;

void main() {
	writeln("filtfiltD");

	// Output of cpp hard coded for comparison
	Matrix y_filter_out_cpp = new Matrix([-0.1, -0.0372095, 0.225333, 0.0873813, 0.14869, 0.209017, 0.468125, 0.325779, 0.381754, 0.435827, 0.687785, 0.537424, 0.584547, 0.628969, 0.670513, 0.709017, 0.744328, 0.776307, 1.00483, 1.02978, 1.05106, 0.868583, 0.882287, 1.09211, 0.898027, 0.9, 1.09803, 0.892115, 0.882287, 0.868583, 1.05106, 0.829776, 1.00483, 0.976307, 0.944328, 0.909017, 0.670513, 0.628969, 0.784547, 0.537424, 0.687785, 0.435827, 0.381754, 0.325779, 0.268125, 0.209017, 0.34869, 0.287381, 0.225333, -0.0372095, -0.1, -0.162791, -0.0253332, -0.0873813, -0.34869, -0.409017, -0.468125, -0.525779, -0.381754, -0.635827, -0.487785, -0.737424, -0.784547, -0.628969, -0.870513, -0.709017, -0.744328, -0.776307, -1.00483, -0.829776, -1.05106, -0.868583, -0.882287, -1.09211, -0.898027, -0.9, -0.898027, -1.09211, -0.882287, -0.868583, -1.05106, -0.829776, -0.804827, -0.776307, -0.744328, -0.709017, -0.870513, -0.828969, -0.784547, -0.737424, -0.487785, -0.435827, -0.581754, -0.525779, -0.468125, -0.209017, -0.14869, -0.287381, -0.225333, 0.0372095, 0.1]);

	// generate input signal
	Matrix X = new Matrix([0:100])*2.0*PI/100.0; // create x-axis ending at pi
	Matrix Y = matrixd.sin(X); // full period sine wave
	Y = Y + matrixd.noise(X, 0.0, 0.1); // add noise vector

	// Fig: unfiltered input 
	plt.plot(Y.toDouble_v, "b-");
	plt.plot(y_filter_out_cpp.toDouble_v, "r-");
	plt.xlabel("Magnitude");
	plt.ylabel("Samples");	
	plt.legend();
	plt.grid();
	plt.savefig("input_signal.png");
	plt.clear();

	Matrix b_coeff = new Matrix([0.1, 0.1]); 
	Matrix a_coeff = new Matrix([0.1, 0.1]); 
	Matrix input_signal = Y; 

	Matrix y_filter_out; Matrix zi = new Matrix([ 0.0 ]);
	//Matrix y_filter_ori = new Matrix();
	//MATLAB output from calling y_filtfilt_ori = filtfilt(b_coeff, a_coeff, input_signal);
}