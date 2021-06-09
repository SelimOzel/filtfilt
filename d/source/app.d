import std.math: PI;
import std.stdio;

import plt = matplotlibd.pyplot;

import matrixd;

void main() {
	writeln("filtfiltD");

	// generate input signal
	Matrix X = new Matrix([0:100])*2.0*PI/100.0; // create x-axis ending at pi
	Matrix Y = matrixd.sin(X); // full period sine wave
	Y = Y + matrixd.noise(X, 0.0, 0.1); // add noise vector

	// Fig: unfiltered input 
	plt.plot(Y.toDouble_v, "b-");
	plt.xlabel("Magnitude");
	plt.ylabel("Samples");	
	plt.legend();
	plt.grid();
	plt.savefig("input_signal.png");
	plt.clear();

	Matrix b_coeff = new Matrix([0.1, 0.1]); // Initialise your feedforward coefficients here
	Matrix a_coeff = new Matrix([0.1, 0.1]); // Initialise your feedback coefficients here

	Matrix input_signal = Y; // Some input data to be filtered

	//Matrix y_filter_ori = new Matrix();
	//MATLAB output from calling y_filtfilt_ori = filtfilt(b_coeff, a_coeff, input_signal);

}