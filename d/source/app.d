import std.stdio;
import matrixd;

void main() {
	writeln("filtfiltD");

	Matrix b_coeff = new Matrix([0.1, 0.1]); // Initialise your feedforward coefficients here
	Matrix a_coeff = new Matrix([0.1, 0.1]); // Initialise your feedback coefficients here

	Matrix input_signal = new Matrix([1.0, 2.0, 2.5, 2.0, 1.0, 0.5, 0.25, 0.5]); // Some input data to be filtered
	//Matrix y_filter_ori = new Matrix();
	//MATLAB output from calling y_filtfilt_ori = filtfilt(b_coeff, a_coeff, input_signal);

}