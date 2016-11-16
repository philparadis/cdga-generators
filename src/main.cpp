#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <Windows.h>

#include "cdga.h"
#include "homology.h"

using namespace std;

// *** Program specifications ***
// Written by Philippe Paradis

// The program "cdga-generator" works as follows (note that /\ represents the symbol lambda below):
// 
// INPUT: The user provides four pieces of data to the program (numbers (3) and (4) can be taken to be 0).
// (1) A graded vector space X: that is, a list of generators each with a unique label and a degree
// (2) A differential X --> /\X: that is, for each generator of X, it's differential must be given as a word in /\X
// (3) A graded vector space T (note again that all labels must be unique and different from those of X).
// (4) A differential T --> /\X (+) (/\X (x) T). Here the the symbol "(+)" stands for direct sum and "(x)" stands
//     for the tensor product.
//
// INPUT: The user also provides the following input:
// (5) The "degree" at which homology shall be computed. It is possible to also specify a range of degrees.
// e.g. writing "5..12" would compute the homology of the space from degree 5 to degree 12.
// (6) The "category" of the space. This parameter is optional. If this is provided, then rather than computing a
// minimal model for /\X, the program will calculate a minimal model for the projection of /\X ---> /\X / /\^{>n} X,
// where 'n' is the category.
//
// NOTATION: Denote by Z the (/\X, d)-differential graded module Z = (/\X (+) (/\X (x) T)).
//
// OUTPUT: The program will then provide partial information about the cdga Z. It will output the following 5 pieces of data:
// (1) An ordered basis of Z in degree n-1
// (2) An ordered basis of Z in degree n
// (3) An ordered basis of Z in degree n+1
// (4) The matrix of the differential Z^{n-1} --> Z^n
// (5) The matrix of the differential Z^n --> Z^{n+1}
//
// LIMITATIONS:
// (1) The vector space X can only have generators in degree 2 or above. (This condition could easily be lifted, if necessary,
// at the cost of performance)
// (2) You cannot give the label "1" to a generator, because this labels stands for the unit in /\X. Also the label "0" is reserved.
// Moreover, the following characters cannot be used in a label: ' ', '+', '-', '(', ')', '*', '^'.

// TODO: Make the program check automatically if d^2 = 0 and if the differential has degree +1?
// TODO: Automate the entire process of computing the rational retraction index. This can be done, although the work required is almost
//       certainly asymptotically exponential with respect the dimension of the model X in each degree. The steps to follow, in order of
//       priority are:
//  (1)  Modify the input and output format such that they are completely standardized. The ultimate goal is that the *output* file produced
//       by the computations performed by cdga-generators, say with input parameter of degree n, should be a valid *input* file that
//       cdga-generators would accept and say perform another computation with input parameter of degree n+1. Hence, the cdga could be
//       computed recursively until the rational index is found. Instead, at the moment, the files must be manually edited between each
//       step or a range of degrees must be used at the same time, which can be very computationally expensive due to exponential growth
//       of work required for the linear algebra step as the degree increases. 
//       For example, we need to switch to symbols that support an arbitrary amount of alphabet letters (rather than 1), otherwise we shall
//       get stuck once Z requires more than 26 generator in a certain degree, among other things.
//  (2)  For each step above (i.e. everytime we increase the degree of T and introduce a new set of generators for that degree, starting
//       from degree 0 and stopping at degree cat_0(X)), compute all possible retractions rho from Z to the minimal model (/\V, d) such
//       and measure the length of all words within rho(T_{cat_0 X}). Find the smallest word length and call it r_d for d = 0, ..., cat_0(X).
//       Finally, the rational retraction index is the largest value r_d over all values of d. For more details and a clearer explanation,
//       you can refer to my thesis "On the Rational Retraction Index" and Chapter 4 (see the first definition).
// TODO: Add a series of unit tests to make sure no errors are introduced in currently working code.

static streambuf* buffer;

static void pause()
{
	cout << endl;
	system("pause");
}

void RunTest1(const string &input_filename, ofstream &output)
{
	FreeCGA cdga;
	Differential diff;
	int degree_start, degree_end;
	int category;
	string output_filename, extension_output_filename;

	if (!ReadInputFromFile(input_filename, output_filename, extension_output_filename, cdga, diff, degree_start, degree_end, category)) {
		cerr << "Failure to read input file '" << input_filename << "'." << endl;
		return;
	}

	if (degree_start <= 1) {
		cerr << "Invaid degree. The degree must be greater or equal to 2." << endl;
		return;
	}
	if (degree_end < degree_start) {
		degree_end = degree_start;
	}

	if (!output_filename.empty()) {
		output.open(output_filename);
		cout << "Redirecting all output to '" << output_filename << "'..." << endl;
		// Make cout redirect to "output"
		buffer = cout.rdbuf();
		cout.rdbuf(output.rdbuf());
	}

	cout << "Successfully parsed input file '" << input_filename << "'..." << endl;
	if (degree_start < degree_end) {
		cout << "Now computing a basis of cocycles for the homology in degrees " << degree_start << " to " << degree_end << " (assuming the category to be " << category << ")..." << endl << endl;
	} else {
		cout << "Now computing a basis of cocycles for the homology in degree " << degree_start << " (assuming the category to be " << category << ")..." << endl << endl;
	}
	
	vector<OrderedBasis> basis;
	cdga.GetDegreeIndexedBasisExtended(basis, degree_end, category+1);
	
	for (int degree = degree_start; degree <= degree_end; degree++) {
		// Here we compute a cocycles basis in the specified degree
		int dim_source = (int)basis[degree].size();
		int dim_target = (int)basis[degree+1].size();
		OrderedLCBasis cocycles_basis, image_basis;
		if (dim_source != 0) {
			// Allocate memory for the matrix
			int **diff_matrix = 0;
		
			if (dim_target != 0) {
				diff_matrix = new int*[dim_target];
				for (int i=0; i<dim_target; i++) {
					diff_matrix[i] = new int[dim_source];
				} 

				// This compute the differential of d : (deg n) ---> (deg n+1) and store it into "diff_matrix"
				diff.ComputeDifferentialMatrix(diff_matrix, basis[degree], basis[degree+1]);
			}

			int dim_prev = (int)basis[degree-1].size();
			if (dim_prev == 0) {
				// The image of d_{n-1} is zero, so it suffices to find the cocycles, because none of them will be boundaries
				FindCocycleBasis(diff_matrix, dim_target, dim_source, cocycles_basis, basis[degree]);
			} else {
				// First, allocate memory for the differential of d : (deg n-1) ---> (deg n)
				int **diff_matrix_prev = new int*[dim_source];
				for (int i=0; i<dim_source; i++) {
					diff_matrix_prev[i] = new int[dim_prev];
				}
			
				// Then, compute the differential of d : (deg n-1) ---> (deg n)
				diff.ComputeDifferentialMatrix(diff_matrix_prev, basis[degree-1], basis[degree]);

				// Finally, find a basis for the homology in degree n
				cerr << "Degree: " << degree << endl;
				FindHomologyBasis(diff_matrix_prev, dim_source, dim_prev, diff_matrix, dim_target, dim_source, cocycles_basis, image_basis, basis[degree]);
			
				// Free up memory
				for (int i=0; i<dim_source; i++) {
					delete [] diff_matrix_prev[i];
				}
				delete [] diff_matrix_prev;
			}
		
			// Free up memory used for the matrix
			if (diff_matrix) {
				for (int i=0; i<dim_target; i++) {
					delete [] diff_matrix[i];
				}
				delete [] diff_matrix;
			}
		}

		int homology_dim = (int)cocycles_basis.size();
		int boundaries_dim = (int)image_basis.size();
		//cout << "The homology has dimension " << homology_dim << ". " << endl << endl;
	
		cout << "HOMOLOGY DEGREE " << degree << " (DIM " << homology_dim << "):" << endl << endl;

		if (homology_dim > 0) {
			OrderedLCBasis::iterator iter;
			//cout << "Here is a basis of cocycles:" << endl << endl;
			for (iter = cocycles_basis.begin(); iter != cocycles_basis.end(); iter++) {
				cout << iter->OutputString() << endl;
			}

			// This only makes senses if homology is being computed for one degree and not a range of degree.
			// So we only output this information for degree == degree_start.
			if (!extension_output_filename.empty() && degree == degree_start) {
				ofstream extension_file(extension_output_filename);
				extension_file << "The homology in degree " << degree << " has dimension " << homology_dim << " and hence can be killed by introducing the following generators:" << endl << endl;
				extension_file << "Extension:" << endl;
				int index = 0;
				int num_digits = (int)ceil(log((double)homology_dim+1.0)/log(10.0));
				for (iter = cocycles_basis.begin(); iter != cocycles_basis.end(); iter++) {
					stringstream ss;
					ss << setfill('0') << setw(num_digits) << ++index;
					extension_file << "�_" + ss.str() << " " << degree-1 << endl;
				}
				index = 0;
				extension_file << endl << "Differential:" << endl;
				for (iter = cocycles_basis.begin(); iter != cocycles_basis.end(); iter++) {
					stringstream ss;
					ss << setfill('0') << setw(num_digits) << ++index;
					extension_file << "d(" << "�_" + ss.str() << ") = " << iter->OutputString() << endl;
				}
				extension_file.close();
			}

			cout << endl;
		}

		/*
		// For now I'm commenting this block out. I don't think it's useful to know the image of d.
		if (boundaries_dim > 0) {
			OrderedLCBasis::iterator iter;
			cout << "Here is a basis for the image of d_" << degree-1 << ":" << endl << endl;
			for (iter = image_basis.begin(); iter != image_basis.end(); iter++) {
				cout << iter->OutputString() << endl;
			}
			cout << endl;
		}*/
	}

	cout << "Here is a basis of the extended cdga from degree 0 up to degree " << degree_end+1 << "." << endl << endl;
	for (int deg=0; deg<=degree_end+1; deg++) {
		int dim = basis[deg].size();
		int i = 0;
		cout << "DEGREE " << deg << " (dim " << dim << "):" << endl;
		OrderedBasis::iterator iter;
		int minLength = -1;
		int maxLength = 0;
		for (iter = basis[deg].begin(); iter != basis[deg].end(); iter++) {
			if (i>0)
				cout << ", ";
			if (iter->GetLength() < minLength || minLength == -1)
				minLength = iter->GetLength();
			if (iter->GetLength() > maxLength)
				maxLength = iter->GetLength();
			cout << iter->OutputString();
			i++;
		}
		cout << endl;
		if (dim > 0) {
			cout << "Min word length: " << minLength << endl;
			cout << "Max word length: " << maxLength << endl;
		}
		cout << endl;
	}

	cout << endl;
}

void RunTest2()
{
	ifstream file("test_file.txt");
	if (!file.is_open()) {
		cerr << "Can't open file." << endl;
		return;
	}
	GradedVectorSpace X;
	X.AddGenerator("a", 2);
	X.AddGenerator("b", 2);
	X.AddGenerator("c", 2);
	X.AddGenerator("x", 3);
	X.AddGenerator("y", 3);
	X.AddGenerator("z", 3);
	X.AddGenerator("attention", 5);
	X.AddGenerator("gamma", 5);
	X.AddGenerator("yo_man_21", 5);
	Word w1, w2, w3;
	file >> w1 >> w2 >> w3;
	cout << w1 << endl
		<< w2 << endl
		<< w3 << endl;

	LinearCombination lc;
	file >> lc;
	cout << "Linear combination is: " << lc << endl;
}

void RunTest3()
{
}

int main(int argc, char **argv)
{
	atexit(pause);

	cout << "Welcome to cdga-generators! Written by Philippe Paradis (June 2011)." << endl;

	ofstream output;
	string input_filename;

	if (argc == 2) {
		input_filename = argv[1];
	} else {
		cout << "Please enter a filename: ";
		getline(cin, input_filename);
	}
	
	DWORD timeBegin = GetTickCount();

	try {
		RunTest1(input_filename, output);
		//RunTest3();
	} catch (logic_error &e) {
		cerr << "An exception has occured: " << e.what() << endl;
	}

	DWORD timeEnd = GetTickCount();
	cout << "Time elapsed: " << (int)(timeEnd - timeBegin) << " milliseconds." << endl;
	cout.rdbuf(buffer);

	return 0;
}
