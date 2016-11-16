#ifndef _CDGA__H
#define _CDGA__H

#include <string>
#include <vector>
#include <map>

using namespace std;

bool isEven(int n);

class FreeCGA;

// A generator is a pair consisting of a label and the degree of the generator
// This label must be unique for each generator
class Generator
{
public:
	Generator()
	{
		label = "Unitialized";
		degree = 0;
	}
	Generator(string _label, int _degree)
	{
		label = _label;
		degree = _degree;
	}
	string label;
	int degree;
};

// A structure needed only by the Word class
class GenPower
{
public:
	GenPower()
	{
		degree = 0;
		power = 0;
	}
	GenPower(int _degree, int _power)
	{
		degree = _degree;
		power = _power;
	}
	int degree;
	int power;
};

// A "word" is a concatenation of generators
class Word
{
public:
	Word();

	void Clear();

	// The following two methods should never be used on the "unit" word because it will keep its unit status
	void AddPowerOfGenerator(const string &label, int deg, int power);
	void AddPowerOfGenerator(const Generator &g, int power); 
	
	// Note: Multiplication has a rather strange behavior.
	// The product obtained is always rearranged in the "canonical lexographical" ordering. However,
	// doing so might introduce a sign change. The sign change will not be accounted for in the object,
	// however, the following methods will return "1" if no sign change occured and "-1" if a sign change
	// occured. It also returns 0 if the resulting product is 0.
	int MultiplyOnLeft(const Generator &g, int power);
	int MultiplyOnLeft(const Word &g);
	int MultiplyOnRight(const Generator &g, int power);
	int MultiplyOnRight(const Word &g);

	string OutputString() const; // Returns a string of the form "a^3 * b^2 * c^6"
	int GetDegree() const; // Returns the degree of the word
	int GetLength() const; // Returns the length of the word
	void SetToUnit(bool _isUnit); // Make this the "unit" element, i.e. the generator "1" in degree 0
	bool IsUnit() const; // Return true if this is the "unit" element

	// Return the first factor of the word as a Generator and return the remaining factors as a Word
	void GetFirstFactor(Generator &first_factor, Word &remaining_factors) const;

	static void ConcatenateWords(Word &result, const Word &w1, const Word &w2);

	// Check equality of two words
	bool operator==(const Word &word) const;

private:
	map<string, GenPower> generators;
	int degree;
	bool isUnit;
};

typedef vector<Word> OrderedBasis;

class Term
{
public:
	int coeff;
	Word word;
};

class LinearCombination
{
public:
	LinearCombination();

	void AddTerm(int coeff, const Word &word);
	void AddTerms(const LinearCombination &lc);
	void ScalarMultiply(int coeff);
	void MultiplyOnLeft(const Generator &g);
	void MultiplyOnLeft(const Word &word);
	void MultiplyOnRight(const Generator &g);
	void MultiplyOnRight(const Word &word);
	void MakeZero();
	
	// This method returns the coordinates of a vector (a linear combination) in terms of a certain ordered basis
	// The pointer "coordinates" must point to an array of the proper size, that is its size must be the number of
	// elements in the basis.
	void GetCoordinates(int *coordinates, const OrderedBasis &basis) const;

	void Simplify(); // Combine all repeated terms in a single term and sum the respective coefficients
	string OutputString() const; // Output a string representing the linear combination (this will not "simplify" the string first)

private:
	vector<Term> terms;
};

istream& operator>>(istream&, LinearCombination&);
ostream& operator<<(ostream&, const LinearCombination&); 

// A graded vector space is created by specifying a basis of generators
class GradedVectorSpace
{
	friend FreeCGA;
public:
	GradedVectorSpace();
	GradedVectorSpace(const vector<Generator> &_basis);

	void AddGenerator(const string &label, int degree);
	static int GetGeneratorDegree(const string &label); // Return the degree of the unique generator with name "label"
	vector<vector<Generator> > GetDegreeIndexedBasis(); // Return a basis ordered by degree

private:
	int maxDegree;
	vector<Generator> even_basis;
	vector<Generator> odd_basis;

	// This list maps unique generator labels to their degree
	static map<string, int> globalGeneratorsList;
};

typedef vector<LinearCombination> OrderedLCBasis;

// The following is the free commutative graded algebra on a graded vector space X
// This objects can be extended to a (/\X, d)-differential module of the form:
//    /\X (+) (/\X (x) T
// where (+) stands for the direct sum, (x) stands for the tensor product, and T is a graded module
// Note: Objects of type vector<OrderedBasis> are simply an array of ordered bases in each degree, indexed by the degree
class FreeCGA
{
public:
	FreeCGA();
	FreeCGA(const GradedVectorSpace &_X);

	void SetGradedVectorSpace(const GradedVectorSpace &_X);

	void GetBasis(vector<OrderedBasis> &basis, int degree, int minLength = 0);

	// This method returns an ordered basis in a given degree on /\^{>=n}X, it is indexed by degree
	// That is, it will return only words of length greater or equal to "minLength"
	void GetDegreeIndexedBasis(vector<OrderedBasis> &basis, int degree, int minLength = 0);
	
	// This method returns an ordered basis in a given degree on /\^{>=n}X (+) (/\^{+}X (x) T), it is indexed by degree
	// That is, it will return only words of length greater or equal to "minLength" in /\X and it will return all words
	// consisting of at least one factor in X and exactly one factor in T
	void GetDegreeIndexedBasisExtended(vector<OrderedBasis> &basis, int degree, int minLength);

	// This method returns an ordered basis in a given degree on /\X, it is indexed by word length
	void GetLengthIndexedBasis(vector<OrderedBasis> &basis, int degree);

	// This adds a generator to the vector space T
	void AddExtensionGenerator(string label, int degree);

	void Test1();
	void Test2();
	void Test3();

private:
	GradedVectorSpace X;
	GradedVectorSpace T;
};

class Differential
{
public:
	Differential();

	// Define the differential of a generator as "a"
	void SetDifferential(const string &generator_label, const LinearCombination &a);
	// Define the differential of a generator to be zero
	void SetDifferentialToZero(const string &generator_label);

	// This method computes the differential from a vector space (/\V)^n ---> (/\V)^{n+1}
	// The argument passed must consist of a basis for (/\V)^n (the source) and a basis for
	// (/\V)^{n+1} (the target). Assuming that dim (/\V)^n = N and dim (/\V)^{n+1} = M, then
	// this method returns the differential as a M x N matrix. Therefore, it is necessary to pass
	// an M x N array as the argument "differential_matrix"
	void ComputeDifferentialMatrix(int **differential_matrix, const OrderedBasis &source, const OrderedBasis &target);
	
	void EvaluateDifferential(LinearCombination &result, const Word & word);

private:

	map<string, LinearCombination> differential;
};



// Use this to read or write a generator to a file
istream& operator>>(istream &stream, Generator &g);
ostream& operator<<(ostream &stream, const Generator &g);

// Use this to read or write a word to a file
istream& operator>>(istream &stream, Word &word);
ostream& operator<<(ostream &stream, const Word &word);

// Use this to read or write a linear combination to a file
// For inptu, we assumed the string is trimmed.
istream& operator>>(istream &stream, LinearCombination &lc);
ostream& operator<<(ostream &stream, const LinearCombination &lc);

// Return 'true' if the file was parsed successfully.
// A FreeCGA and a Differential object will be returned.
bool ReadInputFromFile(const string &filename, string &output_filename, string &extension_output_filename, FreeCGA &cdga, Differential &differential, int &homology_degree_start, int &homology_degree_end, int &category);

#endif