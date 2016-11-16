#include "cdga.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <assert.h>

using namespace std;

bool isEven(int n)
{
	return (n % 2 == 0);
}

int factorial(int n)
{
	int p = 1;
	for (int i=2; i<=n; i++) {
		p = p*i;
	}
	return p;
}

int nChoosek(int n, int k)
{
	if (k > n || k < 0)
		return 0;
	return factorial(n) / (factorial(k) * factorial(n-k));
}

void AddIndexedOrderedBasisElement(vector<OrderedBasis> &basis, Word word)
{
	// Make sure the array has been extended in the required degree
	// If not, add space for 100 more degrees
	if (word.GetDegree() >= (int)basis.size()) {
		basis.resize(basis.size()+100);
	}
	basis[word.GetDegree()].push_back(word);
}

void AddIndexedOrderedBasisElement(vector<vector<Generator> > &basis, Generator gen)
{
	// Make sure the array has been extended in the required degree
	// If not, add space for 100 more degrees
	if (gen.degree >= (int)basis.size()) {
		basis.resize(basis.size()+100);
	}
	basis[gen.degree].push_back(gen);
}

Word::Word()
{
	degree = 0;
	isUnit = false;
}

void Word::Clear()
{
	degree = 0;
	isUnit = false;
	generators.clear();
}

void Word::AddPowerOfGenerator(const string &label, int deg, int power)
{
	assert(!IsUnit());

	if (power == 0) {
		return;
	}

	// Check if g is already part of the word
	if( generators.find(label) == generators.end() ) {
		// Here, g is not part of the word, so we initialize the power of g to "power"
		generators[label] = GenPower(deg, power);
	} else {
		generators[label].power += power;
	}
	degree += deg * power;
}

void Word::AddPowerOfGenerator(const Generator &g, int power)
{
	AddPowerOfGenerator(g.label, g.degree, power);
}

int Word::MultiplyOnLeft(const Generator &g, int power)
{
	assert(!IsUnit());
	assert((power>=1) && (isEven(g.degree) || power == 1));

	if (isEven(g.degree)) {
		// Part of the symmetric algebra, just add the factor to the word
		AddPowerOfGenerator(g, power);
		return 1;
	} else {
		// Part of the exterior algebra, so three things can happen:
		// (1) The word might become zero if the factor is already present
		// (2) We may add the factor
		// (3) We may add the factor and need to multiply by (-1)
		if (generators.find(g.label) != generators.end()) {
			return 0;
		}
		generators[g.label] = GenPower(g.degree, 1);
		map<string, GenPower>::const_iterator iter;
		int position = 0;
		for (iter = generators.begin(); iter != generators.end(); iter++) {
			if (iter->first == g.label)
				break;
			if (!(isEven(iter->second.degree))) {
				position++;
			}
		}
		if (isEven(position)) {
			return 1;
		} else {
			return -1;
		}
	}
}
int Word::MultiplyOnLeft(const Word &word)
{
	int sign = 1;

	map<string, GenPower>::const_reverse_iterator iter;
	for (iter = word.generators.rbegin(); iter != word.generators.rend(); iter++) {
		Generator g(iter->first, iter->second.degree);
		sign *= MultiplyOnLeft(g, iter->second.power);
	}

	return sign;
}

int Word::MultiplyOnRight(const Generator &g, int power)
{
	assert(!IsUnit());
	assert((power>=1) && (isEven(g.degree) || power == 1));

	if (isEven(g.degree)) {
		// Part of the symmetric algebra, just add the factor to the word
		AddPowerOfGenerator(g, power);
		return 1;
	} else {
		// Part of the exterior algebra, so three things can happen:
		// (1) The word might become zero if the factor is already present
		// (2) We may add the factor
		// (3) We may add the factor and need to multiply by (-1)
		if (generators.find(g.label) != generators.end()) {
			return 0;
		}
		generators[g.label] = GenPower(g.degree, 1);
		map<string, GenPower>::const_reverse_iterator iter;
		int position = 0;
		for (iter = generators.rbegin(); iter != generators.rend(); iter++) {
			if (iter->first == g.label)
				break;
			if (!(isEven(iter->second.degree))) {
				position++;
			}
		}
		if (isEven(position)) {
			return 1;
		} else {
			return -1;
		}
	}
}

int Word::MultiplyOnRight(const Word &word)
{
	int sign = 1;
	
	map<string, GenPower>::const_iterator iter;
	for (iter = word.generators.begin(); iter != word.generators.end(); iter++) {
		Generator g(iter->first, iter->second.degree);
		sign *= MultiplyOnRight(g, iter->second.power);
	}

	return sign;
}

int Word::GetDegree() const
{
	return degree;
}

int Word::GetLength() const
{
	if (IsUnit()) {
		return 0;
	}

	map<string, GenPower>::const_iterator iter;
	int len = 0;
	for (iter = generators.begin(); iter != generators.end(); iter++) {
		len += iter->second.power;
	}
	return len;
}

string Word::OutputString() const
{
	map<string, GenPower>::const_iterator iter;
	string output;
	for (iter = generators.begin(); iter != generators.end(); iter++) {
		if (iter->second.power == 1) {
			if (!output.empty())
				output += " * ";
			output += iter->first;
		} else if (iter->second.power > 1) {
			if (!output.empty())
				output += " * ";
			stringstream ss;
			ss << iter->second.power;
			output += iter->first + "^" + ss.str();
		}
	}
	return output;
}

void Word::SetToUnit(bool _isUnit)
{
	if (_isUnit) {
		degree = 0;
		generators.clear();
		generators["1"] = GenPower(0, 1);
		isUnit = true;
	} else {
		isUnit = false;
	}
}

bool Word::IsUnit() const
{
	return isUnit;
}

void Word::GetFirstFactor(Generator &first_factor, Word &remaining_factors) const
{
	map<string, GenPower>::const_iterator iter = generators.begin();
	if (iter == generators.end()) {
		assert(0);
	}
	first_factor.label = iter->first;
	first_factor.degree = iter->second.degree;

	// Now, make "remaining_factors" equal to the remaining factors
	remaining_factors.Clear();
	if (iter->second.power > 1) {
		remaining_factors.AddPowerOfGenerator(iter->first, iter->second.degree, iter->second.power - 1);
	}
	iter++;
	for( ; iter != generators.end(); iter++) {
		remaining_factors.AddPowerOfGenerator(iter->first, iter->second.degree, iter->second.power);
	}
}

void Word::ConcatenateWords(Word &result, const Word &w1, const Word &w2)
{
	// Check if one of the summands is the unit
	if (w1.IsUnit()) {
		result = w2;
		return;
	} else if (w2.IsUnit()) {
		result = w1;
		return;
	}

	result = w2;
	
	map<string, GenPower>::const_iterator iter;
	for (iter = w1.generators.begin(); iter != w1.generators.end(); iter++) {
		result.AddPowerOfGenerator(iter->first, iter->second.degree, iter->second.power);
	}
}

bool Word::operator==(const Word& word) const
{
	if (this->GetLength() != word.GetLength()) {
		return false;
	}
	
	// Note: Since words are always ordered in the "canonical lexicographical" order, it is sufficient
	// to check whether two words have the same factors and the same powers for those factors

	map<string, GenPower>::const_iterator iter;
	for (iter = word.generators.begin(); iter != word.generators.end(); iter++) {
		map<string, GenPower>::const_iterator this_iter;
		this_iter = this->generators.find(iter->first);

		if (this_iter == this->generators.end()) {
			return false;
		}
		
		if (this_iter->second.power != iter->second.power) {
			return false;
		}
	}

	return true;
}

map<string, int> GradedVectorSpace::globalGeneratorsList;

GradedVectorSpace::GradedVectorSpace()
{
	maxDegree = 0;
}

GradedVectorSpace::GradedVectorSpace(const vector<Generator> &_basis)
{
	maxDegree = 0;
	vector<Generator>::const_iterator iter;
	for (iter = _basis.begin(); iter != _basis.end(); iter++) {
		if (isEven(iter->degree)) {
			even_basis.push_back(*iter);
		} else {
			odd_basis.push_back(*iter);
		}
		if (iter->degree > maxDegree) {
			maxDegree = iter->degree;
		}
	}
}

void GradedVectorSpace::AddGenerator(const string &label, int degree)
{
	Generator gen(label, degree);

	// Check if a generator with the same name already exists
	if (globalGeneratorsList.find(label) != globalGeneratorsList.end()) {
		throw logic_error("The generator with label '" + label +"' already exists.");
	}
	globalGeneratorsList[label] = degree;

	if (isEven(degree)) {
		even_basis.push_back(gen);
	} else {
		odd_basis.push_back(gen);
	}
	if (maxDegree < degree)
		maxDegree = degree;
}

int GradedVectorSpace::GetGeneratorDegree(const string &label)
{
	map<string,int>::iterator iter = globalGeneratorsList.find(label);
	if (iter == globalGeneratorsList.end()) {
		throw logic_error("The generator with label '" + label + "' has not been introduced.");
	}
	return iter->second;
}

vector<vector<Generator> > GradedVectorSpace::GetDegreeIndexedBasis()
{
	vector<vector<Generator> > basis(50); // Start with an basis that is empty in degrees 0-49
	vector<Generator>::iterator iter;

	// First iterate over even generators
	for (iter = even_basis.begin(); iter != even_basis.end(); iter++) {
		AddIndexedOrderedBasisElement(basis, *iter);
	}
	// Now iterate over odd generators
	for (iter = odd_basis.begin(); iter != odd_basis.end(); iter++) {
		AddIndexedOrderedBasisElement(basis, *iter);
	}

	return basis;
}

Word CreateWord(const vector<Generator> &basis, int length, int *indices)
{
	Word word;
	for (int i=0; i<length; i++) {
		word.AddPowerOfGenerator(basis[indices[i]], 1);
	}
	return word;
}

// The dimension of the symmetric algebra in length k on a n dimension vector space is (n+k-1) choose k.
void GetOrderedBasisSymmetricAlgebra(OrderedBasis &symm_alg_basis, const vector<Generator> &basis, int wordlength)
{
	int dim = basis.size();
	int k = wordlength;
	if (k <= 0) {
		// Return a basis that only contains the word "1" in degree 0
		OrderedBasis degreeZeroBasis;
		Word word;
		word.SetToUnit(true);
		degreeZeroBasis.push_back(word);
		symm_alg_basis = degreeZeroBasis;
		return;
	}

	// Make sure the OrderedBasis that was passed is empty
	symm_alg_basis.clear();

	// Now we need to iterate through the integers 0 <= i_1 <= ... <= i_k < n to form a basis of Symm^{k}(X).
	// So we create an array of k integers
	int *indices = new int[k];

	// Initialize all indices to 1. That means the first word is going to be the first generator in "basis" times itself, k times.
	for (int i=0; i<k; i++)
		indices[i] = 0;

	// Here we iterate through the integers 0 <= i_1 <= ... <= i_k < n and add each possible word
	int a;
	do
	{
		a = 0;
		while (indices[a] != dim) {
			Word word = CreateWord(basis, k, indices);
			symm_alg_basis.push_back(word);
			indices[a]++;
		}

		while (indices[a] == dim && a < k) {
			a++;
			if (a < k) {
				indices[a]++;
				for (int i=0; i<a; i++) {
					indices[i] = indices[a];
				}
			}
		}
	} while (a < k);

	delete [] indices;
}

void GetOrderedBasisExteriorAlgebra(OrderedBasis &ext_alg_basis, const vector<Generator> &basis, int wordlength)
{
	int dim = basis.size();
	int k = wordlength;

	if (k <= 0) {
		// Return a basis that only contains the word "1" in degree 0
		OrderedBasis degreeZeroBasis;
		Word word;
		word.SetToUnit(true);
		degreeZeroBasis.push_back(word);
		ext_alg_basis = degreeZeroBasis;
		return;
	}

	// The length of words cannot be greater than the dimension of the space, otherwise the exterior algebra is empty
	if (k > dim) {
		assert(0);
	}

	// Make sure the OrderedBasis that was passed is empty
	ext_alg_basis.clear();

	// Now we need to iterate through the integers 0 < i_1 < ... < i_k < n to form a basis of Ext^{k}(X).
	// So we create an array of k integers
	int *indices = new int[k];

	// Initialize the indices
	int a = k-1;
	for (int i=0; i<k; i++) {
		indices[i] = a;
		if (a > 0)
			a--;
	}
			
	// Next, iterate through all indices with 0 <= i_1 < ... < i_k < n
	do {
		// Create a word
		Word word = CreateWord(basis, k, indices);
		ext_alg_basis.push_back(word);

		// Increment the first index that does not surpass the dimension
		int a = 0;
		indices[0]++;
		bool loop = false;
		if (indices[0] == dim) {
			loop = true;
		}
		while (loop && a < k) {
			a++;
			if (a < k) {
				indices[a]++;
				for (int i=a-1; i>=0; i--) {
					indices[i] = indices[i+1] + 1;
				}
			}
			loop = false;
			for (int i=0; i<k; i++) {
				if (indices[i] >= dim)
					loop = true;
			}
		}
	} while (indices[0] < dim && a < k);

	delete [] indices;
}

FreeCGA::FreeCGA()
{
}

FreeCGA::FreeCGA(const GradedVectorSpace &_X)
{
	X = _X;
}

void FreeCGA::SetGradedVectorSpace(const GradedVectorSpace &_X)
{
	X = _X;
}

void FreeCGA::Test1()
{
	OrderedBasis basis;
	cout << "Symmetric algebra of X^even:" << endl;
	
	cout << "X^even has dimension: " << X.even_basis.size() << "." << endl;
	for (int i=1; i<=6; i++) {
		int n = X.even_basis.size();
		cout << "Length " << i << " has dimension " << nChoosek(n+i-1, i) << "." << endl;
	}
		
	system("pause");

	for (int i=1; i<=6; i++) {
		cout << "Length " << i << ": ";
		GetOrderedBasisSymmetricAlgebra(basis, X.even_basis, i);
		for (unsigned int i = 0; i<basis.size(); i++) {
			if (i > 0)
				cout << ", ";
			stringstream ss;
			ss << basis[i].GetDegree();
			cout << basis[i].OutputString(); // << " (deg " << ss.str() << ")";
		}
		cout << " (DIM = " << basis.size() << ")";
		cout << endl;
	}
}

void FreeCGA::Test2()
{
	OrderedBasis basis;
	cout << "Exterior algebra of X^odd:" << endl;
	
	cout << "X^odd has dimension: " << X.odd_basis.size() << "." << endl;
	for (int i=1; i<=6; i++) {
		int n = X.odd_basis.size();
		cout << "Length " << i << " has dimension " << nChoosek(n, i) << "." << endl;
	}

	for (int k=1; k<=6; k++) {
		int dim = X.odd_basis.size();
		int *indices = new int[k];
		cout << "All indices needed for length " << k << ":" << endl;

		if (k <= dim) {
			// Initialize the indices
			int a = k-1;
			for (int i=0; i<k; i++) {
				indices[i] = a;
				if (a > 0)
					a--;
			}
			
			// Next, iterate through all indices with 0 <= i_1 < ... < i_k < n
			do {
				// Print the index
				for (int i=0; i<k; i++) {
					if (i>0)
						cout << ", ";
					cout << indices[i];
				}
				cout << endl;

				// Increment the first index that does not surpass the dimension
				int a = 0;
				indices[a]++;
				while (indices[a] == dim && a < k) {
					a++;
					if (a < k) {
						indices[a]++;
						for (int i=a-1; i>=0; i--) {
							indices[i] = indices[i+1] + 1;
						}
					}
				}
			} while (indices[0] < dim && a < k);


		} else {
			cout << "(empty)" << endl;
		}

		delete [] indices;
	}
}

void FreeCGA::Test3()
{
	OrderedBasis basis;
	cout << "Exterior algebra of X^odd:" << endl;
	
	cout << "X^odd has dimension: " << X.odd_basis.size() << "." << endl;
	int dim = X.odd_basis.size();
	for (int i=1; i<=dim; i++) {
		cout << "Length " << i << " has dimension " << nChoosek(dim, i) << "." << endl;
	}
		
	system("pause");

	for (int i=1; i<=dim; i++) {
		cout << "Length " << i << ": ";
		GetOrderedBasisExteriorAlgebra(basis, X.odd_basis, i);
		for (unsigned int i = 0; i<basis.size(); i++) {
			if (i > 0)
				cout << ", ";
			stringstream ss;
			ss << basis[i].GetDegree();
			cout << basis[i].OutputString(); // << " (deg " << ss.str() << ")";
		}
		cout << " (DIM = " << basis.size() << ")";
		cout << endl;
	}
}

void FreeCGA::GetDegreeIndexedBasis(vector<OrderedBasis> &basis, int degree, int minLength)
{
	// Start with a basis that is empty in each degree. Allocate space for degrees from 0 to 499.
	basis.clear();
	basis.resize(100);

	// Assuming that X = {X^i}, i>=2, it suffices to look at words of length up to ceil(degree+1/2)
	int maxLength = (int)ceil( ((double)(degree+1)) / 2.0 );

	// Recall that /\X = Symm(X^even) (x) Ext(X^odd)
	// Here, "n" represents the length of words in Symm(X^even) and "m" represents the length of words in Ext(X^odd)
	// We iterate over all possible words starting in length 1 and ending in length "maxLength"
	for (int len=minLength; len<=maxLength; len++) {
		int n = len;
		int m = 0;
	
		// As soon as n is negative or m is greater than the dimension of X^odd, we have to stop
		OrderedBasis symm_basis, ext_basis;
		while (n >= 0 && m <= (int)X.odd_basis.size()) {
			GetOrderedBasisSymmetricAlgebra(symm_basis, X.even_basis, n);
			GetOrderedBasisExteriorAlgebra(ext_basis, X.odd_basis, m);
			
			// Add all products of words in symm_basis and ext_basis
			OrderedBasis::iterator iter_symm, iter_ext;
			Word word;
			for (iter_ext = ext_basis.begin() ; iter_ext != ext_basis.end(); iter_ext++) {
				for (iter_symm = symm_basis.begin(); iter_symm != symm_basis.end(); iter_symm++) {
					Word::ConcatenateWords(word, *iter_symm, *iter_ext);
					AddIndexedOrderedBasisElement(basis, word);
				}
			}
			
			// Next, increment the length of m and decrement the length of n
			n--;
			m = len - n;
		}
	}
}

void FreeCGA::GetDegreeIndexedBasisExtended(vector<OrderedBasis> &basis, int degree, int minLength)
{
	// Start with a basis that is empty in each degree. Allocate space for degrees from 0 to 99.
	basis.clear();
	basis.resize(100);

	// First, get a basis for /\^{>=n}X
	// TODO: This is rather unefficient, having the main routine called twice for different word lengths. Could be optimized.
	vector<OrderedBasis> X_basis;
	GetDegreeIndexedBasis(basis, degree, minLength);
	GetDegreeIndexedBasis(X_basis, degree, 1);
	
	// Assuming that X = {X^i}, i>=2, it suffices to look at words of length up to ceil(degree+1/2)
	int maxLength = (int)ceil( ((double)(degree+1)) / 2.0 );

	// Now, compute a basis for /\^{+} X (x) T
	vector<vector<Generator> > T_basis = T.GetDegreeIndexedBasis();
	for (int i=0; i<(int)T_basis.size(); i++) {
		if (!T_basis[i].empty()) {
			vector<Generator>::iterator Titer;
			for (Titer = T_basis[i].begin(); Titer != T_basis[i].end(); Titer++) {
				Word T_word;
				T_word.AddPowerOfGenerator(*Titer, 1);
				// For each element of T, we will consider all possible products with elements of X_basis
				for (int deg=0; deg<(int)X_basis.size(); deg++) {
					OrderedBasis::iterator Xiter;
					Word word;
					for (Xiter = X_basis[deg].begin(); Xiter != X_basis[deg].end(); Xiter++) {
						Word::ConcatenateWords(word, *Xiter, T_word);
						AddIndexedOrderedBasisElement(basis, word);
					}
				}
			}
		}
	}
}

// This adds a generator to the vector space T
void FreeCGA::AddExtensionGenerator(string label, int degree)
{
	T.AddGenerator(label, degree);
}

LinearCombination::LinearCombination()
{
}

void LinearCombination::AddTerm(int coeff, const Word &word)
{
	if (coeff != 0) {
		Term term;
		term.coeff = coeff;
		term.word = word;
		terms.push_back(term);
	}
}

void LinearCombination::AddTerms(const LinearCombination &lc)
{
	vector<Term>::const_iterator iter;
	for (iter = lc.terms.begin(); iter != lc.terms.end(); iter++) {
		AddTerm(iter->coeff, iter->word);
	}
}

void LinearCombination::ScalarMultiply(int coeff)
{
	vector<Term>::iterator iter;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		iter->coeff *= coeff;
	}
}

void LinearCombination::MakeZero()
{
	terms.clear();
}

void LinearCombination::MultiplyOnLeft(const Generator &g)
{
	// TODO
	// OOPS: Big mistake. Did not worry about ordering of factors in a word.
	vector<Term>::iterator iter;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		int sign = iter->word.MultiplyOnLeft(g, 1);
		iter->coeff *= sign;
	}
}

void LinearCombination::MultiplyOnLeft(const Word &word)
{
	vector<Term>::iterator iter;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		int sign = iter->word.MultiplyOnLeft(word);
		iter->coeff *= sign;
	}
}

void LinearCombination::MultiplyOnRight(const Generator &g)
{
	vector<Term>::iterator iter;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		int sign = iter->word.MultiplyOnRight(g, 1);
		iter->coeff *= sign;
	}
}

void LinearCombination::MultiplyOnRight(const Word &word)
{
	vector<Term>::iterator iter;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		int sign = iter->word.MultiplyOnRight(word);
		iter->coeff *= sign;
	}
}

// This method returns the coordinates of a vector (a linear combination) in terms of a certain ordered basis
// The pointer "coordinates" must point to an array of the proper size, that is its size must be the number of
// elements in the basis.
void LinearCombination::GetCoordinates(int *coordinates, const OrderedBasis &basis) const
{
	int dim = (int)basis.size();
	assert(dim > 0);

	// Clear the coordinates first so that they are 0 everywhere
	memset(coordinates, 0, dim*sizeof(int));

	// Return the zero vector if the linear combination is empty
	if (terms.empty()) {
		return;
	}

	vector<Term>::const_iterator iter_term;
	for (iter_term = terms.begin(); iter_term != terms.end(); iter_term++) {
		int i;
		for (i=0; i<dim; i++) {
			if (basis[i] == iter_term->word) {
				coordinates[i] += iter_term->coeff;
				break;
			}
		}
		assert(i<dim); // Assert if the term is not in the basis
	}
}


void LinearCombination::Simplify()
{
	// TODO
}

string LinearCombination::OutputString() const
{
	string output;
	vector<Term>::const_iterator iter;

	// If the linear combination is empty, that means it equals zero
	if (terms.empty()) {
		output = "0";
	}

	int i = 0;
	for (iter = terms.begin(); iter != terms.end(); iter++) {
		if (i == 0 && iter->coeff == -1) {
			output += "-";
		} else if (i > 0 && iter->coeff == -1) {
			output += " - ";
		} else if (i > 0) {
			output += " + ";
		}

		if (iter->coeff != 1 && iter->coeff != -1) {
			stringstream ss;
			ss << iter->coeff;
			output += "(" + ss.str() + ")";
		}

		output += iter->word.OutputString();
		i++;
	}

	return output;
}

/*istream& operator>>(istream&, LinearCombination&)
{
}

ostream& operator<<(ostream&, const LinearCombination&)
{
}*/

Differential::Differential()
{
	SetDifferentialToZero("1");
}

// Define the differential of a generator as "a"
void Differential::SetDifferential(const string &generator_label, const LinearCombination &a)
{
	differential[generator_label] = a;
}

void Differential::SetDifferentialToZero(const string &generator_label)
{
	LinearCombination a;
	differential[generator_label] = a;
}

// Evaluate the differential on a word and store the result in "result"
void Differential::EvaluateDifferential(LinearCombination &result, const Word &word)
{
	// First, make the linear combination zero
	result.MakeZero();

	Generator first_factor;
	Word remaining_factors;
	word.GetFirstFactor(first_factor, remaining_factors);
	
	// If "word" consists of only one factor, we know the result (it was given to us by the user)
	// If "word" has length 2 or greater, we need to apply Leibniz' rule
	if (word.GetLength() <= 1) {
		map<string, LinearCombination>::iterator iter = differential.find(first_factor.label);
		if (iter == differential.end()) {
			throw logic_error("There is no differential defined for generator '" + first_factor.label + "'.");
		}
		result = iter->second;
	}
	if (word.GetLength() > 1) {
		LinearCombination first_term;
		// The first term is the product of the differential of the first factor and of the remaining factors
		map<string, LinearCombination>::iterator iter = differential.find(first_factor.label);
		if (iter == differential.end()) {
			throw logic_error("There is no differential defined for generator '" + first_factor.label + "'.");
		}
		first_term = iter->second;
		first_term.MultiplyOnRight(remaining_factors);

		// Next, we need to add the product "first_factor" and then the differential of "remaining_factors"
		LinearCombination remaining_factors_result;
		EvaluateDifferential(remaining_factors_result, remaining_factors);
		
		if (!isEven(first_factor.degree)) {
			remaining_factors_result.ScalarMultiply(-1);
		}
		remaining_factors_result.MultiplyOnLeft(first_factor);

		result.AddTerms(first_term); // First term in Leibniz' rule
		result.AddTerms(remaining_factors_result); // Second term in Leibniz' rule
	}
}

// The differential matrix must be of size (dim_target) x (dim_source).
// The first index should represent the column number and the second index should represent the row number.

void Differential::ComputeDifferentialMatrix(int **differential_matrix, const OrderedBasis &source, const OrderedBasis &target)
{
	// Compute the differential on each element of source
	int dim_source = (int)source.size();
	int dim_target = (int)target.size();
	assert(dim_source > 0 && dim_target > 0);

	LinearCombination result;

	for (int i=0; i<dim_source; i++) {
		EvaluateDifferential(result, source[i]);
		int *coords = new int[dim_target];
		result.GetCoordinates(coords, target);
		for (int j=0; j<dim_target; j++) {
			differential_matrix[j][i] = coords[j];
		}
		delete [] coords;
	}
}


//// This part of the code deals with input/output from and to files ////

string trim(const string &str)
{
	size_t start = str.find_first_not_of(" \t\n\r");
	if(start == string::npos)
		return "";
	return str.substr(start, str.find_last_not_of(" \t\n\r") - start + 1);
} 

bool caseInsensitiveStringCompare(const string& str1, const string& str2)
{
    string str1Cpy(str1);
    string str2Cpy(str2);
    transform(str1Cpy.begin(), str1Cpy.end(), str1Cpy.begin(), ::tolower);
    transform(str2Cpy.begin(), str2Cpy.end(), str2Cpy.begin(), ::tolower);
    return ( str1Cpy == str2Cpy );
}

// Use this to read or write a generator to a file
istream& operator>>(istream &stream, Generator &g)
{
	stream >> g.label;
	stream >> g.degree;
	return stream;
}

ostream& operator<<(ostream &stream, const Generator &g)
{
	stream << g.label << " " << g.degree;
	return stream;
}

// Use this to read or write a word to a file
istream& operator>>(istream &stream, Word &word)
{
	string factor;
	string line;
	
	word.Clear();
	
	// Only read one line
	getline(stream, line);

	// Make the line into a stream
	stringstream linestream(line);

	while (linestream && !linestream.eof()) {
		getline(linestream, factor, '*');
		factor = trim(factor);
		int pos = factor.find("^");
		if (pos == string::npos) {
			word.AddPowerOfGenerator(factor, GradedVectorSpace::GetGeneratorDegree(factor), 1);
		} else {
			stringstream ss(trim(factor.substr(pos+1)));
			int power;
			if (!(ss >> power)) {
				cerr << "Can't read factor '" << factor << "'." << endl;
				linestream.setstate(ios::failbit);
				return linestream;
			}
			string label = trim(factor.substr(0, pos));
			word.AddPowerOfGenerator(label, GradedVectorSpace::GetGeneratorDegree(label), power);
		}
	}
	return stream;
}

ostream& operator<<(ostream &stream, const Word &word)
{
	return stream << word.OutputString();
}

// Use this to read or write a linear combination to a file
// For inptu, we assumed the string is trimmed.
istream& operator>>(istream &stream, LinearCombination &lc)
{
	string line;

	lc.MakeZero();

	// Only read one line
	string tmp;
	getline(stream, tmp);

	// Remove all whitespace and store into "line"
	for (size_t i=0; i<tmp.size(); i++) {
		if (tmp[i] != ' ' && tmp[i] != '\t')
			line += tmp[i];
	}

	if (line.empty()) {
		stream.setstate(ios::failbit);
		return stream;
	}

	if (line == "0") {
		return stream;
	}

	int term_pos = 0;
	int next_term_pos = 0;

	do {
		int plus_symbol_pos;
		int minus_symbol_pos;

		next_term_pos = term_pos;

		do {
			plus_symbol_pos = line.find("+", next_term_pos+1);
			minus_symbol_pos = line.find("-", next_term_pos+1);

			if (plus_symbol_pos == string::npos && minus_symbol_pos != string::npos) {
				next_term_pos = minus_symbol_pos;
			} else if (plus_symbol_pos != string::npos && minus_symbol_pos == string::npos) {
				next_term_pos = plus_symbol_pos;
			} else {
				next_term_pos = min(minus_symbol_pos, plus_symbol_pos);
			}
			// Keep looping if the minus sign found was the leading coefficient of the term we are looking for
			if (next_term_pos > 0) {
				if (line[next_term_pos - 1] != '(') {
					break;
				}
			} else if (next_term_pos == string::npos) {
				break;
			} else {
				next_term_pos++;
			}
		} while (true);

		string str = line.substr(term_pos, next_term_pos-term_pos);

		int coeff;
		Word word;
		stringstream str_word;

		if (str.compare(0, 1, "-") == 0) {
			coeff = -1;
		} else {
			coeff = 1;
		}

		if (str.compare(0, 1, "(") == 0 || str.compare(1, 1, "(") == 0) {
			// This is a term with the coefficient between parentheses
			size_t begin_pos = 0;
			if (str.compare(1, 1, "(") == 0)
				begin_pos = 1;
			size_t pos = str.find_first_of(")");
			stringstream ss(str.substr(begin_pos+1, pos-(begin_pos+1)));
			int tmp;
			if (!(ss >> tmp)) {
				// Failed to read a coefficient
				cerr << "Could not read a a coefficient ('" << ss.str() << "') in a linear combination." << endl;
				stream.setstate(ios::failbit);
				return stream;
			}
			coeff *= tmp;
			str_word << str.substr(pos+1);
		} else if (str.compare(0, 1, "-") == 0 || str.compare(0, 1, "+") == 0) {
			// This term starts with a + or - symbol
			str_word << str.substr(1);
		} else {
			// This is the first term so it might not have any operator in front of it
			str_word << str;
		}

		if (!(str_word >> word)) {
			cerr << "Could not read a word ('" << str_word.str() << "') in a linear combination." << endl;
			stream.setstate(ios::failbit);
			return stream;
		}

		// If we reached this point, the word was read correctly
		lc.AddTerm(coeff, word);
		term_pos = next_term_pos;
	} while (term_pos != string::npos);

	return stream;
}

ostream& operator<<(ostream &stream, const LinearCombination &lc)
{
	return stream << lc.OutputString();
}


// Return 'true' if the file was parsed successfully.
// A FreeCGA and a Differential object will be returned.
enum INPUT_STAGE
{
	NONE, GENERATORS, EXTENSION, DIFFERENTIAL, OUTPUT
};

bool ReadInputFromFile(const string &filename, string &output_filename, string &extension_output_filename, FreeCGA &cdga, Differential &differential, int &homology_degree_start, int &homology_degree_end, int &category)
{
	string line;
	ifstream file(filename);
	GradedVectorSpace X;
	int line_number = 0;

	category = -1; // By default, the category will be set to -1, so words of length 0 and up (i.e. everything) will be considered
	homology_degree_start = -1;
	homology_degree_end = -1;

	if (file.is_open()) {
		INPUT_STAGE stage = NONE;
		while (file.good()) {
			line_number++;
			getline(file, line);
			line = trim(line);
			stringstream ssline(line);
			if (line.empty() || line.find("#") == 0) {
				// Empty line or comment line, just skip it
			} else if (caseInsensitiveStringCompare(line, "generators:")) {
				stage = GENERATORS;
			} else if (caseInsensitiveStringCompare(line, "extension:")) {
				stage = EXTENSION;
			} else if (caseInsensitiveStringCompare(line, "differential:")) {
				stage = DIFFERENTIAL;
			} else if (caseInsensitiveStringCompare(line, "output:")) {
				stage = OUTPUT;
			} else {
				Generator g;
				switch (stage) {
				case GENERATORS:
					if (ssline >> g) {
						X.AddGenerator(g.label, g.degree);
					} else {
						cerr << "Couldn't read generator on line " << line_number << "." << endl;
						return false;
					}
					break;
				case EXTENSION:
					if (ssline >> g) {
						cdga.AddExtensionGenerator(g.label, g.degree);
					} else {
						cerr << "Couldn't read module extension generator on line " << line_number << "." << endl;
						return false;
					}
					break;
				case DIFFERENTIAL:
					if (line.find("d(") == 0) {
						size_t pos = line.find_first_of(")");
						if (pos != string::npos) {
							string generator_label = trim(line.substr(2, pos-2));
							string str_diff;
							LinearCombination lc;
							pos = line.find_first_of("=");
							str_diff = trim(line.substr(pos+1));
							stringstream ss(str_diff);
							if (ss >> lc) {
								differential.SetDifferential(generator_label, lc);
								break;
							}
						}
					}
					cerr << "Couldn't read differential on line " << line_number << "." << endl;
					return false;
				case OUTPUT:
					if (line.compare(0, strlen("degree"), "degree") == 0) {
						size_t pos = line.find("=");
						size_t range_pos = line.find("..");
						stringstream ss(trim(line.substr(pos+1)));
						if (range_pos != string::npos) {
							stringstream deg_start(trim(line.substr(pos+1, range_pos-(pos+1))));
							stringstream deg_end(trim(line.substr(range_pos+2)));
							if (!(deg_start >> homology_degree_start)) {
								cerr << "'Start of range degree 'for computation of a basis in homology invalid." << endl;
								return false;
							}
							if (!(deg_end >> homology_degree_end)) {
								cerr << "'End of range degree' for computation of a basis in homology invalid." << endl;
								return false;
							}
						} else if (!(ss >> homology_degree_start)) {
							cerr << "Degree for computation of a basis in homology invalid." << endl;
							return false;
						}
					} else if (line.compare(0, strlen("category"), "category") == 0) {
						size_t pos = line.find("=");
						stringstream ss(trim(line.substr(pos+1)));
						if (!(ss >> category)) {
							cerr << "The specified category of the space is invalid." << endl;
							return false;
						}
					} else if (line.compare(0, strlen("filename"), "filename") == 0) {
						size_t pos = line.find("=");
						stringstream ss(trim(line.substr(pos+1)));
						ss >> output_filename;
					} else if (line.compare(0, strlen("extension-output"), "extension-output") == 0) {
						size_t pos = line.find("=");
						stringstream ss(trim(line.substr(pos+1)));
						ss >> extension_output_filename;
					}
					break;
				};
			}
		}
	} else {
		cerr << "Unable to open file '" << filename << "'." << endl;
		return false;
	}
	cdga.SetGradedVectorSpace(X);
	if (homology_degree_start < 0) {
		cerr << "Degree for computation of a basis in homology not specified or invalid." << endl;
		return false;
	}
	return true;
}