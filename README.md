# cdga-generators
This program computes various topological invariants, such as the cohomology of various cdga (commutative differential graded algebra) models in Rational Homotopy Theory. In particular, this program allows tedious computations of the rational retraction index to be automated.

# Purpose
I wrote this program in 2001 during my master's thesis to help me perform computations that would have been extremely tedious, time-consuming and error-prone. In fact, for anything but the simplest topological spaces admitting a minimal Sullivan model (an algebraic counterpart to the topological space), the computations quickly become extremely complex as an intermediate model needs to be computed and new generators are added for each degree d=0 to d=cat(X). So, unless the initial model is very small or the rational LS category is very small, the intermediate model will almost always become huge. This is the reason why a computer program is almost mandatory for most spaces in order to compute this particular rational homotopy invariant (the rational retraction index).

# More details
Much more documentation about the program behaviors, such as the exact input expected and explanations of the input's meaning, as well as a list of all the output produced and explanations of the data produced, limitations of the program, planned features, etc. can be found at the top of the file "main.cpp".

THe program should be considered research code and may contain many undocumented bugs. It is recommended to read the source code to understand exactly what the program is doing. I made an effort to document the entire program so that another person can make sense of it.

# Usage
There are two ways to run the program. The first one is to pass the input file as the first argument, for example:
```
./cdga-generators.exe input_example1.txt
```
Alternatively, you can just run the program with no argument
```
./cdga-generators.exe
```
and you will be prompted for the name of an input file.

# Dependencies
There are only one dependency for compiling cdga-generators, namely NTL, a number theory library which I use exclusively for its LLL algorithm. The website of NTL says "NTL is a high-performance, portable C++ library providing data structures and algorithms for manipulating signed, arbitrary length integers, and for vectors, matrices, and polynomials over the integers and over finite fields.".

See http://www.shoup.net/ntl/ to download NTL.

# Platforms supported
Unfortunately, the cdga-generators currently makes a few Windows system calls and as such will only compile on Windows. It would be very easy to replace those calls with a portable alternative and I would be happy to do so if anyone contacts me to request a Linux or Mac port. I expected the task to be very easy and quick.

# Input file example 1

## Example 1 from p.439 of [FHT]
## This is an elliptic model, hence top class degree: 6+4+2+2 - (1+3) = 10

Generators:
a 2
x 3
u 3
b 4
v 5
w 7

Extension:

Differential:
d(a) = 0
d(x) = 0
d(u) = a^2
d(b) = a * x
d(v) = a * b - u * x
d(w) = b^2 - v * x

Output:
filename = output.txt
extension-output = output-extension.txt
degree = 2..30
category = -1

# Input file example 2

## This is for the example 2 of [GJ]
## Category for this space is: 5
## This space is elliptic
## Top class is in degree: 15 - 2 = 13

Generators:
a 2
b 2
y_1 5
y_2 5
y_3 5

Extension:

Differential:
d(a) = 0
d(b) = 0
d(y_1) = a^3
d(y_2) = b^3
d(y_3) = a^2 * b

Output:
filename = output-gj2-output.txt
extension-output = output-gj2-extension.txt
degree = 2..29
category = 5

# Example 3

## This is for the space of exercise #1, p.405 of [FHT], M = Sp(5) / SU(5)

Generators:
a 6
b 10
x 11
y 15
z 19

Extension:

Differential:
d(a) = 0
d(b) = 0
d(x) = a^2
d(y) = a * b
d(z) = b^2

Output:
filename = output-sp5-su5.txt
extension-output = output-sp5-su5-extension.txt
degree = 2..15
category = -1

# Input file example 4

## This is the formal model (including only generators up to degree 10) of the cohomology
## of the counter-example #3

## Since cup-length = 3, it follows that category = 3
## The conjecture predicts a retraction index of 2.

Generators:
a 2
b 2
c 3
e 5
f 5
g 6
h 6
i 7
j 7
k 7
l 7
m 8
n 8
p 8
q 8
r_1 9
r_2 9
r_3 9
r_4 9
r_5 9
s_1 10
s_2 10
s_3 10
s_4 10
s_5 10
s_6 10
s_7 10
s_8 10
s_9 10

Differential:
d(a) = 0
d(b) = 0
d(c) = a * b
d(e) = a^3
d(f) = b^3
d(g) = -b * e + a^2 * c
d(h) = -a * f + b^2 * c
d(i) = 0
d(j) = 0
d(k) = a * g + c * e
d(l) = b * h + c * f
d(m) = a * i
d(n) = b * j
d(p) = -c * g + b * k
d(q) = -c * h + a * l
d(r_1) = -b * m + c * i
d(r_2) = -a * n + c * j
d(r_3) = -a * p + c * k
d(r_4) = -b * q + c * l
d(r_5) = -e * f - a^2 * h + b^2 * g
d(s_1) = c * m + a * r_1
d(s_2) = c * n + b * r_2
d(s_3) = c * p + b * r_3
d(s_4) = c * q + a * r_4
d(s_5) = -f * g + b * r_5 + a * c * h
d(s_6) = -e * h - a * r_5 + b * c * g
d(s_7) = a^2 * j
d(s_8) = -e * g + a^2 * k
d(s_9) = -f * h + b^2 * l

Extension:
v_1 7
v_2 7
v_3 7
v_4 7
v_5 7
w_1 8
w_2 8
w_3 8
w_4 8
w_5 8
w_6 8
w_7 8
w_8 8

Differential:
d(v_1) = a^4
d(v_2) = a^3 * b
d(v_3) = a^2 * b^2
d(v_4) = a * b^3
d(v_5) = b^4
d(w_1) = -a^3 * c + b * v_1
d(w_2) = -a^3 * c + a * v_2
d(w_3) = -a^2 * b * c + b * v_2
d(w_4) = -a^2 * b * c + a * v_3
d(w_5) = -a * b^2 * c + b * v_3
d(w_6) = -a * b^2 * c + a * v_4
d(w_7) = -b^3 * c + b * v_4
d(w_8) = -b^3 * c + a * v_5

## No perturbation on v_1, ..., v_5

## D(d(w_1)) = D(-a^3 * c + b * v_1) = -a^4 * b + b * a^4 = 0

## rho(d(v_1)) = rho(a^4) = a^4
## rho(v_1) = a * e

## Cannot be avoided, so r_0 <= 2.

Output:
filename = junk.txt
extension-output = extension.txt
degree = 10..10
category = 3
compute = extension

# Input file example 5

## This is for the product of Sp(3) biquotient with itself

## Category: 6
## Top class is in degree: 24
## Retraction index is greater or equal to 4.

Generators:
a 4
b 4
c 4
e 4
x 7
y 11
v 7
w 11

Extension:
s_001 27
s_002 27
s_003 27
s_004 27
s_005 27
s_006 27
s_007 27
s_008 27
s_009 27
s_010 27
s_011 27
s_012 27
s_013 27
s_014 27
s_015 27
s_016 27
s_017 27
s_018 27
s_019 27
s_020 27
s_021 27
s_022 27
s_023 27
s_024 27
s_025 27
s_026 27
s_027 27
s_028 27
s_029 27
s_030 27
s_031 27
s_032 27
s_033 27
s_034 27
s_035 27
s_036 27
s_037 27
s_038 27
s_039 27
s_040 27
s_041 27
s_042 27
s_043 27
s_044 27
s_045 27
s_046 27
s_047 27
s_048 27
s_049 27
s_050 27
s_051 27
s_052 27
s_053 27
s_054 27
s_055 27
s_056 27
s_057 27
s_058 27
s_059 27
s_060 27
s_061 27
s_062 27
s_063 27
s_064 27
s_065 27
s_066 27
s_067 27
s_068 27
s_069 27
s_070 27
s_071 27
s_072 27
s_073 27
s_074 27
s_075 27
s_076 27
s_077 27
s_078 27
s_079 27
s_080 27
s_081 27
s_082 27
s_083 27
s_084 27
s_085 27
s_086 27
s_087 27
s_088 27
s_089 27
s_090 27
s_091 27
s_092 27
s_093 27
s_094 27
s_095 27
s_096 27
s_097 27
s_098 27
s_099 27
s_100 27
s_101 27
s_102 27
s_103 27
s_104 27
s_105 27
s_106 27
s_107 27
s_108 27
s_109 27
s_110 27
s_111 27
s_112 27
s_113 27
s_114 27
s_115 27
s_116 27
s_117 27
s_118 27
s_119 27
s_120 27

Differential:
d(a) = 0
d(b) = 0
d(c) = 0
d(e) = 0
d(x) = a^2 + a * b + b^2
d(y) = 0
#d(y) = a^3
d(v) = c^2 + c * e + e^2
d(w) = 0
#d(w) = c^3

d(s_001) = a^7
d(s_002) = a^6 * b
d(s_003) = a^6 * c
d(s_004) = a^6 * e
d(s_005) = a^5 * b^2
d(s_006) = a^5 * b * c
d(s_007) = a^5 * b * e
d(s_008) = a^5 * c^2
d(s_009) = a^5 * c * e
d(s_010) = a^5 * e^2
d(s_011) = a^4 * b^3
d(s_012) = a^4 * b^2 * c
d(s_013) = a^4 * b^2 * e
d(s_014) = a^4 * b * c^2
d(s_015) = a^4 * b * c * e
d(s_016) = a^4 * b * e^2
d(s_017) = a^4 * c^3
d(s_018) = a^4 * c^2 * e
d(s_019) = a^4 * c * e^2
d(s_020) = a^4 * e^3
d(s_021) = a^3 * b^4
d(s_022) = a^3 * b^3 * c
d(s_023) = a^3 * b^3 * e
d(s_024) = a^3 * b^2 * c^2
d(s_025) = a^3 * b^2 * c * e
d(s_026) = a^3 * b^2 * e^2
d(s_027) = a^3 * b * c^3
d(s_028) = a^3 * b * c^2 * e
d(s_029) = a^3 * b * c * e^2
d(s_030) = a^3 * b * e^3
d(s_031) = a^3 * c^4
d(s_032) = a^3 * c^3 * e
d(s_033) = a^3 * c^2 * e^2
d(s_034) = a^3 * c * e^3
d(s_035) = a^3 * e^4
d(s_036) = a^2 * b^5
d(s_037) = a^2 * b^4 * c
d(s_038) = a^2 * b^4 * e
d(s_039) = a^2 * b^3 * c^2
d(s_040) = a^2 * b^3 * c * e
d(s_041) = a^2 * b^3 * e^2
d(s_042) = a^2 * b^2 * c^3
d(s_043) = a^2 * b^2 * c^2 * e
d(s_044) = a^2 * b^2 * c * e^2
d(s_045) = a^2 * b^2 * e^3
d(s_046) = a^2 * b * c^4
d(s_047) = a^2 * b * c^3 * e
d(s_048) = a^2 * b * c^2 * e^2
d(s_049) = a^2 * b * c * e^3
d(s_050) = a^2 * b * e^4
d(s_051) = a^2 * c^5
d(s_052) = a^2 * c^4 * e
d(s_053) = a^2 * c^3 * e^2
d(s_054) = a^2 * c^2 * e^3
d(s_055) = a^2 * c * e^4
d(s_056) = a^2 * e^5
d(s_057) = a * b^6
d(s_058) = a * b^5 * c
d(s_059) = a * b^5 * e
d(s_060) = a * b^4 * c^2
d(s_061) = a * b^4 * c * e
d(s_062) = a * b^4 * e^2
d(s_063) = a * b^3 * c^3
d(s_064) = a * b^3 * c^2 * e
d(s_065) = a * b^3 * c * e^2
d(s_066) = a * b^3 * e^3
d(s_067) = a * b^2 * c^4
d(s_068) = a * b^2 * c^3 * e
d(s_069) = a * b^2 * c^2 * e^2
d(s_070) = a * b^2 * c * e^3
d(s_071) = a * b^2 * e^4
d(s_072) = a * b * c^5
d(s_073) = a * b * c^4 * e
d(s_074) = a * b * c^3 * e^2
d(s_075) = a * b * c^2 * e^3
d(s_076) = a * b * c * e^4
d(s_077) = a * b * e^5
d(s_078) = a * c^6
d(s_079) = a * c^5 * e
d(s_080) = a * c^4 * e^2
d(s_081) = a * c^3 * e^3
d(s_082) = a * c^2 * e^4
d(s_083) = a * c * e^5
d(s_084) = a * e^6
d(s_085) = b^7
d(s_086) = b^6 * c
d(s_087) = b^6 * e
d(s_088) = b^5 * c^2
d(s_089) = b^5 * c * e
d(s_090) = b^5 * e^2
d(s_091) = b^4 * c^3
d(s_092) = b^4 * c^2 * e
d(s_093) = b^4 * c * e^2
d(s_094) = b^4 * e^3
d(s_095) = b^3 * c^4
d(s_096) = b^3 * c^3 * e
d(s_097) = b^3 * c^2 * e^2
d(s_098) = b^3 * c * e^3
d(s_099) = b^3 * e^4
d(s_100) = b^2 * c^5
d(s_101) = b^2 * c^4 * e
d(s_102) = b^2 * c^3 * e^2
d(s_103) = b^2 * c^2 * e^3
d(s_104) = b^2 * c * e^4
d(s_105) = b^2 * e^5
d(s_106) = b * c^6
d(s_107) = b * c^5 * e
d(s_108) = b * c^4 * e^2
d(s_109) = b * c^3 * e^3
d(s_110) = b * c^2 * e^4
d(s_111) = b * c * e^5
d(s_112) = b * e^6
d(s_113) = c^7
d(s_114) = c^6 * e
d(s_115) = c^5 * e^2
d(s_116) = c^4 * e^3
d(s_117) = c^3 * e^4
d(s_118) = c^2 * e^5
d(s_119) = c * e^6
d(s_120) = e^7

Output:
filename = output-test3.txt
extension-output = output-test3-extension.txt
degree = 29
category = 6
