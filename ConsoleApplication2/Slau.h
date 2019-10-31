#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <locale.h>
#include "finite.h"

using namespace std;

typedef double real;

class Slau
{
public:
	vector<int> ig;
	vector<int> jg;
	vector<real> ggl;
	vector<real> ggu;
	vector<real> di;
	vector<real> f;
	vector<real> y;
	vector<real> x0;
	vector<real> x;
	vector<real> r;
	vector<real> z;
	vector<real> p;
	vector<real>D;
	vector<real>L;
	vector<real>U;
	int sig, sjg, sggl, sggu;
	int n;
	int maxiter;
	Matrix A;
	double eps;
	int iter;
	real nev;
	real normF;

public:
	Slau();
	void mult(vector<real> x);
	real scalarMult(vector<real> a, vector<real> b);
	real calcNev(vector<real> r);
	void bwrd(vector<real> &y);
	void fwrd(vector<real> &y);
	void LOS();
	void LU();
};

