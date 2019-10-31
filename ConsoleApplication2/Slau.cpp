#include "Slau.h"

Slau::Slau()
{
	// открываем файл
	ifstream kuslauf("kuslau.txt");;
	int u = 0;

	// ограничитель
	if (!kuslauf.is_open()) throw "Не удалось открыть файл ";

	//макс кол-во итераций и макс невязка
	kuslauf >> maxiter;
	kuslauf >> eps;
	
	A.K.input();
	n = A.K.v_num;
	di.resize(n);
	f.resize(n);
	y.resize(n);
	x0.resize(n);
	x.resize(n);
	r.resize(n);
	z.resize(n);
	p.resize(n);

	for (int i = 0; i < n; i++)
		x[i] = 0;
}

double Slau::scalarMult(vector<double> a, vector<double> b)
{
	double scalar = 0;
	for (int i = 0; i < n; i++)
	{
		scalar += a[i] * b[i];
	}
	return scalar;
}
//умножение матрицы на вектор
void Slau::mult(vector<double> x)
{
	//по строкам нижнего и столбцам верхнего треугольников
	for (int i = 0; i < n; i++)
	{
		y[i] = di[i] * x[i];
		int i0 = ig[i];
		int i1 = ig[i + 1];
		for (int k = i0; k < i1; k++)
		{
			int j = jg[k];
			y[i] += ggl[k] * x[j];
			y[j] += ggu[k] * x[i];
		}
	}
}

double Slau::calcNev(vector<double> r)
{
	double rnorm = 0, fnorm = 0;
	nev = 0;
	rnorm = scalarMult(r, r);
	fnorm = scalarMult(f, f);
	nev = sqrt(rnorm / fnorm);
	return nev;
}

void Slau::fwrd(vector<real> &yy) {
	for (int k = 0; k < n; k++)
	{
		double s = 0.0;
		int ig_k = ig[k], ig_k_1 = ig[k + 1];
		for (int i = ig_k; i < ig_k_1; i++)
		{
			int j = jg[i];
			s += L[i] * yy[j];
		}
		yy[k] = (yy[k] - s) / D[k];
	}
}

void Slau::bwrd(vector<real> &yy) {
	for (int k = n - 1; k >= 0; k--)
	{
		double yk = yy[k] = yy[k]; // D[k];
		int ig_k = ig[k], ig_k_1 = ig[k + 1];
		for (int j = ig_k; j < ig_k_1; j++)
		{
			int i = jg[j];
			yy[i] -= U[j] * yk;
		}
	}
}

void Slau::LU() {
	A.GetGlobal();
	ggl = A.ggl; ggu = A.ggu; di = A.di;
	L = ggl; D = di; U = ggu; ig = A.ig; jg = A.jg; 
	f = A.global_right;
	normF = scalarMult(f, f);
	for (int i = 0; i < A.K.v_num; i++) {
		int ib = A.ig[i];
		int iend = A.ig[i + 1];
		double diag = 0;
		for (int k = ib, j; k < iend; k++) {
			double suml = 0; double sumu = 0;
			j = A.jg[k]; int k1 = A.ig[j + 1];
			int ii = ib; int jj = A.ig[j];
			while (ii < k && jj < k1)
				if (A.jg[ii] < A.jg[jj])
					++ii;
				else if (A.jg[ii] > A.jg[jj])
					++jj;
				else
				{
					suml += L[ii] * U[jj];
					sumu += L[jj] * U[ii];
					++ii; ++jj;
				}
			L[k] -= suml;
			U[k] = (U[k] - sumu) / D[j];
			diag += L[k] * U[k];
		}
		D[i] -= diag;
	}
}

void Slau::LOS()
{
	LU();
	double rNormNew = 0;
	double nev = 0;
	mult(x);	//Ax
	for (int i = 0; i < n; i++)
		r[i] = f[i] - y[i];
	fwrd(r);
	z = r;
	bwrd(z);
	mult(z);	//Az считает в y
	p = y;		//p = Az
	fwrd(p);
	double alpha = 0;
	double betta = 0;
	double pr = 0, pp = 0;
	vector<real> rr;

	nev = calcNev(r);
	for (int k = 0; k < maxiter && nev > eps; k++)
	{
		pr = pp = 0;
		pr = scalarMult(p, r);
		pp = scalarMult(p, p);
		alpha = pr / pp;
		for (int i = 0; i < n; i++)
		{
			x[i] += alpha*z[i];
			r[i] -= alpha*p[i];
		}
		rr = r;
		bwrd(rr);
		mult(rr);
		fwrd(y);
		pr = scalarMult(p, y);
		betta = -pr / pp;
		for (int i = 0; i < n; i++)
		{
			z[i] = r[i] + betta*z[i];
			p[i] = y[i] + betta*p[i];
		}
		nev = calcNev(r);
		iter = k;
		printf("\r Итерация %d: %.2le", k, nev);
	}
}