#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 
#include <locale.h>
#include "Slau.h"

using namespace std;
#define real double

int main()
{
	setlocale(LC_CTYPE, "Russian");
	ofstream os("result.txt");
	Slau Axb;
	Axb.iter = 0;
	Axb.normF = 0;
	Axb.LOS();
	//����� ���������� � ����
	os.precision(14);
	for (int i = 0; i < Axb.n; i++)
		os << Axb.x[i] << endl;
	os << "\n" << endl;
	os << Axb.nev << endl;
	os << "\n" << endl;
	//system("pause");
}
