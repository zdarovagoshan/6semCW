#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <iterator>
#include <math.h>
struct Vertex {
	double x;
	double y;
	int index;
};
struct Elem {
	int vertex[3];
	int fe;
};
struct KR1 {
	int index;
	int ug;
};
struct KR2 {
	int first_index;
	int second_index;
	int teta;
};
struct KR3 {
	int first_index;
	int second_index;
	int beta;
	int ub;
};
class FE {
public:
	std::vector<Vertex>vertices;
	std::vector<Elem>fe;
	std::vector<KR1>kr1;
	std::vector<KR2>kr2;
	std::vector<KR3>kr3;
	void input();
	void GetLinkList(std::vector<std::vector<int>>&);
	double GetD1(std::vector<std::vector<double>>&, int&);
	double lymda(int&);
	double gamma(int&);
	double function(int&, int&);
	double ug(int&);
	double teta(int&, int);
	double beta(int&);
	double ub(int&, int);
	int GetNumOfGlobal(int&);
	int GetGlobalIndex(int&, int&);
	int v_num;
	int fe_num;
};
class Matrix {
public:
	FE K;
	void GeneratePortret();
	void GetLocal(int&);
	void LocalToGlobal(int&);
	void GetGlobal();
	void allow_boundary_1();
	void allow_boundary_2();
	void allow_boundary_3();
	std::vector<std::vector<double>>local_matrix;
	std::vector<double>local_right;
	std::vector<double>global_right;
	std::vector<double>di;
	std::vector<double>ggl;
	std::vector<double>ggu;
	std::vector<int>ig;
	std::vector<int>jg;
};