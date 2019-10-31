#include "finite.h"
int FE::GetGlobalIndex(int& ii, int& jj) {
	int k1 = fe[ii].vertex[0],
		k2 = fe[ii].vertex[1],
		k3 = fe[ii].vertex[2],
		flag = 0;
	for (int i = 0; i < v_num; i++) {
		int vi = vertices[i].index;
		if (vi == k1 || vi == k2 || vi == k3)
			if (flag++ == jj)return i;
	}
}
int FE::GetNumOfGlobal(int& ii) {
	return 3;
}
void FE::input() {
	std::ifstream F;
	F.open("coords.txt");
	F >> v_num;
	vertices.resize(v_num);
	for (int i = 0; i < v_num; i++) {
		F >> vertices[i].x >> vertices[i].y;
		vertices[i].index = i;
	}
	F.close();
	F.open("nvtr.txt");
	F >> fe_num;
	fe.resize(fe_num);
	for (int i = 0; i < fe_num; i++)
		F >> fe[i].vertex[0] >> fe[i].vertex[1]
		>> fe[i].vertex[2] >> fe[i].fe;
	F.close();
	KR1 k1; KR2 k2; KR3 k3;
	F.open("1ku.txt");
	while (!F.eof()) {
		F >> k1.index >> k1.ug;
		kr1.push_back(k1);
	}
	F.close();
	F.open("2ku.txt");
	while (!F.eof()) {
		F >> k2.first_index >> k2.second_index >> k2.teta;
		kr2.push_back(k2);
	}
	F.close();
	F.open("3ku.txt");
	while (!F.eof()) {
		F >> k3.first_index >> k3.second_index
			>> k3.beta >> k3.ub;
		kr3.push_back(k3);
	}
	F.close();
}
void FE::GetLinkList(std::vector<std::vector<int>>& L) {
	for (int k = 0; k < fe_num; k++)
		for (int i = 0; i < GetNumOfGlobal(k) - 1; i++) {
			int ii = GetGlobalIndex(k, i);
			for (int j = i + 1; j < GetNumOfGlobal(k); j++) {
				int ind1 = ii;
				int ind2 = GetGlobalIndex(k, j);
				if (ind2 < ind1)
					std::swap(ind1, ind2);
				if (L[ind2].end() == std::find(L[ind2].begin(), L[ind2].end(), ind1))
					L[ind2].push_back(ind1);
			}
		}
}
void Matrix::allow_boundary_1() {
	for (int i = 0; i < K.kr1.size(); i++) {
		int ii = K.kr1[i].index;
		global_right[ii] = K.ug(i);
		di[ii] = 1;
		for (int j = ig[ii]; j < ig[ii + 1]; j++)
			ggl[j] = 0;
		for (int j = 1; j < K.v_num; j++)
			for (int k = ig[j]; k < ig[j + 1]; k++)
				if (jg[k] == ii)ggu[k] = 0;

	}
}
void Matrix::allow_boundary_2() {
	for (int i = 0; i < K.kr2.size(); i++) {
		double h = sqrt(pow(K.vertices[K.kr2[i].second_index].x -
			K.vertices[K.kr2[i].first_index].x, 2) +
			pow(K.vertices[K.kr2[i].second_index].y -
				K.vertices[K.kr2[i].first_index].y, 2));
		global_right[K.kr2[i].first_index] += K.teta(i, 1)*h / 6;
		global_right[K.kr2[i].second_index] += K.teta(i, 2)*h / 6;
	}
}
void Matrix::allow_boundary_3() {
	for (int i = 0; i < K.kr3.size(); i++) {
		double h = sqrt(pow(K.vertices[K.kr3[i].second_index].x -
			K.vertices[K.kr3[i].first_index].x, 2) +
			pow(K.vertices[K.kr3[i].second_index].y -
				K.vertices[K.kr3[i].first_index].y, 2));
		di[K.kr3[i].first_index] += K.beta(i)*h / 3;
		di[K.kr3[i].second_index] += K.beta(i)*h / 3;
		for (int j = ig[K.kr3[i].second_index]; j<ig[K.kr3[i].second_index + 1]; j++)
			if (jg[j] == K.kr3[i].first_index) {
				ggu[j] += K.beta(i)*h / 6;
				ggl[j] += K.beta(i)*h / 6;
				break;
			}
		global_right[K.kr3[i].first_index] += K.beta(i)*K.ub(i, 1)*h / 6;
		global_right[K.kr3[i].second_index] += K.beta(i)*K.ub(i, 2)*h / 6;
	}
}
double FE::GetD1(std::vector<std::vector<double>>& D, int& k) {
	std::vector<int>vert;
	for (int i = 0; i < 3; i++)
		vert.push_back(GetGlobalIndex(k, i));
	double det = (vertices[vert[1]].x - vertices[vert[0]].x)*
		(vertices[vert[2]].y - vertices[vert[0]].y) -
		(vertices[vert[2]].x - vertices[vert[0]].x)*
		(vertices[vert[1]].y - vertices[vert[0]].y);
	D[0][0] = (vertices[vert[1]].y - vertices[vert[2]].y) / det;
	D[1][0] = (vertices[vert[2]].y - vertices[vert[0]].y) / det;
	D[2][0] = (vertices[vert[0]].y - vertices[vert[1]].y) / det;
	D[0][1] = (vertices[vert[2]].x - vertices[vert[1]].x) / det;
	D[1][1] = (vertices[vert[0]].x - vertices[vert[2]].x) / det;
	D[2][1] = (vertices[vert[1]].x - vertices[vert[0]].x) / det;
	vert.clear();
	return det;
}
double FE::lymda(int& k) {
	switch (fe[k].fe) {
	
	case 1:return 3 * 10;
	case 2:return 3 * 1;
	/*
	case 1:return 3 * 1;
	case 2:return 3 * 1;
	*/
	}
}
double FE::gamma(int& k) {
	switch (fe[k].fe) {
	
	case 1:return 0;
	case 2:return 0;
	
	/*case 1:return 5;
	case 2:return 0;
	*/
	}
}
double FE::function(int& k, int& i) {
	switch (fe[k].fe) {
	
	case 1:return -20 * 4;
	case 2:return 0;
	

	
	/*case 1:if (i == 0)return 2 * (5 * vertices[fe[k].vertex[0]].x + 30 * vertices[fe[k].vertex[0]].y - 10)
		+ (5 * vertices[fe[k].vertex[1]].x + 30 * vertices[fe[k].vertex[1]].y - 10)
		+ (5 * vertices[fe[k].vertex[2]].x + 30 * vertices[fe[k].vertex[2]].y - 10);
		   else if (i == 1)return 2 * (5 * vertices[fe[k].vertex[1]].x + 30 * vertices[fe[k].vertex[1]].y - 10)
			   + (5 * vertices[fe[k].vertex[0]].x + 30 * vertices[fe[k].vertex[0]].y - 10)
		+ (5 * vertices[fe[k].vertex[2]].x + 30 * vertices[fe[k].vertex[2]].y - 10);
		   else
			   return 2 * (5 * vertices[fe[k].vertex[2]].x + 30 * vertices[fe[k].vertex[2]].y - 10)
		+ (5 * vertices[fe[k].vertex[1]].x + 30 * vertices[fe[k].vertex[1]].y - 10)
		+ (5 * vertices[fe[k].vertex[0]].x + 30 * vertices[fe[k].vertex[0]].y - 10);
	case 2:return 0;
	*/
	}
}
double FE::ug(int& k) {
	return pow(vertices[kr1[k].index].y, 2);
	//return 6 * vertices[kr1[k].index].y + 2;
}
double FE::teta(int& k, int i) {
	switch (kr2[k].teta) {
	
	case 1:return 20 * 3;
	case 2:return 0;
	

	/*case 1:return -6 * 3;
	case 2:return -1 * 3;
	case 3:return 6 * 3;*/
	}
}
double FE::beta(int& k) {
	return 2;

	/*return 10;*/
}
double FE::ub(int& k, int i) {
	
	if (i == 1)return 2 * (20 * vertices[kr3[k].first_index].y - 27) +
		20 * vertices[kr3[k].second_index].y - 27;//-20y+27;
	else return 20 * vertices[kr3[k].first_index].y - 27 +
				2 * (20 * vertices[kr3[k].second_index].y - 27);
	

	//if (i == 1)return 2 * (6 * vertices[kr3[k].first_index].y + 2.1) +
	//	(6 * vertices[kr3[k].second_index].y + 2.1);//6y+ 2.1;
	//else return (6 * vertices[kr3[k].first_index].y + 2.1) +
	//		   2 * (6 * vertices[kr3[k].second_index].y + 2.1);
}
void Matrix::GeneratePortret() {
	std::vector<std::vector<int>>L(K.v_num);
	K.GetLinkList(L);
	ig.push_back(0);
	for (int i = 0; i < K.v_num; i++) {
		std::sort(L[i].begin(), L[i].end());
		jg.resize(jg.size() + L[i].size());
		std::copy(L[i].begin(), L[i].end(), jg.begin() + ig[i]);
		ig.push_back(ig[i] + L[i].size());
	}
	global_right.resize(K.v_num);
	di.resize(K.v_num);
	ggl.resize(ig[K.v_num]);
	ggu.resize(ig[K.v_num]);
	L.clear();
}
void Matrix::GetLocal(int& k) {
	std::vector<std::vector<double>>D(3, std::vector<double>(2));
	double det = fabs(K.GetD1(D, k)), lyambda, y, G, M;
	lyambda = K.lymda(k);
	y = K.gamma(k);
	for (int i = 0; i < K.GetNumOfGlobal(k); i++) {
		for (int j = 0; j < K.GetNumOfGlobal(k); j++) {
			G = (D[i][0] * D[j][0] + D[i][1] * D[j][1])*lyambda*det / 6;
			M = y*det / 24;
			if (i == j)M *= 2;
			local_matrix[i][j] = G + M;
		}
		local_right[i] = K.function(k, i)*det / 24;
	}
	D.clear();
}
void Matrix::LocalToGlobal(int& k) {
	std::vector<int>L;
	for (int i = 0; i < 3; i++)
		L.push_back(K.GetGlobalIndex(k, i));
	for (int i = 0; i < K.GetNumOfGlobal(k); i++) {
		for (int j = 0; j < i; j++) {
			int ii = ig[L[i]];
			int iend = ig[L[i] + 1] - 1;
			for (int c = ii; c <= iend; c++)
				if (jg[c] == L[j]) { ii = c; break; }
			ggl[ii] += local_matrix[i][j];
			ggu[ii] += local_matrix[j][i];
			++ii;
		}
		di[L[i]] += local_matrix[i][i];
		global_right[L[i]] += local_right[i];
	}
	L.clear();
}
void Matrix::GetGlobal() {
	GeneratePortret();
	local_matrix.resize(3, std::vector<double>(3));
	local_right.resize(3);
	for (int k = 0; k < K.fe_num; k++) {
		GetLocal(k);
		LocalToGlobal(k);
	}
	allow_boundary_2();
	allow_boundary_3();
	allow_boundary_1();
}