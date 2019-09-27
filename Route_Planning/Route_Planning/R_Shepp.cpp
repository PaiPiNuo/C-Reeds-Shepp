#include <cmath>
#include <vector>
#include "R_Shepp.h"
#include "main.h"
using namespace std;

PATH  FindRSPath(float x, float y, float phi) {
	float Lmin = INF;
	x = x / rmin;
	y = y / rmin;
	vector<PATH> path;
	PATH path_CSC= CSC(x, y, phi);
	PATH path_CCC = CCC(x, y, phi);
	PATH path_CCCC = CCCC(x, y, phi);
	//PATH path_CCSC = CCSC(x, y, phi);
	//PATH path_CCSCC = CCSCC(x, y, phi);
	path.push_back(path_CSC);
	path.push_back(path_CCC);
	path.push_back(path_CCCC);
	//path.push_back(path_CCSC);
	//path.push_back(path_CCSCC);
	int ind_Record = 0;
	for (int i = 0; i < path.size();i++) {
		if (true==path[i].isok) {
			if (Lmin > path[i].total_Len) {
				ind_Record = i;
				Lmin = path[i].total_Len;
			}
		}
	}
	return path[ind_Record];
}

PATH CSC(float x, float y, float phi) {
	float Lmin = INF,L;
	PATH tmpPath;
	tmpPath.pathType[0] = RS_NOP;
	tmpPath.pathType[1] = RS_NOP;
	tmpPath.pathType[2] = RS_NOP;
	tmpPath.pathType[3] = RS_NOP;
	tmpPath.pathType[4] = RS_NOP;
	tmpPath.t = 0;
	tmpPath.u = 0;
	tmpPath.v = 0;
	tmpPath.w = 0;
	tmpPath.x = 0;

	pathLength pLength = LpSpLp(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[14][0];
			tmpPath.pathType[1] = Path_Type[14][1];
			tmpPath.pathType[2] = Path_Type[14][2];
			tmpPath.pathType[3] = Path_Type[14][3];
			tmpPath.pathType[4] = Path_Type[14][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}

	pLength = LpSpLp(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[14][0];
			tmpPath.pathType[1] = Path_Type[14][1];
			tmpPath.pathType[2] = Path_Type[14][2];
			tmpPath.pathType[3] = Path_Type[14][3];
			tmpPath.pathType[4] = Path_Type[14][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpSpLp(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[15][0];
			tmpPath.pathType[1] = Path_Type[15][1];
			tmpPath.pathType[2] = Path_Type[15][2];
			tmpPath.pathType[3] = Path_Type[15][3];
			tmpPath.pathType[4] = Path_Type[15][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpSpLp(-x, -y, phi);//timeflp + reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[15][0];
			tmpPath.pathType[1] = Path_Type[15][1];
			tmpPath.pathType[2] = Path_Type[15][2];
			tmpPath.pathType[3] = Path_Type[15][3];
			tmpPath.pathType[4] = Path_Type[15][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	//------------------------
	pLength = LpSpRp(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[12][0];
			tmpPath.pathType[1] = Path_Type[12][1];
			tmpPath.pathType[2] = Path_Type[12][2];
			tmpPath.pathType[3] = Path_Type[12][3];
			tmpPath.pathType[4] = Path_Type[12][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpSpRp(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[12][0];
			tmpPath.pathType[1] = Path_Type[12][1];
			tmpPath.pathType[2] = Path_Type[12][2];
			tmpPath.pathType[3] = Path_Type[12][3];
			tmpPath.pathType[4] = Path_Type[12][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpSpRp(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[13][0];
			tmpPath.pathType[1] = Path_Type[13][1];
			tmpPath.pathType[2] = Path_Type[13][2];
			tmpPath.pathType[3] = Path_Type[13][3];
			tmpPath.pathType[4] = Path_Type[13][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpSpRp(-x, -y, phi);//timeflip + reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[13][0];
			tmpPath.pathType[1] = Path_Type[13][1];
			tmpPath.pathType[2] = Path_Type[13][2];
			tmpPath.pathType[3] = Path_Type[13][3];
			tmpPath.pathType[4] = Path_Type[13][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	if (INF == Lmin)
		tmpPath.isok = false;
	else
		tmpPath.isok = true;
	tmpPath.total_Len = abs(tmpPath.t)+ abs(tmpPath.u)+ abs(tmpPath.v)+ abs(tmpPath.w)+ abs(tmpPath.x);
	return tmpPath;
}

PATH CCC(float x, float y, float phi) {
	float Lmin = INF, L;
	PATH tmpPath;
	tmpPath.pathType[0] = RS_NOP;
	tmpPath.pathType[1] = RS_NOP;
	tmpPath.pathType[2] = RS_NOP;
	tmpPath.pathType[3] = RS_NOP;
	tmpPath.pathType[4] = RS_NOP;
	tmpPath.t = 0;
	tmpPath.u = 0;
	tmpPath.v = 0;
	tmpPath.w = 0;
	tmpPath.x = 0;

	pathLength pLength = LpRmL(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[0][0];
			tmpPath.pathType[1] = Path_Type[0][1];
			tmpPath.pathType[2] = Path_Type[0][2];
			tmpPath.pathType[3] = Path_Type[0][3];
			tmpPath.pathType[4] = Path_Type[0][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[0][0];
			tmpPath.pathType[1] = Path_Type[0][1];
			tmpPath.pathType[2] = Path_Type[0][2];
			tmpPath.pathType[3] = Path_Type[0][3];
			tmpPath.pathType[4] = Path_Type[0][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[1][0];
			tmpPath.pathType[1] = Path_Type[1][1];
			tmpPath.pathType[2] = Path_Type[1][2];
			tmpPath.pathType[3] = Path_Type[1][3];
			tmpPath.pathType[4] = Path_Type[1][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(-x, -y, phi);//timeflip+reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[1][0];
			tmpPath.pathType[1] = Path_Type[1][1];
			tmpPath.pathType[2] = Path_Type[1][2];
			tmpPath.pathType[3] = Path_Type[1][3];
			tmpPath.pathType[4] = Path_Type[1][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.v;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	// backwards
	float xb, yb;
	xb= x*cos(phi) + y*sin(phi);
	yb = x*sin(phi) - y*cos(phi);
	pLength = LpRmL(xb, yb, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[0][0];
			tmpPath.pathType[1] = Path_Type[0][1];
			tmpPath.pathType[2] = Path_Type[0][2];
			tmpPath.pathType[3] = Path_Type[0][3];
			tmpPath.pathType[4] = Path_Type[0][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.t;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(-xb, yb, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[0][0];
			tmpPath.pathType[1] = Path_Type[0][1];
			tmpPath.pathType[2] = Path_Type[0][2];
			tmpPath.pathType[3] = Path_Type[0][3];
			tmpPath.pathType[4] = Path_Type[0][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.t;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(xb, -yb, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[1][0];
			tmpPath.pathType[1] = Path_Type[1][1];
			tmpPath.pathType[2] = Path_Type[1][2];
			tmpPath.pathType[3] = Path_Type[1][3];
			tmpPath.pathType[4] = Path_Type[1][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.t;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmL(-xb, -yb, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[1][0];
			tmpPath.pathType[1] = Path_Type[1][1];
			tmpPath.pathType[2] = Path_Type[1][2];
			tmpPath.pathType[3] = Path_Type[1][3];
			tmpPath.pathType[4] = Path_Type[1][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.t;
			tmpPath.w = 0;
			tmpPath.x = 0;
		}
	}
	if (INF == Lmin)
		tmpPath.isok = false;
	else
		tmpPath.isok = true;
	tmpPath.total_Len = abs(tmpPath.t) + abs(tmpPath.u) + abs(tmpPath.v) + abs(tmpPath.w) + abs(tmpPath.x);
	return tmpPath;
}

PATH CCCC(float x, float y, float phi) {
	float Lmin = INF, L;
	PATH tmpPath;
	tmpPath.pathType[0] = RS_NOP;
	tmpPath.pathType[1] = RS_NOP;
	tmpPath.pathType[2] = RS_NOP;
	tmpPath.pathType[3] = RS_NOP;
	tmpPath.pathType[4] = RS_NOP;
	tmpPath.t = 0;
	tmpPath.u = 0;
	tmpPath.v = 0;
	tmpPath.w = 0;
	tmpPath.x = 0;
	pathLength pLength = LpRupLumRm(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + 2*abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[2][0];
			tmpPath.pathType[1] = Path_Type[2][1];
			tmpPath.pathType[2] = Path_Type[2][2];
			tmpPath.pathType[3] = Path_Type[2][3];
			tmpPath.pathType[4] = Path_Type[2][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = -pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRupLumRm(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[2][0];
			tmpPath.pathType[1] = Path_Type[2][1];
			tmpPath.pathType[2] = Path_Type[2][2];
			tmpPath.pathType[3] = Path_Type[2][3];
			tmpPath.pathType[4] = Path_Type[2][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRupLumRm(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[3][0];
			tmpPath.pathType[1] = Path_Type[3][1];
			tmpPath.pathType[2] = Path_Type[3][2];
			tmpPath.pathType[3] = Path_Type[3][3];
			tmpPath.pathType[4] = Path_Type[3][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = -pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRupLumRm(-x, -y, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[3][0];
			tmpPath.pathType[1] = Path_Type[3][1];
			tmpPath.pathType[2] = Path_Type[3][2];
			tmpPath.pathType[3] = Path_Type[3][3];
			tmpPath.pathType[4] = Path_Type[3][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	
	pLength = LpRumLumRp(-x, -y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[2][0];
			tmpPath.pathType[1] = Path_Type[2][1];
			tmpPath.pathType[2] = Path_Type[2][2];
			tmpPath.pathType[3] = Path_Type[2][3];
			tmpPath.pathType[4] = Path_Type[2][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRumLumRp(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[2][0];
			tmpPath.pathType[1] = Path_Type[2][1];
			tmpPath.pathType[2] = Path_Type[2][2];
			tmpPath.pathType[3] = Path_Type[2][3];
			tmpPath.pathType[4] = Path_Type[2][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRumLumRp(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[3][0];
			tmpPath.pathType[1] = Path_Type[3][1];
			tmpPath.pathType[2] = Path_Type[3][2];
			tmpPath.pathType[3] = Path_Type[3][3];
			tmpPath.pathType[4] = Path_Type[3][4];
			tmpPath.t = pLength.t;
			tmpPath.u = pLength.u;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRumLumRp(-x, -y, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + 2 * abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[3][0];
			tmpPath.pathType[1] = Path_Type[3][1];
			tmpPath.pathType[2] = Path_Type[3][2];
			tmpPath.pathType[3] = Path_Type[3][3];
			tmpPath.pathType[4] = Path_Type[3][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = -pLength.u;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	if (INF == Lmin)
		tmpPath.isok = false;
	else
		tmpPath.isok = true;
	tmpPath.total_Len = abs(tmpPath.t) + abs(tmpPath.u) + abs(tmpPath.v) + abs(tmpPath.w) + abs(tmpPath.x);
	return tmpPath;
}
PATH CCSC(float x, float y, float phi) {
	float Lmin = INF, L;
	PATH tmpPath;
	tmpPath.pathType[0] = RS_NOP;
	tmpPath.pathType[1] = RS_NOP;
	tmpPath.pathType[2] = RS_NOP;
	tmpPath.pathType[3] = RS_NOP;
	tmpPath.pathType[4] = RS_NOP;
	tmpPath.t = 0;
	tmpPath.u = 0;
	tmpPath.v = 0;
	tmpPath.w = 0;
	tmpPath.x = 0;
	pathLength pLength = LpRmSmLm(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[4][0];
			tmpPath.pathType[1] = Path_Type[4][1];
			tmpPath.pathType[2] = Path_Type[4][2];
			tmpPath.pathType[3] = Path_Type[4][3];
			tmpPath.pathType[4] = Path_Type[4][4];
			tmpPath.t = pLength.t;
			tmpPath.u = -PI/2;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[4][0];
			tmpPath.pathType[1] = Path_Type[4][1];
			tmpPath.pathType[2] = Path_Type[4][2];
			tmpPath.pathType[3] = Path_Type[4][3];
			tmpPath.pathType[4] = Path_Type[4][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[5][0];
			tmpPath.pathType[1] = Path_Type[5][1];
			tmpPath.pathType[2] = Path_Type[5][2];
			tmpPath.pathType[3] = Path_Type[5][3];
			tmpPath.pathType[4] = Path_Type[5][4];
			tmpPath.t = pLength.t;
			tmpPath.u = -PI / 2;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(-x, -y, phi);//reflect+timeflip 
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[5][0];
			tmpPath.pathType[1] = Path_Type[5][1];
			tmpPath.pathType[2] = Path_Type[5][2];
			tmpPath.pathType[3] = Path_Type[5][3];
			tmpPath.pathType[4] = Path_Type[5][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[8][0];
			tmpPath.pathType[1] = Path_Type[8][1];
			tmpPath.pathType[2] = Path_Type[8][2];
			tmpPath.pathType[3] = Path_Type[8][3];
			tmpPath.pathType[4] = Path_Type[8][4];
			tmpPath.t = pLength.t;
			tmpPath.u =- PI / 2;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(-x, y, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[8][0];
			tmpPath.pathType[1] = Path_Type[8][1];
			tmpPath.pathType[2] = Path_Type[8][2];
			tmpPath.pathType[3] = Path_Type[8][3];
			tmpPath.pathType[4] = Path_Type[8][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(x, -y, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[9][0];
			tmpPath.pathType[1] = Path_Type[9][1];
			tmpPath.pathType[2] = Path_Type[9][2];
			tmpPath.pathType[3] = Path_Type[9][3];
			tmpPath.pathType[4] = Path_Type[9][4];
			tmpPath.t = pLength.t;
			tmpPath.u = -PI / 2;
			tmpPath.v = pLength.u;
			tmpPath.w = pLength.v;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(-x, -y, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[9][0];
			tmpPath.pathType[1] = Path_Type[9][1];
			tmpPath.pathType[2] = Path_Type[9][2];
			tmpPath.pathType[3] = Path_Type[9][3];
			tmpPath.pathType[4] = Path_Type[9][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = -pLength.v;
			tmpPath.x = 0;
		}
	}
	//----------------backwards------------
	float xb = x*cos(phi) + y*sin(phi);
	float yb = x*sin(phi) - y*cos(phi);
	pLength = LpRmSmLm(xb, yb, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[6][0];
			tmpPath.pathType[1] = Path_Type[6][1];
			tmpPath.pathType[2] = Path_Type[6][2];
			tmpPath.pathType[3] = Path_Type[6][3];
			tmpPath.pathType[4] = Path_Type[6][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = -PI/2;
			tmpPath.w = pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(-xb, yb, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[6][0];
			tmpPath.pathType[1] = Path_Type[6][1];
			tmpPath.pathType[2] = Path_Type[6][2];
			tmpPath.pathType[3] = Path_Type[6][3];
			tmpPath.pathType[4] = Path_Type[6][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = PI / 2;
			tmpPath.w = -pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(xb, -yb, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[7][0];
			tmpPath.pathType[1] = Path_Type[7][1];
			tmpPath.pathType[2] = Path_Type[7][2];
			tmpPath.pathType[3] = Path_Type[7][3];
			tmpPath.pathType[4] = Path_Type[7][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = -PI / 2;
			tmpPath.w = pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmLm(-xb, -yb, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[7][0];
			tmpPath.pathType[1] = Path_Type[7][1];
			tmpPath.pathType[2] = Path_Type[7][2];
			tmpPath.pathType[3] = Path_Type[7][3];
			tmpPath.pathType[4] = Path_Type[7][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = PI / 2;
			tmpPath.w = -pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(xb, yb, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[10][0];
			tmpPath.pathType[1] = Path_Type[10][1];
			tmpPath.pathType[2] = Path_Type[10][2];
			tmpPath.pathType[3] = Path_Type[10][3];
			tmpPath.pathType[4] = Path_Type[10][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = -PI / 2;
			tmpPath.w = pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(-xb, yb, -phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[10][0];
			tmpPath.pathType[1] = Path_Type[10][1];
			tmpPath.pathType[2] = Path_Type[10][2];
			tmpPath.pathType[3] = Path_Type[10][3];
			tmpPath.pathType[4] = Path_Type[10][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = PI / 2;
			tmpPath.w = -pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(xb, -yb, -phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[11][0];
			tmpPath.pathType[1] = Path_Type[11][1];
			tmpPath.pathType[2] = Path_Type[11][2];
			tmpPath.pathType[3] = Path_Type[11][3];
			tmpPath.pathType[4] = Path_Type[11][4];
			tmpPath.t = pLength.v;
			tmpPath.u = pLength.u;
			tmpPath.v = -PI / 2;
			tmpPath.w = pLength.t;
			tmpPath.x = 0;
		}
	}
	pLength = LpRmSmRm(-xb, -yb, phi);//reflect+timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[11][0];
			tmpPath.pathType[1] = Path_Type[11][1];
			tmpPath.pathType[2] = Path_Type[11][2];
			tmpPath.pathType[3] = Path_Type[11][3];
			tmpPath.pathType[4] = Path_Type[11][4];
			tmpPath.t = -pLength.v;
			tmpPath.u = -pLength.u;
			tmpPath.v = PI / 2;
			tmpPath.w = -pLength.t;
			tmpPath.x = 0;
		}
	}
	if (INF == Lmin)
		tmpPath.isok = false;
	else
		tmpPath.isok = true;
	tmpPath.total_Len = abs(tmpPath.t) + abs(tmpPath.u) + abs(tmpPath.v) + abs(tmpPath.w) + abs(tmpPath.x);
	return tmpPath;
}
PATH CCSCC(float x, float y, float phi) {
	float Lmin = INF, L;
	PATH tmpPath;
	tmpPath.pathType[0] = RS_NOP;
	tmpPath.pathType[1] = RS_NOP;
	tmpPath.pathType[2] = RS_NOP;
	tmpPath.pathType[3] = RS_NOP;
	tmpPath.pathType[4] = RS_NOP;
	tmpPath.t = 0;
	tmpPath.u = 0;
	tmpPath.v = 0;
	tmpPath.w = 0;
	tmpPath.x = 0;
	pathLength pLength = LpRmSLmRp(x, y, phi);
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[16][0];
			tmpPath.pathType[1] = Path_Type[16][1];
			tmpPath.pathType[2] = Path_Type[16][2];
			tmpPath.pathType[3] = Path_Type[16][3];
			tmpPath.pathType[4] = Path_Type[16][4];
			tmpPath.t = pLength.t;
			tmpPath.u = -PI / 2;
			tmpPath.v = pLength.u;
			tmpPath.w = -PI / 2;
			tmpPath.x = pLength.v;
		}
	}
	//-------------------**Maybe there is some problems,Be Careful------
	pLength = LpRmSLmRp(x, y, phi);//timeflip
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[16][0];
			tmpPath.pathType[1] = Path_Type[16][1];
			tmpPath.pathType[2] = Path_Type[16][2];
			tmpPath.pathType[3] = Path_Type[16][3];
			tmpPath.pathType[4] = Path_Type[16][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = PI / 2;
			tmpPath.x = -pLength.v;
		}
	}
	pLength = LpRmSLmRp(x, y, phi);//reflect
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[17][0];
			tmpPath.pathType[1] = Path_Type[17][1];
			tmpPath.pathType[2] = Path_Type[17][2];
			tmpPath.pathType[3] = Path_Type[17][3];
			tmpPath.pathType[4] = Path_Type[17][4];
			tmpPath.t = pLength.t;
			tmpPath.u = -PI / 2;
			tmpPath.v = pLength.u;
			tmpPath.w = -PI / 2;
			tmpPath.x = pLength.v;
		}
	}
	pLength = LpRmSLmRp(x, y, phi);//reflect+timeflip 
	if (pLength.isok) {
		L = abs(pLength.t) + abs(pLength.u) + abs(pLength.v);
		if (Lmin > L) {
			Lmin = L;
			tmpPath.pathType[0] = Path_Type[17][0];
			tmpPath.pathType[1] = Path_Type[17][1];
			tmpPath.pathType[2] = Path_Type[17][2];
			tmpPath.pathType[3] = Path_Type[17][3];
			tmpPath.pathType[4] = Path_Type[17][4];
			tmpPath.t = -pLength.t;
			tmpPath.u = PI / 2;
			tmpPath.v = -pLength.u;
			tmpPath.w = PI / 2;
			tmpPath.x = -pLength.v;
		}
	}
	if (INF == Lmin)
		tmpPath.isok = false;
	else
		tmpPath.isok = true;
	tmpPath.total_Len = abs(tmpPath.t) + abs(tmpPath.u) + abs(tmpPath.v) + abs(tmpPath.w) + abs(tmpPath.x);
	return tmpPath;
}
//L+S+L+ formula 8.1
pathLength LpSpLp(float x, float y, float phi) {
	pathLength pathLen;
	float theta, uLen;
	cart2pol(x - sinf(phi), y - 1 + cosf(phi), uLen, theta);
	pathLen.t = theta;
	pathLen.u = uLen;
	if (pathLen.t >= 0) {
		pathLen.v = mod2pi(phi - pathLen.t);
		if (pathLen.v >= 0) {
			pathLen.isok = true;
			return pathLen;
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
//formula 8.2
pathLength LpSpRp(float x, float y, float phi) {
	pathLength pathLen;
	float theta1, u1,tmpTheta;
	cart2pol(x + sin(phi), y - 1 - cos(phi), u1, theta1);
	if (u1*u1>=4) {
		pathLen.u = sqrt(u1*u1 - 4);
		tmpTheta= atan2(2, pathLen.u);	
		pathLen.t = mod2pi(theta1 + tmpTheta);
		pathLen.v = mod2pi(pathLen.t- phi);
		if (pathLen.t>=0&& pathLen.v>=0) {
			pathLen.isok = true;
			return pathLen;
		}	
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
//L+R-L+  and  L+R-L-     formula 8.3/8.4
pathLength LpRmL(float x, float y, float phi) {
	pathLength pathLen;
	float u,theta;
	cart2pol(x-sin(phi),y-1+cos(phi),u, theta);
	if (u<=4) {
		pathLen.u = -2 * asin(u/4);
		pathLen.t = mod2pi(theta+ pathLen.u/2+PI);
		pathLen.v = mod2pi(phi- pathLen.t+ pathLen.u);
		if (pathLen.t>=0 && pathLen.u<=0) {
			pathLen.isok = true;
			return pathLen;
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
pathLength LpRupLumRm(float x, float y, float phi) {
	pathLength pathLen;
	float xi = x + sin(phi);
	float eta= y - 1 - cos(phi);
	float rho = (2 + sqrt(xi * xi  + eta *eta)) / 4;
	if (rho <= 1) {
		pathLen.u = acos(rho);
		tauOmega(pathLen.u,-pathLen.u,xi,eta,phi, pathLen.t, pathLen.v);
		if (pathLen.t>=0&& pathLen.v<=0) {
			pathLen.isok = true;
			return pathLen;
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
pathLength LpRumLumRp(float x, float y, float phi) {
	pathLength pathLen;
	float xi= x + sin(phi);
	float eta= y - 1 - cos(phi);
	float rho= (20 - xi *xi - eta *eta) / 16;
	if (rho>=0&&rho<=1) {
		pathLen.u = -acos(rho);
		if (pathLen.u>=(PI/2)) {
			tauOmega(pathLen.u, pathLen.u,xi,eta,phi, pathLen.t, pathLen.v);
			if (pathLen.t>=0&& pathLen.v>=0) {
				pathLen.isok = true;
				return pathLen;
			}
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
pathLength LpRmSmLm(float x, float y, float phi) {
	pathLength pathLen;
	float xi = x - sin(phi);
	float eta= y - 1 + cos(phi);
	float theta, rho;
	cart2pol(xi,eta, rho, theta);
	if (rho >= 2) {
		float r = sqrt(rho*rho-4);
		pathLen.u = 2 - r;
		pathLen.t = mod2pi(theta+atan2(r,-2));
		pathLen.v = mod2pi(phi-PI/2- pathLen.t);
		if (pathLen.t>=0&& pathLen.u<=0&& pathLen.v<=0) {
			pathLen.isok = true;
			return pathLen;
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;
	return pathLen;
}
pathLength LpRmSmRm(float x, float y, float phi) {
	pathLength pathLen;
	float xi= x + sin(phi);
	float eta= y - 1 - cos(phi);
	float theta, rho;
	cart2pol(-eta,xi, rho, theta);
	if (rho>=2) {
		pathLen.t = theta;
		pathLen.u = 2 - rho;
		pathLen.v = mod2pi(pathLen.t+PI/2-phi);
		if (pathLen.t>=0&& pathLen.u<=0&& pathLen.v<=0) {
			pathLen.isok = true;
			return pathLen;
		}
	}
	pathLen.isok = true;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;

	return pathLen;
}
pathLength LpRmSLmRp(float x, float y, float phi) {
	pathLength pathLen;
	float xi= x + sin(phi);
	float eta= y - 1 - cos(phi);
	float rho, theta;
	cart2pol(xi, eta,rho,theta);
	if (rho>=2) {
		pathLen.u = 4 - sqrt(rho*rho-4);
		if (pathLen.u<0) {
			pathLen.t= mod2pi(atan2((4-pathLen.u)*xi-2*eta,-2*xi+(pathLen.u-4)*eta));
			pathLen.v= mod2pi(pathLen.t - phi);
			if (pathLen.t>=0&& pathLen.v>=0) {
				pathLen.isok = true;
				return pathLen;
			}
		}
	}
	pathLen.isok = false;
	pathLen.t = 0;
	pathLen.u = 0;
	pathLen.v = 0;

	return pathLen;
}
void cart2pol(float x, float y,float &length,float &angle) {
	length = sqrt(x*x+y*y);
	angle = atan2(y, x);
}
float mod2pi(float x) {  // x去周期到 -pi,pi之间
	float v;
	v = fmod(x, 2 * PI);  //取余 结果符号跟x
	if (v < -PI)
		v = v + 2 * PI;
	else if (v > PI)
		v = v - 2 * PI;
	return v;
}
void tauOmega(float u,float v,float xi,float eta,float phi,float &tau,float &omega) {
	float delta = mod2pi(u - v);
	float A= sin(u) - sin(delta);
	float B= cos(u) - cos(delta) - 1;
	float t1 = atan2(eta*A - xi*B, xi*A + eta*B);
	float t2= 2 * (cos(delta) - 2 * cos(v) - 2 * cos(u)) + 3;
	if (t2 < 0)
		tau = mod2pi(t1 + PI);
	else
		tau = mod2pi(t1);
	omega = mod2pi(tau-u+v-phi);
}