#pragma once

enum PathDir { RS_NOP, RS_LEFT, RS_STRAIGHT, RS_RIGHT };

typedef struct PATH {
	int pathType[5];
	float t = 0;
	float u = 0;
	float v = 0;
	float w = 0;
	float x = 0;
	float total_Len=0;
	bool isok = false;
};

typedef struct pathLength {
	float t;
	float u;
	float v;
	float w;
	bool isok;
};


static int Path_Type[18][5]= {
	{ RS_LEFT,  RS_RIGHT, RS_LEFT,  RS_NOP,  RS_NOP },       //00-LRL
	{ RS_RIGHT, RS_LEFT,  RS_RIGHT, RS_NOP,  RS_NOP } ,      //01 RLR
	{ RS_LEFT,  RS_RIGHT, RS_LEFT,  RS_RIGHT,RS_NOP } ,     //02 LRLR
	{ RS_RIGHT, RS_LEFT,  RS_RIGHT, RS_LEFT, RS_NOP },      //03 RLRL
	{ RS_LEFT,  RS_RIGHT, RS_STRAIGHT, RS_LEFT,  RS_NOP } , //04 LRSL
	{ RS_RIGHT, RS_LEFT,  RS_STRAIGHT, RS_RIGHT, RS_NOP } , //05 RLSR
	{ RS_LEFT,  RS_STRAIGHT, RS_RIGHT, RS_LEFT,  RS_NOP },  //06 LSRL
	{ RS_RIGHT, RS_STRAIGHT, RS_LEFT,  RS_RIGHT, RS_NOP } , //07 RSLR
	{ RS_LEFT,  RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP } , //08 LRSR
	{ RS_RIGHT, RS_LEFT,  RS_STRAIGHT, RS_LEFT,  RS_NOP },  //09 RLSL
	{ RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT,  RS_NOP } , //10 RSRL
	{ RS_LEFT,  RS_STRAIGHT, RS_LEFT,  RS_RIGHT, RS_NOP } , //11 LSLR
	{ RS_LEFT,  RS_STRAIGHT, RS_RIGHT, RS_NOP,   RS_NOP },  //12 LSR
	{ RS_RIGHT, RS_STRAIGHT, RS_LEFT,  RS_NOP,   RS_NOP },  //13 RSL
	{ RS_LEFT,  RS_STRAIGHT, RS_LEFT,  RS_NOP,   RS_NOP } , //14 LSL
	{ RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP,   RS_NOP } , //15 RSR
	{ RS_LEFT,  RS_RIGHT, RS_STRAIGHT, RS_LEFT,  RS_RIGHT },//16 LRSLR
	{ RS_RIGHT, RS_LEFT,  RS_STRAIGHT, RS_RIGHT, RS_LEFT }//17 RLSRL
};

PATH  FindRSPath(float x, float y, float phi);
PATH CSC(float x, float y, float phi);
PATH CCC(float x, float y, float phi);
PATH CCCC(float x, float y, float phi);
PATH CCSC(float x, float y, float phi);
PATH CCSCC(float x, float y, float phi);
pathLength LpSpLp(float x, float y, float phi);
pathLength LpSpRp(float x, float y, float phi);
pathLength LpRmL(float x, float y, float phi);
pathLength LpRupLumRm(float x, float y, float phi);
pathLength LpRumLumRp(float x, float y, float phi);
pathLength LpRmSmLm(float x, float y, float phi);
pathLength LpRmSmRm(float x, float y, float phi);
pathLength LpRmSLmRp(float x, float y, float phi);
void cart2pol(float x, float y, float &length, float &angle);
float mod2pi(float x);
void tauOmega(float u, float v, float xi, float eta, float phi, float &tau, float &omega);

