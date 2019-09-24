////
//// Created by luoyy on 19-7-5.
//#pragma once
//#include "common.h"
//#include "xpilot/core/geometry/Pose2D.h"
//#include "xpilot/core/geometry/Point2D.h"
//#include "xpilot/core/shapes/line_string/LineString.h"
//#include "algorithm/astar/AstarNode.h"
//
//
//using namespace xpilot::core;
///**************数据类型定义***************/
//struct polar{
//    float distance;
//    float angel;
//};
//
//struct TO{
//    float tau;
//    float omega;
//};
///**************其他函数声明***************/
//polar rect_to_polar(float x, float y); //笛卡尔坐标系转极坐标系
//float mod2pi(float x); //令处于-pi,pi 之间
//TO tauOmega(float u,float v,float xi,float eta,float theta);
//Pose2Dd global2local(const Pose2Dd& startpose,const Pose2Dd& endpose); //将endpose在全局中的坐标转换到以startpose为原点/方向的局部坐标
//class PathReedsShepp;
//void plot_path(const Pose2Dd& startpose,const PathReedsShepp& path_find); //画出最后的reedsshepp曲线
//void plot_start_traj(int start_num, std::vector<float> &start_path_x,std::vector<float> &start_path_y);//画起点轨迹线
//
///**************路径类声明***************/
//class PathReedsShepp {
//public:
//    PathReedsShepp(){ cost=100000;if_all_front= false;};
//    std::string PathType[5]={"RS_NOP","RS_NOP" ,"RS_NOP" ,"RS_NOP","RS_NOP" };
//    float pathlength[5]={0};
//    float totallength=pathlength[0]+pathlength[1]+pathlength[2]+pathlength[3]+pathlength[4];
//    float cost;
//    //将reedsshepp曲线离散化至lineString类型
//    void reedsshepp_path_to_linestring(AstarNode curnode,LineString& reedsshepp_point);
//
//    bool if_all_front;
//private:
//
//};
//
//
///**************ReedsShepp类声明***************/
//class ReedsShepp {
//
//public:
//    ReedsShepp(Pose2Dd Startpose, Pose2Dd Endpose):startpose(Startpose),endpose(Endpose){};
//    std::string PathTypeList[18][5] = {
//            {"RS_LEFT",  "RS_RIGHT", "RS_LEFT",  "RS_NOP",  "RS_NOP" },     //0 LRLOO
//            {"RS_RIGHT", "RS_LEFT",  "RS_RIGHT", "RS_NOP",  "RS_NOP"} ,     //1 RLROO
//            {"RS_LEFT",  "RS_RIGHT", "RS_LEFT",  "RS_RIGHT","RS_NOP"} ,     //2 LRLRO
//            {"RS_RIGHT", "RS_LEFT",  "RS_RIGHT", "RS_LEFT", "RS_NOP"},      //3 RLRLO
//            {"RS_LEFT",  "RS_RIGHT", "RS_STRAIGHT", "RS_LEFT",  "RS_NOP"} , //4 LRSLO
//            {"RS_RIGHT", "RS_LEFT",  "RS_STRAIGHT", "RS_RIGHT", "RS_NOP"} , //5 RLSRO
//            {"RS_LEFT",  "RS_STRAIGHT", "RS_RIGHT", "RS_LEFT",  "RS_NOP"},  //6 LSRLO
//            {"RS_RIGHT", "RS_STRAIGHT", "RS_LEFT",  "RS_RIGHT", "RS_NOP"} , //7 RSLRO
//            {"RS_LEFT",  "RS_RIGHT", "RS_STRAIGHT", "RS_RIGHT", "RS_NOP"} , //8 LRSRO
//            {"RS_RIGHT", "RS_LEFT",  "RS_STRAIGHT", "RS_LEFT",  "RS_NOP" }, //9 RLSLO
//            {"RS_RIGHT", "RS_STRAIGHT", "RS_RIGHT", "RS_LEFT",  "RS_NOP"} , //10 RSRLO
//            {"RS_LEFT",  "RS_STRAIGHT", "RS_LEFT",  "RS_RIGHT", "RS_NOP"} , //11 LSLRO
//            {"RS_LEFT",  "RS_STRAIGHT", "RS_RIGHT", "RS_NOP",   "RS_NOP" }, //12 LSROO
//            {"RS_RIGHT", "RS_STRAIGHT", "RS_LEFT",  "RS_NOP",   "RS_NOP" }, //13 RSLOO
//            {"RS_LEFT",  "RS_STRAIGHT", "RS_LEFT",  "RS_NOP",   "RS_NOP"} , //14 LSLOO
//            {"RS_RIGHT", "RS_STRAIGHT", "RS_RIGHT", "RS_NOP",   "RS_NOP"} , //15 RSROO
//            {"RS_LEFT",  "RS_RIGHT", "RS_STRAIGHT", "RS_LEFT",  "RS_RIGHT"},//16 LRSLR
//            {"RS_RIGHT", "RS_LEFT",  "RS_STRAIGHT", "RS_RIGHT", "RS_LEFT"}};//17 RLSRL
//    float calcu_cost_CCCC(float temp_segm[5],int type_index);
//    float calcu_cost_CSC(float temp_segm[5], int type_index,int symbol);
//    PathReedsShepp FindRPath();
//    PathReedsShepp path; //最后返回的最短路径
//    void CSC(const Pose2Dd& new_endpose, bool &isok, PathReedsShepp &paths);
//    void CCC(const Pose2Dd& new_endpose, bool &isok, PathReedsShepp &paths);
//    void CCCC(const Pose2Dd& new_endpose, bool &isok, PathReedsShepp &paths);
//    void CCSC(const Pose2Dd& new_endpose, bool &isok, PathReedsShepp &paths);
//    void CCSCC(const Pose2Dd& new_endpose, bool &isok, PathReedsShepp &paths);
//
//    static bool LpRmL(float x,float y,float theta, float temp_segm[3]);
//    static bool LpRmSLmRp(float x,float y,float theta, float temp_segm[3]);
//    static bool LpRmSmLm(float x,float y,float theta, float temp_segm[3]);
//    static bool LpRmSmRm(float x,float y,float theta, float temp_segm[3]);
//    static bool LpRumLumRp(float x,float y,float theta, float temp_segm[3]);
//    static bool LpRupLumRm(float x,float y,float theta, float temp_segm[3]);
//    static bool LpSpLp(float x,float y,float theta, float temp_segm[3]);
//    static bool LpSpRp(float x,float y,float theta, float temp_segm[3]);
//
//private:
//    Pose2Dd startpose ,endpose;
//
//
//};
//
//const std::string Type[4]={"RS_LEFT","RS_RIGHT","RS_STRAIGHT","RS_NOP"};
//
//
//
