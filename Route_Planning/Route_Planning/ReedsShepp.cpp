//
//// Created by luoyy on 19-7-5.
//// Ctrl + / 行注释
//// Ctrl + Shift + / 块注释
////文档快速对齐 Ctrl+Shift+Alt+L
//#include "ReedsShepp.h"
//#include "main.h"
//#include <cmath>
//#include <vector>
////#include "matplotlibcpp.h"
//
//// namespace plt = matplotlibcpp;
//using namespace xpilot::core;
//
///**************PathReedsShepp函数实现***************/
//void PathReedsShepp::reedsshepp_path_to_linestring(AstarNode curnode, LineString &reedsshepp_point) {
//
//    Pose2Dd curpose;
//    curpose.x() = curnode.getPose().x();
//    curpose.y() = curnode.getPose().y();
//    curpose.theta() = curnode.getPose().theta();
//    // std::cout<<curpose<<std::endl;
//
//
//    double dl = 0.0;                         //第i段路径长度，有正负
//    double theta;                            //当前yaw
//    double dvec[3];                          //第i段路径末点的全局位姿变化量 delta_x,delta_y,delta_theta
//    double dtheta;                           //第i段弧度（计算中应该考虑上符号）
//    double cenx, ceny;                       //第i段圆心坐标
//
//    for (int i = 0; i < 5; ++i) {
//        if (PathType[i] == Type[2]) {  //若是 RS_STRAIGHT
//            theta = curpose.theta();
//            dl = pathlength[i];
//            dvec[0] = dl * cos(theta);
//            dvec[1] = dl * sin(theta);
//            dvec[2] = 0;
//            int n = std::fabs(this->pathlength[i])/0.5;  // plot间隔
//            if(!n)continue;
//            std::vector<double> dx(n), dy(n), t(n);  //第i段图像坐标,细分弧度
//            for (int j = 0; j < n; ++j) {  // 0,99
//                dx[j] = curpose.x() + dvec[0] * j / (n - 1);
//                dy[j] = curpose.y() + dvec[1] * j / (n - 1);
//                reedsshepp_point.push_back(Vertex2D(dx[j], dy[j]));
//            }
//            curpose.x() = curpose.x() + dvec[0];
//            curpose.y() = curpose.y() + dvec[1];
//            curpose.theta() = curpose.theta() + dvec[2];
//        } else if (PathType[i] == Type[0]) {  //若是 RS_LEFT
//            theta = curpose.theta();
//            dtheta = pathlength[i] / rmin;
//            cenx = curpose.x() - rmin * sin(theta);
//            ceny = curpose.y() + rmin * cos(theta);
//            int n = std::fabs(this->pathlength[i])/0.5;  // plot间隔
//            if(!n)continue;
//            std::vector<double> dx(n), dy(n), t(n);  //第i段图像坐标,细分弧度
//            for (int j = 0; j < n; ++j) {
//                t[j] = theta - M_PI_2 + dtheta * j / (n - 1);
//                dx[j] = cenx + rmin * cosf(t[j]);
//                dy[j] = ceny + rmin * sinf(t[j]);
//                reedsshepp_point.push_back(Vertex2D(dx[j], dy[j]));
//            }
//            theta = theta + dtheta;
//            curpose.x() = dx[n - 1];
//            curpose.y() = dy[n - 1];
//            curpose.theta() = theta;
//            dl = dtheta;
//        } else if (PathType[i] == Type[1]) {  //若是 RS_RIGHT
//            theta = curpose.theta();
//            dtheta = -pathlength[i] / rmin;  //右转 yaw是减小的，带负号
//            cenx = curpose.x() + rmin * sin(theta);
//            ceny = curpose.y() - rmin * cos(theta);
//            int n = std::fabs(this->pathlength[i])/0.5;  // plot间隔
//            if(!n)continue;
//            std::vector<double> dx(n), dy(n), t(n);  //第i段图像坐标,细分弧度
//            for (int j = 0; j < n; ++j) {
//                t[j] = theta + M_PI_2 + dtheta * j / (n - 1);
//                dx[j] = cenx + rmin * cos(t[j]);
//                dy[j] = ceny + rmin * sin(t[j]);
//                reedsshepp_point.push_back(Vertex2D(dx[j], dy[j]));
//            }
//            theta = theta + dtheta;
//            curpose.x() = dx[n - 1];
//            curpose.y() = dy[n - 1];
//            curpose.theta() = theta;
//            dl = -dtheta;
//        } else {
//        }  // do none
//    }
//}
//
///**************ReedsShepp类函数实现***************/
//
//PathReedsShepp ReedsShepp::FindRPath() {
//    // std::cout<<startpose.y()<<std::endl;
//
//    PathReedsShepp paths[5], path_temp;
//    bool isok[5] = {false};
//    float Lmin = 100000;
//    //缩放最小转弯半径到1
//    Pose2Dd new_endpose(endpose.x() / rmin, endpose.y() / rmin, endpose.theta());
//
//    CSC(new_endpose, isok[0], paths[0]);
//    CCC(new_endpose, isok[1], paths[1]);
//    CCCC(new_endpose, isok[2], paths[2]);
//    CCSC(new_endpose, isok[3], paths[3]);
//    CCSCC(new_endpose, isok[4], paths[4]);
//
//    for (int(i) = 0; (i) < 5; ++(i)) {
//        if (!isok[i]) continue;  //若找到路径,把该路径存入path.temp
//        for (int j = 0; j < 5; ++j) {
//            path_temp.PathType[j] = paths[i].PathType[j];
//            path_temp.pathlength[j] = paths[i].pathlength[j];
//        }
//        path_temp.totallength = paths[i].totallength;
//        path_temp.cost = paths[i].cost;
//
//        if (Lmin - paths[i].cost > precision_float) {  //若该路径是当前最短路径，把该路径赋值给path
//            Lmin = paths[i].cost;
//            for (int j = 0; j < 5; ++j) {
//                path.PathType[j] = path_temp.PathType[j];
//                path.pathlength[j] = path_temp.pathlength[j] * rmin;
//            }
//            path.totallength = path_temp.totallength * rmin;
//            path.cost = path_temp.cost * rmin;
//            if (path.pathlength[0] > 0 && path.pathlength[1] > 0 && path.pathlength[2] > 0 && path.pathlength[3] > 0 &&
//                path.pathlength[4] > 0)
//                path.if_all_front = true;
//        }
//    }
//
//    return path;
//}
//float ReedsShepp::calcu_cost_CSC(float *temp_segm, int type_index, int symbol) {
//    if (symbol > 0) {
//        return (std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][0])) +
//                std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][1])) +
//                std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][2])) +
//                std::fabs(temp_segm[3] * Parameter::judge_penalty_back(temp_segm[3]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][3])) +
//                std::fabs(temp_segm[4] * Parameter::judge_penalty_back(temp_segm[4]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][4])));
//    } else {
//        return (std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][0])) +
//                std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][1])) +
//                std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][2])) +
//                std::fabs(temp_segm[3] * Parameter::judge_penalty_back(-temp_segm[3]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][3])) +
//                std::fabs(temp_segm[4] * Parameter::judge_penalty_back(-temp_segm[4]) *
//                          Parameter::judge_penalty_turn(PathTypeList[type_index][4])));
//    }
//}
//
//void ReedsShepp::CSC(const Pose2Dd &new_endpose, bool &isoki, PathReedsShepp &pathi) {
//    float Lmin = 100000;
//    float temp_segm[5] = {0};  // t,u,v,w,x
//    float x = new_endpose.x(), y = new_endpose.y(), theta = new_endpose.theta();
//    float temp_cost = 100000;
//    //该类型路径第一种情况
//    isoki = LpSpLp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 14, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[14][i];  // LSLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpSpLp(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 14, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[14][i];  // LSLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpSpLp(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 15, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[15][i];  // LSLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpSpLp(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 15, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[15][i];  // LSLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // new condition
//    isoki = LpSpRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 12, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[12][i];  // LSLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timrflip
//    isoki = LpSpRp(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 12, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[12][i];  // LSLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpSpRp(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 13, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[13][i];  // LSLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timflip+reflect
//    isoki = LpSpRp(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 13, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[13][i];  // LSLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // std::cout<<"temp_cost: "<<temp_cost<<std::endl;
//    isoki = Lmin != 100000;
//}
//
//void ReedsShepp::CCC(const Pose2Dd &new_endpose, bool &isoki, PathReedsShepp &pathi) {
//    float Lmin = 100000;
//    float temp_segm[5] = {0};  // t,u,v,w,x
//    float x = new_endpose.x(), y = new_endpose.y(), theta = new_endpose.theta();
//    float temp_cost = 100000;
//    //该类型路径第一种情况
//    isoki = LpRmL(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 0, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[0][i];  // LRLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmL(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 0, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[0][i];  // LRLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmL(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 1, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[1][i];  // LRLOO
//                pathi.pathlength[i] = temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmL(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 1, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[1][i];  // LRLOO
//                pathi.pathlength[i] = -temp_segm[i];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // backwards
//    float xb, yb;
//    xb = x * cosf(theta) + y * sinf(theta);
//    yb = x * sinf(theta) - y * cosf(theta);
//    isoki = LpRmL(xb, yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 0, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[0][i];  // LRLOO
//                pathi.pathlength[i] = temp_segm[i];
//                if (i == 0) pathi.pathlength[i] = temp_segm[2];
//                if (i == 1) pathi.pathlength[i] = temp_segm[1];
//                if (i == 2) pathi.pathlength[i] = temp_segm[0];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmL(-xb, yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 0, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[0][i];  // LRLOO
//                pathi.pathlength[i] = -temp_segm[i];
//                if (i == 0) pathi.pathlength[i] = -temp_segm[2];
//                if (i == 1) pathi.pathlength[i] = -temp_segm[1];
//                if (i == 2) pathi.pathlength[i] = -temp_segm[0];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmL(xb, -yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 1, 1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[1][i];  // LRLOO
//                pathi.pathlength[i] = temp_segm[i];
//                if (i == 0) pathi.pathlength[i] = temp_segm[2];
//                if (i == 1) pathi.pathlength[i] = temp_segm[1];
//                if (i == 2) pathi.pathlength[i] = temp_segm[0];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmL(-xb, -yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = calcu_cost_CSC(temp_segm, 1, -1);
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) +
//                                std::fabs(temp_segm[3]) + std::fabs(temp_segm[4]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[1][i];  // LRLOO
//                pathi.pathlength[i] = -temp_segm[i];
//                if (i == 0) pathi.pathlength[i] = -temp_segm[2];
//                if (i == 1) pathi.pathlength[i] = -temp_segm[1];
//                if (i == 2) pathi.pathlength[i] = -temp_segm[0];
//            }
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // std::cout<<"temp_cost: "<<temp_cost<<std::endl;
//    isoki = 100000 != Lmin;
//}
//
////*Parameter::judge_penalty_back(temp_segm[i]) * Parameter::judge_penalty_turn(PathTypeList[type_index][i])
//
//void ReedsShepp::CCCC(const Pose2Dd &new_endpose, bool &isoki, PathReedsShepp &pathi) {
//    float Lmin = 100000;
//    float temp_segm[3] = {0};  // t,u,v
//    float x = new_endpose.x(), y = new_endpose.y(), theta = new_endpose.theta();
//    float temp_cost = 100000;
//    //该类型路径第一种情况
//    isoki = LpRupLumRm(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[2][i];
//            }                                     // LRLRO
//            pathi.pathlength[0] = temp_segm[0];   // t
//            pathi.pathlength[1] = temp_segm[1];   // u
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = temp_segm[2];   // v
//
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRupLumRm(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[2][i];
//            }                                     // LRLROO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = temp_segm[1];   // u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRupLumRm(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[3][i];
//            }                                     // LRLROO
//            pathi.pathlength[0] = temp_segm[0];   //-t
//            pathi.pathlength[1] = temp_segm[1];   //-u
//            pathi.pathlength[2] = -temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];   //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRupLumRm(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[3][i];
//            }                                     // LRLROO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = temp_segm[1];   // u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // new condition
//    isoki = LpRumLumRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[2][i];
//            }                                    // LRLROO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRumLumRp(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[2][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[2][i];
//            }                                     // LRLROO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRumLumRp(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[3][i];
//            }                                    // LRLRO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRumLumRp(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[3][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + 2 * std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]);
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[3][i];
//            }                                     // LRLROO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // std::cout<<"temp_cost: "<<temp_cost<<std::endl;
//    isoki = Lmin != 100000;
//}
//
//void ReedsShepp::CCSC(const Pose2Dd &new_endpose, bool &isoki, PathReedsShepp &pathi) {
//    float Lmin = 100000;
//    float temp_segm[3] = {0};  // t,u,v
//    float x = new_endpose.x(), y = new_endpose.y(), theta = new_endpose.theta();
//    float temp_cost = 100000;
//    //该类型路径第一种情况
//    isoki = LpRmSmLm(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[4][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmSmLm(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[4][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[4][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmSmLm(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[5][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmSmLm(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[5][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[5][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // new
//    isoki = LpRmSmRm(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[8][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmSmRm(-x, y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[8][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[8][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmSmRm(x, -y, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[9][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmSmRm(-x, -y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][2])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[9][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[9][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // backwards
//    float xb, yb;
//    xb = x * cosf(theta) + y * sinf(theta);
//    yb = x * sinf(theta) - y * cosf(theta);
//    isoki = LpRmSmLm(xb, yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[6][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[2];  // v
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = -M_PI_2;       //-pi/2
//            pathi.pathlength[3] = temp_segm[0];  // t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmSmLm(-xb, yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[6][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[6][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[2];  //-v
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = M_PI_2;         // pi/2
//            pathi.pathlength[3] = -temp_segm[0];  //-t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmSmLm(xb, -yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[7][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[2];  // v
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = -M_PI_2;       //-pi/2
//            pathi.pathlength[3] = temp_segm[0];  // t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmSmLm(-xb, -yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[7][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[7][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[2];  //-v
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = M_PI_2;         // pi/2
//            pathi.pathlength[3] = -temp_segm[0];  //-t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // new
//    isoki = LpRmSmRm(xb, yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            Lmin = pathi.totallength;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[10][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[2];  // v
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = -M_PI_2;       //-pi/2
//            pathi.pathlength[3] = temp_segm[0];  // t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmSmRm(-xb, yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[10][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[10][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[2];  //-v
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = M_PI_2;         // pi/2
//            pathi.pathlength[3] = -temp_segm[0];  //-t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmSmRm(xb, -yb, -theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[11][i];
//            }                                    // LRSLOO
//            pathi.pathlength[0] = temp_segm[2];  // v
//            pathi.pathlength[1] = temp_segm[1];  // u
//            pathi.pathlength[2] = -M_PI_2;       //-pi/2
//            pathi.pathlength[3] = temp_segm[0];  // t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmSmRm(-xb, -yb, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][0])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][1])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][2])) +
//                    std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[11][3]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI_2;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[11][i];
//            }                                     // LRSLOO
//            pathi.pathlength[0] = -temp_segm[2];  //-v
//            pathi.pathlength[1] = -temp_segm[1];  //-u
//            pathi.pathlength[2] = M_PI_2;         // pi/2
//            pathi.pathlength[3] = -temp_segm[0];  //-t
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // std::cout<<"temp_cost: "<<temp_cost<<std::endl;
//    isoki = Lmin != 100000;
//}
//
//void ReedsShepp::CCSCC(const Pose2Dd &new_endpose, bool &isoki, PathReedsShepp &pathi) {
//    float Lmin = 100000;
//    float temp_segm[3] = {0};  // t,u,v
//    float x = new_endpose.x(), y = new_endpose.y(), theta = new_endpose.theta();
//    float temp_cost = 100000;
//    //该类型路径第一种情况
//    isoki = LpRmSLmRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][2])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][3])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][4]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[16][i];
//            }                                    // LRSLR
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = -M_PI_2;       //-pi/2
//            pathi.pathlength[4] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip
//    isoki = LpRmSLmRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][2])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][3])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[16][4]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[16][i];
//            }                                     // LRSLR
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = M_PI_2;         // pi/2
//            pathi.pathlength[4] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // reflect
//    isoki = LpRmSLmRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][2])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(-M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][3])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][4]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[17][i];
//            }                                    // LRSLR
//            pathi.pathlength[0] = temp_segm[0];  // t
//            pathi.pathlength[1] = -M_PI_2;       //-pi/2
//            pathi.pathlength[2] = temp_segm[1];  // u
//            pathi.pathlength[3] = -M_PI_2;       //-pi/2
//            pathi.pathlength[4] = temp_segm[2];  // v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // timeflip+reflect
//    isoki = LpRmSLmRp(x, y, theta, temp_segm);
//    if (isoki) {
//        temp_cost = std::fabs(temp_segm[0] * Parameter::judge_penalty_back(-temp_segm[0]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][0])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][1])) +
//                    std::fabs(temp_segm[1] * Parameter::judge_penalty_back(-temp_segm[1]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][2])) +
//                    std::fabs(M_PI_2 * Parameter::judge_penalty_back(M_PI_2) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][3])) +
//                    std::fabs(temp_segm[2] * Parameter::judge_penalty_back(-temp_segm[2]) *
//                              Parameter::judge_penalty_turn(PathTypeList[17][4]));
//        //若当前找到的路径长度比之前最短路径短，更新pathi
//        if ((Lmin - temp_cost) > precision_float) {
//            pathi.totallength = std::fabs(temp_segm[0]) + std::fabs(temp_segm[1]) + std::fabs(temp_segm[2]) + M_PI;
//            for (int i = 0; i < 5; ++i) {
//                pathi.PathType[i] = PathTypeList[17][i];
//            }                                     // LRSLR
//            pathi.pathlength[0] = -temp_segm[0];  //-t
//            pathi.pathlength[1] = M_PI_2;         // pi/2
//            pathi.pathlength[2] = -temp_segm[1];  //-u
//            pathi.pathlength[3] = M_PI_2;         // pi/2
//            pathi.pathlength[4] = -temp_segm[2];  //-v
//            pathi.cost = temp_cost;
//            Lmin = pathi.cost;
//        }
//    }
//    // std::cout<<"temp_cost: "<<temp_cost<<std::endl;
//    isoki = Lmin != 100000;
//}
//
//bool ReedsShepp::LpRmL(float x, float y, float theta, float temp_segm[3]) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float t = 0, u = 0, v = 0;
//    float xi, eta;  //坐标转换
//    xi = x - sinf(theta);
//    eta = y - 1 + cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(xi, eta);
//    if (pl.distance < 4) {
//        u = -2 * asinf(pl.distance / 4);
//        t = mod2pi(pl.angel + u / 2 + M_PI);
//        v = mod2pi(theta - t + u);
//        if (t >= 0 && u <= 0) {  //有解
//            isoki = true;
//            temp_segm[0] = t;
//            temp_segm[1] = u;
//            temp_segm[2] = v;
//        }
//    }
//    return isoki;
//}
//bool ReedsShepp::LpRmSLmRp(float x, float y, float theta, float temp_segm[3]) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float t = 0, u = 0, v = 0;
//    float xi, eta;  //坐标转换
//    xi = x + sinf(theta);
//    eta = y - 1 - cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(xi, eta);
//    if (pl.distance >= 2) {
//        u = 4 - sqrt(pow(pl.distance, 2) - 4);
//        if (u <= 0) {
//            t = mod2pi(std::atan2((4 - u) * xi - 2 * eta, -2 * xi + (u - 4) * eta));
//            v = mod2pi(t - theta);
//            if (t >= 0 && v >= 0) {
//                isoki = true;
//                temp_segm[0] = t;
//                temp_segm[1] = u;
//                temp_segm[2] = v;
//            }
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpRmSmLm(float x, float y, float theta, float temp_segm[3]) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float t = 0, u = 0, v = 0, r;
//    float xi, eta;  //坐标转换
//    xi = x - sinf(theta);
//    eta = y - 1 + cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(xi, eta);
//    if (pl.distance >= 2) {
//        r = sqrt(pow(pl.distance, 2) - 4);
//        u = 2 - r;
//        t = mod2pi(pl.angel + atan2(r, -2));
//        v = mod2pi(theta - M_PI_2 - t);
//        if (t >= 0 && u <= 0 && v <= 0) {
//            isoki = true;
//            temp_segm[0] = t;
//            temp_segm[1] = u;
//            temp_segm[2] = v;
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpRmSmRm(float x, float y, float theta, float *temp_segm) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float t = 0, u = 0, v = 0;
//    float xi, eta;  //坐标转换
//    xi = x + sinf(theta);
//    eta = y - 1 - cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(-eta, xi);
//    if (pl.distance >= 2) {
//        t = pl.angel;
//        u = 2 - pl.distance;
//        v = mod2pi(t + M_PI_2 - theta);
//        if (t >= 0 && u <= 0 && v <= 0) {
//            isoki = true;
//            temp_segm[0] = t;
//            temp_segm[1] = u;
//            temp_segm[2] = v;
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpRumLumRp(float x, float y, float theta, float *temp_segm) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float u = 0, rho;
//    float xi, eta;  //坐标转换
//    TO to{};
//    xi = x + sinf(theta);
//    eta = y - 1 - cosf(theta);
//    rho = (20 - pow(xi, 2) - pow(eta, 2)) / 16;
//    if (rho >= 0 && rho <= 1) {
//        u = -acosf(rho);
//        if (u >= -M_PI_2) {
//            to = tauOmega(u, u, xi, eta, theta);
//            if (to.tau >= 0 && to.omega >= 0) {
//                isoki = true;
//                temp_segm[0] = to.tau;
//                temp_segm[1] = u;
//                temp_segm[2] = to.omega;
//            }
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpRupLumRm(float x, float y, float theta, float *temp_segm) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float u = 0, rho;
//    float xi, eta;  //坐标转换
//    TO to{};
//    xi = x + sinf(theta);
//    eta = y - 1 - cosf(theta);
//    rho = (2 + sqrt(pow(xi, 2) + pow(eta, 2))) / 4;
//    if (rho <= 1) {
//        u = acosf(rho);
//        to = tauOmega(u, -u, xi, eta, theta);
//        if (to.tau >= 0 && to.omega <= 0) {
//            isoki = true;
//            temp_segm[0] = to.tau;
//            temp_segm[1] = u;
//            temp_segm[2] = to.omega;
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpSpLp(float x, float y, float theta, float *temp_segm) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float v = 0;
//    float xi, eta;  //坐标转换
//    xi = x - sinf(theta);
//    eta = y - 1 + cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(xi, eta);
//    if (pl.angel >= 0) {
//        v = mod2pi(theta - pl.angel);
//        if (v >= 0) {
//            isoki = true;
//            temp_segm[0] = pl.angel;
//            temp_segm[1] = pl.distance;
//            temp_segm[2] = v;
//        }
//    }
//    return isoki;
//}
//
//bool ReedsShepp::LpSpRp(float x, float y, float theta, float *temp_segm) {
//    //初始化
//    bool isoki = false;
//    temp_segm[0] = 0;
//    temp_segm[1] = 0;
//    temp_segm[2] = 0;
//    float t = 0, u = 0, v = 0, theta_tem;
//    float xi, eta;  //坐标转换
//    xi = x + sinf(theta);
//    eta = y - 1 - cosf(theta);
//    polar pl{};
//    pl = rect_to_polar(xi, eta);
//    if (pow(pl.distance, 2) >= 4) {
//        u = sqrt(pow(pl.distance, 2) - 4);
//        theta_tem = atan2(2, u);
//        t = mod2pi(pl.angel + theta_tem);
//        v = mod2pi(t - theta);
//        if (t >= 0 && v >= 0) {
//            isoki = true;
//            temp_segm[0] = t;
//            temp_segm[1] = u;
//            temp_segm[2] = v;
//        }
//    }
//    return isoki;
//}
//
///**************其它函数实现***************/
//polar rect_to_polar(float x, float y) {  //直角坐标转极坐标实现
//    polar pl{};
//    pl.distance = std::sqrt(x * x + y * y);
//    pl.angel = std::atan2(y, x);  // atan2（y,x）求的是y/x的反正切，其返回值为[-pi,+pi]之间的一个数。
//    return pl;
//}
//
//float mod2pi(float x) {  // x去周期到 -pi,pi之间
//    float v;
//    v = fmod(x, 2 * M_PI);  //取余 结果符号跟x
//    if (v < -M_PI)
//        v = v + 2 * M_PI;
//    else if (v > M_PI)
//        v = v - 2 * M_PI;
//    return v;
//}
//
//TO tauOmega(float u, float v, float xi, float eta, float theta) {
//    TO to{};
//    float delta = mod2pi(u - v);
//    float A = sinf(u) - sinf(delta);
//    float B = cosf(u) - cosf(delta) - 1;
//    float t1 = std::atan2(eta * A - xi * B, xi * A + eta * B);
//    float t2 = 2 * (cosf(delta) - cosf(v) - cosf(u)) + 3;
//    if (t2 < 0)
//        to.tau = mod2pi(t1 + M_PI);
//    else
//        to.tau = mod2pi(t1);
//    to.omega = mod2pi(to.tau - u + v - theta);
//    return to;
//}
//
//Pose2Dd global2local(const Pose2Dd &startpose, const Pose2Dd &endpose) {
//    Pose2Dd transEnd;
//    transEnd.x() = cosf(startpose.theta()) * (endpose.x() - startpose.x()) +
//                   sinf(startpose.theta()) * (endpose.y() - startpose.y());
//    transEnd.y() = -sinf(startpose.theta()) * (endpose.x() - startpose.x()) +
//                   cosf(startpose.theta()) * (endpose.y() - startpose.y());
//    transEnd.theta() = endpose.theta() - startpose.theta();
//    return transEnd;
//}
//void plot_start_traj(int start_num, std::vector<float> &start_path_x, std::vector<float> &start_path_y) {
//    //画起点轨迹
//    float s = 0, angle = 0, t;
//
//    for (int j = 0; j < start_num; ++j) {
//        t = j * 0.05;
//        s = 10 * (t + 1);
//        angle = t * 360;
//        start_path_x[j] = s * cosf(angle);
//        start_path_y[j] = s * sinf(angle);
//    }
//    // plt::plot(start_path_x,start_path_y,"bo");
//    // plt::plot(start_path_x,start_path_y,"b");
//    // plt::show();
//}
//
//void plot_sdf(std::vector<std::vector<float>> x, std::vector<std::vector<float>> y, std::vector<std::vector<float>> z) {
//    // matplotlibcpp::plot_surface(x,y,z);
//    // matplotlibcpp::show();
//}