#include<iostream>
#include "main.h"
#include "R_Shepp.h"

using namespace std;

int main() {
	float x =1 ;
	float y = 1;
	float theta = PI;
	PATH finalPath;
	finalPath = FindRSPath(x, y, theta);

	cout << "HelloWorld!"<< endl;
	system("pause");
	return 0;
}