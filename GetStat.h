#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include<fstream>
#include <chemfiles.hpp>
#include<iomanip>
#include<vector>

#include"./Static/Init.h"
#include"./Static/RDF.h"
#include"./Static/SSF.h"
#include"./Static/PSP.h"
#include"./Static/ANG.h"
#include"./Static/LB2.h"
#include"./Static/MID.h"


using namespace std;
using namespace chemfiles;
using std::string;
using std::stringstream;

void Properity(int Switch,int *NPart,int *NBead,int *NPoly,double *Box,double *Vol, vector<vector<double>>& RX, vector<vector<double>>& PX);

int main(int argc, char* argv[]);