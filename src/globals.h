#ifndef GLOBALS_H
#define GLOBALS_H
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <math/math.h>
#include <unistd.h>
#include <algorithm>
#include <thread>
#include <locale.h>
//#include<omp.h>
#ifdef Q_OS_UNIX
#endif

namespace glob{

  struct SPACEGROUP{
    std::string name = "";
    std::vector<std::string> opStrings;
    //std::vector<std::array< std::array<float, 4>,3>> opMatrices;
    unsigned int tableNumber = 0;
    unsigned int numOperations;
 };

  void initGlobalProperties();

}
#endif // GLOBALS_H

