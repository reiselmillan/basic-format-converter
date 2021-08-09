#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <array>
#include <vector>
#include <algorithm>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>



namespace props{
  // symmetry related ------------------------
  struct SPACEGROUP{
    std::string name = "";
    std::vector<std::string> opStrings;
    //std::vector<std::array< std::array<float, 4>,3>> opMatrices;
    unsigned int tableNumber = 0;
    unsigned int numOperations;
  };

  struct LatticeParams{
    LatticeParams(float ai, float bi, float ci, float alphai, float betai, float gammai)
    {
      a = ai; b = bi; c = ci; alpha = alphai; beta = betai; gamma = gammai;
    }
    LatticeParams(){}
    float a = 0;
    float b = 0;
    float c = 0;
    float alpha = 0;
    float beta = 0;
    float gamma = 0;
  };

  extern std::array<SPACEGROUP, 230> spaceGroups;
  extern std::map<std::string, std::array<std::array<float, 4>, 3>> symStrMatOperations;
  void setSpaceGroups();
  // symmetry related ------------------------
  
  const short NUMELEM = 100;
  extern bool init;
  //extern std::map<std::string, float> radios;
  extern float cylinders[NUMELEM];
  extern float masas[NUMELEM];
  extern std::array<glm::vec4, NUMELEM> colors;
  extern char atomicSymbols[NUMELEM][3];
  extern float radios[NUMELEM];
  extern float spheres[NUMELEM];
  	                       
  short atomicNumber(const char* symbol);
  void initProps();
  void getUnitCellVectorsComponents(LatticeParams &lp, glm::mat3 &cell);
  void getLatticeParameters(LatticeParams &lp, glm::mat3 &cell);

}; //end namespace prop


class Style
{
  public:
    Style(char* name);
    Style(char *fileName, std::string name);
    std::string name;
    float cylinders[props::NUMELEM];
    float masas[props::NUMELEM];
    std::array<glm::vec4, props::NUMELEM> colors;
    float radios[props::NUMELEM];
    float spheres[props::NUMELEM];
    bool drawWithLine[props::NUMELEM];

};

extern std::vector<Style> styles;


#endif // PROPERTIES_H
