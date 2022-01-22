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
  
  const short NUMELEM = 101;
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
    Style(std::string filePath);
    std::string name;
    float cylinders[props::NUMELEM] = { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, //first and second rows
                                      0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                      0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
                                      0.2, 0.2, 0.2, 0.2 };

    float masas[props::NUMELEM] = { 0.0, 1.007970, 0.0, 6.941, 0.0, 0.0, 12.0, 0.0, 15.99940, 
                                    0.0, 0.0, 0.0, 0.0, 26.98153, 28.08600, 30.97376, 0.0, 0.0, 
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                                    0.0};

    std::array<glm::vec4, props::NUMELEM> colors;

    float radios[props::NUMELEM] = {0.5, 0.5, 0.5, 1.0, 1.0, 0.9, 0.9, 0.8, 1.0, 0.7, 0.7,
                                    1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9,
                                    0.9,0.9,0.9,0.9,1.5,1.5,1.5,1.5,1.5,1.5,1.5};

    float spheres[props::NUMELEM] = { 0.3, 0.25, 0.25, 0.4, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, //first and second rows
                                      0.6, 0.6, 0.6, 0.3, 0.6, 0.6, 0.5, 0.5, 0.6, 0.6, 0.6, 
                                      0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
                                      0.6, 0.6, 0.6, 0.6 };
    bool drawWithLine[props::NUMELEM];

    bool defaultStyle = false;

    //funcs ------------------------
    void setProperties(short atNum, float sphere, float radius, float cyl, glm::vec4 &color);
    bool operator==(const Style& other);
    bool operator!=(const Style& other);

};

extern std::vector<Style> styles;


#endif // PROPERTIES_H
