#include "properties.h"
#include "space_gourps.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "math/math.h"

using props::NUMELEM;
using std::string;
using std::getline;

// var vars vars

//symmetry related
std::array<props::SPACEGROUP, 230> props::spaceGroups;
std::map<std::string, std::array<std::array<float, 4>, 3>> props::symStrMatOperations;

bool props::init = false;
float props::cylinders[NUMELEM];
std::array<glm::vec4, NUMELEM> props::colors;


float props::radios[NUMELEM] = {0.5, 0.5, 0.5, 1.0, 1.0, 0.9, 0.9, 0.8, 1.0, 0.7, 0.7,
                        1.1, 1.1, 1.1, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9,0.9,0.9,0.9,0.9,0.9,1.5,1.5,1.5,1.5,1.5,1.5,1.5};
float props::spheres[NUMELEM] = { 0.3, 0.25, 0.25, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, //first and second rows
                            0.6, 0.6, 0.6, 0.3, 0.6, 0.6, 0.5, 0.5,
          0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6
};

float props::masas[NUMELEM] = {
  0.0, 1.007970, 0.0, 6.941, 0.0, 0.0, 12.0, 0.0, 15.99940, 0.0, 0.0, 0.0, 0.0, 26.98153, 28.08600,
  30.97376, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

char props::atomicSymbols[NUMELEM][3] = {
    "*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn","Fr", "Ra"
};

std::vector<Style> styles;

// FUNCS FUNCS FUNCS
bool readSettings(){
  return false;
}


short props::atomicNumber(const char* symbol){
    for(short i = 0; i < NUMELEM; i++){
      if(strcmp(symbol, atomicSymbols[i]) == 0){
        return i;
      }
    }
    return -1;
  }
  
void props::initProps(){
  if(!readSettings()){
    colors[0] = glm::vec4(0.1, 0.5, 0.9, 1.0);
    colors[1] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[2] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[3] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[4] = glm::vec4(0.4, 0.4, 0.0, 1.0);
    colors[5] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[6] = glm::vec4(0.4, 0.4, 0.4, 1.0);
    colors[7] = glm::vec4(0.0, 0.0, 0.7, 1.0);
    colors[8] = glm::vec4(0.7, 0.0, 0.0, 1.0);
    colors[9] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[10] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[11] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[12] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[13] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[14] = glm::vec4(0.7, 0.5, 0.0, 1.0);
    colors[15] = glm::vec4(0.56, 0.56, 0.0, 1.0);
    colors[16] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[17] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[18] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[19] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[20] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[21] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[22] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[23] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[24] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[25] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[26] = glm::vec4(0.14,0.45,0.19, 1.0);
    colors[27] = glm::vec4(0.5, 0.1, 0.8, 1.0);
    colors[28] = glm::vec4(0.14,0.45,0.19, 1.0);
    colors[29] = glm::vec4(0.0, 0.47, 0.0, 1.0);
    colors[30] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[31] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[32] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[33] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    colors[34] = glm::vec4(0.7, 0.7, 0.7, 1.0);
    
    //setup the styles
    Style st((char*)"style_0");
    styles.push_back(st);
  }
  props::setSpaceGroups();
  props::init = true;
}

void props::getLatticeParameters(LatticeParams &lp, glm::mat3 &cell){
  lp.a = sqrt(cell[0][0] * cell[0][0] + 
              cell[0][1] * cell[0][1] +
              cell[0][2] * cell[0][2]);
  lp.b = sqrt(cell[1][0] * cell[1][0] + 
              cell[1][1] * cell[1][1] +
              cell[1][2] * cell[1][2]);
  lp.c = sqrt(cell[2][0] * cell[2][0] + 
              cell[2][1] * cell[2][1] +
              cell[2][2] * cell[2][2]);

  glm::vec3 zero(0,0,0);
  lp.alpha = math::calcAngle(cell[1], zero, cell[2]);
  lp.beta = math::calcAngle(cell[0], zero, cell[2]);
  lp.gamma = math::calcAngle(cell[0], zero, cell[1]);
}

void props::getUnitCellVectorsComponents(props::LatticeParams &lp, glm::mat3 &cell){
  float PI = 3.14159265;

  cell[0][0] = lp.a;
  cell[0][1] = 0.0;
  cell[0][2] = 0.0;
  float ang_rad = PI/180*lp.gamma;  //gamma to radians
  float cosang = cos(ang_rad);
  float sinang = sin(ang_rad);
  cell[1][0] = lp.b * cosang;
  cell[1][1] = lp.b * sinang;
  cell[1][2] = 0.0;

  //z vector components
  ang_rad = PI/180.0*lp.beta;  //
  cosang = cos(ang_rad);
  cell[2][0] = (cosang*lp.a*lp.c)/cell[0][0];
  ang_rad = PI/180.0*lp.alpha;  //
  cosang = cos(ang_rad);
  cell[2][1] = (cosang*lp.b*lp.c - cell[1][0]*cell[2][0])/cell[1][1];
  cell[2][2] = sqrt(lp.c * lp.c - cell[2][0]*cell[2][0] - cell[2][1]*cell[2][1]);
}

void props::setSpaceGroups(){
  // FILE *file;
  // char buffer[1024];
  int tableNumber, numOperations;
  char groupName[1024];
  char operationString[20];
  //std::array<SPACEGROUP, 230> spaceGroups;
  props::SPACEGROUP sp;
  float c1, c2, c3, c4;
  std::array< std::array<float, 4>,3> opMati; //temporarily stores the operation matrix
  unsigned int count = 0;
  //std::string spaceGroupFile = "space_groups";

  std::istringstream f(space_groups_str);
  std::string line;

  // file = fopen(spaceGroupFile.c_str(), "r");
  // while(fgets(buffer, 1024, file) != NULL){
  while (std::getline(f, line)) {
    // sscanf(buffer, "%i  %i  %s", &tableNumber, &numOperations, groupName);
    sscanf(line.c_str(), "%i  %i  %s", &tableNumber, &numOperations, groupName);
    sp.tableNumber = tableNumber;
    sp.numOperations = numOperations;
    sp.name = groupName;
    
    //read all the operations
    for(int i=0; i< numOperations; ++i){
        // read the string
        // fgets(buffer, 1024, file);
        getline(f, line);
        // sscanf(buffer, "%s", operationString);
        sscanf(line.c_str(), "%s", operationString);
        sp.opStrings.push_back(operationString);
        
        //read the matrix
        for(int j = 0; j < 3; ++j){
          // fgets(buffer, 1024, file);
          getline(f, line);
          sscanf(line.c_str(), "%f  %f   %f   %f", &c1, &c2, &c3, &c4);
          opMati[j][0] = c1; opMati[j][1] = c2; opMati[j][2] = c3; opMati[j][3] = c4; 
        }
        props::symStrMatOperations[operationString] = opMati;
    }
    props::spaceGroups[count] = sp;
    sp.opStrings.clear();
    //sp.opMatrices.clear();
    count += 1;
  }
  return;
}


Style::Style(char* fileName, std::string name){
  this->name = name;
  FILE *file;
  char buffer[1024], symbol[3];
  float radius, sphere, cyl, color1, color2, color3, color4;

  file = fopen(fileName, "r");
  if (file == NULL){
    printf("local setting elements file not found\n");
  }

  while(fgets(buffer, 1024, file) != NULL){
    if(strncmp(buffer, "#", 1) == 0) continue;

    int n = sscanf(buffer, "%s %f %f %f %f %f %f %f", symbol, &sphere, &radius,
	                                          &cyl, &color1, &color2, &color3, &color4);
	  if(n != 8) continue;
    short atNum = props::atomicNumber(symbol);
    if(atNum == -1)continue;
    
    spheres[atNum] = sphere;
    radios[atNum] = radius;
    cylinders[atNum] = cyl;
    colors[atNum][0] = color1;
    colors[atNum][1] = color2;
    colors[atNum][2] = color3;
    colors[atNum][3] = color4;
  }

  fclose(file);
}


Style::Style(char* name){
  this->name = name;
  for(short i = 0; i < NUMELEM; i++){
    spheres[i] = props::spheres[i];
    radios[i] = props::radios[i];
    cylinders[i] = props::cylinders[i];
    colors[i][0] = props::colors[i][0];
    colors[i][1] = props::colors[i][1];
    colors[i][2] = props::colors[i][2];
    colors[i][3] = props::colors[i][3];
  }
}
