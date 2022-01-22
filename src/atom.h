// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "properties.h"

class BaseAtom
{
  public:
    BaseAtom(){};

    glm::vec3 coor = glm::vec3(0.0f, 0.0f, 0.0f);
    short atomicNum = 0;

    //functions
    const char* symbol(){return props::atomicSymbols[atomicNum];}
    void setSymbol(const char* s){atomicNum = props::atomicNumber(s);}

    float sphere();
    void setSphere(const float s);
    void setRadius(const float r);

    void setStyle(unsigned int id);
    void setStyle(std::string style); 
    Style& style(){return styles[styleIndex];}

    float radius();
    float cylRadius();
    void setCylRadius(float s);
    glm::vec4& color();

  private:
    short styleIndex = 0;


};

class Atom: public BaseAtom
{
  public:
    Atom(){};

    char atomType[10] = "         ";
    const char* getAtomType(){ return atomType; }  //only the python wrapper
    void setAtomType(const char* at){strcpy(atomType, at);} //only for the python wrapper

    char residue[10] =  "RES      ";
    void setResidue(const char* res){strcpy(residue, res);} //only for the python wrapper
    float charge = 0;
    float chemShift = 0.0;
    bool asymUnit = true;

    short molType = -1;
    short molIndex = -1;

    std::array<bool, 3> fixed = {false, false, false};

    std::array<short, 10> bondInds = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

    bool selected = false;
    bool hidden = false;
    bool locked = false;

    //functions -----------------
    //Indices of the bonds of the corresponding atom
    unsigned short numBonds();
    void clearBonds();
    void addBond(short i);

    bool operator==(const Atom &at)const{
      if(coor.x == at.coor.x && coor.y == at.coor.y && coor.z == at.coor.z && atomicNum == at.atomicNum){
          return true;
      }else{
          return false;
      }
    }

  private:


};

class PeriodicAtom : public BaseAtom
{
  public:
    PeriodicAtom():BaseAtom(){};
    int oIndex = -1;

};