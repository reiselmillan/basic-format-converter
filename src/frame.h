#ifndef FRAME_H
#define FRAME_H

#include "atom.h"

struct Bond
{
  int at1 = -1;
  int at2 = -1;
  float cylRadius = 0.2;
  float length = -1;

  bool operator==(const Bond &bd)const{
    if(at1 == bd.at1 && at2 == bd.at2){
        return true;
    }else{
        return false;
    }
  }

};

struct Angle
{
  int at1 = -1;
  int at2 = -1;
  int at3 = -1;
};

struct Molecule
{
  char name[10];
  int numMols = 0;
};

class Frame
{
  public:
    Frame(){};

    //variables
    //check whether the constraints should be writen in output file
    bool cartesian = true;
    bool cellIsSet = false;

    unsigned int step = 0;
    //unsigned int tableNumber = 1;

    float energy = 0.0;
    float time = 0.0;       // simulation time
    std::string comment = "";

    glm::vec3 centroid = {0.0, 0.0, 0.0};
    glm::mat3 cell;

    //functions
    Atom& addAtom(const char* symbol, const float x, const float y, const float z); //It might return a reference to the Atom instead of the index
    Atom& addAtom(const char* symbol, glm::vec3 &coor);
    Atom& addAtom(short atomicNumber, glm::vec3 &coor);
    Atom& addAtom(Atom& atom);
    unsigned int numAtoms(){return atoms_.size();}

    Atom& atom(unsigned int i);

    std::vector<Atom> deleteSelected();

    void removeAtom(unsigned int index);
    void removeDuplicates();

    void setUnitCell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz );
    void setUnitCell(glm::mat3 &cell);
    void setUnitCell(props::LatticeParams &lp);

  private:

    std::vector<Atom> atoms_;

};

#endif // FRAME_H
