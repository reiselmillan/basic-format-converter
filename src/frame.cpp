#include "frame.h"

using std::string;
using std::vector;
using std::array;
using std::cout; 
using std::endl;


Atom& Frame::addAtom(const char* symbol, const float x, const float y, const float z){
  Atom at;
  at.atomicNum = props::atomicNumber(symbol);
  at.coor = glm::vec3(x, y, z);
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(const char* symbol, glm::vec3 &coor){
  Atom at;
  at.atomicNum = props::atomicNumber(symbol);
  at.coor = coor;
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(short atomicNumber, glm::vec3 &coor){
  Atom at;
  at.atomicNum = atomicNumber;
  at.coor = coor;
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(Atom& atom){
  Atom at;
  at.atomicNum = atom.atomicNum;
  at.coor = atom.coor;
  strncpy(at.atomType, atom.atomType, 10);
  at.clearBonds();
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::atom(unsigned int i){
  return atoms_[i];
}

vector<Atom> Frame::deleteSelected(){
  std::vector<Atom> tempatoms;
  std::vector<Atom> deletedatoms;
  tempatoms.reserve(this->atoms_.size());
  for(auto &at:atoms_){
    if(at.selected){
      deletedatoms.emplace_back(at);
    }
    else{
      tempatoms.emplace_back(at);
    }
  }

  atoms_ = tempatoms; 
  return deletedatoms;
}

void Frame::removeAtom(unsigned int index){
  atoms_.erase(atoms_.begin() + index);
}

void Frame::setUnitCell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz ){

  this->cell[0][0] = xx;
  this->cell[0][1] = xy;
  this->cell[0][2] = xz;
  this->cell[1][0] = yx;
  this->cell[1][1] = yy;
  this->cell[1][2] = yz;
  this->cell[2][0] = zx;
  this->cell[2][1] = zy;
  this->cell[2][2] = zz;
  return;
}

void Frame::setUnitCell(glm::mat3 &cell) {
  this->cell = cell; 
  this->cellIsSet = true;  
}

void Frame::setUnitCell(props::LatticeParams &lp){
  props::getUnitCellVectorsComponents(lp, this->cell);
  this->cellIsSet;
}