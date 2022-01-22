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
  strcpy(at.atomType, symbol);
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(const char* symbol, glm::vec3 &coor){
  Atom at;
  at.atomicNum = props::atomicNumber(symbol);
  at.coor = coor;
  strcpy(at.atomType, symbol);
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(short atomicNumber, glm::vec3 &coor){
  Atom at;
  at.atomicNum = atomicNumber;
  at.coor = coor;
  strcpy(at.atomType, props::atomicSymbols[atomicNumber]);
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(short atomicNumber, const float x, const float y, const float z){
  Atom at;
  at.atomicNum = atomicNumber;
  at.coor = glm::vec3(x, y, z);
  strcpy(at.atomType, props::atomicSymbols[atomicNumber]);
  atoms_.emplace_back(at);
  return atoms_.back();
}

Atom& Frame::addAtom(Atom& atom){
  Atom at;
  at.atomicNum = atom.atomicNum;
  at.coor = atom.coor;
  strncpy(at.atomType, atom.atomType, 10);
  strncpy(at.residue, atom.residue, 10);
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

std::vector<std::string> Frame::orderByAtomType(){
  std::vector<Atom> tempatoms;
  tempatoms.reserve(this->atoms_.size());
  std::vector<std::string> atomtypes;
  //get the atom types
  for(auto &at:atoms_){
    if(std::find(atomtypes.begin(), atomtypes.end(), at.atomType) == atomtypes.end()){ //not been added
      atomtypes.push_back(at.atomType);
    }
  }
  for(auto &s:atomtypes){
    for(auto &at:atoms_){
      if(s.compare(at.atomType) == 0){
        tempatoms.push_back(at);
      }
    }
  }
  atoms_ = tempatoms;
  return atomtypes;
}


std::vector<AtomStatis> Frame::orderByAtomicNum(){
  std::vector<Atom> tempatoms;
  tempatoms.reserve(this->atoms_.size());
  std::vector<short> atomic_nums;
  std::vector<AtomStatis> statis;
  //get the atom types
  for(auto &at:atoms_){
    if(std::find(atomic_nums.begin(), atomic_nums.end(), at.atomicNum) == atomic_nums.end()){ //not been added
      atomic_nums.push_back(at.atomicNum);
      AtomStatis as; //create an AtomStatic instance , the number will be updated in the following loop
      as.at = at;
      statis.push_back(as);
    }
  }
  for(unsigned int i = 0; i < atomic_nums.size(); i++){
    for(auto &at:atoms_){
      if(atomic_nums[i] == at.atomicNum){
        tempatoms.push_back(at);
        statis[i].num += 1;
      }
    }
  }
  atoms_ = tempatoms;
  return statis;
}

std::vector<AtomStatis> Frame::orderByAtomicNum(const std::vector<short> &atomic_nums){
  std::vector<Atom> tempatoms;
  tempatoms.reserve(this->atoms_.size());
  std::vector<AtomStatis> statis;

  for(unsigned int i=0; i<atomic_nums.size(); i++){
    AtomStatis as;
    statis.push_back(as);
  }

  for(unsigned int i = 0; i < atomic_nums.size(); i++){
    for(auto &at:atoms_){
      if(atomic_nums[i] == at.atomicNum){
        tempatoms.push_back(at);
        if(statis[i].num == 0){
          statis[i].at = at;
        }
        statis[i].num += 1;
      }
    }
  }
  atoms_ = tempatoms;
  return statis;
}


void Frame::removeAtom(unsigned int index){
  atoms_.erase(atoms_.begin() + index);
}

void Frame::reserve(int n){
  atoms_.reserve(n);
}

void Frame::setUnitCell(glm::mat3 &cell) {
  this->cell_ = cell; 
  this->cellIsSet = true;  
}

void Frame::setUnitCell(props::LatticeParams &lp){
  props::getUnitCellVectorsComponents(lp, this->cell_);
  this->cellIsSet = true;
}

void Frame::setUnitCell(Cell &cell){
  this->cell_ = cell.vectorsComponets();
  this->cellIsSet = true;
}

float Frame::totalCharge(){
  float charge = 0;
  for(auto &at:atoms_){
    charge += at.charge;
  }
  return charge;
}

bool Frame::replicateAtoms(int x, int y, int z){
  if(! cellIsSet) return false;

  for(unsigned int i = 0; i < atoms_.size(); i++){
    PeriodicAtom PAtom;
    Atom& at = atoms_[i];
    PAtom.coor = at.coor + glm::vec3(x * cell().a(), 0.0f, 0.0f);
    // PAtom.coor = glm::vec3{5.0f, 5.0f, 5.0f};
    PAtom.oIndex = i;
    periodicAtoms_.push_back(PAtom);
  }

  return true;
}