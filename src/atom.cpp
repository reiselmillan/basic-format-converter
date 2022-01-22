#include "atom.h"

float BaseAtom::sphere(){
  // if(!props::init){
  //   props::initProps();
  // }
  return styles[this->styleIndex].spheres[atomicNum];

}

void BaseAtom::setSphere(const float s){
  styles[this->styleIndex].spheres[atomicNum] = s;
}

void BaseAtom::setRadius(const float s){
  styles[this->styleIndex].radios[atomicNum] = s;
}


float BaseAtom::radius(){
  // if(!props::init){
  //   props::initProps();
  // }
  return styles[this->styleIndex].radios[atomicNum];

}

float BaseAtom::cylRadius(){
  // if(!props::init){
  //   props::initProps();
  // }
  return styles[this->styleIndex].cylinders[atomicNum];

}

void BaseAtom::setCylRadius(float s){
  styles[this->styleIndex].cylinders[atomicNum] = s;
}


glm::vec4& BaseAtom::color(){
  // if(!props::init){
  //   props::initProps();
  // }
  return styles[this->styleIndex].colors[atomicNum];
}

unsigned short Atom::numBonds() {
  unsigned short n = 0;
  for(auto &bi:bondInds){
    if(bi != -1){
      n += 1;
    }
  }
  return n;
}

void Atom::clearBonds(){
  for(unsigned short i=0; i < 10; i++){
    bondInds[i] = -1;
  }
}

void Atom::addBond(short index){
  for(unsigned short i=0; i < 10; i++){
    if(bondInds[i] == -1){
      bondInds[i] = index;
      return;
    }
  }
}

void BaseAtom::setStyle(unsigned int index){
  if(index >= styles.size()){
    return;
  }
  this->styleIndex = index;
}


void BaseAtom::setStyle(std::string style){
  for(unsigned int i = 0; i<styles.size(); i++){
    if(styles[i].name.compare(style) == 0){
      this->styleIndex = i;
      return;
    }
  }
}

