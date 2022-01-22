#include "trajectory.h"
#include <math/math.h>
#include <chrono>

using std::string;
using std::vector;
using std::array;
using std::cout; 
using std::endl;



Trajectory::Trajectory()
{
  this->initCellLines();
}

void Trajectory::addBond(const unsigned int i, const unsigned int j, const float &bdist, bool updateAtoms){
  Bond b;
  if(atom(i).cylRadius() >= atom(j).cylRadius()){
    b.cylRadius = atom(j).cylRadius();
  }
  else{
    b.cylRadius = atom(i).cylRadius();
  }
  b.at1 = i;
  b.at2 = j;
  b.length = bdist;
  bonds_.push_back(b);
  if(updateAtoms){
    this->atom(i).addBond(bonds_.size() - 1);
    this->atom(j).addBond(bonds_.size() - 1);
  }
}


void Trajectory::addBond(const unsigned int i, const unsigned int j, const float &bdist, glm::vec3& offset){
  Bond b;
  if(atom(i).cylRadius() >= atom(j).cylRadius()){
    b.cylRadius = atom(j).cylRadius();
  }
  else{
    b.cylRadius = atom(i).cylRadius();
  }
  b.at1 = i;
  b.at2 = j;
  b.length = bdist;
  b.offset = offset;
  bonds_.push_back(b);
}

void Trajectory::addBond(const unsigned int i, const unsigned int j, 
                          const float &bdist, const float &cylRadius, bool updateAtoms){
  // I am not using it right now
  Bond b;
  b.at1 = i;
  b.at2 = j;
  b.cylRadius = cylRadius;
  b.length = bdist;
  bonds_.push_back(b);
  
  if(updateAtoms){
    this->atom(i).addBond(bonds_.size() - 1);
    this->atom(j).addBond(bonds_.size() - 1);
  }
}

void Trajectory::addBondedAtomsGeometrically(unsigned int index, short n, short atomicNumber, float blength){       
  Atom &at = this->atom(index);
  for(short i = 0; i < n; i++){
	  glm::vec3 aveVector(0.0f, 0.0f, 0.0f);
	  //iterate over all bonds to sum the vector of bonds to determine average
	  for(auto &aindex:this->getBondedIndices(index)){  // iterate over the atomic neighbors of atom index
      Atom &at2 = this->atom(aindex);
      //get position vectors with respect to atom with index index
      glm::vec3 relvec = at2.coor - at.coor; 
      aveVector += relvec;
	  }

    Frame& fr = currentFrame();

	  if(aveVector[0] == 0.0 && aveVector[1] == 0.0 && aveVector[2] == 0.0){
      // add the atom displaced blength along the x axis
      glm::vec3 newcoor = at.coor + glm::vec3(blength, 0.0f, 0.0f);
	    Atom& at = fr.addAtom(atomicNumber, newcoor); 
      strcpy(at.atomType, props::atomicSymbols[atomicNumber]);
		  //add bond first atom
		  this->addBond(index, this->numAtoms() - 1, blength); //add bond with the indices of each atom
	  }
	  else{
		  //normalize vector
      aveVector = glm::normalize(aveVector);
		  // change the direction of the vector to opposite direction and scale to blength
      aveVector *= -blength;
		  // translate to position relative to atom at
	    aveVector += at.coor;

      //add atom and bond it to atom at
      Atom& at = fr.addAtom(atomicNumber, aveVector);  
      strcpy(at.atomType, props::atomicSymbols[atomicNumber]); 
		  //add bond first atom
		  this->addBond(index, this->numAtoms() - 1, blength); //add bond with the indices of each atom

	  } 
	} // end for n
  // this->generateBonds();

}

void Trajectory::addHydrogen(const unsigned int index, const float blength){   
  if(cfi_ == -1) return;    
  Atom &at = this->atom(index);
  Frame& fr = currentFrame();
  //define atoms to add
  unsigned int atomsToAdd = 0;
  //int geometry = 0; //tetrahedro will use it later
  if(at.atomicNum == 8 || at.atomicNum == 16){
    atomsToAdd = 2;
  }
  else if(at.atomicNum == 6 || at.atomicNum == 14){
    atomsToAdd = 4;
  }

  std::vector<unsigned int> ni = this->getBondedIndices(index);

  for(unsigned int i = (unsigned int)ni.size(); i < atomsToAdd; i++){
    ni = this->getBondedIndices(index); // update the neighbors
    printf("neighbors %i %i\n", (int)ni.size(), i);
    
    if(i == 0){  //first hydrogen atom to add
      glm::vec3 relvec = glm::normalize(glm::vec3{1.0f, 1.0f, 1.0f});
      glm::vec3 ncoor = at.coor + relvec * blength;
      Atom& at = fr.addAtom(1, ncoor); 
      strncpy(at.atomType, "H", 2);
      //add bond with the indices of each atom and update the atoms bonds indices
      this->addBond(index, this->numAtoms() - 1, blength, true); 
    }
    else if(i == 1){  //2nd hydrogen atom to add
      Atom &at2 = this->atom(ni[0]);
      //get position vectors with respect to atom with index index
      glm::vec3 relvec = at2.coor - at.coor; 
      relvec = glm::normalize(relvec);
      glm::vec3 rotVec = glm::cross(relvec, glm::vec3{1.0f, 0.0f, 0.0f});
      glm::mat3 RotationMatrix;
      math::rotationMatrix3(rotVec, 1.9094, RotationMatrix); //1.9094 is 109.4 in radians
      relvec = RotationMatrix * relvec;
      relvec *= blength;
      relvec += at.coor;
      Atom& addedAt = fr.addAtom(1, relvec); 
      strncpy(addedAt.atomType, "H", 2);
      //add bond with the indices of each atom and update the atoms bonds indices
      this->addBond(index, this->numAtoms() - 1, blength, true); 
    }
    else if(i == 2){ // have to add the 3rd
      Atom& atn2 = atom(ni[0]);
      Atom& atn3 = atom(ni[1]);
      glm::vec3 relvec1 = (atn2.coor - at.coor);
      glm::vec3 relvec2 = (atn3.coor - at.coor);
      glm::vec3 aveVec =  relvec1 + relvec2 + glm::vec3{0.01, 0.01, 0.01}; //average vector of the two bonds
      //normalize it
      aveVec = glm::normalize(aveVec);
      glm::vec3 planeVec = glm::cross(relvec1, relvec2);
      planeVec = glm::normalize(planeVec);
      glm::vec3 rotVec = glm::cross(planeVec, aveVec);

      glm::mat3 RotationMatrix;
      math::rotationMatrix3(rotVec, 1.9094, RotationMatrix); //1.9094 is 109.4 in radians
      aveVec = RotationMatrix * aveVec;
      aveVec *= blength;
      aveVec += at.coor;
      Atom& addedAt = fr.addAtom(1, aveVec); 
      strncpy(addedAt.atomType, "H", 2);
      //add bond with the indices of each atom and update the atoms bonds indices
      this->addBond(index, this->numAtoms() - 1, blength, true); 

    }
    else if(i == 3){ // have to add the 3rd
      glm::vec3 aveVec;
      for(auto &nIndex:ni){
        Atom& atn = atom(nIndex);
        aveVec += atn.coor - at.coor;
      }
      //normalize it
      aveVec = glm::normalize(aveVec);
      aveVec *= -blength;
      aveVec += at.coor;
      Atom& addedAt = fr.addAtom(1, aveVec); 
      strncpy(addedAt.atomType, "H", 2);
      //add bond with the indices of each atom and update the atoms bonds indices
      this->addBond(index, this->numAtoms() - 1, blength, true); 

    }
  }
  

}

void Trajectory::addFrame(){
  Frame fr;
  frames_.emplace_back(fr);
  if(cfi_ == -1) cfi_ = 0;
  //this->cfi_ = (int)frames_.size() - 1;
}

unsigned int Trajectory::addFrame(Frame &fr){
  frames_.emplace_back(fr);
  if(cfi_ == -1) cfi_ = 0;
  return frames_.size() -1;
  //this->cfi_ = (int)frames_.size() - 1;
}


bool Trajectory::atBondDistance(const unsigned int i, const unsigned int j, float &bdist){
  float sradio = atom(i).radius() + atom(j).radius();

  float sum2 = glm::distance(atom(i).coor, atom(j).coor);

  // float sradio2 = sradio * sradio;
  // glm::vec3 diff = (atom(i).coor - atom(j).coor);
  // glm::vec3 diff2 = diff * diff;
  // float sum2 = diff2[0] + diff2[1] + diff2[2];
  if( sum2 > sradio){
    return false;
  }
  else{
    // bdist = sqrt(sum2);
    bdist = sum2;
  }
     
  return true;
}


bool Trajectory::atBondDistancePBC(const unsigned int i, const unsigned int j, vector<glm::vec3> &trans){
  float sradio = atom(i).radius() + atom(j).radius();

  //go periodic
  for(auto &t:trans){
    glm::vec3 ncoor = t + atom(j).coor;
    float sum2 = glm::distance(atom(i).coor, ncoor);
    if(sum2 <= sradio){
      addBond(i, j, sum2, t);
      // bdist = sum2;
      return true;
    }
  } //end for periodic
   
  return false;
}

//return angle bewteen 3 atoms
float Trajectory::angle(const unsigned int i, const unsigned int j, const unsigned int k){
  /**
   * Esto es un comentario
   */

  glm::vec3 v1 = frames_[cfi_].atom(i).coor - frames_[cfi_].atom(j).coor;
  glm::vec3 v2 = frames_[cfi_].atom(k).coor - frames_[cfi_].atom(j).coor;

  float angle = acos(glm::dot( glm::normalize(v1), glm::normalize(v2)));
  return angle*180/3.14159265359;
}

//calculate distance between two atom
float Trajectory::distance(const unsigned int i, const unsigned int j){
  //auto start = std::chrono::high_resolution_clock::now();
	
  float d = glm::distance(frames_[cfi_].atom(i).coor, frames_[cfi_].atom(j).coor);

  //auto stop = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds> (stop-start);
  //std::cout<<duration.count()<<std::endl;

  return d;
}

void Trajectory::bringInsideCellAllFrames(){
  int prevCfi = cfi_;  //store the current frame index
  for(unsigned int i = 0; i < this->numFrames(); i++){
    this->setCurrentFrameIndex(i);
    this->bringInsideCell();
  }
  cfi_ = prevCfi; //restore the initial current frame
}

void Trajectory::bringInsideCell(){
  if(cfi_ == -1) { printf("Traj has no Frames\n"); return;} //if no frames return
  if(!this->currentFrame().cellIsSet){
    printf("Currrent Frame does not have a defined unit cell\n");
    return;
  }
  //tengo que verificar si hay celda o no
  //hacer un warning cuando alguno es mayor de 2
  this->cartesianToFractional();
  for(unsigned int i=0; i < this->numAtoms(); i++){
    float zero = 0.0;
    float one = 1.0;
    glm::vec3 &coor = this->atom(i).coor;
    for(int j=0;j<3;j++){
      while(coor[j] > one){
        coor[j] -= one;
      }
      while(coor[j] < zero){
        coor[j] += one;
      }
      if(coor[j] < 0.00001 && coor[j] > -0.00001) coor[j] = 0.0f;
    }
  }
  this->fractionalToCartesian();
}

glm::vec3 Trajectory::calcCentroid(){
  this->centroid_ = glm::vec3(0.0,0.0,0.0);
  for(unsigned int i=0; i < this->numAtoms(); i++){
    centroid_ += this->atom(i).coor;
  }
  this->centroid_ /= this->numAtoms();
  return centroid_;
}

void Trajectory::cartToFracAllFrames(){
  int prevCfi = cfi_;
  for(unsigned int i = 0; i<this->frames_.size(); i++){
    cfi_ = i;
    this->cartesianToFractional();
  }
  cfi_ = prevCfi;
}

void Trajectory::cartesianToFractional(){
    if(!cartesian()) return;
    if(cfi_ == -1) return;

    glm::mat3 matrix;
    float suma[3];
    setCartFracMatrix(matrix);

    // multiplico las coordenadas ****************************************
    for(unsigned int i=0; i<this->numAtoms(); i++){
      Atom &at = this->atom(i);
	
      suma[0] = 0.0;
      suma[1] = 0.0;
      suma[2] = 0.0;
      for(unsigned int row=0;row<3;row++){
          for(unsigned int rowcol=0;rowcol<3;rowcol++){
              suma[row] += at.coor[rowcol]*matrix[row][rowcol];
          }
      }
      at.coor[0] = suma[0];
      at.coor[1] = suma[1];
      at.coor[2] = suma[2];
    }
    // multiplico las coordenadas ****************************************

    this->setCartesian(false);
}

void Trajectory::changeSphereSelected(const float s){ 
  for(unsigned int i=0; i < this->numAtoms(); i++){
    if(this->atom(i).selected){
      this->atom(i).setSphere(s);
    }
  }
}

bool Trajectory::checkEqualCoors(glm::vec3 &a1, glm::vec3 &a2, float umb){
    for(unsigned int i=0;i<3;i++){
        if(abs(a1[i] - a2[i]) > umb){
            return false;
        }
    }
    return true;
}

bool Trajectory::checkEqualReplicated(glm::vec3 &a1, glm::vec3 &a2, float umb){
  for(unsigned int i=0;i<3;i++){
	//if (a1[i] == 0.0f) a1[i] = 0.0f;
	//if (a2[i] == 0.0f) a2[i] = 0.0f;
    if(abs(a1[i] - a2[i]) > 1.0f - umb && abs(a1[i] - a2[i]) < 1.0f + umb){
      return true;
    }
  }
  return false;
}

short Trajectory::checkRingClosed(RingSearch rs, vector<unsigned int>&inds, vector<unsigned int> checked){
  printf("searching in level %i , atom Index %i\n", rs.level, rs.atIndex);
  checked.push_back(rs.atIndex);
  //Caution !! this is a recursive function
  Atom &at = this->atom(rs.atIndex);
  if(at.locked) return 2;
  //checked.push_back(atIndex); //added to checked, before entering sublevels, otherwise the algorithm returns to the first point

  //temp neighbors and code, to decide later which will be return based on the codes
  std::vector<std::vector<unsigned int>> ringCandidates;
  std::vector<short> codes;

  //start search bonds
  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    //get bond reference
    Bond &b = bonds_[bi];
    int neighbor = this->getAtomNeighborIndex(rs.atIndex, b);
    if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings

    //not going back the same path to target index. This will exclude the target in the first level
    //printf("neigh %i    previndex %i\n", neighbor, rs.prevAtIndex);
    if(neighbor == rs.prevAtIndex) continue; 
    
    if(rs.level == rs.ringSize){ // matched level, correct ring size, no more recursive calls
      if(neighbor == rs.atIndexTarget){ //will always be true for level 1, but excluded in the previous conditional
        inds.push_back(neighbor);
        return 0; // useless to continue the loop, Found my guy.
      }
       continue; //do not call recursive but continue with other neighbors
    }

    if(neighbor == rs.atIndexTarget){ //found smaller ring
      std::vector<unsigned int> tempCandidates;
      codes.push_back(1);
      tempCandidates.push_back(neighbor);
      ringCandidates.push_back(tempCandidates);
      continue;
    }

    //cannot be before because conflicts with checking the target index
    if(std::find(checked.begin(), checked.end(), neighbor) != checked.end()) continue; 
    
    //default call recursive
    std::vector<unsigned int> tempCandidates;
    RingSearch rs2(rs);
    rs2.level = rs.level + 1;
    rs2.atIndex = neighbor;
    rs2.prevAtIndex = rs.atIndex;

    short ringSizeCode = checkRingClosed(rs2,  tempCandidates, checked);
    //add candidates if match right of smaller ring size
    if(ringSizeCode == 0 || ringSizeCode == 1){  //right or smaller rings // I should remove small rings
      codes.push_back(ringSizeCode);
      tempCandidates.push_back(neighbor);
      ringCandidates.push_back(tempCandidates);
    }
  } // end for loop

  //check for smaller ring sizes. It there is one smaller ring the function return 2
  short returnCode = 2;
  for(unsigned int i=0; i<codes.size(); i++){
    // if(codes[i] == 1){ //there is one smaller ring, return false
    //   printf("smaller ring found in level %i\n", rs.level+1);
    //   return 1;
    // }
    if(codes[i] == 0){ //right ring found, store the vector index to add to final inds reference
      for(auto &ind:ringCandidates[i]){
        inds.push_back(ind);
      }

      returnCode = 0; //there is at least one right size ring
    }
  }

  return returnCode;
/*
  printf("CHECKING LEVEL: %i index: %i  %s\n", level, atIndex, at.symbol());
  if(level == 1){  //must check match
    for(auto &bi:at.bondInds){
      if(bi == -1) break;
      //get bond reference
      Bond &b = bonds_[bi];
      int neighbor = this->getAtomNeighborIndex(atIndex, b);
      if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings
      //check if any of the atom indices of the bond matches the target
      if(neighbor == (int)atIndexTarget){
        //do not have to append it, it is the initial one, appended outside
        printf("found level:%i  , %i  %i\n", level, neighbor, atIndexTarget);
        return true; //return true
      }
    }
  }
  else{  //level is not zero, should not match
    //decrease one level and search recursively
    level -= 1;
    for(auto &bi:at.bondInds){
      if(bi == -1) break;
      //get bond reference
      Bond &b = bonds_[bi];
      //get the neighbor
      int neighbor = this->getAtomNeighborIndex(atIndex, b);
      if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings
      
      // exclude if already analysed, avoid returning to initial atom through smaller ring
      if(std::find(checked.begin(), checked.end(), neighbor) != checked.end()) continue; 

      bool foundCloseRing = checkRingClosed(level, neighbor, atIndexTarget, checked, inds);
      //add neighbor index and return true if found match with target
      if(foundCloseRing){
        inds.push_back(neighbor);
        return true;
      }
    } // end for bond indices
  }
  //return false, not found any match
  return false;
  */
} // end recursive function


float Trajectory::energy(){
  if(cfi_ == -1) return 0;
  return frames_[cfi_].energy;
}

void Trajectory::findRing(short ringSize, unsigned int atIndex, vector<unsigned int>&inds){
  RingSearch rs(atIndex, atIndex, ringSize);
  vector<unsigned int> checked;
  this->checkRingClosed(rs, inds, checked);
  printf("ring vector size %i\n", (int)inds.size());
  return;
}

void Trajectory::findMoleculesFromScratch(){
  molecules_.clear(); //from scratch
  //first store all the previously unlocked atoms and unlock the locked ones
  //all atoms must be unlocked before searching for the molecules
  vector<unsigned int > unlockedInds;
  unlockedInds.reserve(numAtoms());
  for(unsigned int i = 0; i < numAtoms(); i++){
    if(!this->atom(i).locked) unlockedInds.push_back(i);
    else this->atom(i).locked = false;
  }
  //first store all the previously locked atoms
  int counterName = 0;
  for(unsigned int i=0; i < numAtoms(); i++){
    if(atom(i).locked) continue;
    //pass the index i to start looking for the molecule 
    //and lock all the atoms found. Those will be skipped in next iterations 
    //and will not increase the number of molecules
    Molecule mol;
    findMoleculeRecursively(i, mol.inds, true);
    //printf("molecule n has %i atoms\n", mol.inds.size());
    string name = "m_";
    name += std::to_string(counterName);
    strncpy(mol.name, name.c_str(), 10);
    molecules_.push_back(mol); 
    counterName += 1;
  }

  //restore the previously unlocked atoms
  for(auto i:unlockedInds){
    this->atom(i).locked = false;
  }

}

void Trajectory::findMoleculeRecursively(unsigned int atIndex, vector<unsigned int>&inds, bool lock){
  inds.push_back(atIndex);
  Atom &at = this->atom(atIndex);
  if(lock) {at.locked = true;}

  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    //get bond reference
    Bond &b = bonds_[bi];
    //get the neighbor
    int neighbor = this->getAtomNeighborIndex(atIndex, b);
    if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings
    //check if neighbor has already been added to vector inds
    auto it = std::find(inds.begin(), inds.end(), neighbor);
    if(it == inds.end()){ //only add if it is not found in vector inds
      // inds.push_back(neighbor);
      this->findMoleculeRecursively(neighbor, inds, lock);
    }
  } // end for bond indices
}

void Trajectory::fractToCartAllFrames(){
  int prevCfi = cfi_;
  for(unsigned int i = 0; i<this->frames_.size(); i++){
    cfi_ = i;
    this->fractionalToCartesian();
  }
  cfi_ = prevCfi;
}

void Trajectory::fractionalToCartesian(){
    if(cartesian()) return;
    if(cfi_ == -1) return;

    glm::mat3 matrix;
    float suma[3];
    setFracCartMatrix(matrix);
    // multiplico las coordenadas ****************************************
    for(unsigned int i=0; i<this->numAtoms(); i++){
      Atom &at = this->atom(i);
      suma[0] = 0.0;
      suma[1] = 0.0;
      suma[2] = 0.0;
      for(unsigned int row=0;row<3;row++){
        for(unsigned int rowcol=0;rowcol<3;rowcol++){
          suma[row] += at.coor[rowcol]*matrix[row][rowcol];
        }
      }
      at.coor[0] = suma[0];
      at.coor[1] = suma[1];
      at.coor[2] = suma[2];
    }
    // multiplico las coordenadas ****************************************
    // set cartesian true, to know the coordinumberOfAtoms_es are cartesian and not direc
    this->setCartesian(true);
}

//return reference to frame object
Frame& Trajectory::frame(const unsigned int i){
  if(frames_.size() == 0) {
    throw std::invalid_argument( "Trajectory has no frames\n" );
  }
  if(i >= frames_.size() || (int)i < 0){
    return frames_[0];
  }
  else{
    return frames_[i];
  }
}


void Trajectory::generateBonds(){
  /*
  generates bond from scratch, using the information of the atomic properties, that is,
  cylinder radius and bond radius
  */
  // float totclear = 0, totadd = 0, totabd = 0;
  vector<Bond> tempbonds = bonds_;
  this->bonds_.clear();
  
  //auto start0 = std::chrono::high_resolution_clock::now();

  //clear all atom bonds
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    this->atom(i).clearBonds();
  }
  //clear all atom bonds

  //take care of locked atoms first
  for(auto &bd:tempbonds){
    Atom& at1 = atom(bd.at1);
    Atom& at2 = atom(bd.at2);
    if(at1.locked && at2.locked){
      bonds_.push_back(bd);
    }
  }

  for(unsigned int i = 0; i < numAtoms(); i++){
    //printf("loop first level\n");
    // this->atom(i).clearBonds();  //clear the bonds of atom i
    for(unsigned int j = 0; j < i; j++){
      if(this->atom(i).locked && this->atom(j).locked){
        continue;
      }
      float length;
      bool abd = atBondDistance(i, j, length);
      if(abd){
        addBond(i, j, length);
      }
    }

  }
  this->updateAtomsBondsIndices();
  //auto  stop0 = std::chrono::high_resolution_clock::now();
  //auto duration0 = std::chrono::duration_cast<std::chrono::microseconds> (stop0-start0);
	//std::cout<<"total bonds duration:"<<duration0.count()*0.001<<"milisec "<<std::endl;
}


void Trajectory::generateNewBonds(std::vector<unsigned int> &inds){
  /*
  generates new bonds from scratch. Old bonds are not deleted
  */
  for(unsigned int i = 0; i < inds.size(); i++){
    for(unsigned int j = 0; j < i; j++){
      if(this->atom(inds[i]).locked && this->atom(inds[j]).locked){
        continue;
      }
      float length;
      bool abd = atBondDistance(inds[i], inds[j], length);
      if(abd){
        addBond(inds[i], inds[j], length);
        this->atom(inds[i]).addBond(this->bonds_.size()-1);       //add bond to atom
        this->atom(inds[j]).addBond(this->bonds_.size()-1);       //add bond to atom
      }
    }
  }
}


void Trajectory::generateBondsPBC(){
  //create PBC translations
  std::vector<glm::vec3> trans;
  Cell cell(frames_[cfi_].cell_);
  float a = cell.a();
  float b = cell.b();
  float c = cell.c();
  
  trans.push_back(glm::vec3(a, 0, 0));
  trans.push_back(glm::vec3(0, b, 0));
  trans.push_back(glm::vec3(0, 0, c));
  trans.push_back(glm::vec3(a, b, 0));
  trans.push_back(glm::vec3(a, 0, c));
  trans.push_back(glm::vec3(0, b, c));
  trans.push_back(glm::vec3(a, b, c));

  trans.push_back(glm::vec3(-a, 0, 0));
  trans.push_back(glm::vec3(0, -b, 0));
  trans.push_back(glm::vec3(0, 0, -c));
  trans.push_back(glm::vec3(-a, -b, 0));
  trans.push_back(glm::vec3(-a, 0, -c));
  trans.push_back(glm::vec3(0, -b, -c));
  trans.push_back(glm::vec3(-a, -b, -c));
 
  trans.push_back(glm::vec3(a, -b, 0));
  trans.push_back(glm::vec3(-a, b, 0));

  trans.push_back(glm::vec3(a, 0, -c));
  trans.push_back(glm::vec3(-a, 0, c));

  trans.push_back(glm::vec3(0, b, -c));
  trans.push_back(glm::vec3(0, -b, c));

  trans.push_back(glm::vec3(-a, b, c));
  trans.push_back(glm::vec3(a, -b, c));
  trans.push_back(glm::vec3(a, b, -c));

  trans.push_back(glm::vec3(-a, -b, c));
  trans.push_back(glm::vec3(-a, b, -c));
  trans.push_back(glm::vec3(a, -b, -c));

  /*
  generates bond from scratch, using the information of the atomic properties, that is,
  cylinder radius and bond radius
  */
  // float totclear = 0, totadd = 0, totabd = 0;
  vector<Bond> tempbonds = bonds_;
  this->bonds_.clear();
  
  //auto start0 = std::chrono::high_resolution_clock::now();

  //clear all atom bonds
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    this->atom(i).clearBonds();
  }
  //clear all atom bonds

  //take care of locked atoms first
  for(auto &bd:tempbonds){
    Atom& at1 = atom(bd.at1);
    Atom& at2 = atom(bd.at2);
    if(at1.locked && at2.locked){
      bonds_.push_back(bd);
    }
  }

  for(unsigned int i = 0; i < numAtoms(); i++){
    //printf("loop first level\n");
    // this->atom(i).clearBonds();  //clear the bonds of atom i
    for(unsigned int j = 0; j < i; j++){
      if(this->atom(i).locked && this->atom(j).locked){
        continue;
      }
      float length;
      bool abd = atBondDistance(i, j, length);
      if(abd){
        addBond(i, j, length);
      }
      else{
        abd = atBondDistancePBC(i, j, trans);
        // if(abd){
        //   addBond(i, j, length);
        // }
      }
    }
  }
  this->updateAtomsBondsIndices();
  //auto  stop0 = std::chrono::high_resolution_clock::now();
  //auto duration0 = std::chrono::duration_cast<std::chrono::microseconds> (stop0-start0);
	//std::cout<<"total bonds duration:"<<duration0.count()*0.001<<"milisec "<<std::endl;
}

int Trajectory::getAtomNeighborIndex(int i, Bond &bond){
  if(bond.at1 == i){
    return bond.at2;
  }
  else if(bond.at2 == i){
    return bond.at1;
  }
  else{
    return -1;
  }
}

Atom& Trajectory::getAtomNeighbor(const Atom& at, Bond &bond){
  Atom& at1 = this->atom(bond.at1);
  if(at == at1){
    return this->atom(bond.at2);
  }
  else {
    return at1;
  }
}

std::vector<unsigned int> Trajectory::getBondedIndices(unsigned int atIndex){
  std::vector<unsigned int> inds;
  if(atIndex > this->numAtoms()){
    printf("TRYING TO ACCESS ATOM INDEX GREATER THAN NUMMBER OF PROPERTIES\n");
    return inds;
  }
  Atom& at = this->atom(atIndex); //get reference to property atIndex
  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    Bond &b = bonds_[bi];
    int neighbor = this->getAtomNeighborIndex(atIndex, b);
    if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings
    
    inds.push_back((unsigned int)neighbor);
  }
  return inds;
}

std::vector<Atom> Trajectory::getBondedAtoms(unsigned int atIndex){
  std::vector<Atom> ats;
  if(atIndex > this->numAtoms()){
    printf("TRYING TO ACCESS ATOM INDEX GREATER THAN NUMMBER OF PROPERTIES\n");
    return ats;
  }

  Atom& at = this->atom(atIndex); //get reference to property atIndex
  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    Bond &b = bonds_[bi];
    int neighbor = this->getAtomNeighborIndex(atIndex, b);
    if(neighbor == -1) continue;   //if by mistake indices are wrong, check. No neighbor do nothings
    
    ats.push_back(this->atom(neighbor));
  }
  return ats;
}

std::vector<Atom> Trajectory::getBondedAtoms(const Atom& at){
  std::vector<Atom> ats;

  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    Bond &b = bonds_[bi];
    Atom& at2 = this->getAtomNeighbor(at, b);    
    ats.push_back(at2);
  }
  return ats;
}

std::vector<Bond> Trajectory::getBonds(const unsigned int i){
  std::vector<Bond> atombonds;
  for(auto &bi:atom(i).bondInds){
    if(bi == -1) break;
    atombonds.push_back(bonds_[bi]);
  }

  return atombonds;
}

std::vector<Bond> Trajectory::getBonds(const Atom& at){
  std::vector<Bond> atombonds;
  for(auto &bi:at.bondInds){
    if(bi == -1) break;
    atombonds.push_back(bonds_[bi]);
  }

  return atombonds;
}

std::vector<Bond> Trajectory::getBondFromMolecule(Molecule &mol){
  std::vector<Bond> bonds;
  // bonds.reserve(mol.atomInds.size()*2); //might not be enough
  // for(auto &ai:mol.atomInds){
  //   //the info for bonds is stored in atomic properties
  //   //now iterate over the bond indices stored in atomic property ai
  //   for(auto &bi:this->atom(ai).bondInds){
  //     if(bi == -1) break;
  //     Bond &b = bonds_[bi] ;  // make a reference to Bond object bi
  //     //check whether or not the bond bi has alreay be added to the vector bonds
  //     if(std::find(bonds.begin(), bonds.end(), b) == bonds.end()){
  //       //add it
  //       bonds.push_back(b);
  //     }
  //   }
  // }
  return bonds;
}

void Trajectory::hideSelectedCurrentFrame(){
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    Atom& at = this->atom(i);
    if(at.selected){
      // at.selected = false;
      at.hidden = true;
      at.selected = false;
    }
  }
}
void Trajectory::hideSelectedAllFrames(){
  for(auto &fr:frames_){
    for(unsigned int i = 0; i < fr.numAtoms(); i++){
      Atom& at = fr.atom(i);
      if(at.selected){
        // at.selected = false;
        at.hidden = true;
        at.selected = false;
      }
    }  
  }
}

std::string Trajectory::info(){
  string inf;
  inf += "Frames: " + std::to_string(numFrames()) + "\n";
  inf += "Atoms : " + std::to_string(numAtoms()) + "\n";
  inf += "Bonds :   " + std::to_string(numBonds()) + "\n";
  findMoleculesFromScratch();
  inf += "Molecules :   " + std::to_string(molecules_.size()) + "\n";
  return inf;
}

void Trajectory::initCellLines(){
  //here I have to read file
  for(unsigned int i=0; i<24; i++){
      cellLines_[i][0] = 0.0f;
      cellLines_[i][1] = 0.0f;
      cellLines_[i][2] = 0.0f;
  }
    
  cellLinesColors_[0] = {1.0f, 0.0f, 0.0f, 1.0f};
  cellLinesColors_[1] = {0.0f, 1.0f, 0.0f, 1.0f};
  cellLinesColors_[2] = {0.0f, 0.0f, 1.0f, 1.0f};
  for(unsigned int i = 3; i < 12; i++){
      cellLinesColors_[i] = {0.5f, 0.5f, 0.5f, 1.0f};
  }
  
  for(unsigned int i = 0; i < 12; i++){
      cellLinesWidths_[i] = 2;
  }

}

int Trajectory::isRepeated(glm::vec3 &coor, Frame &fr){
  //check whether the given coordinumberOfAtoms_es overlaps with an existing atom
  float umb;
  if(fr.cartesian) umb = 0.01f;
  else{
    umb = 0.0001f;
  }
  for(unsigned int i = 0; i < fr.numAtoms(); i++){
    if(checkEqualCoors(coor, fr.atom(i).coor, umb)){
      return i;		
    }
    if(!fr.cartesian && checkEqualReplicated(coor, fr.atom(i).coor, umb)){
      return i;	
    }
  }
  return -1;
}

int Trajectory::isRepeated(glm::vec3& coor){
  if(cfi_ == -1) return -1;  //return false for repetition
  //check whether the given coordinumberOfAtoms_es overlaps with an existing atom
  float umb;
  if(cartesian()) umb = 0.001f;
  else umb = 0.0001f;
  for(unsigned int i=0; i<this->numAtoms(); i++){
    if(checkEqualCoors(coor, this->atom(i).coor, umb)){
      return i;		
    }
    if(!cartesian() && checkEqualReplicated(coor, this->atom(i).coor, umb)){
      return i;	
    }
  }
  return -1;
}


void Trajectory::lockSelectedCurrentFrame(){
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    Atom& at = this->atom(i);
    if(at.selected){
      // at.selected = false;
      at.locked = true;
    }
  }
}
void Trajectory::lockSelectedAllFrames(){
  for(auto &fr:frames_){
    for(unsigned int i = 0; i < fr.numAtoms(); i++){
      Atom& at = fr.atom(i);
      if(at.selected){
        // at.selected = false;
        at.locked = true;
      }
    }  
  }
}

int Trajectory::nextFrame(int step){
  if(cfi_ + step >= (int)frames_.size()){
    return -1;  
  }
  cfi_ += step;
  this->updateBondLengths();
  return cfi_;
}

void Trajectory::removeDuplicated(){
  if(cfi_ == -1) return;
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    
  }
}

void Trajectory::rotateCoors(const float x, const float y, const float z, const float angle){
  if(cfi_ == -1) return;

  glm::vec3 rotVec(x, y, z);
  glm::mat3 rotMat;
  math::rotationMatrix3(rotVec, angle, rotMat); 

  for(unsigned int i = 0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    at.coor = rotMat * at.coor;
  }
}

// vector<unsigned int> Trajectory::recursiveBondSearch(int bonds, const vector<unsigned int> atomsList, bool periodic){
//     vector<unsigned int> partialList; 
//     vector<unsigned int> totalList;

//     if(bonds > 0){
//       bonds--;
//       for(unsigned int a = 0; a < atomsList.size(); a++){
//         vector<unsigned int> neighbors;

// 	    //iterate over the IDs of the neighbors of atom a
//         for(auto &nb:atoms_[atomsList[a]].neighborsIds()){
//           //check whether it has already been add to totallist
//           auto it=find(totalList.begin(), totalList.end(), nb);
//           if(it == totalList.end()){
// 		        // add the bonded atom id to the totallist
//             totalList.push_back(nb);
//           }
//                 //check whether it has already been add to totallist
		
// 		//append to neighbors list
// 		neighbors.push_back(nb);

// 		/* need to fix this later
// 		if(periodic){
//         	    //iterate over the periodic neighbors of atom a
//                     for(auto &nb:atoms_[atomsList[a]].periodicNeighbors){
//                         //check whether it has already been add to totallist
//                         auto it=find(totalList.begin(), totalList.end(), nb.id);
//                         if(it == totalList.end()){
//         		    // add the bonded atom id to the totallist
//                             totalList.push_back(nb.id);
//                         }
//                         //check whether it has already been add to totallist
		
// 		        //append to neighbors list
// 		        neighbors.push_back(nb.id);
// 		    }
// 		}
// 		*/
//             }
//             //partialList = recursiveBondSearch(bonds, atoms_[a].neighbors);
//             //partialList = recursiveBondSearch(bonds, atoms_[atomsList[a]].neighborsIds());
//             partialList = recursiveBondSearch(bonds, neighbors);
//             for(unsigned int p = 0;p < partialList.size();p++){
//                 //check whether it has already been add to totallist
//                 auto it=find(totalList.begin(), totalList.end(), partialList[p]);
//                 if(it == totalList.end()){
//                     totalList.push_back(partialList[p]);
//                 }
//             }
//         }
//     }

//     return totalList;
// }
    
//set the current frame
bool Trajectory::setCurrentFrameIndex(const int index){
  if(index >= (int)frames_.size() || index < 0){
    // the value if greater than or lower than allowed, do nothing
    return false;
  }

  cfi_ = index;
  this->updateBondLengths();
  return true;
}


void Trajectory::setCartFracMatrix(glm::mat3 &matrix){
  if(cfi_ == -1) return;
  glm::mat3 &cell = this->frames_[cfi_].cell_; 

  matrix = math::cartToFracMatrix(cell);
}

bool Trajectory::setDistance(Atom &at1, Atom &at2, float d){
  glm::vec3 vec = at1.coor - at2.coor;
  vec = glm::normalize(vec);
  at1.coor = at2.coor + vec * d;
  return true;
}

bool Trajectory::setDistance(unsigned int i, unsigned int j, float d){
  Atom& at1 = this->atom(i);
  Atom& at2 = this->atom(j);
  glm::vec3 vec = at1.coor - at2.coor;
  vec = glm::normalize(vec);
  at1.coor = at2.coor + vec * d;
  return true;
}

void Trajectory::setFracCartMatrix(glm::mat3  &matrix){
  if(cfi_ == -1) return;
  glm::mat3 &cell = this->frames_[cfi_].cell_; 

  matrix = math::fracToCartMatrix(cell);
}

void Trajectory::setSymbolSelected(const char* s){
  if(cfi_ == -1) return;
  for(unsigned int i=0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    if(at.selected){
      at.atomicNum = props::atomicNumber(s);
    }
  } 
}

void Trajectory::setAtomTypeSelected(const char* s){
  if(cfi_ == -1) return;
  for(unsigned int i=0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    if(at.selected){
      strncpy(at.atomType, s, 9);
    }
  } 
}

void Trajectory::setChargeSelected(const float charge){
  if(cfi_ == -1) return;
  for(unsigned int i=0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    if(at.selected){
      at.charge = charge;
    }
  } 
}

void Trajectory::setAtomSphereSelected(const float sphere){
  if(cfi_ == -1) return;
  for(unsigned int i=0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    if(at.selected){
      at.setSphere(sphere);
    }
  } 
}

void Trajectory::setUnitCell(float cellVecComps[9]){
    //set the lattice parametes a, b, c , alpha, beta, gamma from the vectors components
    //cellVectComps stands for unit cell vector components
    if(cfi_ == -1) return;
    for(int i=0; i<3; i++){
      for(int j = 0; j <3; j++){
        this->frames_[cfi_].cell_[i][j] = cellVecComps[3*i + j] * this->volScale;
      }
    }
    
    this->setUnitCellLines();
    return;
}

void Trajectory::setUnitCell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz ){
  if(cfi_ == -1) return;

  this->frames_[cfi_].cell_[0][0] = xx * this->volScale;
  this->frames_[cfi_].cell_[0][1] = xy * this->volScale;
  this->frames_[cfi_].cell_[0][2] = xz * this->volScale;
  this->frames_[cfi_].cell_[1][0] = yx * this->volScale;
  this->frames_[cfi_].cell_[1][1] = yy * this->volScale;
  this->frames_[cfi_].cell_[1][2] = yz * this->volScale;
  this->frames_[cfi_].cell_[2][0] = zx * this->volScale;
  this->frames_[cfi_].cell_[2][1] = zy * this->volScale;
  this->frames_[cfi_].cell_[2][2] = zz * this->volScale;
  
  this->setUnitCellLines();
  return;
}
    

void Trajectory::setUnitCellLines(){
  if (cfi_ == -1) return;
  if(!this->frames_[cfi_].cellIsSet){ 
    this->cellIsSet_ = false;
    return;
  }
  glm::mat3 &cell = this->frames_[cfi_].cell_; 
//    // de cero a a
  cellLines_[0] = {0.0,0.0,0.0};
  cellLines_[1] = {cell[0][0], cell[0][1], cell[0][2]};
  // de cero a b
  cellLines_[2] = {0.0,0.0,0.0};
  cellLines_[3] = {cell[1][0], cell[1][1], cell[1][2]};
  // de cero a c
  cellLines_[4] = {0.0,0.0,0.0};
  cellLines_[5] = {cell[2][0], cell[2][1], cell[2][2]};
  // de a a a+b
  cellLines_[6] = {cell[0][0], cell[0][1], cell[0][2]};
  cellLines_[7] = {cell[0][0]+cell[1][0], cell[0][1]+cell[1][1], cell[0][2]+cell[1][2]};
  // de  b a a+b
  cellLines_[8] = {cell[1][0], cell[1][1], cell[1][2]};
  cellLines_[9] = {cell[0][0]+cell[1][0], cell[0][1]+cell[1][1], cell[0][2]+cell[1][2]};
  // de a a a+c
  cellLines_[10] = {cell[0][0], cell[0][1], cell[0][2]};
  cellLines_[11] = {cell[0][0]+cell[2][0], cell[0][1]+cell[2][1], cell[0][2]+cell[2][2]};
  // de c a a+c
  cellLines_[12] = {cell[2][0], cell[2][1], cell[2][2]};
  cellLines_[13] = {cell[0][0]+cell[2][0], cell[0][1]+cell[2][1], cell[0][2]+cell[2][2]};
  // de  b a b+c
  cellLines_[14] = {cell[1][0], cell[1][1], cell[1][2]};
  cellLines_[15] = {cell[1][0]+cell[2][0], cell[1][1]+cell[2][1], cell[1][2]+cell[2][2]};
  // de  c a b+c
  cellLines_[16] = {cell[2][0], cell[2][1], cell[2][2]};
  cellLines_[17] = {cell[1][0]+cell[2][0], cell[1][1]+cell[2][1], cell[1][2]+cell[2][2]};
  //de b+c a a+b+c
  cellLines_[18] = {cell[1][0]+cell[2][0], cell[1][1]+cell[2][1], cell[1][2]+cell[2][2]};
  cellLines_[19] = {cell[0][0]+cell[1][0]+cell[2][0], cell[0][1]+cell[1][1]+cell[2][1],cell[0][2]+cell[1][2]+cell[2][2]};
  //de a+c a a+b+c
  cellLines_[20] = {cell[0][0]+cell[2][0], cell[0][1]+cell[2][1], cell[0][2]+cell[2][2]};
  cellLines_[21] = {cell[0][0]+cell[1][0]+cell[2][0], cell[0][1]+cell[1][1]+cell[2][1],cell[0][2]+cell[1][2]+cell[2][2]};
  //de a+b a a+b+c
  cellLines_[22] = {cell[0][0]+cell[1][0], cell[0][1]+cell[1][1], cell[0][2]+cell[1][2]};
  cellLines_[23] = {cell[0][0]+cell[1][0]+cell[2][0], cell[0][1]+cell[1][1]+cell[2][1],cell[0][2]+cell[1][2]+cell[2][2]};

  this->cellIsSet_ = true;

}

void Trajectory::setUnitCellVectorsComponents(props::LatticeParams &lp){
  float PI = 3.14159265;
  if(cfi_ == -1) return;
  glm::mat3 cell;
  
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

  this->frames_[cfi_].setUnitCell(cell);

  setUnitCellLines();
  this->cellIsSet_ = true;
}

void Trajectory::toCentroid(){
  this->calcCentroid();
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    at.coor -= centroid_;
  }
}

void Trajectory::translateCoors(const glm::vec3 &transVec){
  if(cfi_ == -1) return;

  for(unsigned int i = 0; i < this->numAtoms(); i++){
    Atom &at = this->atom(i);
    at.coor += transVec;
  }
}

void Trajectory::unlockObjectsCurrentFrame(){
  for(unsigned int i = 0; i < this->numAtoms(); i++){
    this->atom(i).locked = false;
  }
}

void Trajectory::unlockObjectsAllFrames(){
  for(auto &fr:frames_){
    for(unsigned int i = 0; i < fr.numAtoms(); i++){
      fr.atom(i).locked = false;
    }  
  }
}

void Trajectory::updateBondLengths(){
  for(unsigned int i=0; i < this->numBonds(); i++){
    Bond &b = this->bond(i);
    b.length = this->distance(b.at1, b.at2);
  }
}

void Trajectory::updateBondRadius(){
  for(unsigned int i=0; i < this->numBonds(); i++){
    Bond &b = this->bond(i);
    if(atom(b.at1).cylRadius() >= atom(b.at2).cylRadius()){
      b.cylRadius = atom(b.at2).cylRadius();
    }
    else{
      b.cylRadius = atom(b.at1).cylRadius();
    }
  }
}

void Trajectory::updateAtomsBondsIndices(){
  for(unsigned int i = 0; i < bonds_.size(); i++){
    Bond &bd = this->bonds_[i];
    this->atom(bd.at1).addBond(i);
    this->atom(bd.at2).addBond(i);
  }
}

/*
void Trajectory::toCentroid(int num){
    if(num == -1){
        for(auto &st:structures_){
	    st.toCentroid();
	}
    }else{
        structures_[num].toCentroid();
    }
}

void Trajectory::recalcCentroid(int num){
    if(num == -1){
        for(auto &st:structures_){
	    st.recalcCentroid();
	}
    }else{
        structures_[num].recalcCentroid();
    }
}

void Trajectory::generateBonds(int index){
    if(index == -1){
        for(auto &st:structures_){
	    st.generateBonds();
	}
    }else{
        structures_[index].generateBonds();
    }
}

void Trajectory::generateBondsFromFirstFrame(){
  structures_[0].generateBonds();
  Structure &st0 = structures_[0];

  for(unsigned int i=1; i < structures_.size(); i++){
    for(unsigned int j = 0; j < st0.totNumberOfBonds(); j++){
      if(st0.bond(j).deleted){
        continue;
      }
     
      //add the atom
      //this has to be better programed, need to check first if atoms 1 and 2 exists
      //this catch is not gonna work
      try{
        structures_[i].addBond(st0.bond(j).at1, st0.bond(j).at2, st0.bond(j).cylRadius);
      }catch(...){
        printf("problem generating bonds\n");
      }
    }    
  }
}

unsigned int Trajectory::setNextFrame(){

    if(numStructures_ == 0){
        return -1;
    }
    if(currentFrameIndex_ == (int)numStructures_ - 1) {
        currentFrameIndex_ = 0;
	return 0;
    }

    currentFrameIndex_ ++;
    return currentFrameIndex_;
}


void Trajectory::selectAllAtomsCurrentFrame(){
    structures_[currentFrameIndex_].selectAll();
}

void Trajectory::selectAllAtoms(int frame){
    if(frame == -1){  //only for selected trajectories
        for(auto &st:structures_){
	    if(st.deleted) continue;
	    if(!st.selected) continue;
	    st.selectAll();
	}
    }
    else{  //for specific index
	for(auto &st:structures_){
	    if(st.deleted) continue;
	    if(frame == st.index){
	        st.selectAll();
	    }
	}
    }
}


void Trajectory::selectAllFrames(){
    for(auto &st:structures_){
        st.selected = true;
    }
}

void Trajectory::selectByElementsCurrentFrame(vector<string> symbols){
    structures_[currentFrameIndex_].selectByElements(symbols);
}

void Trajectory::selectByAtomTypeCurrentFrame(vector<string> symbols){
    structures_[currentFrameIndex_].selectByAtomType(symbols);
}

void Trajectory::selectByResidueCurrentFrame(vector<string> symbols){
    structures_[currentFrameIndex_].selectByResidue(symbols);
}

void Trajectory::selectByNumberNeighborsCurrentFrame(std::vector<unsigned short> numNeighbors, unsigned short n, bool pbc){
    structures_[currentFrameIndex_].selectByNumberNeighbors(numNeighbors, n, pbc);
    return;
}

void Trajectory::selectByNBondsDistanceCurrentFrame(int bonds){
    structures_[currentFrameIndex_].selectByNBondsDistance(bonds);
}

void Trajectory::info(){
    structures_[currentFrameIndex_].info();
}

void Trajectory::invertSelectionCurrentFrame(){
    structures_[currentFrameIndex_].invertSelection();
}

void Trajectory::invertSelectionElementsSelectedCurrentFrame(){
    structures_[currentFrameIndex_].invertSelectionElementsSelected();
}

void Trajectory::unselectAllCurrentFrame(){
    structures_[currentFrameIndex_].unselectAll();
}

void Trajectory::hideCurrentFrame(){
    structures_[currentFrameIndex_].hide();
}

void Trajectory::deleteSelectedCurrentFrame(){
    //deletes the selected atoms, bonds and other objects of the current frame
    structures_[currentFrameIndex_].deleteSelected();

    // registers the changes
    TrajChange tc;
    tc.action = ACTION_DEL_OBJECTS;
    tc.structIds.push_back(currentFrameIndex_);
    beforeChanges.push_back(tc);
}

void Trajectory::changeElementCurrentFrame(std::string stringElement){
    structures_[currentFrameIndex_].changeElementSelected(stringElement);

    //I have to track this and put in beforeCHanges

    return;    
}

void Trajectory::setVolScaleCurrentFrame(float scale){
    structures_[currentFrameIndex_].setVolScale(scale);
}

void Trajectory::setUnitCellCurrentFrame(const array<float, 9> &cellParams){
    structures_[currentFrameIndex_].setUnitCell(cellParams);
}

void Trajectory::setUnitCellCurrentFrame(const array<float, 6> &cellParams){
    structures_[currentFrameIndex_].setUnitCell(cellParams);
}

void Trajectory::setUnitCellCurrentFrame(const array<array<float, 3>, 3> &cellParams){
    structures_[currentFrameIndex_].setUnitCell(cellParams);
}

void Trajectory::changeAtomColorCurrentFrame(const float r, 
		                             const float g, 
					     const float b, 
					     const float a){


    structures_[currentFrameIndex_].changeAtomColorSelected(r, g, b, a);
}


void Trajectory::undo(){
    if(beforeChanges.size() == 0) return;

    TrajChange tc = beforeChanges.back();
    afterChanges.push_back(tc);
    
    if(tc.action == ACTION_DEL_OBJECTS){
	for(auto &id:tc.structIds){
	    structures_[id].undo();
	}
    }
    else if(tc.action == ACTION_REPLICATE){
	for(auto &id:tc.structIds){
	    structures_[id].undo();
	}
    }

    //delete change from before changes
    beforeChanges.pop_back();
}

void normalModeAnalysisCurrentFrame(int itemnum, float temp, float pressure, float vol){
    

}

void Trajectory::replicateCurrentFrame(const int xi, const int xf, const int yi, const int yf, const int zi, const int zf){
  structures_[currentFrameIndex_].replicateStructure(xi, xf, yi, yf, zi, zf);
    
  // registers the changes
  TrajChange tc;
  tc.action = ACTION_REPLICATE;
  tc.structIds.push_back(currentFrameIndex_);
  beforeChanges.push_back(tc);

}
*/
