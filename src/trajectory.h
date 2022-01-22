#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "frame.h"
///#include "properties.h"   dara lugar a multiple definitions porque no use extern en el header

// Tracks the changes of on trajectories objects or
// on the structures
/*
struct TrajChange{
    // action indicates whether the performed change 
    unsigned int action;

    // structsIds indicates which structures where modified
    std::vector<unsigned int> structIds;
};
*/

struct RingSearch{

  RingSearch(int a, int b, short d){
    atIndex = a; atIndexTarget = b; ringSize = d;
  }
  RingSearch(const RingSearch& rs){
    atIndex = rs.atIndex;
    atIndexTarget = rs.atIndexTarget;
    prevAtIndex = rs.prevAtIndex;
    ringSize = rs.ringSize;
    level = rs.level;
  }

  int atIndex;
  int atIndexTarget;
  int prevAtIndex = -1; //garantees that the loop will not get stuck in the first level
  short ringSize;
  short level = 1; //always starts with 1
};

struct Freq{
  float vib = 0;
  std::vector<glm::vec3> disp;

  bool operator==(const Freq& other) {return false;}
  bool operator!=(const Freq& other) {return true;}
};

class Trajectory
{
  public:
    Trajectory();

    //variables
    bool selected = false;

    float volScale = 1.0f;

    int tableNumber = 1;

    std::string spaceGroup;
    std::string name;
    std::string dir;
    std::vector<Freq> freqs;
    
    //functions
    void addBond(const unsigned int i, const unsigned int j, const float &bdist, bool updateAtoms = false);
    void addBond(const unsigned int i, const unsigned int j, const float &bdist, glm::vec3 &offset);
    void addBond(const unsigned int i, const unsigned int j, const float &bdist, const float &cylRadius, bool updateAtoms = false);
    void addBondedAtomsGeometrically(unsigned int index, short n, short atomicNumber, float blength);
    void addHydrogen(const unsigned int index, const float blength);
    void addFrame();
    unsigned int addFrame(Frame &fr);

    bool atBondDistance(const unsigned int i, const unsigned int j, float &bdist);
    bool atBondDistancePBC(const unsigned int i, const unsigned int j, std::vector<glm::vec3>& trans);

    Atom& atom(unsigned int i){return frames_[cfi_].atom(i);}
    
    //return angle bewteen 3 atoms
    float angle(const unsigned int i, const unsigned int j, const unsigned int k);

    //get bond
    Bond& bond(unsigned int i){return bonds_[i];}   //check the size

    void bringInsideCell();
    void bringInsideCellAllFrames();

    glm::vec3 calcCentroid();

    bool cartesian(){ if(cfi_ == -1) return true;  return currentFrame().cartesian;}
    void cartesianToFractional(); 
    void cartToFracAllFrames();

    glm::vec3 centroid(){return centroid_;}

    glm::mat3& cell(){return frames_[cfi_].cell_;}
    glm::vec3& cellLine(unsigned int i){return cellLines_[i];}
    glm::vec4& cellLineColor(unsigned int i){return cellLinesColors_[i];}
    int cellLineWidth(unsigned int i){return cellLinesWidths_[i];}
    bool cellIsSet(){return cellIsSet_;}
    
    void changeSphereSelected(const float s);
    bool checkEqualCoors(glm::vec3 &a1, glm::vec3 &a2, float umb);
    bool checkEqualReplicated(glm::vec3 &a1, glm::vec3 &a2, float umb);
    // bool checkRingClosed(int level, int atIndex, int atIndexTarget, int excPrevIndex, std::vector<unsigned int>&inds);
    short checkRingClosed(const RingSearch rc, std::vector<unsigned int>&inds, std::vector<unsigned int> checked);
  
    Bond createBond(const unsigned int i, const unsigned int j, const float &bdist);
    //return reference to frame object 
    Frame& currentFrame(){return frames_[cfi_];} //CAUTION COULD HAVE SEGFAULT
    
    //return reference to frame object
    int currentFrameIndex(){return cfi_;}

    //calculate distance between two atom
    float distance(const unsigned int i, const unsigned int j);

    float energy();

    void findMoleculesFromScratch();
    void findRing(short ringSize, unsigned int atIndex, std::vector<unsigned int>&inds);
    void findMoleculeRecursively(unsigned int atIndex, std::vector<unsigned int>&inds, bool lock = false);
    
    void fractionalToCartesian();
    void fractToCartAllFrames();

    //return reference to frame object
    Frame& frame(const unsigned int i);
    std::vector<Frame> frames(){return frames_;}

    //generate bonds for the current frame
    void generateBonds();
    void generateBondsPBC();
    void generateNewBonds(std::vector<unsigned int> &inds);

    std::vector<Bond> getBondFromMolecule(Molecule &mol);

    int getAtomNeighborIndex(int i, Bond &bond);
    Atom& getAtomNeighbor(const Atom& at, Bond &bond);
    std::vector<unsigned int> getBondedIndices(unsigned int atIndex);
    std::vector<Atom> getBondedAtoms(unsigned int atIndex);
    std::vector<Atom> getBondedAtoms(const Atom& at);
    std::vector<Bond> getBonds(const Atom& at);
    std::vector<Bond> getBonds(const unsigned int i);

    void hideSelectedCurrentFrame();
    void hideSelectedAllFrames();
    std::string info();
    int isRepeated(glm::vec3 &coor, Frame &fr);
    int isRepeated(glm::vec3 &coor);
    void lockSelectedCurrentFrame();
    void lockSelectedAllFrames();
    
    std::vector<Molecule> molecules(){return molecules_;}

    //advanece one frame forward
    int nextFrame(int step = 1);

    unsigned int numAtoms(){return frames_[cfi_].numAtoms();}
    unsigned int numBonds(){return bonds_.size();}
    unsigned int numFrames(){return frames_.size();}

    std::vector<std::vector<unsigned int>> orderByMolecule();

    void removeDuplicated();

    void rotateCoors(const float x, const float y, const float z, const float angle);

    bool setCurrentFrameIndex(const int index);
    void setCartesian(bool iscart){currentFrame().cartesian = iscart;}
    bool setDistance(Atom &at1, Atom &at2, float d);
    bool setDistance(unsigned int i, unsigned int j, float d);
    
    void setFracCartMatrix(glm::mat3 &matrix);
    void setCartFracMatrix(glm::mat3 &matrix);

    void setSymbolSelected(const char* s);
    void setAtomTypeSelected(const char* s);
    void setChargeSelected(const float charge);  
    void setAtomSphereSelected(const float sphere);

    void setUnitCell(float cellVecComps[9]);
    void setUnitCell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz );
    void setUnitCellLines();
    void setUnitCellVectorsComponents(props::LatticeParams &lp);

    void toCentroid();
    void translateCoors(const glm::vec3 &transVec);

    float totalCharge(){return frames_[cfi_].totalCharge();}

    void unlockObjectsCurrentFrame();
    void unlockObjectsAllFrames();

    void updateBondLengths();
    void updateBondRadius();
    void updateAtomsBondsIndices();


  private:
    //variables
    bool cellIsSet_ = false;

    std::vector<Bond> bonds_; 
    std::vector <Frame> frames_;
    std::vector<Molecule> molecules_;
    int cfi_ = -1;  //current Frame Index

    std::array<glm::vec3, 24> cellLines_;
    std::array<glm::vec4, 12> cellLinesColors_;
    std::array<int, 12> cellLinesWidths_;
    
    glm::vec3 centroid_ = glm::vec3(0.0,0.0,0.0);

    //std::vector<TrajChange> beforeChanges, afterChanges;

    //action constants
    unsigned int ACTION_DEL_FRAMES = 0;
    unsigned int ACTION_DEL_OBJECTS = 1; //delete atoms or bonds of structures
    unsigned int ACTION_REPLICATE = 2; //delete atoms or bonds of structures

    // funcs funcs funcs
    void initCellLines();

};


#endif
