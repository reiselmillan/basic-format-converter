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
    
    //functions
    void addBond(const unsigned int i, const unsigned int j, const float &bdist, bool updateAtoms = false);
    void addBond(const unsigned int i, const unsigned int j, const float &bdist, const float &cylRadius, bool updateAtoms = false);
    void addBondedAtomsGeometrically(unsigned int index, short n, short atomicNumber, float blength);
    void addHydrogen(const unsigned int index, const float blength);
    void addFrame();
    unsigned int addFrame(Frame &fr);

    void addMoleculeFromSelected(const std::string name);

    bool atBondDistance(const unsigned int i, const unsigned int j, float &bdist);

    Atom& atom(unsigned int i){return frames_[cfi_].atom(i);}
    
    //return angle bewteen 3 atoms
    float angle(const unsigned int i, const unsigned int j, const unsigned int k);

    //get bond
    Bond& bond(unsigned int i){return bonds_[i];}   //check the size

    void bringInsideCell();

    glm::vec3 calcCentroid();

    bool cartesian(){ if(cfi_ == -1) return true;  return currentFrame().cartesian;}
    void cartesianToFractional(); 
    void cartToFracAllFrames();

    glm::vec3 centroid(){return centroid_;}

    glm::mat3& cell(){return frames_[cfi_].cell;}
    glm::vec3& cellLine(unsigned int i){return cellLines_[i];}
    glm::vec4& cellLineColor(unsigned int i){return cellLinesColors_[i];}
    int cellLineWidth(unsigned int i){return cellLinesWidths_[i];}
    bool cellIsSet(){return cellIsSet_;}
    
    void changeSphereSelected(const float s);
    bool checkEqualCoors(glm::vec3 &a1, glm::vec3 &a2, float umb);
    bool checkEqualReplicated(glm::vec3 &a1, glm::vec3 &a2, float umb);
    // bool checkRingClosed(int level, int atIndex, int atIndexTarget, int excPrevIndex, std::vector<unsigned int>&inds);
    bool checkRingClosed(int level, int atIndex, int atIndexTarget, std::vector<unsigned int>&checked, std::vector<unsigned int>&inds);
  
    Bond createBond(const unsigned int i, const unsigned int j, const float &bdist);
    //return reference to frame object
    Frame& currentFrame(){return frames_[cfi_];} //CAUTION COULD HAVE SEGFAULT
    
    //return reference to frame object
    int currentFrameIndex(){return cfi_;}

    //calculate distance between two atom
    float distance(const unsigned int i, const unsigned int j);

    void findRing(short ringSize, unsigned int atIndex, std::vector<unsigned int>&inds);
    void findMoleculeRecursively(unsigned int atIndex, std::vector<unsigned int>&inds);
    
    void fractionalToCartesian();
    void fractToCartAllFrames();

    //return reference to frame object
    Frame& frame(const unsigned int i);

    //generate bonds for the current frame
    void generateBonds();
    void generateNewBonds(std::vector<unsigned int> &inds);

    std::vector<Bond> getBondFromMolecule(Molecule &mol);

    int getAtomNeighborIndex(int i, Bond &bond);
    std::vector<unsigned int> getAtomNeighborsIndices(unsigned int atIndex);

    void hideSelectedCurrentFrame();
    void hideSelectedAllFrames();
    int isRepeated(glm::vec3 &coor, Frame &fr);
    int isRepeated(glm::vec3 &coor);
    void lockSelectedCurrentFrame();
    void lockSelectedAllFrames();
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
