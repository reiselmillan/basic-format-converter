#ifndef FRAME_H
#define FRAME_H

#include "atom.h"
#include "cell.h"

struct AtomStatis{
  Atom at;
  unsigned int num = 0;
};

struct Bond
{
  int at1 = -1;
  int at2 = -1;
  float cylRadius = 0.2;
  float length = -1;
  glm::vec3 offset = glm::vec3(0.0f, 0.0f, 0.0f);


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
  std::vector<unsigned int> inds;
  
  //just for python wrapper
  const char * getName(){return name;}
  void setName(const char* n){
    strcpy(name, n);
  }

  bool operator==(const Molecule& other) {return false;}
  bool operator!=(const Molecule& other) {return true;}

};

class Frame
{
  public:
    Frame(){};
    virtual ~Frame(){};

    // class EvenIterator {
    //   public:
    //     EvenIterator(std::vector<Atom>::iterator it, std::vector<Atom>::iterator end) : it(it), end(end) {
    //       while (true) {
    //         // if (isEven(*it)) {
    //         //   break;
    //         // } else if (it == end) {
    //         //   break;
    //         // }

    //         if(it == end){
    //           break;
    //         }
    //         it++;
    //       }
    //     }

    //     bool operator != (const EvenIterator& evenIt) {
    //       return evenIt.it != this->it;
    //     }

    //     Atom operator * () {
    //       return *it;
    //     }

    //     EvenIterator operator ++ () {
    //       while (true) {
    //         it++;
    //         if (isEven(*it)) {
    //           return EvenIterator(it, end);
    //         } else if (it == end) {
    //           return EvenIterator(it, end);
    //         }
    //       }
    //     }
    //   private:
    //     std::vector<Atom>::iterator it;
    //     std::vector<Atom>::iterator end;
    //   };

    //   static bool isEven(Atom number) {
    //     return true;
    //     // return number % 2 == 0;
    //   }

    //   // void add(int number) {
    //   //   atoms_.push_back(number);
    //   // }

    //   EvenIterator begin() {
    //     return EvenIterator(atoms_.begin(), atoms_.end());
    //   }

    //   EvenIterator end() {
    //     return EvenIterator(atoms_.end(), atoms_.end());
    //   }


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
    glm::mat3 cell_;

    //functions
    Atom& addAtom(const char* symbol, const float x, const float y, const float z); //It might return a reference to the Atom instead of the index
    Atom& addAtom(const char* symbol, glm::vec3 &coor);
    Atom& addAtom(short atomicNumber, glm::vec3 &coor);
    Atom& addAtom(short atomicNumber, const float x, const float y, const float z);
    Atom& addAtom(Atom& atom);
    unsigned int numAtoms(){return atoms_.size();}

    Atom& atom(unsigned int i);

    Cell cell(){return Cell(cell_);}

    std::vector<Atom> deleteSelected();

    std::vector<std::string> orderByAtomType();
    std::vector<AtomStatis> orderByAtomicNum();
    std::vector<AtomStatis> orderByAtomicNum(const std::vector<short> &anvec);

    void removeAtom(unsigned int index);
    void removeDuplicates();
    bool replicateAtoms(int x = 1, int y = 1, int z = 1);
    void reserve(int n);

    void setUnitCell(glm::mat3 &cell);
    void setUnitCell(props::LatticeParams &lp);
    void setUnitCell(Cell &cell);

    float totalCharge();
    std::vector<PeriodicAtom> periodicAtoms_;

  private:

    std::vector<Atom> atoms_;

};

#endif // FRAME_H
