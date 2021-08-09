#ifndef BOND_H
#define BOND_H
//#include<cylinder.h>

class Bond  //: public Cylinder
{
public:
    Bond(){}
    //unsigned int index = -1;
    int oId = -1;
    unsigned int id;
    unsigned int at1, at2;
    float cylRadius = 0.1;
    bool deleted = false;
    bool selected = false;
    bool pbc = false;
    bool hidden = false;
    bool locked = false;
    bool drawWithLine = false;
    unsigned int order = 1;
    //translation along the cell vectors retore the neighborhood
    std::array<int, 3> pdir;
};


#endif // ATOM_H
