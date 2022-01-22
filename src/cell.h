#ifndef CELL_H
#define CELL_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Cell
{
  public:
    Cell(){};
    Cell(glm::mat3 &cell);
    Cell(const float a, const float b, const float c, const float alpha = 90.0f, const float beta = 90.0f, const float gamma = 90.0f);
    Cell(glm::vec3 &v1, glm::vec3 &v2, glm::vec3 &v3);
    Cell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz );
    glm::mat3 vectorsComponets(){ return vecComps_; }

    //glm::mat3 vecComps_; //vectors components

    float a();
    float b();
    float c();
    float alpha();
    float beta();
    float gamma();

    void scale(const float volScale);

  protected:
    glm::mat3 vecComps_; //vectors components
};

#endif