#ifndef MATH_H
#define MATH_H

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <array>
#include <vector>
#include <math.h>
#include "../frame.h"


namespace math{

  extern float PI;
  
  glm::mat3 fracToCartMatrix(props::LatticeParams &lp);
  glm::mat3 fracToCartMatrix(glm::mat3 &cell);
  glm::mat3 cartToFracMatrix(props::LatticeParams &lp);
  glm::mat3 cartToFracMatrix(glm::mat3 &cell);
  void bringInsideCell(glm::vec3 &coor);

  void rotationMatrix3(glm::vec3 &rotVec, const float angle, glm::mat3 &rotMat);
  void matMult4x1(const std::array<std::array<float, 4>, 3> &mat, std::array<float, 3> &vec);
  void matMult4x1(const std::array<std::array<float, 4>, 3> &mat, glm::vec3 &vec);
  void rotateVector(float refcoors[3], float rotmat[3][3]);
  void rotateMatrix(float mat[3][3], float rotmat[3][3]);
  void rotateVector(glm::vec3 &refcoors, float rotmat[3][3]);
  void rotateVector(float ang, std::array<float, 3> &rotvec, std::array<float, 3> &refcoors);
  void VectorNormalisation(float v[3]);
  float calcAngle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);
  float calcAngleRadians(const std::array<float,3> &p1, const std::array<float,3> &p2, const std::array<float,3> &p3);
  float calcCosAngle(std::array<float,3>v1, std::array<float,3>v2);
  float calcDihedral(std::array<float,3>p1, std::array<float,3>p2, std::array<float,3>p3, std::array<float,3>p4);
  void CrossProduct(float w[3], float v[3], float cross[3]);
  void CrossProduct(std::array<float, 3> &p1, std::array<float, 3> &p2, std::array<float, 3> &p3, std::array<float, 3> &cross);
  float DotProduct(float v[3], float w[3]);
  float calcDistance(std::array<float,3> &point1, std::array<float,3> &point2);
  float calcSquaredDistance(std::array<float,3> &point1, std::array<float,3> &point2);
  float calcDistance(float x1, float y1, float z1, float x2, float y2, float z2);
  float calcAverage(std::vector<float> list);
  float distance(const Atom& at1, const Atom& at2);
  //float average(vector<int> list);
  bool matrixIsIdentity(std::array<std::array<float, 4>, 3> &rmat);
  bool setDistance(Atom &at1, Atom &at2, float d);
};
//inline functions

#endif //MATH_H
