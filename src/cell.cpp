#include "cell.h"
#include "math/math.h"


Cell::Cell(glm::mat3 &cell){
  this->vecComps_ = cell;
}

Cell::Cell(const float a, const float b, const float c, const float alpha, const float beta, const float gamma){
  float PI = 3.14159265;

  vecComps_[0][0] = a;
  vecComps_[0][1] = 0.0;
  vecComps_[0][2] = 0.0;
  float ang_rad = PI/180*gamma;  //gamma to radians
  float cosang = cos(ang_rad);
  float sinang = sin(ang_rad);
  vecComps_[1][0] = b * cosang;
  vecComps_[1][1] = b * sinang;
  vecComps_[1][2] = 0.0;

  //z vector components
  ang_rad = PI/180.0*beta;  //
  cosang = cos(ang_rad);
  vecComps_[2][0] = (cosang*a*c)/vecComps_[0][0];
  ang_rad = PI/180.0*alpha;  //
  cosang = cos(ang_rad);
  vecComps_[2][1] = (cosang*b*c - vecComps_[1][0]*vecComps_[2][0])/vecComps_[1][1];
  vecComps_[2][2] = sqrt(c * c - vecComps_[2][0]*vecComps_[2][0] - vecComps_[2][1]*vecComps_[2][1]);
}

Cell::Cell(glm::vec3 &v1, glm::vec3 &v2, glm::vec3 &v3){
  this->vecComps_[0] = v1;
  this->vecComps_[1] = v2;
  this->vecComps_[3] = v3;
}

Cell::Cell(float xx, float xy, float xz, float yx, float yy, float yz, float zx, float zy, float zz ){
  this->vecComps_[0][0] = xx;
  this->vecComps_[0][1] = xy;
  this->vecComps_[0][2] = xz;
  this->vecComps_[1][0] = yx;
  this->vecComps_[1][1] = yy;
  this->vecComps_[1][2] = yz;
  this->vecComps_[2][0] = zx;
  this->vecComps_[2][1] = zy;
  this->vecComps_[2][2] = zz;
  return;
}

float Cell::a(){
  return sqrt(vecComps_[0][0] * vecComps_[0][0] + 
            vecComps_[0][1] * vecComps_[0][1] +
              vecComps_[0][2] * vecComps_[0][2]);
}

float Cell::b(){
    return sqrt(vecComps_[1][0] * vecComps_[1][0] + 
            vecComps_[1][1] * vecComps_[1][1] +
              vecComps_[1][2] * vecComps_[1][2]);
}

float Cell::c(){
  return sqrt(vecComps_[2][0] * vecComps_[2][0] + 
            vecComps_[2][1] * vecComps_[2][1] +
              vecComps_[2][2] * vecComps_[2][2]);

}

float Cell::alpha(){
  glm::vec3 zero(0,0,0);
  return math::calcAngle(vecComps_[1], zero, vecComps_[2]);
}


float Cell::beta(){
  glm::vec3 zero(0,0,0);
  return math::calcAngle(vecComps_[0], zero, vecComps_[2]);
}


float Cell::gamma(){
  glm::vec3 zero(0,0,0);
  return math::calcAngle(vecComps_[0], zero, vecComps_[1]);
}

void Cell::scale(const float volScale){
  this->vecComps_[0][0] *= volScale;
  this->vecComps_[0][1] *= volScale;
  this->vecComps_[0][2] *= volScale;
  this->vecComps_[1][0] *= volScale;
  this->vecComps_[1][1] *= volScale;
  this->vecComps_[1][2] *= volScale;
  this->vecComps_[2][0] *= volScale;
  this->vecComps_[2][1] *= volScale;
  this->vecComps_[2][2] *= volScale;

}