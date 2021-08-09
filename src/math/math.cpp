#include "math.h"
#include <chrono>
#include <iostream>

using std::array;
using std::vector;
using std::pow;
using math::PI;


float math::PI = 3.1415926;

glm::mat3 math::fracToCartMatrix(props::LatticeParams &lp){
  glm::mat3 matrix;
  float cos_alpha = cos(lp.alpha * PI / 180.0f);
  float cos_beta = cos(lp.beta * PI / 180.0f);
  float cos_gamma = cos(lp.gamma * PI / 180.0f);
  float sin_gamma = sqrt(1-pow(cos_gamma,2));
  float v = sqrt( 1 - pow(cos_alpha,2) - pow(cos_beta,2) -  pow(cos_gamma,2) + 2 * cos_alpha * cos_beta * cos_gamma );

  matrix[0] = {lp.a , lp.b * cos_gamma, lp.c * cos_beta };
  matrix[1] = {0, lp.b * sin_gamma, lp.c * (cos_alpha - cos_beta * cos_gamma)/(sin_gamma)};
  matrix[2] = {0,0, lp.c * v/sin_gamma};

  return matrix;

}

glm::mat3 math::fracToCartMatrix(glm::mat3 &cell){
  glm::mat3 matrix;

  float a,b,c, v, sin_gamma;
  float cos_alpha, cos_beta, cos_gamma;

  a = sqrt(pow(cell[0][0],2) + pow(cell[0][1],2) + pow(cell[0][2],2));
  b = sqrt(pow(cell[1][0],2) + pow(cell[1][1],2) + pow(cell[1][2],2));
  c = sqrt(pow(cell[2][0],2) + pow(cell[2][1],2) + pow(cell[2][2],2));

  cos_alpha = (cell[1][0]*cell[2][0] + cell[1][1]*cell[2][1]+cell[1][2]*cell[2][2])/(b*c);
  cos_beta = (cell[0][0]*cell[2][0] + cell[0][1]*cell[2][1]+cell[0][2]*cell[2][2])/(a*c);
  cos_gamma = (cell[0][0]*cell[1][0] + cell[0][1]*cell[1][1]+cell[0][2]*cell[1][2])/(a*b);
  sin_gamma = sqrt(1-pow(cos_gamma,2));
  v = sqrt( 1 - pow(cos_alpha,2) - pow(cos_beta,2) -  pow(cos_gamma,2) + 2 * cos_alpha * cos_beta * cos_gamma );

  matrix[0] = {a , b * cos_gamma, c * cos_beta };
  matrix[1] = {0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma)/(sin_gamma)};
  matrix[2] = {0,0, c * v/sin_gamma};

  return matrix;

}

glm::mat3 math::cartToFracMatrix(glm::mat3 &cell){
  glm::mat3 matrix;
  
  float a,b,c, v, sin_gamma;
  float cos_alpha, cos_beta, cos_gamma;

  a = sqrt(pow(cell[0][0],2) + pow(cell[0][1],2) + pow(cell[0][2],2));
  b = sqrt(pow(cell[1][0],2) + pow(cell[1][1],2) + pow(cell[1][2],2));
  c = sqrt(pow(cell[2][0],2) + pow(cell[2][1],2) + pow(cell[2][2],2));

  cos_alpha = (cell[1][0]*cell[2][0] + cell[1][1]*cell[2][1]+cell[1][2]*cell[2][2])/(b*c);
  cos_beta = (cell[0][0]*cell[2][0] + cell[0][1]*cell[2][1]+cell[0][2]*cell[2][2])/(a*c);
  cos_gamma = (cell[0][0]*cell[1][0] + cell[0][1]*cell[1][1]+cell[0][2]*cell[1][2])/(a*b);
  sin_gamma = sqrt(1-pow(cos_gamma,2));
  v = sqrt( 1 - pow(cos_alpha,2) - pow(cos_beta,2) -  pow(cos_gamma,2) + 2 * cos_alpha * cos_beta * cos_gamma );

  matrix[0] = {1/a, -cos_gamma/(a*sin_gamma), (cos_alpha * cos_gamma - cos_beta)/(a*v*sin_gamma)};
  matrix[1] = { 0, 1/(b * sin_gamma), (cos_beta * cos_gamma - cos_alpha)/(b*v*sin_gamma)};
  matrix[2] = {0, 0, sin_gamma/(c*v)};

  return matrix;
}


glm::mat3 math::cartToFracMatrix(props::LatticeParams &lp){
  glm::mat3 matrix;
  float cos_alpha = cos(lp.alpha * PI / 180.0f);
  float cos_beta = cos(lp.beta * PI / 180.0f);
  float cos_gamma = cos(lp.gamma * PI / 180.0f);
  float sin_gamma = sqrt(1-pow(cos_gamma,2));
  float v = sqrt( 1 - pow(cos_alpha,2) - pow(cos_beta,2) -  pow(cos_gamma,2) + 2 * cos_alpha * cos_beta * cos_gamma );

  matrix[0] = {1/lp.a, -cos_gamma/(lp.a*sin_gamma), (cos_alpha * cos_gamma - cos_beta)/(lp.a * v * sin_gamma)};
  matrix[1] = { 0, 1/(lp.b * sin_gamma), (cos_beta * cos_gamma - cos_alpha)/(lp.b*v*sin_gamma)};
  matrix[2] = {0, 0, sin_gamma/(lp.c*v)};

  return matrix;

}

void math::bringInsideCell(glm::vec3 &coor){
  //hacer un warning cuando alguno es mayor de 2

  for(int j=0;j<3;j++){
    while(coor[j] > 1.0f){
      coor[j] -= 1.0f;
    }
    while(coor[j] < 0.0f){
      coor[j] += 1.0f;
    }
    if(coor[j] < 0.00001 && coor[j] > -0.00001) coor[j] = 0.0f;
  }
}

void math::rotationMatrix3(glm::vec3 &rotVec, const float ang, glm::mat3 &rotMatrix){
    glm::vec3 unitVec = glm::normalize(rotVec);
    float x = unitVec.x;
    float y = unitVec.y;
    float z = unitVec.z;
    //make the rotation matrix with updated vector
    rotMatrix[0][0] = cos(ang) + pow(x,2)*(1-cos(ang));
    rotMatrix[0][1] = x*y*(1-cos(ang))-z*sin(ang);
    rotMatrix[0][2] = x*z*(1-cos(ang)) + y*sin(ang);
    rotMatrix[1][0] = x*y*(1-cos(ang)) + z*sin(ang);
    rotMatrix[1][1] = cos(ang)+pow(y,2)*(1-cos(ang));
    rotMatrix[1][2] = y*z*(1-cos(ang))-x*sin(ang);
    rotMatrix[2][0] = z*x*(1-cos(ang))-y*sin(ang);
    rotMatrix[2][1] = z*y*(1-cos(ang))+x*sin(ang);
    rotMatrix[2][2] = cos(ang)+pow(z,2)*(1-cos(ang));
}

void math::matMult4x1(const std::array<std::array<float, 4>, 3> &mat, std::array<float, 3> &vec){
    float suma[3] = {0.0,0.0,0.0};
    for(unsigned int row=0;row<3;row++){
       for(unsigned int col=0;col<3;col++){
          suma[row] += vec[col]*mat[row][col];
       }
    }
    for(unsigned int i=0;i<3;i++){
        vec[i] = suma[i]+mat[i][3];
    }
}

void math::matMult4x1(const std::array<std::array<float, 4>, 3> &mat, glm::vec3 &vec){
    float suma[3] = {0.0,0.0,0.0};
    for(unsigned int row=0;row<3;row++){
       for(unsigned int col=0;col<3;col++){
          suma[row] += vec[col]*mat[row][col];
       }
    }
    for(unsigned int i=0;i<3;i++){
        vec[i] = suma[i]+mat[i][3];
    }
}

void math::rotateVector(const float ang, array<float, 3> &rotvec, std::array<float, 3> &refcoors){
    float x = rotvec[0];
    float y = rotvec[1];
    float z = rotvec[2];
    float mod = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    float rotMatrix[3][3];
//
    if(mod == 0.0) return;
    //ang *= mod;
//
    x /= mod;
    y /= mod;
    z /= mod;
    //make the rotation matrix with updated vector
    rotMatrix[0][0] = cos(ang) + pow(x,2)*(1-cos(ang));
    rotMatrix[0][1] = x*y*(1-cos(ang))-z*sin(ang);
    rotMatrix[0][2] = x*z*(1-cos(ang)) + y*sin(ang);
    rotMatrix[1][0] = x*y*(1-cos(ang)) + z*sin(ang);
    rotMatrix[1][1] = cos(ang)+pow(y,2)*(1-cos(ang));
    rotMatrix[1][2] = y*z*(1-cos(ang))-x*sin(ang);
    rotMatrix[2][0] = z*x*(1-cos(ang))-y*sin(ang);
    rotMatrix[2][1] = z*y*(1-cos(ang))+x*sin(ang);
    rotMatrix[2][2] = cos(ang)+pow(z,2)*(1-cos(ang));
    
    float suma[3] = {0.0,0.0,0.0};
    for(unsigned int row=0;row<3;row++){
       for(unsigned int rowcol=0;rowcol<3;rowcol++){
          suma[row] += refcoors[rowcol]*rotMatrix[row][rowcol];
       }
    }
    for(unsigned int i=0;i<3;i++){
        refcoors[i] = suma[i];
    }

}

void math::rotateVector(float refcoors[3], float rotmat[3][3]){
    float suma[3] = {0.0,0.0,0.0};
    for(unsigned int row=0;row<3;row++){
       for(unsigned int rowcol=0;rowcol<3;rowcol++){
          suma[row] += refcoors[rowcol]*rotmat[row][rowcol];
       }
    }
    for(unsigned int i=0;i<3;i++){
        refcoors[i] = suma[i];
    }
}

void math::rotateMatrix(float mat[3][3], float rotmat[3][3]){
    /* multiplico dos matrices 3x3, este es diferente a las otras implementaciones
   porque en las otras multiplico un vector para hacerlo como las otras tendria que
transponer los vectores filas por columnas*/

    float tempmatmul[3][3] ={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    for(int row=0;row<3;row++){
            for(int col=0;col<3;col++){
                for(int rowcol=0;rowcol<3;rowcol++){
                    tempmatmul[row][col] += rotmat[row][rowcol]* mat[rowcol][col];
                }
            }
        }
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                mat[i][j] = tempmatmul[i][j];
            }
        }
}

void math::rotateVector(glm::vec3 &refcoors, float rotmat[3][3]){
    float suma[3] = {0.0,0.0,0.0};
    for(unsigned int row=0;row<3;row++){
       for(unsigned int rowcol=0;rowcol<3;rowcol++){
          suma[row] += refcoors[rowcol]*rotmat[row][rowcol];
       }
    }
    for(unsigned int i=0;i<3;i++){
        refcoors[i] = suma[i];
    }
}


float math::calcDistance(std::array<float,3> &point1, std::array<float,3> &point2){
    //return sqrt(pow((point1[0]-point2[0]), 2) + pow((point1[1]-point2[1]), 2) + pow((point1[2]-point2[2]), 2));
    //auto t1 = std::chrono::high_resolution_clock::now();
    float d = sqrt((point1[0]-point2[0])*(point1[0]-point2[0]) + 
		    (point1[1]-point2[1])*(point1[1]-point2[1]) + 
	    (point1[2]-point2[2])*(point1[2]-point2[2])); 
    //auto t2 = std::chrono::high_resolution_clock::now();

    //auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
    //std::cout<<duration<<std::endl;    

    return d;
}

float math::calcSquaredDistance(std::array<float,3> &point1, std::array<float,3> &point2){
    return (point1[0]-point2[0])*(point1[0]-point2[0]) + 
		    (point1[1]-point2[1])*(point1[1]-point2[1]) + 
	    (point1[2]-point2[2])*(point1[2]-point2[2]);
}

float math::calcDistance(float x1, float y1, float z1, float x2, float y2, float z2){
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}


float math::calcAngle(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3){
  glm::vec3 v1 = p1 - p2;
  glm::vec3 v2 = p3 - p2;

  float angle = acos(glm::dot( glm::normalize(v1), glm::normalize(v2)));
  return angle*180/3.14159265359;

}

float math::calcAngleRadians(const std::array<float, 3> &p1, const std::array<float, 3> &p2, const std::array<float, 3> &p3){
    float v1[3],v2[3];
    float n1,n2, angle;

    v1[0] = p2[0] - p1[0];
    v1[1] = p2[1] - p1[1];
    v1[2] = p2[2] - p1[2];

    v2[0] = p2[0] - p3[0];
    v2[1] = p2[1] - p3[1];
    v2[2] = p2[2] - p3[2];

    n1 = sqrt(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2));
    n2 = sqrt(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2));

    angle = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(n1*n2);
    angle = acos(angle);
    return angle;

}

float math::calcCosAngle(std::array<float, 3> v1, std::array<float, 3> v2){
    float n1,n2, cosangle;

    n1 = sqrt(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2));
    n2 = sqrt(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2));

    cosangle = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])/(n1*n2);
//    angle = acos(angle);
//    angle = angle*180/PI;

    return cosangle;

}

float math::calcAverage(std::vector<float> numbers){
    if(numbers.size()==0) return (float)0.0;

    float suma=0.0, ave;
    for(auto &n:numbers){
        suma += n;
    }
    ave = suma/numbers.size();
    return ave;
}

float math::calcDihedral(std::array<float, 3> p1, std::array<float, 3> p2, std::array<float, 3> p3, std::array<float, 3> p4){
    float b_a[3], b_c[3], c_d[3];
    float n1[3];
    float n2[3];
    float m[3];
    float x, y, ang;

    //b1
    b_a[0] = -(p1[0] - p2[0]);
    b_a[1] = -(p1[1] - p2[1]);
    b_a[2] = -(p1[2] - p2[2]);

    //b2
    b_c[0] = p2[0] - p3[0];
    b_c[1] = p2[1] - p3[1];
    b_c[2] = p2[2] - p3[2];

    //b3
    c_d[0] = p4[0] - p3[0];
    c_d[1] = p4[1] - p3[1];
    c_d[2] = p4[2] - p3[2];

    //normalize std::vectors ba, bc, cd
//    b_c =
    VectorNormalisation(b_c);
    VectorNormalisation(b_a);
    VectorNormalisation(c_d);

    CrossProduct(b_a, b_c, n1);
    CrossProduct(b_c, c_d, n2);
    CrossProduct(n1, b_c, m);

    x = DotProduct(n1, n2);
    y = DotProduct(m, n2);

    ang = 180.0 / PI * atan2(y, x);
    return ang;
}

void math::VectorNormalisation(float v[3]){
    float norm = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
    v[0] /= norm; v[1] /= norm; v[2] /= norm;
}

void math::CrossProduct(float w[3], float v[3], float cross[3]){
    cross[0] = w[1] * v[2] - w[2] * v[1];
    cross[1] = w[2] * v[0] - w[0] * v[2];
    cross[2] = w[0] * v[1] - w[1] * v[0];
}

void math::CrossProduct(array<float, 3> &p1, array<float, 3> &p2, array<float, 3> &p3, array<float, 3> &cross){
    float v[3], w[3]; 
    v[0] = p2[0] - p1[0];
    v[1] = p2[1] - p1[1];
    v[2] = p2[2] - p1[2];

    w[0] = p2[0] - p3[0];
    w[1] = p2[1] - p3[1];
    w[2] = p2[2] - p3[2];

    cross[0] = w[1] * v[2] - w[2] * v[1];
    cross[1] = w[2] * v[0] - w[0] * v[2];
    cross[2] = w[0] * v[1] - w[1] * v[0];
}

float math::DotProduct(float v[3], float w[3]) {
    return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]);
}

float math::distance(const Atom& at1, const Atom& at2){
    return glm::distance(at1.coor, at2.coor);
}

bool math::matrixIsIdentity(std::array<std::array<float, 4>, 3> &rmat){
   for(int i=0; i<3; i++){
       for(int j = 0; j < 4; j ++){
           if(i != j && rmat[i][j] != 0.0f ){
	       return false;
	   }
	   if(i == j && rmat[i][j] != 1.0f){
	       return false;
	   }
       }
   }
  return true; 
}

bool math::setDistance(Atom &at1, Atom &at2, float d){
  glm::vec3 vec = at1.coor - at2.coor;
  vec = glm::normalize(vec);
  at1.coor = at2.coor + vec * d;
  return true;
}