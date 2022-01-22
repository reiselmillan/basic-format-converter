#ifndef RW_H
#define RW_H

#include "globals.h"
#include "math/math.h"
#include "trajectory.h"  // trajectory already includes Atom and Struct

namespace rw{
  const static int XYZ_FORMAT = 0;
  const static int CIF_FORMAT = 1;
  const static int VASP_FORMAT = 2;
  const static int XTL_FORMAT = 3;  
  const static int GULP_FORMAT = 4;
  const static int DL_POLY_FORMAT = 5;


  struct FrameChooser{
    FrameChooser(){
      first = 1;
      last = -1;
      stride = 1;
    }
    FrameChooser(int f, int l, int s){
      first = f; last = l; stride = s;
    }

    int first = 1;
    int last = -1;
    int stride = 1;
  };

  struct Scene{
    bool defaultParams = true;
    bool unitCellMustBeDrawn = false;
    int width = 1200;
    int height = 1200;
    //std::string cmd = "povray +Itemp.pov +O" + outputfile + " +W"+w+" +H"+h+" +V -D +FN +Q9 -P -UD +UL +UV +A +AM2";
    float roughness = 0.002f;
    float specular = 0.6f;
    float diffuse = 1.0f;
    float ambient = 0.4f;
    float metallic = 0.1f;
    float fov = 45;
    std::string shadowless = "";
    std::string orthographic = "orthographic";
    glm::vec4 background{1.0f, 1.0f, 1.0f, 1.0f};
    glm::vec3 cameraLocation{0.0f, 0.0f, 1.0f};
    glm::vec3 cameraUp{0.0f, 1.0f, 0.0f};             // aspect ratio
    glm::vec3 cameraRight{1.0f, 0.0f, 0.0f};          // aspect ratio
    glm::vec3 sky{0.0f, 1.0f, 0.0f};
    glm::vec3 cameraLookAt{0.0f, 0.0f, 0.0f};
    glm::vec4 rangeView{-15.0f,15.0f,-15.0f,15.0f};

    glm::vec4 ambienLightColor{1.0f, 1.0f, 1.0f, 1.0f};
    glm::vec3 lightLocation{0.0f, 0.0f, 50.0f};
    glm::vec4 lightColor{1.0f, 1.0f, 1.0f, 1.0f};
    float lightFadeDistance = 500.0f;
    float lightFadePower = 0.0f;
    glm::vec3 lightPointAt{0.0, 0.0, 0.0};

  };


    std::vector<std::string> split(const std::string &line);
    std::vector<std::string> split(const std::string &s, char delim);
    void trimCString(char *atomtype);

        
    unsigned int loadFile(const char *fileName, Trajectory &traj);
    unsigned int loadFile(const char *fileName, Trajectory &traj, FrameChooser &fc);
    unsigned int loadFile(const char *fileName, Trajectory &traj, FrameChooser &fc, Scene &sc);
    Trajectory loadFile(const char *fileName);
    Trajectory loadFile(const char *fileName, FrameChooser &fc);

    unsigned int writeFile(const char *fileName, Trajectory &traj, const FrameChooser &f);
    unsigned int writeFile(const char *fileName, Trajectory &traj);
    unsigned int writeFile(const char *fileName, Trajectory &traj, const char* mode);
    unsigned int writeFile(const char *fileName, Trajectory &traj, const FrameChooser &f, const char* mode);

    unsigned int loadPDB(const char *fileName, Trajectory &traj);
    unsigned int loadCIF(const char *fileName, Trajectory &traj);
    unsigned int loadCP2K_Input(const char *fileName, Trajectory &traj);
    unsigned int loadOUTCAR(const char *fileName, Trajectory &traj);
    unsigned int loadPOSCAR(const char *fileName, Trajectory &traj);
    unsigned int loadXDATCAR(const char *fileName, Trajectory &traj);
    unsigned int loadXYZ(const char *fileName, Trajectory &traj);
    unsigned int loadXYZ(const char *fileName, Trajectory &traj, FrameChooser &fc);
    unsigned int loadXTL(const char *fileName, Trajectory &traj); //for the moment I will leave it without traj
    unsigned int loadDLPOLY_HISTORY(const char *fileName, Trajectory &traj, FrameChooser &fc);
    unsigned int loadCAR(const char *fileName, Trajectory &traj, FrameChooser &fc);
    unsigned int loadDLPOLY_CONFIG(const char *fileName, Trajectory &traj);
    unsigned int loadGULP_Input(const char *fileName, Trajectory &traj);

    unsigned int loadTigre(const char *fileName, Trajectory &traj, const FrameChooser &fc, Scene &sc);
    
    unsigned int writeCIF(const char *fileName, Trajectory &traj);
    unsigned int writePOSCAR(const char *fileName, Trajectory &traj);
    unsigned int writeXYZ(const char *fileName, Trajectory &traj, const FrameChooser &fc = FrameChooser(), const char* mode="w");
    unsigned int writeXTL(const char *fileName, Trajectory &traj);
    unsigned int writeXDATCAR(const char *fileName, Trajectory &traj, const FrameChooser &fc = FrameChooser(), const char* mode="w");
    unsigned int writePDB(const char * fileName, Trajectory &traj, const FrameChooser &fc);
    unsigned int writeDLPOLY_CONFIG(const char *filaName, Trajectory &traj);
    unsigned int writeGULP_Input(const char *filaName, Trajectory &traj);
    unsigned int writeImage(const char *fileName, Trajectory &traj, Scene &sc);
    unsigned int writeGIF(const char *fileName, Trajectory &traj, Scene &sc, rw::FrameChooser &fc);
    unsigned int writeTigre(const char *fileName, Trajectory &traj, Scene &pr, const FrameChooser &fc);
    //

};
#endif // RW_H
