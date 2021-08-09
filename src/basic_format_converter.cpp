#include <stdlib.h>
#include <rw.h>

int main(int argc, char** argv){
  if(argc < 3) {
    printf("Please provide:\n the file to be converted as first argument\n"
              " the output file name as second argument\n Formats are guess from file extension\n");
    return 0;
  }

  Trajectory  traj;
  rw::FrameChooser fr{1, 1, 1};
  int out = rw::loadFile(argv[1], traj, fr);
  if(out != 0){
    printf("Input file could not be read properly\n");
    return 0;
  }

  out = rw::writeFile(argv[2], traj, fr);
  if(out != 0){
    printf("Output file was not written properly\n");
    return 0;
  }

  return 0;

}