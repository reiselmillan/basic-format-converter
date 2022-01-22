#include <rw.h>
#include <sstream>
#include <iomanip>



using std::string;
using std::vector;
using std::array;
using rw::loadFile;
using rw::loadXYZ;
using rw::loadPOSCAR;
using rw::split;
using rw::trimCString;

void parseXYZCommentLine(char* comment, Frame &fr){
  int step = 0;
  float time = 0.0f, energy = 0.0f;
  char s1[10], s2[10], s3[10], s4[10], s5[10], s6[10];
  int n = sscanf(comment, "%s %s %i, %s %s %f, %s %s %f", s1, s2, &step, s3, s4, &time, s5, s6, &energy);
  if(n == 9 && strncmp(s1, "i", 1) == 0){
    //printf("output of parsing comment %i ,  %s: %i  time: %f  energy: %f\n", n, s1, step, time, energy);
    fr.step = step;
    fr.time = time;
    fr.energy = energy;
  }
  return;
}

void setTrajName(Trajectory& traj, const char* fileName){
   string fullDir = fileName;
   size_t found = fullDir.find_last_of("/\\");
   string justName = fullDir.substr(found+1, string::npos);
   traj.name = justName;
   traj.dir = fullDir;
}

std::vector<std::string> rw::split(const std::string &line){
    //split the string
    //boost::split(lsplitted, line, boost::is_any_of("\t "));
    std::vector<std::string> lsplitted;
    std::string buf;
    std::stringstream ss(line);       // Insert the string into a stream
    
    while (ss >> buf)  lsplitted.push_back(buf);
    return lsplitted;
}

std::vector<std::string> rw::split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        if(item.compare("") != 0){
                tokens.push_back(item);
        }
    }
    return tokens;
}

unsigned int cifOperationMatFromStr(std::array<std::array<float, 4>, 3> &rmat, 
		               std::string &line){
    
    line.erase(remove_if(line.begin(), line.end(), isspace), line.end());	
    std::vector<std::string> eq_pos = split(line, ',');  //split string line such like this  +x,-y,1/2-z

    unsigned int output = 0;

    if (eq_pos.size()!=3) {output = 1; return output;}

    const char plus = '+';
    const char minus = '-';
    const char x = 'x';
    const char y = 'y';
    const char z = 'z';
    int index;
    std::size_t found;

    //initialize matrix with zeros
    for(int i=0;i<4;i++){
        for(int j=0;j<3;j++){
            rmat[i][j] = 0.0;
        }
    }

    for(unsigned int i=0;i<3;i++){
    //parsing x---------------------------------
        // finding index of x      
	found = eq_pos[i].find(x);
	if(found == std::string::npos) index = -1;
	else index = found;
        // finding index of x      

        if(index == -1){
            rmat[i][0] = 0.0;
        }
        else if(index==0){
            rmat[i][0]=1.0;
        }
        else{
            if(eq_pos[i][index-1]==plus){
                rmat[i][0] = 1.0;
            }
            else if(eq_pos[i][index-1]==minus){
                rmat[i][0]= -1.0;
            }else{
                 rmat[i][1]= (float)0.0;
                //devolver una matriz en perticular para mostrar luego un mesaje de error
                std::cout<<"this file might be corrupted, weird character found before x\n";  
                output = 1; 
            }
        }
        //parsing y---------------------------------
        // finding index of y  
	found = eq_pos[i].find(y);
	if(found == std::string::npos) index = -1;
	else index = found;
        // finding index of y      
        if(index == -1){
            rmat[i][1] = 0.0;
        }
        else if(index==0){
            rmat[i][1]=1.0;
        }
        else{
            if(eq_pos[i][index-1]==plus){
                rmat[i][1] = 1.0;
            }
            else if(eq_pos[i][index-1]==minus){
                rmat[i][1]= -1.0;
            }else{
                rmat[i][1]= (float)0.0;
                //devolver una matriz en perticular para mostrar luego un mesaje de error
                std::cout<<"this file might be corrupted, weird character found before y\n";  
                output = 1; 
            }
        }
        //parsing z---------------------------------
        // finding index of z
	found = eq_pos[i].find(z);
	if(found == std::string::npos) index = -1;
	else index = found;
        // finding index of z      
        if(index == -1){
            rmat[i][2] = 0.0;
        }
        else if(index==0){
            rmat[i][2]=1.0;
        }
        else{
            if(eq_pos[i][index-1]==plus){
                rmat[i][2] = 1.0;
            }
            else if(eq_pos[i][index-1]==minus){
                rmat[i][2]= -1.0;
            }else{
                 rmat[i][1]= (float)0.0;
                 //devolver una matriz en perticular para mostrar luego un mesaje de error
                std::cout<<"this file might be corrupted, weird character found before z\n";  
                 output = 1; 
            }
        }

        // parse translation
        // finding index of /
	found = eq_pos[i].find('/');
	if(found == std::string::npos) index = -1;
	else index = found;
        // finding index of /      
        if(index == -1){
            rmat[i][3] = 0.0;
        }else if(index == 0){
            //devolver una matriz en perticular para mostrar luego un mesaje de error
            std::cout<<"this file might be corrupted, / character found in 0 position of string\n";  
            output = 1; 
        }else if(index == 1){    
	   //check if / is the second character, then 1/2 would be the starting of the string,  ojo: might get an error
            rmat[i][3] = (float)(eq_pos[i].at(index-1)-'0')/(float)(eq_pos[i].at(index+1)-'0');
        }else{
            if(eq_pos[i].at(index-2)==plus){
                rmat[i][3] = (float)(eq_pos[i].at(index-1)-'0')/(float)(eq_pos[i].at(index+1)-'0');
            }else if(eq_pos[i].at(index-2)==minus){
                rmat[i][3] = -(float)(eq_pos[i].at(index-1)-'0')/(float)(eq_pos[i].at(index+1)-'0');
            }else{
                rmat[i][3] = (float)(eq_pos[i].at(index-1)-'0')/(float)(eq_pos[i].at(index+1)-'0');
            }
        }
    }

    return output;
}


std::string parseCifElementSymbol(const std::string &symbol){
    //fails if first character is a whitespace
    std::string finalStr;
    size_t len = symbol.length();
    for(size_t i=0; i<len; i++){
      if(isalpha(symbol[i])){
	  finalStr.push_back(symbol[i]); 
	}
      else{return finalStr;}
    }
    return finalStr;
}

void rw::trimCString(char * s) {
    char * p = s;
    int l = strlen(p);

    while(isspace(p[l - 1])) p[--l] = 0;
    while(* p && isspace(* p)) ++p, --l;

    memmove(s, p, l + 1);
}

//function to split the lines from a cif file that contains single or double quotes
std::vector<std::string> splitCifString(const std::string &str){
    std::vector <std::string> vs;
    std::string temp;
    
    for(size_t i=0; i<str.length(); i++){
        char c = str[i];
        if( c == ' ' || c == '\t' || c == '\r' ){
    	    if (temp.empty()) continue;
            vs.push_back(temp);
    	    temp.clear();
        }else if(c == '\'' ){
            i++;
            while( str[i] != '\'' ){
    	    temp.push_back(str[i]); i++; 
    	}
            vs.push_back(temp);
    	temp.clear();
        }else if(c == '"' ){
            i++;
            while( str[i] != '"' ){
	        temp.push_back(str[i]); i++; 
	    }
            vs.push_back(temp);
	    temp.clear();
        }else if(c == '(' ){
            i++;
            while( str[i] != ')' ){
	        i++; 
            }
	    if(!temp.empty()) vs.push_back(temp);
	    temp.clear();
        }else{
    	//std::string ts (1, c);
            //temp.append(ts);
    	temp.push_back(c);
        }
    }
    vs.push_back(temp);
    return vs;    
}
// --------------------- functions to read and write files ---------------

unsigned int rw::loadFile(const char *fileName, Trajectory &traj){
    rw::FrameChooser fc;
    unsigned int out;
    out = loadFile(fileName, traj, fc);
    return out;
}

unsigned int rw::loadFile(const char *fileName, Trajectory &traj, FrameChooser &fc){
    string sfileName = fileName;
    unsigned int out = 3;

    if(sfileName.substr(sfileName.find_last_of(".") + 1) == "xyz"){
       out = rw::loadXYZ(fileName, traj, fc);
	    //may check something
    } 
    if(sfileName.substr(sfileName.find_last_of(".") + 1) == "xtl"){
       out = rw::loadXTL(fileName, traj);
	    //may check something
    }     
    else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "cif"){
        out = rw::loadCIF(fileName, traj);
	    //may check something
    }
    else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "gin"){
        out = rw::loadGULP_Input(fileName, traj);
	    //may check something
    }
//    else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "pdb"){
//        out = loadPDB(fileName, traj);
// 	//may check something
// 	return out;  
//    } 
    else if (sfileName.find("OUTCAR") != string::npos){
       out = rw::loadOUTCAR(fileName, traj);
    }
    else if (sfileName.find("POSCAR") != string::npos){
        out = rw::loadPOSCAR(fileName, traj);
    }
    else if (sfileName.find("CONTCAR") != string::npos){
        out = rw::loadPOSCAR(fileName, traj);
    }
    else if (sfileName.find("XDATCAR") != string::npos){
        out = rw::loadXDATCAR(fileName, traj);
    }
    else if (sfileName.find("HISTORY") != string::npos){
        out = rw::loadDLPOLY_HISTORY(fileName, traj, fc);
    }
    else if (sfileName.find("CONFIG") != string::npos){
        out = rw::loadDLPOLY_CONFIG(fileName, traj);
    }
    else if (sfileName.find("REVCON") != string::npos){
        out = rw::loadDLPOLY_CONFIG(fileName, traj);
    }
    setTrajName(traj, fileName);
    return out;
}
//
unsigned int rw::writeFile(const char *fileName, Trajectory &traj, const FrameChooser &fc){
  string sfileName = fileName;
  unsigned int out;

  if(sfileName.substr(sfileName.find_last_of(".") + 1) == "xyz"){
    out = rw::writeXYZ(fileName, traj, fc);
    //may check something
    return out;  
  }
  else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "xtl"){
    out = rw::writeXTL(fileName, traj);
    //may check something
    return out;  
  }
  else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "gin"){
    out = rw::writeGULP_Input(fileName, traj);
	//may check something
	return out;  
  }
  else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "cif"){
    out = rw::writeCIF(fileName, traj);
	return out;
  }
//    else if(sfileName.substr(sfileName.find_last_of(".") + 1) == "pdb"){
//        out = writePDB(fileName, traj);
// 	return out;
//    }
   else if(sfileName.find("POSCAR") != string::npos || sfileName.find("CONTCAR") != string::npos){
     out = rw::writePOSCAR(fileName, traj);
	 return out;
   } 
   else if(sfileName.find("CONFIG") != string::npos ){
      out = rw::writeDLPOLY_CONFIG(fileName, traj);
	  return out;
    }

   return 0;
}

unsigned int rw::loadDLPOLY_HISTORY(const char *fileName, Trajectory &traj, rw::FrameChooser &fc){
  FILE *file;
  file = fopen(fileName, "r");

  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }

  if(!props::init){  //init global properties 
    props::initProps();
  }

  Frame fr;
  char buffer[1024], atomType[10];
  float x, y, z, time, charge;
  int numAtoms, numFrames, frameCount = 1, cycleToRead = 0;
  int step;
  glm::mat3 cell;
  std::string symbol;

  //reading comment
  if(fgets(buffer, 1024, file) == NULL) return 1; 
  //reading header
  if(fgets(buffer, 1024, file) == NULL) return 1;
  sscanf(buffer, "%*s %*s %i %i", &numAtoms, &numFrames);

  while (fgets(buffer, 1024, file) != NULL){
    sscanf(buffer, "%*s %i %i %*s %*s %*s %f", &step, &numAtoms, &time);   //store the number of atoms to read
    if(fc.last != -1 && frameCount > fc.last){ break;  }

    if(frameCount == fc.first + cycleToRead * fc.stride ){ // frame to read
      //printf("reading frame %i %i\n", fc.last , frameCount);
        //reading cell ------------------------------------------
      if(fgets(buffer, 1024, file) == NULL) break;  //first cell line
      sscanf(buffer, "%f %f %f ", &cell[0][0], &cell[0][1], &cell[0][2]);
      if(fgets(buffer, 1024, file) == NULL) break;  //first cell line        
      sscanf(buffer, "%f %f %f ", &cell[1][0], &cell[1][1], &cell[1][2]);
      if(fgets(buffer, 1024, file) == NULL) break;  //first cell line
      sscanf(buffer, "%f %f %f ", &cell[2][0], &cell[2][1], &cell[2][2]);

      fr.setUnitCell(cell);
      //reading cell

      //reading atom i
      for(int i = 0; i < numAtoms; i++){
        if(fgets(buffer, 1024, file) == NULL){ break;   }
        sscanf(buffer, "%s %*s %*s %f ", atomType, &charge);
        symbol = parseCifElementSymbol(atomType);
        
        // read coordinates
        if(fgets(buffer, 1024, file) == NULL){ break;   }
        sscanf(buffer, "%f %f %f", &x, &y, &z);
        //add the atom to the frame
        Atom& lat = fr.addAtom(symbol.c_str(), x, y, z);
        strncpy(lat.atomType, atomType, 10);
        lat.charge = charge;
      }// end for atoms
      fr.step = step;
      fr.time = time;
      traj.addFrame(fr);
      fr = Frame();
      cycleToRead += 1;
      frameCount += 1;
    }
    else{ // pass this frame
        if(fgets(buffer, 1024, file) == NULL) break;  //first cell line
        if(fgets(buffer, 1024, file) == NULL) break;  //second cell line
        if(fgets(buffer, 1024, file) == NULL) break;  //third cell line
        //start to read atoms
        for(int i = 0; i < numAtoms; i++){
            if(fgets(buffer, 1024, file) == NULL) break; 
            if(fgets(buffer, 1024, file) == NULL) break; 
        }
        frameCount += 1;
    }
  } // end while
  printf("num frames: %i   %i\n", traj.numFrames(), traj.numAtoms());
  traj.setUnitCellLines();
  fclose(file);
  return 0;
}

unsigned int rw::loadGULP_Input(const char* fileName, Trajectory &traj){
  FILE *file;
  file = fopen(fileName, "r");
  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }  

  Frame fr;
  char buffer[1024], atomType[1024];
  float x, y, z, charge;
  props::LatticeParams lp;
  std::string symbol;
  bool readingCoors = false;
  
  while (fgets(buffer, 1024, file) != NULL){
    //check if it's the cell line
    if(strncmp(buffer, "cell", 4) == 0) {
      fgets(buffer, 1024, file);
      sscanf(buffer, "%f %f %f %f %f %f", &lp.a, &lp.b, &lp.c, &lp.alpha, &lp.beta, &lp.gamma);
      fr.setUnitCell(lp);
      continue;
    }
    //check if it's the coordinates line
    if(strncmp(buffer, "cartesian", 9) == 0) {
      traj.setCartesian(true);
      readingCoors = true;
      continue;
    }    
    if(strncmp(buffer, "fractional", 10) == 0) {
      traj.setCartesian(true);
      readingCoors = true;
      continue;
    }
    if(readingCoors){
      //reading atom i
      int n = sscanf(buffer, "%s %*s %f %f %f %f", atomType, &x, &y, &z, &charge);
      if(n != 5){
        continue; // The info to create a new atom is not complete or not atom line
      }
      symbol = parseCifElementSymbol(atomType);
      Atom &at = fr.addAtom(symbol.c_str(), x, y, z);
      strcpy(at.atomType, atomType);
      at.charge = charge;
    }
  } // end while
  traj.addFrame(fr);
  traj.setUnitCellLines();
  traj.fractionalToCartesian();
  fclose(file);
  return 0;
}

unsigned int rw::loadDLPOLY_CONFIG(const char *fileName, Trajectory &traj){
  FILE *file;
  file = fopen(fileName, "r");
  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }

  if(!props::init){  //init global properties 
    props::initProps();
  }

  Frame fr;
  char buffer[1024], atomType[1024];
  float x, y, z;
  glm::mat3 cell;
  std::string symbol;
  int levcfg = 0;

  //reading comment
  if(fgets(buffer, 1024, file) == NULL) return 1;
  //reading header
  if(fgets(buffer, 1024, file) == NULL) return 1;
  sscanf(buffer, "%i ", &levcfg);
  //reading cell ------------------------------------------
  if(fgets(buffer, 1024, file) == NULL) return 1;
  sscanf(buffer, "%f %f %f ", &cell[0][0], &cell[0][1], &cell[0][2]);
  if(fgets(buffer, 1024, file) == NULL) return 1;
  sscanf(buffer, "%f %f %f ", &cell[1][0], &cell[1][1], &cell[1][2]);
  if(fgets(buffer, 1024, file) == NULL) return 1;
  sscanf(buffer, "%f %f %f ", &cell[2][0], &cell[2][1], &cell[2][2]);
  fr.setUnitCell(cell);
  //reading cell

  while (fgets(buffer, 1024, file) != NULL){
    //reading atom i
    sscanf(buffer, "%s ", atomType);
    symbol = parseCifElementSymbol(atomType);
        
    // read coordinates
    if(fgets(buffer, 1024, file) == NULL) break;
      sscanf(buffer, "%f %f %f", &x, &y, &z);

    //add the atom to the frame
    Atom& lat = fr.addAtom(symbol.c_str(), x, y, z);
    strcpy(lat.atomType, atomType);

    if(levcfg == 1){ // read velocities
        if(fgets(buffer, 1024, file) == NULL) break;
    }else if(levcfg == 2){ // read velocities and forces
        if(fgets(buffer, 1024, file) == NULL) break;
        if(fgets(buffer, 1024, file) == NULL) break;
    }

  } // end while
  traj.addFrame(fr);
  traj.setUnitCellLines();
  fclose(file);
  return 0;
}

unsigned int rw::loadOUTCAR(const char *fileName, Trajectory &traj){
  FILE *file;
  file = fopen(fileName, "r");

  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }

   Frame fr;
   char buffer[1024], symbol[3];
   vector<string> vaspSymbols;
   vector<string> lsplitted;
   vector<unsigned int> ionsPerType; //summed, not the actual number per type
   int ibrion = -10;
   unsigned int nAtoms = 0;
   string stringLine;
   float x, y, z;
   glm::mat3 cell;

   while(fgets(buffer, 1024, file) != NULL){
        if(strstr(buffer, "POTCAR:") != NULL){
            sscanf(buffer, "%*s  %*s  %s", symbol);
            vaspSymbols.emplace_back(symbol);
	    }
        else if(strstr(buffer, "IBRION") != NULL){
            sscanf(buffer, "%*s  %*s  %i", &ibrion);
        }
        else if(strstr(buffer, "ions per type =") != NULL){
            stringLine = buffer;
            lsplitted = split(stringLine);
            for(unsigned int i=4; i<lsplitted.size(); ++i){
              nAtoms += stoi(lsplitted[i]);
              ionsPerType.emplace_back(nAtoms);
            }
        }
        else if(strstr(buffer, "direct lattice vectors") != NULL){
          if(fgets(buffer, 1024, file) == NULL) break;
          sscanf(buffer, "%f %f %f", &cell[0][0], &cell[0][1], &cell[0][2]);
          if(fgets(buffer, 1024, file) == NULL) break;
          sscanf(buffer, "%f %f %f", &cell[1][0], &cell[1][1], &cell[1][2]);
          if(fgets(buffer, 1024, file) == NULL) break;
          sscanf(buffer, "%f %f %f", &cell[2][0], &cell[2][1], &cell[2][2]);
          fr.setUnitCell(cell);
        }
        else if(strstr(buffer, "POSITION") != NULL){
            if(fgets(buffer, 1024, file) == NULL) break; // pass the ------ line
            unsigned int indexIonPerType = 0;
            for(unsigned int i = 0; i < nAtoms; ++i){
               fgets(buffer, 1024, file);
               sscanf(buffer, "%f %f %f", &x, &y, &z);
               // update the position of the label
               if(i >= ionsPerType[indexIonPerType]){
                   indexIonPerType ++;
               }
               // update the position of the label
               Atom& at = fr.addAtom(vaspSymbols[indexIonPerType].c_str(), x, y, z);
               strcpy(at.atomType, vaspSymbols[indexIonPerType].c_str());
            }

            //skip a bunch of lines to get to the energy
            for(unsigned int i = 0; i < 12; ++i) fgets(buffer, 1024, file);

            if(fgets(buffer, 1024, file) == NULL) break; //get the energy line
            sscanf(buffer, "%*s %*s %*s %*s %*s %*s %f", &fr.energy);

            //end the reading of atoms
            // traj.addFrame(fr);
            traj.addFrame(fr);
            fr = Frame();
        }//end the reading of atoms
	 else if(strstr(buffer, "  (absolute, valence and core)") != NULL){
	     // check whether traj has at least one struct
	     if(traj.numFrames() == 0) continue;
	     float nmr;
	     for(unsigned int i = 0; i < nAtoms; ++i){
            fgets(buffer, 1024, file);
            sscanf(buffer, "%*s %*s %*s %*s %f", &nmr);
            traj.atom(i).chemShift = nmr;
	     }
	 }

   }
   fclose(file);
   traj.setUnitCellLines();
   return 0;
}

unsigned int rw::loadXYZ(const char *fileName, Trajectory &traj){
  FILE *file; // Open for reading
  Frame fr;
  
  float x, y, z;
  unsigned int count = 0, n;
  char buffer[1024], symbol[3];
  
  file = fopen(fileName, "r");
  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }

  fgets(buffer, 1024, file);  //read first line to take number of atoms
  sscanf(buffer, "%i", &n);
          
  fgets(buffer, 1024, file);  //comment line
  fr.comment = buffer;
  	
  while (fgets(buffer, 1024, file) != NULL){
    if(count < n){
      //read coordinates
      sscanf(buffer, "%s %f %f %f", symbol, &x, &y, &z);
      Atom& at = fr.addAtom(symbol, x, y, z);
      strcpy(at.atomType, symbol);
      count ++;
    }
  		
    //check whether count reached number of atoms
    if(count == n){
      //reset to start reading the next strucure
      count = 0;
      traj.addFrame(fr);
      fr = Frame();
                      
      fgets(buffer, 1024, file);  //comment number of atoms
      sscanf(buffer, "%i", &n);
      fgets(buffer, 1024, file);  //comment line
      fr.comment = buffer;
    }
  } // end while
   
   fclose(file);
   return 0;
}

unsigned int rw::loadXYZ(const char *fileName, Trajectory &traj, rw::FrameChooser &fc){
  FILE *file; // Open for reading
  file = fopen(fileName, "r");
  if( file == NULL ){
    printf("Impossible to open the file named  %s !\n", fileName);
    return 1;
  }

  Frame fr;
  float x, y, z;
  //a frame will be read if frameCount == fc.first + cycleToRead*fc.stride
  unsigned int frameCount = 1, cycleToRead = 0, numAtoms;
  char buffer[1024], symbol[3];
  
  while (fgets(buffer, 1024, file) != NULL){
    sscanf(buffer, "%i", &numAtoms);   //store the number of atoms to read

    if(frameCount == fc.first + cycleToRead * fc.stride){ // frame to read
        if(fgets(buffer, 1024, file) != NULL){
            //comment line
            fr.comment = buffer;
            parseXYZCommentLine(buffer, fr);
        }  

        //start to read atoms
        for(unsigned int i = 0; i < numAtoms; i++){
            if(fgets(buffer, 1024, file) != NULL){
                sscanf(buffer, "%s %f %f %f", symbol, &x, &y, &z);
                Atom& at = fr.addAtom(symbol, x, y, z);
                strcpy(at.atomType, symbol);
            }
        }
        
        traj.addFrame(fr);
        fr = Frame();
        cycleToRead += 1;
        frameCount += 1;
    }
    else{
        if(fgets(buffer, 1024, file) == NULL) break;  //comment line
        // read atoms but pass, do not store this frames
        for(unsigned int i = 0; i < numAtoms; i++){
          if(fgets(buffer, 1024, file) == NULL) break; 
        }
        frameCount += 1;
    }
    if(fc.last != -1 && frameCount > (unsigned int)fc.last){
      break;
    }
  } // end while
   
   fclose(file);
   return 0;
}


unsigned int rw::loadXTL(const char *fileName, Trajectory &traj){
   Frame fr;
   std::ifstream file;
   std::string line;
   std::vector<std::string> lsplitted;
   std::vector<std::string> names;
   std::vector<float> charges;
   bool readingCoors = false, foundSpaceGroup = false;
   int nameCol = -1 , xCol = -1, yCol = -1, zCol = -1, chargeCol = -1;
   std::vector<glm::vec3> coors;
   std::vector<int> multiplicities;
   glm::vec3 tempcoor;
   props::LatticeParams lp;
   std::vector<glm::mat4> symmetryMatrices;

   file.open(fileName, std::ios::in);
   if(!file.is_open()){
     printf("file %s can not be read\n", fileName);
     return 1;
   }
    // setup initial properties in case they were not initialized
    if(!props::init){
      props::initProps();
    }

  while(getline(file, line)){
    if(line[0] == '#' || line[0] == '\r'){
      continue;
    }
    if(line.empty()){
      continue;
    }
    lsplitted = splitCifString(line);

    // detect space group number --------------------------------------------------------
    if(lsplitted[0].compare("SYMMETRY") == 0 && lsplitted[1].compare("NUMBER") == 0){
      traj.tableNumber = std::stoi(lsplitted[2]);
	    traj.spaceGroup = props::spaceGroups[traj.tableNumber-1].name;
	    foundSpaceGroup = true;
      continue;
	  }
       
    // detect unit cell lattice parameters
    if(lsplitted[0].compare("CELL") == 0){
      getline(file, line); //get the next line
      lsplitted = splitCifString(line); //split the cell line
      lp.a = std::stof(lsplitted[0]);
      lp.b = std::stof(lsplitted[1]);
      lp.c = std::stof(lsplitted[2]);
      lp.alpha = std::stof(lsplitted[3]);
      lp.beta = std::stof(lsplitted[4]);
      lp.gamma = std::stof(lsplitted[5]);
      continue;
    }

    if(lsplitted[0].compare("ATOMS") == 0){
      // go the next line to find the columns to read
      getline(file, line); //get the next line
      lsplitted = splitCifString(line); //split the cell line
      // determine in which column is the info
      for(unsigned int i = 0; i < lsplitted.size(); i++){
          if(lsplitted[i].compare("NAME") == 0) { nameCol = i; }
          else if(lsplitted[i].compare("X") == 0) { xCol = i; }
          else if(lsplitted[i].compare("Y") == 0) { yCol = i; }
          else if(lsplitted[i].compare("Z") == 0) { zCol = i; }
          else if(lsplitted[i].compare("CHARGE") == 0) { chargeCol = i; }
      }
      readingCoors = true;
      continue;
    }
    if(lsplitted[0].compare("EOF") == 0){
      break;
    }
    // get symmetry matrices
    if(lsplitted[0].compare("SYM") == 0 && lsplitted[1].compare("MAT") == 0){
      glm::mat4 matrix;
      matrix[0][0] = stof(lsplitted[2]);
      matrix[0][1] = stof(lsplitted[3]);
      matrix[0][2] = stof(lsplitted[4]);
      matrix[0][3] = stof(lsplitted[5]);
      matrix[1][0] = stof(lsplitted[6]);
      matrix[1][1] = stof(lsplitted[7]);
      matrix[1][2] = stof(lsplitted[8]);
      matrix[1][3] = stof(lsplitted[9]);
      matrix[2][0] = stof(lsplitted[10]);
      matrix[2][1] = stof(lsplitted[11]);
      matrix[2][2] = stof(lsplitted[12]);
      matrix[2][3] = stof(lsplitted[13]);
    //   matrix[3][0] = stof(lsplitted[14]);
    //   matrix[3][1] = stof(lsplitted[15]);
    //   matrix[3][2] = stof(lsplitted[16]);
    //   matrix[3][3] = stof(lsplitted[17]);
    }

    if(readingCoors){
      if(xCol > -1  && yCol > -1 && zCol > -1){
        tempcoor[0] = std::stof(lsplitted[xCol]);
        tempcoor[1] = std::stof(lsplitted[yCol]);
        tempcoor[2] = std::stof(lsplitted[zCol]);
        coors.push_back(tempcoor);
      }
      if(nameCol > -1){
        names.push_back(lsplitted[nameCol]);
      }
      if(chargeCol > -1){
        charges.push_back(std::stof(lsplitted[chargeCol]));
      }
      continue;
    }

   } //close the while
   
   //prepare atomic symbols
  std::vector<std::string> symbols;
  if(names.size() > 0){
    for(auto &s:names){
      symbols.push_back(parseCifElementSymbol(s));
    }
  }
  else{printf("NO LABELS COULD BE READ FROM CIF FILE\n");}

  if(symbols.size() != coors.size()){ 
    printf("READ LABELS AND COORS DO NOT MATCH: %i  %i \n",
			    (int)symbols.size(), (int) coors.size());
  }

  //get the matrix to convert from fractional to cartesian
  // the coordinates will be converted one by one as being read
  glm::mat3 fcMatrix = math::fracToCartMatrix(lp);
   //add atoms of asymetric unit
  for(unsigned int i=0; i < coors.size(); ++i){
	  if(foundSpaceGroup){
      for(auto &os:props::spaceGroups[traj.tableNumber-1].opStrings){
          tempcoor = coors[i]; 
          //take the matrix value with the key se
          math::matMult4x1(props::symStrMatOperations[os], tempcoor); 
          math::bringInsideCell(tempcoor);
          tempcoor = tempcoor * fcMatrix; // convert to cartesian
          int rep = traj.isRepeated(tempcoor, fr);
          if(rep != -1) continue;
          Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
          
          //check identity operation to assign asym unit
          if(!math::matrixIsIdentity(props::symStrMatOperations[os])){
              lat.asymUnit = false;
          }
          if(names.size() == symbols.size()){
            strcpy(lat.atomType, names[i].c_str());
          }
          if(chargeCol != -1 &&  charges.size() == symbols.size()){
            lat.charge = charges[i];
          }
      } //end for space groups
	  } //end first if
	// else if(symEquivalentSites.size() > 0){
    //     for(auto &os:symEquivalentSites){
    //         tempcoor = coors[i];
    //         //take the matrix value with the key os(operation string)
    //         std::array<std::array<float, 4>, 3> rmat;
    //         unsigned int out = cifOperationMatFromStr(rmat, os);
    //         if(out != 0) 
    //             printf("THE INFORMATION OF THE SYMMETRY OPERATIONS MIGHT NOT BE RIGHT !!");
    //         math::matMult4x1(rmat, tempcoor); 
    //         traj.bringInsideCell(tempcoor);
    //         int rep = traj.isRepeated(tempcoor, fr);
    //         if(rep != -1)  continue;
                
    //         Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
            
    //         if(atomSiteLabels.size() == symbols.size()){
    //             strcpy(lat.atomType, atomSiteLabels[i].c_str());
    //         }

    //         //define space group before doing this
    //         //check identity operation to assign asym unit
    //         //if(!matrixIsIdentity(rmat)){
    //         //    st.lastAtom().asymUnit = false;
    //         //}
	//     } //end for eq site
    //    }
	  else{
	    tempcoor = coors[i];
	    //take the matrix value with the key os(operation string)
      tempcoor = tempcoor * fcMatrix; // convert to cartesian
	    Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
      strcpy(lat.atomType, names[i].c_str());
      if(chargeCol != -1 && charges.size() == symbols.size()){
        lat.charge = charges[i];
      }	
	  }
  }// end for coors

  if(!foundSpaceGroup){ //SHOULD INCLUDE A CHECK FOR SYMMETRY SITES
    printf("!!!!!!!!!!!!  WARNING !!!!!!!!!\n");	
    printf("IMPORTANT FEATURES ARE MISSING IN THIS FILE\n");	
    printf("PLEASE CHECK THAT THE STRUCUTRE MAKES SENSE \n");	
    printf("!!!!!!!!!!!!  WARNING !!!!!!!!!\n");
  }

  fr.setUnitCell(lp);
  traj.addFrame(fr);
  traj.setUnitCellLines();
  return 0;
}

//
unsigned int rw::loadXDATCAR(const char *fileName, Trajectory &traj){
  Frame fr;

  FILE *file;
  string format = "%s";
  vector<string> symbols, lsplitted;
  vector <unsigned int> numElements;
  char buffer[1024], symbol[3];
  float x, y, z, scaleVol;
  glm::mat3 cell;
  
  file = fopen(fileName, "r");

  fgets(buffer, 1024, file);
  //read comment
  fr.comment = buffer;
  //read scale vol  -------
  fgets(buffer, 1024, file);
  sscanf(buffer, "%f", &scaleVol);
  traj.volScale = scaleVol;
  //read scale vol  -------
   
  //read cell----------------------
  for(unsigned int i=0; i<3; ++i){
    fgets(buffer, 1024, file);
	  int nread = sscanf(buffer, "%f %f %f", &cell[i][0], &cell[i][1], &cell[i][2]);
    if(nread != 3){
      printf("problem reading cell in XDATCAR file\n");
    }
    cell[i][0] *=  traj.volScale;
    cell[i][1] *=  traj.volScale;
    cell[i][2] *=  traj.volScale;
  } // cell is set below, because it is set for all frames
  
  //read cell----------------------
  //read atom symbols
  fgets(buffer, 1024, file);
  string buffCopy = buffer;
  symbols = rw::split(buffCopy);

  //read num elements
  fgets(buffer, 1024, file);
  buffCopy = buffer;
  lsplitted = rw::split(buffCopy);
  for(auto &s:lsplitted){
    numElements.push_back(stoi(s));
  }


  char *lineRead;
  lineRead = fgets(buffer, 1024, file); //read empty line or direct/cartesian line
  //DO SOMETHING, IF DIRECT OR WHATEVER
  while(lineRead != NULL){
    for(unsigned int i = 0; i < symbols.size(); ++i){
	    for(unsigned int j = 0; j < numElements[i]; ++j){
	    
        lineRead = fgets(buffer, 1024, file); 
	        
		    if(lineRead == NULL) break;
	        
		    int nread = sscanf(buffer, "%f %f %f", &x, &y, &z);
        if (nread < 3){
		      printf("error reading coordinates from XDATCAR\n");
		    }
		      Atom& at = fr.addAtom(symbols[i].c_str(), x, y, z);
          strcpy(at.atomType, symbols[i].c_str());
	      }
	      if(lineRead == NULL) break;
	    }

      fr.cell = cell; 
      fr.cartesian = false;
      traj.addFrame(fr);
	    fr = Frame();
      lineRead = fgets(buffer, 1024, file); //read empty line or direct/cartesian line
   }
   traj.setCurrentFrameIndex(0);
	 traj.fractToCartAllFrames();
   traj.setUnitCellLines();
   fclose(file);

  //printf("num frames: %i, frame 1 num at: %i\n", traj.numFrames(), traj.numAtoms());

   return 0;
}
//
unsigned int rw::loadPOSCAR(const char *fileName, Trajectory &traj){

   std::ifstream file; // Open for reading
   std::vector <std::string> lsplitted;
   std::vector <std::string> lines;
   std::vector <std::string> symbols;
   std::vector <unsigned int> numElements;
   std::vector <float> temp;
   std::string linea, firstchar;
   glm::mat3 cell;
   Frame fr;
   int initcoorpos;
   float x, y, z;
   bool vaspOld = false, selDynamics = false;

   file.open(fileName, std::ios::in);
   if(!file.is_open()){
      printf("file %s can not be read\n", fileName);
      return 1;
   }

  while(getline(file, linea)){
    lines.push_back(linea);
  }
  file.close();
  //check size of lines to see if it's a bad file if lines.size() < k not good
  fr.comment = lines[0];
  traj.volScale = atof(lines[1].c_str());

  for(unsigned int i = 2; i < 5; i++){
    lsplitted = split(lines[i]);
    cell[(i-2)][0] = atof(lsplitted[0].c_str()) * traj.volScale;
    cell[(i-2)][1] = atof(lsplitted[1].c_str()) * traj.volScale;
    cell[(i-2)][2] = atof(lsplitted[2].c_str()) * traj.volScale;
  }
  fr.setUnitCell(cell);

  lsplitted = split(lines[5]);
  if(atoi(lsplitted[0].c_str()) == 0){
    vaspOld = false;
    for(unsigned int i=0;i<lsplitted.size();i++){
      symbols.push_back(lsplitted[i]);
    }
  }else{
    vaspOld = true;
    for(unsigned int i = 0; i < lsplitted.size(); i++){
      numElements.push_back(atoi(lsplitted[i].c_str()));
    }
    lsplitted = split(fr.comment);
    for(unsigned int i=0;i<lsplitted.size();i++){
      symbols.push_back(lsplitted[i]);
    }
  }

  if(!vaspOld){
    lsplitted = split(lines[6]);
    for(unsigned int i=0;i<lsplitted.size();i++){
      numElements.push_back(atoi(lsplitted[i].c_str()));
    }
    lsplitted = split(lines[7]);
    firstchar = lsplitted[0].substr(0,1);
    if(firstchar.compare("s") == 0 || firstchar.compare("S") == 0){
      //printf("este es sel dynamics");
      selDynamics = true;
      lsplitted = split(lines[8]);
      initcoorpos = 9;
      firstchar = lsplitted[0].substr(0,1);
    }else{
      selDynamics = false;
      initcoorpos = 8;
    }
    // I am reading the line 8 if seldynamics was true. See line 1096, lsplitted = split(lines[8]);
    // or reading the line 7 of POSCAR if seldynamics is false
    if(firstchar.compare("d") == 0 || firstchar.compare("D") == 0){
      fr.cartesian = false;
    }
    else{
      fr.cartesian = true;
    }
  }else{
    lsplitted = split(lines[6]);
    firstchar = lsplitted[0].substr(0,1);
    if(firstchar.compare("s") == 0 || firstchar.compare("S") == 0){
        selDynamics = true;
        lsplitted = split(lines[7]);
        firstchar = lsplitted[0].substr(0,1);
        initcoorpos = 8;
    }else{
        selDynamics = false;
        initcoorpos = 7;
    }

    if(firstchar.compare("d") == 0 || firstchar.compare("D") == 0){
      fr.cartesian = false;
    }
    else{
      fr.cartesian = true;
    }
  }
  if(symbols.size() != numElements.size()){
       printf("Atoms and labels on line 5 and 6 of POSCAR are not the same\n");
       return 2;
  }

  for(unsigned int i = 0; i < numElements.size(); i++){
    for(unsigned int j = 0; j < numElements[i]; j++){
      lsplitted = split(lines[initcoorpos]);
      initcoorpos++;
      x = atof(lsplitted[0].c_str());
      y = atof(lsplitted[1].c_str());
      z = atof(lsplitted[2].c_str());
	    Atom& at = fr.addAtom(symbols[i].c_str(), x, y, z);
      strcpy(at.atomType, symbols[i].c_str());
      if(selDynamics){
        for(int k=3;k<6;k++){
          if(lsplitted[k].compare("F") == 0){
            at.fixed[k-3] = true;
          }else{
            at.fixed[k-3] = false;
          }
        }
      }
          //MM += masas[vasp.labels[i]];
    }
  }

  file.close();
  traj.addFrame(fr);
  traj.setUnitCellLines();

  traj.fractionalToCartesian();
  return 0;
}

unsigned int loadPDB(const char *fileName, Trajectory &traj){
  //  FILE *file;
  char buffer[1024], symbol[3]; //spaceGroup[10];
  char vector[9], x[9], y[9], z[9]; //for unit cell and coordinates
  float fx, fy, fz;
  //int index;
  int counter = 0;
  FILE* file = fopen(fileName, "r");
  Frame fr;

  bool foundBondInfo = false;

  while(fgets(buffer, 1024, file) != NULL){
    if(strncmp(buffer, "MODEL ", 6) == 0) {
      continue;
	  }

	//read unit cell    
    if(strncmp(buffer, "CRYST1", 6) == 0) {  //"CRYST1 xxxxx xxxxx xxxxx xxxxx xxxxxxx"
      float temp_a, temp_b, temp_c, temp_alpha, temp_beta, temp_gamma;
	    strncpy(vector, buffer+6, 9);
	    vector[9] = '\0';
	    temp_a = atof(vector);
	    strncpy(vector, buffer+15, 9);
	    vector[9] = '\0';
	    temp_b = atof(vector);
	    strncpy(vector, buffer+24, 9);
	    vector[9] = '\0';
	    temp_c = atof(vector);
	    strncpy(vector, buffer+33, 7);
	    vector[7] = '\0';
	    temp_alpha = atof(vector);
	    strncpy(vector, buffer+40, 7);
	    vector[7] = '\0';
	    temp_beta = atof(vector);
	    strncpy(vector, buffer+47, 7);
	    vector[7] = '\0';
	    temp_gamma = atof(vector);

      Cell cell(temp_a, temp_b, temp_c, temp_alpha, temp_beta, temp_gamma);
      fr.setUnitCell(cell);

	    continue;
	}

	// read atom
	if(strncmp(buffer, "ATOM  ", 6) == 0 || strncmp(buffer, "HETATM", 6)==0){
    strncpy(x, buffer+30, 8);
    x[8] = '\0';
    fx = atof(x);
    strncpy(y, buffer+38, 8);
    y[8] = '\0';
    fy = atof(y);
    strncpy(z, buffer+46, 8);
    z[8] = '\0';
    fz = atof(z);
    strncpy(symbol, buffer+76, 2);
    symbol[2] = '\0';
	  Atom &at = fr.addAtom(symbol, fx, fy, fz);

    //atom type 
    char atomtype[5], residue[4];
    strncpy(atomtype, buffer+12, 4);
    atomtype[4] = '\0';
    strncpy(at.atomType, atomtype, 4);
    
    // read residue
    strncpy(residue, buffer+17, 3);
    residue[3] = '\0';
    strncpy(at.residue, residue, 3);
	  continue;
	}

	// if(strncmp(buffer, "CONECT", 6) == 0 ){
  //   foundBondInfo = true;
  //   char atom1[5], bonded[5];
  //   strncpy(atom1, buffer+7, 5);
  //   atom1[5] = '\0';
  //   int atom1Index = atoi(atom1) - 1;

  //   //starting adding bonds
  //   //loop over the 4 fields
  //   array<int, 4> initColumn = {12, 17, 22, 27};
  //   for(auto column:initColumn){
  //     strncpy(bonded, buffer + column, 5);
  //     bonded[5] = '\0';
  //     int bondedIndex = atoi(bonded) - 1; //if zero probably because of empty string
  //     if(bondedIndex == -1)
  //       continue;
  //     }
          
  //   continue;
	// }
	
	if(strncmp(buffer, "ENDMDL", 6) == 0 || strncmp(buffer, "END", 3) == 0){
	    // it is used at the end of the while loop to check if structure has been added to traj
	    counter ++; 
      traj.addFrame(fr);
	    fr = Frame();
	    continue;
	  }
  } // end while
  
   fclose(file);
   return 0;
}



unsigned int rw::loadCIF(const char *fileName, Trajectory &traj){
   Frame fr;
   std::ifstream file;
   std::string line;
   std::vector<std::string> loopHeaders;
   std::vector<std::string> lsplitted;
   std::vector<std::string> symEquivalentSites, labels, atomSiteLabels, atomSiteTypeSymbols;
   std::vector<std::string> atomTypeSymbols;
   bool readingLoop = false;
   bool readingLoopHeaders = false;
   bool foundSpaceGroup = false;
   int loopHeaderX = -1, loopHeaderY = -1, loopHeaderZ = -1, symetricPos = -1;
   int atomSymMult = -1, atomSiteTypeSymbol = -1, atomSiteLabel = -1;
   int atomTypeSymbol = -1;
   std::vector<glm::vec3> coors;
   std::vector<int> multiplicities;
   glm::vec3 tempcoor;
   props::LatticeParams lp;

   file.open(fileName, std::ios::in);
   if(!file.is_open()){
     printf("file %s can not be read\n", fileName);
     return 1;
   }
    // setup initial properties in case they were not initialized
    if(!props::init){
      props::initProps();
    }

   while(getline(file, line)){
     if(line[0] == '#' || line[0] == '\r'){
           continue;
       }
       if(line.empty()){
           readingLoop = false;
           continue;
       }

       lsplitted = splitCifString(line);

       // detect space group number --------------------------------------------------------
       if(lsplitted[0].compare("_symmetry_Int_Tables_number") == 0){
         traj.tableNumber = std::stoi(lsplitted[1]);
	     traj.spaceGroup = props::spaceGroups[traj.tableNumber-1].name;
	     foundSpaceGroup = true;
         //sym_equivalent_sites = spaceGroups[tableNumber-1].symOperations;
         continue;
	 }
       
	 if(lsplitted[0].compare("_space_group_IT_number") == 0){
       traj.tableNumber = std::stoi(lsplitted[1]);
	   traj.spaceGroup = props::spaceGroups[traj.tableNumber-1].name;
	   foundSpaceGroup = true;
        //sym_equivalent_sites = spaceGroups[tableNumber-1].symOperations;
        continue;
     }
       // detect unit cell lattice parameters
       if(lsplitted[0].compare("_cell_angle_alpha") == 0){
           lp.alpha = std::stof(lsplitted[1]);
           continue;
       }
       if(lsplitted[0].compare("_cell_angle_beta") == 0){
           lp.beta = std::stof(lsplitted[1]);
           continue;
       }
       if(lsplitted[0].compare("_cell_angle_gamma") == 0){
           lp.gamma = std::stof(lsplitted[1]);
           continue;
       }
       if(lsplitted[0].compare("_cell_length_a") == 0){
           lp.a = std::stof(lsplitted[1]);
           continue;
       }
       if(lsplitted[0].compare("_cell_length_b") == 0){
           lp.b = std::stof(lsplitted[1]);
           continue;
       }
       if(lsplitted[0].compare("_cell_length_c") == 0){
           lp.c = std::stof(lsplitted[1]);
           continue;
       }

       //check if I am reading some other header outside the loop
       if(readingLoop && !readingLoopHeaders && lsplitted[0][0] == '_'){
           readingLoop = false;
           //I shoudl read something here
           continue;
       }
	
       if(lsplitted[0].compare("loop_") == 0){
           readingLoop = true;
           readingLoopHeaders = true;
           loopHeaders.clear();
           loopHeaderX = -1;
           loopHeaderY = -1;
           loopHeaderZ = -1;
           symetricPos = -1;
           atomSymMult = -1;
           atomSiteTypeSymbol = -1;
           atomSiteLabel = -1;
	       atomTypeSymbol = -1;
           continue;
       }

       if(readingLoop){
           if(lsplitted[0][0] == '_'){
               readingLoopHeaders = true;
               loopHeaders.push_back(lsplitted[0]);
           }
           else{
		if(readingLoopHeaders){
            readingLoopHeaders = false;
            //find out the indices of interest in headers list ---------------------
            for(unsigned int i=0;i<loopHeaders.size();i++){
                if(loopHeaders[i].compare("_symmetry_equiv_pos_as_xyz")==0){
                    symetricPos = i;
                }
                else if(loopHeaders[i].compare("_space_group_symop_operation_xyz")==0){
                    // this is just another keyword in the cif file for the same thing
                    symetricPos = i;
                }
                else if(loopHeaders[i].compare("_atom_site_fract_x") == 0){
                    loopHeaderX = i;
                }
                else if(loopHeaders[i].compare("_atom_site_fract_y") == 0){
                    loopHeaderY = i;
                }
                else if(loopHeaders[i].compare("_atom_site_fract_z") == 0){
                    loopHeaderZ = i;
                }
                else if(loopHeaders[i].compare("_atom_site_symmetry_multiplicity") == 0){
                    atomSymMult = i;
                }
                else if(loopHeaders[i].compare("_atom_site_type_symbol") == 0){
                    atomSiteTypeSymbol = i;
                }
                else if(loopHeaders[i].compare("_atom_site_label") == 0){
                    atomSiteLabel = i;
                }
                else if(loopHeaders[i].compare("_atom_type_symbol") == 0){
                    atomTypeSymbol = i;
                }
            }//end for, find out the indices of interest ------------------------------------
		}// end if readingLoopHeaders
           } //end else setting readLoopHeader to false and find out the indices of interest ------------------------------------

	    if(readingLoopHeaders) continue;
           if(loopHeaderX >-1){
               tempcoor[0] = stof(lsplitted[loopHeaderX]);
               tempcoor[1] = stof(lsplitted[loopHeaderY]);
               tempcoor[2] = stof(lsplitted[loopHeaderZ]);
               coors.push_back(tempcoor);
           }
           if(atomSiteTypeSymbol > -1){
               atomSiteTypeSymbols.push_back(lsplitted[atomSiteTypeSymbol]);
           }
           if(atomSymMult >-1){
               multiplicities.push_back(stoi(lsplitted[atomSymMult]));
           }
           if(symetricPos > -1){
               symEquivalentSites.push_back(lsplitted[symetricPos]);
           }
           if(atomSiteLabel > -1){
              atomSiteLabels.push_back(lsplitted[atomSiteLabel]);
           }
           if(atomTypeSymbol > -1){
              atomTypeSymbols.push_back(lsplitted[atomTypeSymbol]);
           }
       }// end if readingLoop

   } //close the while
   
  //prepare atomic symbols
  std::vector<std::string> symbols;
  if(atomTypeSymbols.size() > 0){
    for(auto &s:atomTypeSymbols){
	    symbols.push_back(parseCifElementSymbol(s));
	  }
  }
  else if(atomSiteTypeSymbols.size()>0){
    for(auto &s:atomSiteTypeSymbols){
	    symbols.push_back(parseCifElementSymbol(s));
	  }
  }
  else if(atomSiteLabels.size()>0){
    for(auto &s:atomSiteLabels){
	    symbols.push_back(parseCifElementSymbol(s));
	  }
  }
  else{printf("NO LABELS COULD BE READ FROM CIF FILE\n");}

  if(symbols.size() != coors.size()){ 
	    printf("READ LABELS AND COORS DO NOT MATCH: %i  %i \n",
			    (int)symbols.size(), (int) coors.size());
  }

  //get the matrix to convert from fractional to cartesian
  // the coordinates will be converted one by one as being read
  glm::mat3 fcMatrix = math::fracToCartMatrix(lp);
  //add atoms of asymetric unit
  for(unsigned int i=0; i < coors.size(); ++i){
    if(foundSpaceGroup){
      for(auto &os:props::spaceGroups[traj.tableNumber-1].opStrings){
        tempcoor = coors[i];
        
        //take the matrix value with the key se
        math::matMult4x1(props::symStrMatOperations[os], tempcoor); 
        math::bringInsideCell(tempcoor);
        tempcoor = tempcoor * fcMatrix; // convert to cartesian
        int rep = traj.isRepeated(tempcoor, fr);
        if(rep != -1) continue;
        Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
        
        //check identity operation to assign asym unit
        if(!math::matrixIsIdentity(props::symStrMatOperations[os])){
          lat.asymUnit = false;
        }
        if(atomSiteLabels.size() == symbols.size()){
          strcpy(lat.atomType, atomSiteLabels[i].c_str());
        }
      } //end for space groups
    } //end first if
    else if(symEquivalentSites.size() > 0){
      for(auto &os:symEquivalentSites){
          tempcoor = coors[i];
          //take the matrix value with the key os(operation string)
          std::array<std::array<float, 4>, 3> rmat;
          unsigned int out = cifOperationMatFromStr(rmat, os);
          if(out != 0) 
              printf("THE INFORMATION OF THE SYMMETRY OPERATIONS MIGHT NOT BE RIGHT !!");
          math::matMult4x1(rmat, tempcoor); 
          math::bringInsideCell(tempcoor);
          tempcoor = tempcoor * fcMatrix; // convert to cartesian
          int rep = traj.isRepeated(tempcoor, fr);
          if(rep != -1) continue;
                
          Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
          
          if(atomSiteLabels.size() == symbols.size()){
            strcpy(lat.atomType, atomSiteLabels[i].c_str());
          }

          //define space group before doing this
          //check identity operation to assign asym unit
          //if(!matrixIsIdentity(rmat)){
          //    st.lastAtom().asymUnit = false;
          //}
      } //end for eq site
    }
    else{
        tempcoor = coors[i];
        //take the matrix value with the key os(operation string)
        Atom& lat = fr.addAtom(symbols[i].c_str(), tempcoor);
        strcpy(lat.atomType, symbols[i].c_str());
    }
  }// end for coors

  if(!foundSpaceGroup && symEquivalentSites.size() == 0){
    printf("!!!!!!!!!!!!  WARNING !!!!!!!!!\n");	
    printf("IMPORTANT FEATURES ARE MISSING IN THIS FILE\n");	
    printf("PLEASE CHECK THAT THE STRUCUTRE MAKES SENSE \n");	
    printf("!!!!!!!!!!!!  WARNING !!!!!!!!!\n");	
  }

  fr.setUnitCell(lp);
  traj.addFrame(fr);
  traj.setUnitCellLines();
  printf("traj num atoms %i \n", traj.numAtoms());
  return 0;
}



unsigned int rw::writeCIF(const char * fileName, Trajectory &traj){
   traj.cartesianToFractional();

   FILE *file;
   
   file = fopen(fileName, "w");
   fprintf(file, "#*******************************************#\n");
   fprintf(file, "#File created by Tigre\n");
   fprintf(file, "#*******************************************#\n");

   props::LatticeParams lp;
   props::getLatticeParameters(lp, traj.cell());
   fprintf(file, "_cell_length_a               % 10.4f\n", lp.a);
   fprintf(file, "_cell_length_b               % 10.4f\n", lp.b);
   fprintf(file, "_cell_length_c               % 10.4f\n", lp.c);
   fprintf(file, "_cell_angle_alpha               % 10.4f\n", lp.alpha);
   fprintf(file, "_cell_angle_beta               % 10.4f\n", lp.beta);
   fprintf(file, "_cell_angle_gamma               % 10.4f\n", lp.gamma);
   fprintf(file, "\n");
   fprintf(file, "_symmetry_Int_Tables_number      %i \n", traj.tableNumber);
   fprintf(file, "_symmetry_space_group_name_H-M      %s \n", traj.spaceGroup.c_str());
   //fprintf(file, "_symmetry_Int_Tables_number      1 \n");
   fprintf(file, "\n");
   fprintf(file, "loop_\n");
   fprintf(file, "_symmetry_equiv_pos_as_xyz     \n");
   //fprintf(file, "'x,y,z'     \n");
   for(auto s:props::spaceGroups[traj.tableNumber-1].opStrings){
     fprintf(file, " '%s'  \n", s.c_str());
   }

   fprintf(file, "\n");
   fprintf(file, "loop_\n");
   fprintf(file, "_atom_site_label\n");
   fprintf(file, "_atom_site_type_symbol\n");
   fprintf(file, "_atom_site_fract_x\n");
   fprintf(file, "_atom_site_fract_y\n");
   fprintf(file, "_atom_site_fract_z\n");

   for(unsigned int i=0; i<traj.numAtoms(); i++){
    Atom& at = traj.atom(i);
	if(!at.asymUnit) continue;
    if(strcmp(at.atomType, "") == 0) strcpy(at.atomType, at.symbol()); //esto no funciona

    fprintf(file, "%-10s    %s    % 10.6f   % 10.6f   % 10.6f \n", 
               at.atomType,
               at.symbol(),
               at.coor.x,
               at.coor.y,
               at.coor.z);
   }

   fprintf(file, "\n");
   fclose(file);
   traj.fractionalToCartesian();
   return 0;
}


unsigned int rw::writeXTL(const char * fileName, Trajectory &traj){
   traj.cartesianToFractional();

   FILE *file;
   
   file = fopen(fileName, "w");
   props::LatticeParams lp;
   props::getLatticeParameters(lp, traj.cell());
   fprintf(file, "CELL\n");
   fprintf(file, "% 10.4f % 10.4f % 10.4f % 10.4f % 10.4f % 10.4f \n", lp.a, lp.b, lp.c,
      lp.alpha, lp.beta, lp.gamma);

   fprintf(file, "SYMMETRY  NUMBER %i  LABEL %s\n", traj.tableNumber, traj.spaceGroup.c_str());
   fprintf(file, "ATOMS\n");
   fprintf(file, "NAME         X          Y          Z\n");

   for(unsigned int i=0; i<traj.numAtoms(); i++){
    Atom& at = traj.atom(i);
	if(!at.asymUnit) continue;
    if(strcmp(at.atomType, "") == 0) strcpy(at.atomType, at.symbol()); //esto no funciona

    fprintf(file, "%s     % 10.4f   % 10.4f   % 10.4f \n", 
               at.atomType,
               at.coor.x,
               at.coor.y,
               at.coor.z);
   }

   fprintf(file, "EOF\n");
   fclose(file);
   traj.fractionalToCartesian();
   return 0;
}


//
//unsigned int writePDB(const char * fileName, Trajectory &traj){
//    FILE *file;
//    file = fopen(fileName, "w");
//    
//    for(unsigned int i = 0; i < traj.numStructures(); i++){
//        Structure &st = traj.structureById(i);
//
//	st.shiftBackToOriginalCentroid();
//
//        //can fail for long cell unit vectors ~100
//	fprintf(file, "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f\n", st.a, st.b, st.c, st.alpha, st.beta, st.gamma);
//        fprintf(file, "MODEL     %d\n", i+1);
//        //printf("REMARK 250 time=%.3f picoseconds\n", timeInPs);
//        //printf("REMARK \n");
//	for(unsigned int j = 0; j < st.totNumberOfAtoms(); j++){
//            if(st.atom(j).deleted) continue; 
//	    // Hetatm label
//            fprintf(file, "HETATM");
//            fprintf(file, "%5d", j+1);  //serial number
//	    fprintf(file, " ");  //empty character
//	    fprintf(file, "%4s", st.atom(j).atomType.c_str()); //atom name
//	    fprintf(file, " ");   //character alternate location 
//	    fprintf(file, "UNK"); //residue name
//            fprintf(file, " A");  //chain identifier
//	    fprintf(file, "   1"); //residue sequence number
//            fprintf(file, "    "); //code for insertion of residues
//            fprintf(file, "%8.3f%8.3f%8.3f", 
//			    st.atom(j).x(), st.atom(j).y(), st.atom(j).z());
//            fprintf(file, "  1.00"); //occupancy (6.2)
//	    fprintf(file, "  0.00"); //temperature factor 6.2
//	    fprintf(file, "      "); //empty  no especification
//            fprintf(file, "    "); //segment identifier 4 left justified
//            //atom symbol , right justified, in C is right justified by default,
//            fprintf(file, "%2s", st.atom(j).cymbol());  
//	    fprintf(file, "%2d\n", st.atom(j).charge); //charge      
//        }
//        fprintf(file, "ENDMDL\n"); // end of trajectory frame
//
//	st.shiftToZeroCentroid();
//    }
//
//    fclose(file);
//    return 0;
//}


unsigned int rw::writePOSCAR(const char *fileName, Trajectory &traj){
   FILE *outfile;
   outfile = fopen(fileName, "w");
   int lpos;
   bool found = false;
   std::vector <std::vector<int>> index_group;
   //temporary vector that stores the indices that go inside index groups
   std::vector<int> tempvector, atoms; 
   std::vector<std::string> labels;

   //reviso el numero de atoms y labels
   for(unsigned int i=0; i < traj.numAtoms(); ++i){
       found = false;
       tempvector.clear();
       for(unsigned int j = 0; j < labels.size(); j++){
         if(labels[j].compare(props::atomicSymbols[traj.atom(i).atomicNum]) == 0){
             found = true;
             lpos = j;
             break;
         }
       }
       if(!found){
           labels.push_back(props::atomicSymbols[traj.atom(i).atomicNum]);
           atoms.push_back(1);
           tempvector.push_back(i);
           index_group.push_back(tempvector);
       }else{
           atoms[lpos]++;
           index_group[lpos].push_back(i);
       }
   }

   //comment line
   for(unsigned int i=0;i < labels.size(); i++){
       fprintf(outfile, "   %s  ", labels[i].c_str());
   }
   fprintf(outfile, " \n  %f \n", traj.volScale);

   // write cell 
   glm::mat3 cell = traj.cell();
   for(unsigned int i=0; i<3; i++){
       for(unsigned int j=0; j<3; j++){
           fprintf(outfile, "  % 20.16f    ", cell[i][j]/traj.volScale );
       }
       fprintf(outfile, " \n");
   }
   // write labels
  
   for(unsigned int i=0; i < labels.size(); i++){
     fprintf(outfile, "   %s  ", labels[i].c_str());
   }
   fprintf(outfile, " \n");
   
   //write atoms
   for(unsigned int i=0;i < atoms.size();i++){
       fprintf(outfile, "   %i  ", atoms[i]);
   }
   fprintf(outfile, " \n");

//    if(st.outConstraints){
//         fprintf(outfile, "Selective dynamics \n");
//    }
   if(traj.cartesian()){
        fprintf(outfile, "Cartesian \n");
   }else{
        fprintf(outfile, "Direct \n");
   }

   for(unsigned int i=0;i < index_group.size(); i++){
       for(unsigned int k=0; k < index_group[i].size(); k++){
           unsigned int index = index_group[i][k];
           for(unsigned int j=0; j<3; j++){
               fprintf(outfile, " % 20.16f", traj.atom(index).coor[j]);
           }
           //write selective dynamics flags
        //    if(st.outConstraints){
        //        for(unsigned int j=0;j<3;j++){
        //            if(st.atom(index).fixed[j]){
        //                fprintf(outfile, "   F   ");
        //            }else{
        //                fprintf(outfile, "   T   ");
        //            }
        //        }
        //    }
           fprintf(outfile, " \n");
       }
   }

   fclose(outfile);
   return 0;
}


unsigned int rw::writeXYZ(const char *fileName, Trajectory &traj, const FrameChooser &fc){
  FILE *file;
  file = fopen(fileName, "w");
  if( file == NULL ){
    printf("Impossible to open the file named  %s for writing!\n", fileName);
    return 1;
  }
  int first = fc.first - 1;
  int last;
  if(fc.last < 0 || (unsigned int)fc.last >= traj.numFrames()){
      last = traj.numFrames();
  }
  else{
      last = fc.last;
  }
  for(int i = first; i < last; i+= fc.stride){
      Frame &fr = traj.frame(i);
      fprintf(file, "%i\n", fr.numAtoms());
      fprintf(file, "\n"); //maybe put comment
      //write symbol and coordinates
      for(unsigned int j = 0; j < fr.numAtoms(); j++){
        Atom &at = fr.atom(j);
        fprintf(file, "%s   %16.8f %16.8f %16.8f\n", 
          props::atomicSymbols[at.atomicNum], at.coor.x, at.coor.y, at.coor.z);
      }
  }
   fclose(file);
   return 0;
}


unsigned int rw::writeDLPOLY_CONFIG(const char *fileName, Trajectory &traj){
  FILE *file, *fieldfile;
  string fieldFileName = fileName;
  fieldFileName += ".field";
  fieldfile = fopen(fieldFileName.c_str(), "w");  //fix this later, file must be in the same direcoty
  file = fopen(fileName, "w");
  if( file == NULL ){
    printf("Impossible to open the file named  %s for writing!\n", fileName);
    return 1;
  }
  //write CONFIG FILE  
  fprintf(file, "%s  file created with tigre\n", traj.name.c_str());
  fprintf(file, "         0         3\n");
  glm::mat3 &cell = traj.cell();  //write cell
  for(int i=0; i<3; i++){
    fprintf(file, "%20.10f %20.10f %20.10f\n", cell[i][0], cell[i][1], cell[i][2]);
  } //write cell

  for(unsigned int i = 0; i<traj.numAtoms(); i++){
    Atom &at = traj.atom(i);
    fprintf(file, "%s          %i\n", at.atomType, i+1);
    fprintf(file, "%20.10f %20.10f %20.10f\n", at.coor.x, at.coor.y, at.coor.z);
  }

  //write field file
  if( fieldfile == NULL ){
    printf("Impossible to open the file named FIELD for writing!\n");
    return 1;
  }

  // i have to write it with right order, by molecules to match the FIELD FILE
  fprintf(fieldfile, "NUMMOLS%14d\n", 1);
  fprintf(fieldfile, "ATOMS%17d\n", (int)traj.numAtoms());
  for(unsigned int i = 0; i < traj.numAtoms(); i++){
    Atom &at = traj.atom(i);
    fprintf(fieldfile, "%s    %f     %f       %i     %i\n", at.atomType, props::masas[at.atomicNum], at.charge, 1, 0);
  }
  fprintf(fieldfile, "BONDS%17d\n", (int)traj.numBonds());
  for(unsigned int i = 0; i< traj.numBonds(); i++){
    Bond& b = traj.bond(i);
    fprintf(fieldfile, "harm     %i    %i     \n", b.at1+1, b.at2+1);
  }
  fprintf(fieldfile, "FINISH\n");

  fclose(file);
  fclose(fieldfile);
  return 0;
}



unsigned int rw::writeGULP_Input(const char *fileName, Trajectory &traj){
  FILE *file;
  file = fopen(fileName, "w");
  if( file == NULL ){
    printf("Impossible to open the file named  %s for writing!\n", fileName);
    return 1;
  }

  //if(traj.numAtoms() > )  make a var to check whether atomicPropertes have been setup
  fprintf(file, "single\n");
  fprintf(file, "cell\n");
  props::LatticeParams lp;
  props::getLatticeParameters(lp, traj.cell());
  fprintf(file, "%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f\n", lp.a, lp.b, lp.c, lp.alpha, lp.beta, lp.gamma);
  fprintf(file, "cartesian\n");
  for(unsigned int i = 0; i < traj.numAtoms(); i++){
    Atom &at = traj.atom(i);
    fprintf(file, "%-5s core %12.6f %12.6f %12.6f %12.6f\n", at.atomType, 
                                                            at.coor.x, 
                                                            at.coor.y, 
                                                            at.coor.z, 
                                                            at.charge);
  }
  fclose(file);
  return 0;
}

unsigned int rw::writeImage(const char *fileName, Trajectory &traj, rw::Povray &pr){
  setlocale(LC_NUMERIC, "en_EN");
  FILE *file;
  string fileNameTemp = "temp.pov";
  file = fopen(fileNameTemp.c_str(), "w");
  if( file == NULL ){
    printf("Impossible to open the file named  %s for writing!\n", fileName);
    return 1;
  }
  fprintf(file, "global_settings {\n"
                    "        ambient_light rgb <%f, %f, %f>\n"
                    "       max_trace_level 15\n"
                    "}\n\n", pr.ambienLightColor[0], pr.ambienLightColor[1], pr.ambienLightColor[2]);
  
  fprintf(file, "background { color rgb <%f, %f, %f> }\n\n", pr.background[0], pr.background[1], pr.background[2]);
  fprintf(file, "camera {\n       %s\n       angle %f\n", pr.orthographic.c_str(), pr.fov);
  fprintf(file, "       location <%f, %f, %f>\n", pr.cameraLocation[0], pr.cameraLocation[1], pr.cameraLocation[2]);
  fprintf(file, "       sky <%f, %f, %f>\n", pr.sky[0], pr.sky[1], pr.sky[2]);
  fprintf(file, "       up <%f, %f, %f>\n", pr.cameraUp[0], pr.cameraUp[1], pr.cameraUp[2]);
  fprintf(file, "       right <%f, %f, %f>\n", pr.cameraRight[0], pr.cameraRight[1], pr.cameraRight[2]);
  fprintf(file, "       look_at <%f, %f, %f> }\n", pr.cameraLookAt[0], pr.cameraLookAt[1], pr.cameraLookAt[2]);

  fprintf(file, "light_source {\n");
  fprintf(file, "        <%f, %f, %f>\n", pr.lightLocation[0], pr.lightLocation[1], pr.lightLocation[2]);
  fprintf(file, "        color rgb <%f, %f, %f> %s \n", pr.lightColor[0], pr.lightColor[1], pr.lightColor[2], pr.shadowless.c_str());
  //                  "        fade_distance %f\n"
  //                  "        fade_power %f\n"
  fprintf(file,       "        parallel\n");
  fprintf(file, "        point_at <%f, %f, %f>\n}\n", pr.lightPointAt[0], pr.lightPointAt[1], pr.lightPointAt[2]);
  fprintf(file,  "#default {\n        finish {ambient %f diffuse %f specular %f roughness %f metallic %f}\n}\n",
                                pr.ambient, pr.diffuse, pr.specular, pr.roughness, pr.metallic);


  //fprintf(file,  "union {\n");
  //start with bonds
  for(unsigned int i=0; i < traj.numBonds(); i++){
    Bond &bd = traj.bond(i);
    Atom &at1 = traj.atom(bd.at1);
    Atom &at2 = traj.atom(bd.at2);
    if(at1.hidden || at2.hidden){
      continue;
    }
    glm::vec3 dirVec = at2.coor - at1.coor;
    dirVec = glm::normalize(dirVec);
    //glm::vec3 fpoint = at2.coor - dirVec * at2.sphere(); //ends just where the spheres is supposed to be drawn
    glm::vec3 ipoint = at1.coor + dirVec * at1.sphere(); //ends just where the spheres is supposed to be drawn
    float distance = bd.length - at1.sphere() - at2.sphere(); //distance from initial sphere to final sphere
    //glm::vec3 mp = at1.coor + dirVec*(bd.length*0.5f); //midpoint between two atomic coordinates
    glm::vec3 mp = ipoint + dirVec*(distance * 0.5f); //midpoint between two atomic coordinates

    //first bond
    fprintf(file, "cylinder {\n        <%f, %f, %f>,  <%f, %f, %f>, %f  \n", 
                    at1.coor[0], at1.coor[1], at1.coor[2], mp[0], mp[1], mp[2], bd.cylRadius);

    fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at1.color()[0],at1.color()[1],at1.color()[2],1-at1.color()[3]);
    
    //second bond
    fprintf(file, "cylinder {\n        <%f, %f, %f>,  <%f, %f, %f>, %f  \n", 
                    mp[0], mp[1], mp[2], at2.coor[0], at2.coor[1], at2.coor[2], bd.cylRadius);

    fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at2.color()[0],at2.color()[1],at2.color()[2],1-at2.color()[3]);
   
  }

  //start with atoms
  for(unsigned int i = 0; i < traj.numAtoms(); i++){
      Atom &at = traj.atom(i);
      if(at.hidden){
        continue;
      }
      fprintf(file, "sphere {\n        <%f, %f, %f>, %f  \n", 
                                at.coor[0], at.coor[1], at.coor[2], at.sphere());

    fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at.color()[0],at.color()[1],at.color()[2],1-at.color()[3]);
    
  }

  //fprintf(file, "}\nmerge {\n}");
  fclose(file);
  string f(fileName);
  string cmd = "povray -Itemp.pov -O" + f + " +W"+ std::to_string(pr.width) + 
                " +H" + std::to_string(pr.height) + "+V  +FN +Q9 -P -UD +UL +UV +A +AM2 -d";
  system(cmd.c_str());
  //system("rm ./temp.pov");
  return 0;
}


unsigned int rw::writeGIF(const char *fileName, Trajectory &traj, rw::Povray &pr, rw::FrameChooser& fc){
  setlocale(LC_NUMERIC, "en_EN");
  FILE *file;
  string fileNameTemp = "temp.pov";
  file = fopen(fileNameTemp.c_str(), "w");
  if( file == NULL ){
    printf("Impossible to open the file named  %s for writing!\n", fileName);
    return 1;
  }
  printf("frame chooser: %i %i %i\n", fc.first, fc.last, fc.stride);

  for(unsigned int frameInd = fc.first; frameInd < fc.last; frameInd+=fc.stride){
    bool setf = traj.setCurrentFrameIndex(frameInd);
    printf("creating movie for frame %i   set? %i\n", frameInd, setf);
    traj.generateBonds();
    fprintf(file, "global_settings {\n"
                        "        ambient_light rgb <%f, %f, %f>\n"
                        "       max_trace_level 15\n"
                        "}\n\n", pr.ambienLightColor[0], pr.ambienLightColor[1], pr.ambienLightColor[2]);
    
    fprintf(file, "background { color rgb <%f, %f, %f> }\n\n", pr.background[0], pr.background[1], pr.background[2]);
    fprintf(file, "camera {\n       %s\n       angle %f\n", pr.orthographic.c_str(), pr.fov);
    fprintf(file, "       location <%f, %f, %f>\n", pr.cameraLocation[0], pr.cameraLocation[1], pr.cameraLocation[2]);
    fprintf(file, "       sky <%f, %f, %f>\n", pr.sky[0], pr.sky[1], pr.sky[2]);
    fprintf(file, "       up <%f, %f, %f>\n", pr.cameraUp[0], pr.cameraUp[1], pr.cameraUp[2]);
    fprintf(file, "       right <%f, %f, %f>\n", pr.cameraRight[0], pr.cameraRight[1], pr.cameraRight[2]);
    fprintf(file, "       look_at <%f, %f, %f> }\n", pr.cameraLookAt[0], pr.cameraLookAt[1], pr.cameraLookAt[2]);

    fprintf(file, "light_source {\n");
    fprintf(file, "        <%f, %f, %f>\n", pr.lightLocation[0], pr.lightLocation[1], pr.lightLocation[2]);
    fprintf(file, "        color rgb <%f, %f, %f> %s \n", pr.lightColor[0], pr.lightColor[1], pr.lightColor[2], pr.shadowless.c_str());
    //                  "        fade_distance %f\n"
    //                  "        fade_power %f\n"
    fprintf(file,       "        parallel\n");
    fprintf(file, "        point_at <%f, %f, %f>\n}\n", pr.lightPointAt[0], pr.lightPointAt[1], pr.lightPointAt[2]);
    fprintf(file,  "#default {\n        finish {ambient %f diffuse %f specular %f roughness %f metallic %f}\n}\n",
                                    pr.ambient, pr.diffuse, pr.specular, pr.roughness, pr.metallic);


    //fprintf(file,  "union {\n");
    //start with bonds
    for(unsigned int i=0; i < traj.numBonds(); i++){
        Bond &bd = traj.bond(i);
        Atom &at1 = traj.atom(bd.at1);
        Atom &at2 = traj.atom(bd.at2);
        if(at1.hidden || at2.hidden){
        continue;
        }
        glm::vec3 dirVec = at2.coor - at1.coor;
        glm::vec3 mp = at1.coor + dirVec*(bd.length*0.25f); //midpoint between two atomic coordinates

        //first bond
        fprintf(file, "cylinder {\n        <%f, %f, %f>,  <%f, %f, %f>, %f  \n", 
                        at1.coor[0], at1.coor[1], at1.coor[2], mp[0], mp[1], mp[2], bd.cylRadius);

        fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at1.color()[0],at1.color()[1],at1.color()[2],1-at1.color()[3]);
        
        //second bond
        fprintf(file, "cylinder {\n        <%f, %f, %f>,  <%f, %f, %f>, %f  \n", 
                        mp[0], mp[1], mp[2], at2.coor[0], at2.coor[1], at2.coor[2], bd.cylRadius);

        fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at2.color()[0],at2.color()[1],at2.color()[2],1-at2.color()[3]);
    
    }

    //start with atoms
    for(unsigned int i = 0; i < traj.numAtoms(); i++){
        Atom &at = traj.atom(i);
        if(at.hidden){
            continue;
        }
        fprintf(file, "sphere {\n        <%f, %f, %f>, %f  \n", 
                                    at.coor[0], at.coor[1], at.coor[2], at.sphere());

        fprintf(file, "        pigment { rgbt <%f, %f, %f, %f> }\n}\n", at.color()[0],at.color()[1],at.color()[2],1-at.color()[3]);
        
    }

    //fprintf(file, "}\nmerge {\n}");
    fclose(file);
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << frameInd;
    string f = "_tempImg-" + ss.str() + ".png";
    string cmd = "povray -Itemp.pov -O" + f + " +W"+ std::to_string(pr.width) + 
                    " +H" + std::to_string(pr.height) + "+V  +FN +Q9 -P -UD +UL +UV +A +AM2";
    system(cmd.c_str());
  }//end for temp.pov
  //system("rm ./temp.pov");
  string f = fileName;
  string gifcmd = "conver -delay 5 _tempImg*png  " + f + ".gif";
  system(gifcmd.c_str());
  return 0;
}


//// --------------------- functions to read and write files ---------------
