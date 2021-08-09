// foo.cpp
// g++ -o foo.so foo.cpp -std=c++11 -fPIC -shared -Wall -Wextra `python2.7-config --includes --libs` -lboost_python

#include <boost/python.hpp>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <globals.h>
#include <rw.h>



BOOST_PYTHON_MODULE(pytigre) {
  using namespace boost::python;

  props::initProps();

  class_<std::array<float, 3>>("arrFloat3")
    .def("__iter__", iterator<std::array<float, 3> >())
;
  class_<std::array<int, 3>>("arrInt3")
    .def("__iter__", iterator<std::array<int, 3> >())
;
  class_<std::array<float, 6>>("arrFloat6")
    .def("__iter__", iterator<std::array<float, 6> >())
;
  class_<std::vector<unsigned int>>("vecInt")
    .def("__iter__", iterator<std::vector<unsigned int> >())
;

  class_<std::array<short, 10>>("bondIndices")
    .def("__iter__", iterator<std::array<short, 10>>())
;
  
  class_<glm::vec3>("vec3", init<>())
    .def_readwrite("x", &glm::vec3::x)
    .def_readwrite("y", &glm::vec3::y)
    .def_readwrite("z", &glm::vec3::z)
    //.def("__iter__", iterator<std::vector<unsigned int> >())
;
  
  class_<glm::mat3>("mat3", init<>())

;

  class_<props::LatticeParams>("LatticeParams", 
            init<float, float, float, float, float, float>(args("a", "b", "c", "alpha", "beta", "gamma")))
    .def(init<>())
    .def_readwrite("a", &props::LatticeParams::a)
    .def_readwrite("b", &props::LatticeParams::b)
    .def_readwrite("c", &props::LatticeParams::c)
    .def_readwrite("alpha", &props::LatticeParams::alpha)
    .def_readwrite("beta", &props::LatticeParams::beta)
    .def_readwrite("gamma", &props::LatticeParams::gamma)
;

  class_<rw::FrameChooser>("FrameChooser", init<int, int, int>(args("first","last", "stride")))
    .def(init<>())
    .def_readwrite("first", &rw::FrameChooser::first)
    .def_readwrite("last", &rw::FrameChooser::last)
    .def_readwrite("stride", &rw::FrameChooser::stride)
  ;

  class_<Bond>("Bond", init<>())
    .add_property("at1", &Bond::at1)
    .add_property("at2", &Bond::at2)
    .add_property("length", &Bond::length)
  ;

  class_<Atom>("Atom", init<>())
	  .def_readwrite("coor", &Atom::coor)
    .def_readwrite("atomicNum", &Atom::atomicNum)
    .add_property("symbol", &Atom::symbol, &Atom::setSymbol)
    .def_readwrite("residue", &Atom::residue)
    .add_property("atomType", &Atom::getAtomType, &Atom::setAtomType)
    .def_readwrite("charge", &Atom::charge)
    .def_readwrite("chemShift", &Atom::chemShift)
    .def_readwrite("asymUnit", &Atom::asymUnit)
    .def_readwrite("selected", &Atom::selected)
    .def_readwrite("locked", &Atom::locked)
    .def_readwrite("hidden", &Atom::hidden)
    .def_readwrite("bondInds", &Atom::bondInds)
//	  .def_readwrite("color", &Atom::color)
	  .def_readwrite("fixed", &Atom::fixed)
    .add_property("sphere", &Atom::sphere, &Atom::setSphere)
    .add_property("radius", &Atom::radius, &Atom::setRadius)
    .add_property("numBonds", &Atom::numBonds)
	  
  ;

  class_<Frame>("Frame",  init<>())
    .def <Atom& (Frame::*)(const char*, const float, const float, const float)> ("addAtom", &Frame::addAtom,boost::python::return_internal_reference<>())
    .def <Atom& (Frame::*)(short, glm::vec3&)> ("addAtom", &Frame::addAtom,boost::python::return_internal_reference<>())

    .def("atom", &Frame::atom, boost::python::return_internal_reference<>())
    .def("numAtoms", &Frame::numAtoms)
    
    .def_readwrite("energy", &Frame::energy)
    .def_readwrite("time", &Frame::time)
    .def_readwrite("step", &Frame::step)
    .def_readwrite("comment", &Frame::comment)
    .def_readwrite("centroid", &Frame::centroid)
    .def_readwrite("cell", &Frame::cell)

    .def <void (Frame::*)(float, float, float, float, float, float, float, float, float)>("setUnitCell", &Frame::setUnitCell)
    
;

  class_<Trajectory>("Trajectory", init<>())
    .def ("addBondedAtomsGeometrically", &Trajectory::addBondedAtomsGeometrically)
    .def <void (Trajectory::*)()> ("addFrame", &Trajectory::addFrame)
    .def <unsigned int (Trajectory::*)(Frame&)> ("addFrame", &Trajectory::addFrame)
    .def("atom", &Trajectory::atom, boost::python::return_internal_reference<>())
    .def <float (Trajectory::*)(const unsigned int, const unsigned int, const unsigned int)> ("angle", &Trajectory::angle)
    .def("bond", &Trajectory::bond, boost::python::return_internal_reference<>())
    .def<void (Trajectory::*)(glm::vec3&)>("bringInsideCell", &Trajectory::bringInsideCell)
    .def <glm::vec3 (Trajectory::*)()> ("centroid", &Trajectory::centroid)
    .def("cartesianToFractional", &Trajectory::cartesianToFractional)
    .def("fractionalToCartesian", &Trajectory::fractionalToCartesian)
    .def("currentFrame", &Trajectory::currentFrame, boost::python::return_internal_reference<>())
    .add_property("currentFrameIndex", &Trajectory::currentFrameIndex)
    .def <float (Trajectory::*)(const unsigned int, const unsigned int)> ("distance", &Trajectory::distance)
    .def("generateBonds", &Trajectory::generateBonds)
    .def("getAtomNeighborsIndices", &Trajectory::getAtomNeighborsIndices)
    .def("frame", &Trajectory::frame, boost::python::return_internal_reference<>())
    .add_property("name", &Trajectory::name)
    .def("nextFrame", &Trajectory::nextFrame)
    .add_property("numAtoms", &Trajectory::numAtoms)
    .add_property("numBonds", &Trajectory::numBonds)
    .add_property("numFrames", &Trajectory::numFrames)
    .def<void (Trajectory::*)(const float, const float, const float, const float)>("rotateCoors", &Trajectory::rotateCoors)
    .def <bool (Trajectory::*)(const int)> ("setCurrentFrameIndex", &Trajectory::setCurrentFrameIndex)
    .def <bool (Trajectory::*)(Atom&, Atom&, float)> ("setDistance", &Trajectory::setDistance)
    .def <bool (Trajectory::*)(unsigned int i, unsigned int j, float)> ("setDistance", &Trajectory::setDistance)
    .def <void (Trajectory::*)(props::LatticeParams&)>("setUnitCell", &Trajectory::setUnitCellVectorsComponents)
    .def <void (Trajectory::*)(float, float, float, float, float, float, float, float, float)>("setUnitCell", &Trajectory::setUnitCell)
    .def ("toCentroid", &Trajectory::toCentroid)
    .def <void (Trajectory::*)(const glm::vec3 &)>("translateCoors", &Trajectory::translateCoors)
    ;
  
  class_<std::vector<Frame>>("vecFrame")
    .def("__iter__", iterator<std::vector<Frame> >())
;
	  
  // read write ------------------------------------------------
  boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadXYZ", &rw::loadXYZ);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadOUTCAR", &loadOUTCAR);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadCIF", &loadCIF);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadPDB", &loadPDB);
  boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadPOSCAR", &rw::loadPOSCAR);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadXDATCAR", &loadXDATCAR);
  boost::python::def <unsigned int  (const char*, Trajectory&)> ("loadFile", &rw::loadFile);
  boost::python::def <unsigned int  (const char*, Trajectory&, rw::FrameChooser&)> ("loadFile", &rw::loadFile);
  
  boost::python::def <unsigned int  (const char*, Trajectory&, const rw::FrameChooser&)> ("writeFile", &rw::writeFile);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("writePOSCAR", &rw::writePOSCAR);
  // boost::python::def <unsigned int  (const char*, Trajectory&)> ("writeGULP_INPUT", &rw::writeGULP_INPUT);
  // boost::python::def <unsigned int  (const char*, Structure&)> ("writeXYZ", &writeXYZ);



  // math  -------------------------------------------------
 // boost::python::def<float (std::array<float, 3>&, std::array<float, 3>&)>("calcDistance", &calcDistance);
  //boost::python::def<float (const std::array<float, 3>&, const std::array<float, 3>&, const std::array<float, 3>&)>("calcAngle", &calcAngle);
  boost::python::def<float (const Atom&, const Atom&)>("distance", &math::distance);
  boost::python::def<float (const glm::vec3&, const glm::vec3&)>("distance", glm::distance);
  boost::python::def <bool (Atom&, Atom&, float)> ("setDistance", &math::setDistance);


 // globalss ---------------------------------------------
  boost::python::def<short (const char*)> ("atomicNumber", &props::atomicNumber);
 // boost::python::def<void (char*)> ("setInitProperties", &setInitProperties);


}
