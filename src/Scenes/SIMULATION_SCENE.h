#ifndef SIMULATION_SCENE_H
#define SIMULATION_SCENE_H

#ifndef GL_DISABLED
#include "util/DRAW_GL.h"
#endif

#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/STRAND_NET_MESH.h"
#include "Geometry/STRAND_MESH.h"
#include "util/TIMER.h"

namespace HOBAK {

class SIMULATION_SCENE {
public:

  // initialize the scene
  SIMULATION_SCENE() {
    _pauseFrame = -2;
    _exitOnPause = true;
    _arrowCounter = -1;
    _leftArrow = false;
    _rightArrow = false;
    _drawFeature = false;
    _frameNumber = 0;
    _normalizedVertices = false;
    _sceneName = std::string("default");
    _initialA = MATRIX3::Identity();
    _initialTranslation = VECTOR3::Zero();

    _gravity.setZero();

    _movieInterval = -1;
    _autoplay = true;
  };

  virtual ~SIMULATION_SCENE()
  {

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      delete _kinematicShapes[x];
  };

  // TODO: Build the actual scene. You have to implement this!
  virtual bool buildScene() = 0;

  // TODO: Describe the scene build built. You have to do this!
  virtual void printSceneDescription() = 0;

  const VECTOR3& eye() const         { return _eye; };
  const VECTOR3& lookAt() const      { return _lookAt; };
  const VECTOR3& up() const          { return _up; };
  const VECTOR3& worldCenter() const { return _worldCenter; };
  const VECTOR3& gravity() const  { return _gravity; };
  const int& pauseFrame() const   { return _pauseFrame; };
  const int& frameNumber() const  { return _frameNumber; };
  const bool& drawFeature() const { return _drawFeature; };
  const int& arrowCounter() const { return _arrowCounter; };
  const vector<KINEMATIC_SHAPE*>& kinematicShapes() const { return _kinematicShapes; };
  const bool& normalizedVertices() const { return _normalizedVertices; };
  const bool& exitOnPause() const { return _exitOnPause; };

  const MATRIX3& initialA() const           { return _initialA; };
  const VECTOR3& initialTranslation() const { return _initialTranslation; };

  REAL& E()              { return _E; };
  const REAL& E() const  { return _E; };
  //REAL& nu()             { return _nu; };
  //const REAL& nu() const { return _nu; };
  REAL& G()             { return _G; };
  const REAL& G() const { return _G; };

  VECTOR3& eye()         { return _eye; };
  VECTOR3& lookAt()      { return _lookAt; };
  VECTOR3& up()          { return _up; };
  VECTOR3& worldCenter() { return _worldCenter; };
  VECTOR3& gravity()     { return _gravity; };
  int& pauseFrame()      { return _pauseFrame; };
  int& frameNumber()     { return _frameNumber; };
  bool& drawFeature()    { return _drawFeature; };
  int& arrowCounter()    { return _arrowCounter; };
  bool& leftArrow()      { return _leftArrow; };
  bool& rightArrow()     { return _rightArrow; };
  bool& normalizedVertices() { return _normalizedVertices; };
  vector<KINEMATIC_SHAPE*>& kinematicShapes() { return _kinematicShapes; };

  const string sceneName() const { return _sceneName; };
  const string movieName() const { return _sceneName + std::string(".mov"); };
  const string jsonName() const  { return _sceneName + std::string(".json"); };
  const bool autoplay() const { return _autoplay; };


  // simulation loop
  virtual void stepSimulation(const bool verbose = true) {

    _frameNumber++;
  };

  // should we write a movie at this frame?
  const bool writeMovie() const { 
    if (_movieInterval == -1) return false;

    if ((_frameNumber % _movieInterval) == 0)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " Frame number: " << _frameNumber << " movie interval: " << _movieInterval << endl;
      cout << " mod: " << _frameNumber % _movieInterval << endl;
      return true;
    }
    return false;
  };

  virtual bool writeFrameToFile() { return true; };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
  };
#endif

protected:
  // set the positions to previous timestep, in case the user wants to
  // look at that instead of the current step
  void setToPreviousTimestep()
  {
  }

  // restore positions from previous timestep, in case the user just drew
  // the previous timestep, but now we want the state to be consistent
  // when drawing the next frame
  void restoreToCurrentTimestep()
  {
  }

  // scene geometry
  vector<KINEMATIC_SHAPE*> _kinematicShapes;

  // solver and materials

  // simulation parameters
  VECTOR3 _gravity;

  // initial rotation-scale and translation of tet mesh
  MATRIX3 _initialA;
  VECTOR3 _initialTranslation;

  // drawing parameters
  int _pauseFrame;

  // counter that can be incremented and decremented by the arrow keys
  int _arrowCounter;

  // bools that can be toggled by the arrow keys
  bool _leftArrow;
  bool _rightArrow;

  // flag for whether or not to draw some user-specific feature
  bool _drawFeature;

  // what frame are we on?
  int _frameNumber;

  // should we save this frame?
  bool _writeToFile = false;


  // did we normalize the vertices when we read them in?
  bool _normalizedVertices;

  // camera parameters
  VECTOR3 _eye;
  VECTOR3 _lookAt;
  VECTOR3 _up;
  VECTOR3 _worldCenter;

  // scene name, used to determine the JSON and MOV filenames
  std::string _sceneName;

  // The path to save the files to
  std::string _filePath;

  // should the scene automatically start simulating?
  bool _autoplay;

  // should we exit when we hit the pause frame?
  bool _exitOnPause;
 
  // what interval should we write out a movie at?
  int _movieInterval;
 
  // Young's modulus and Poisson's ratio. There aren't well-integrated into
  // the volumetric files yet
  REAL _E;
  //REAL _nu;     // let's prefer shear modulus to Poisson's ratio
  REAL _G;
};

}

#endif
