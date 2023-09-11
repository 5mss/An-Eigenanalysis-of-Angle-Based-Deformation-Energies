#include <iostream>
#include <random>
#include "util/FFMPEG_MOVIE.h"


// #include "Scenes/Strand/STRAND_SCENE.h"
// #include "Scenes/Strand/STRAND_COLLISION_SCENE.h"
// #include "Scenes/Strand/STRAIGHT_STRAND_SCENE.h"
// #include "Scenes/Strand/MULTI_STRAND_SCENE.h"
// #include "Scenes/Strand/STRAND_L_SCENE.h"
// #include "Scenes/Strand/STRAND_DVT_SCENE.h"
// #include "Scenes/Strand/STRAND_LINE_SCENE.h"
// #include "Scenes/Strand/STRAND_WAVE_SCENE.h"

#include "Scenes/Strand/NET_SCENE_CORNER.h"
#include "util/DRAW_GL.h"

#include <snapshot.h>

using namespace HOBAK;
using namespace std;

GLVU glvu;
FFMPEG_MOVIE movie;

bool animate = false;
//bool animate = true;
bool singleStep = false;
//bool singleStep = true;
int timestep = -1;

// animation should pause on this frame;
int pauseFrame = -1;

// scene being simulated
SIMULATION_SCENE* scene = NULL;
void finiteDiffForces();
void finiteDiffTwistHessian();

///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    scene->drawScene();
    if (animate)
      movie.addFrameGL();
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate) 
  {
    scene->stepSimulation();
    //recordFrameToJSON(*scene);
    if (singleStep) 
    {
      animate = false;
      singleStep = false;
    }
    //finiteDiffForces();
    //finiteDiffTwistHessian();
    timestep++;
  }
  
  if (timestep == scene->pauseFrame() && animate)
  {
    cout << " Hit pause frame specific in scene file " << endl;
    TIMER::printTimings();
    // movie.streamWriteMovie(scene->movieName().c_str());
    //writeSceneJSON(scene->jsonName().c_str());
    cout << " Wrote out movie: " << scene->movieName().c_str() << endl;
    animate = false;
    exit(0);
  }

  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 27:  // escape key
  case 'Q':
    TIMER::printTimings();
    exit(0);
    break;
  case 'q':
    TIMER::printTimings();
    // movie.streamWriteMovie(scene->movieName().c_str());
    //writeSceneJSON(scene->jsonName().c_str());
    cout << " Wrote out movie: " << scene->movieName().c_str() << endl;
    exit(0);
    break;
  case 'a':
    animate = !animate;
    break;
  case 'd':
    scene->drawFeature() = !scene->drawFeature();
    break;
  case 'm':
    TIMER::printTimings();
    movie.streamWriteMovie(scene->movieName().c_str());
    cout << " Wrote out movie: " << scene->movieName().c_str() << endl;
    return;
    break;
  case ' ':
    animate = true;
    singleStep = true;
    TIMER::printTimings();
    break;
  case 'v': {
      Camera *camera = glvu.GetCurrentCam();
      glvuVec3f eye;
      glvuVec3f ref;
      glvuVec3f up;
      camera->GetLookAtParams(&eye, &ref, &up);
      cout << "    _eye    = VECTOR3(" << eye[0] << ", " << eye[1] << ", " << eye[2] << ");" << endl;
      cout << "    _lookAt = VECTOR3(" << ref[0] << ", " << ref[1] << ", " << ref[2] << ");" << endl;
      cout << "    _up     = VECTOR3(" << up[0] << ", " << up[1] << ", " << up[2] << ");" << endl;
    } 
    break;
  default:
    break;
  }

  glvu.Keyboard(key, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key) {
  case GLUT_KEY_LEFT:
    scene->leftArrow() = true;
    break;
  case GLUT_KEY_RIGHT:
    scene->rightArrow() = true;
    break;
  case GLUT_KEY_UP:
    scene->arrowCounter()++;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_DOWN:
    scene->arrowCounter()--;
    cout << "Arrow counter: " << scene->arrowCounter() << endl;
    break;
  case GLUT_KEY_HOME:
    break;
  case GLUT_KEY_END:
    break;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y) 
{ 
  glvu.Mouse(button, state, x, y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y) 
{ 
  glvu.Motion(x, y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  const string title = string("Simulating scene: ") + scene->sceneName();
  char buffer[title.length() + 1];
  strcpy(buffer, title.c_str());
  glvu.Init(buffer, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE, 0, 0, 800, 800);
  //glvu.Init(buffer, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE, 0, 0, 1000, 1000);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  GLfloat lightZeroPosition[] = { 10.0, -4.0, -10.0, 1.0 };
  GLfloat lightZeroColor[] = { 0.2, 0.2, 0.1, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);

  GLfloat lightPosition1[] = { -10.0, 4.0, 10.0, 1.0 };
  GLfloat lightColor1[] = { 0.2, 0.2, 0.1, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);

  GLfloat lightPosition2[] = { 0.0, 0, -10.0, 1.0 };
  GLfloat lightColor2[] = { 0.2, 0.2, 0.1, 1.0 };
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
  glLightfv(GL_LIGHT2, GL_POSITION, lightPosition2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, lightColor2);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  glEnable(GL_COLOR_MATERIAL);
  glShadeModel(GL_SMOOTH);
  glClearColor(1, 1, 1, 0);
  //glClearColor(0, 0, 0, 0);

  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  glEnable(GL_MULTISAMPLE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glvuVec3f ModelMin(-10, -10, -10), ModelMax(10, 10, 10);
  glvuVec3f Eye, LookAtCntr, Up, WorldCenter;

  const VECTOR3& eye    = scene->eye();
  const VECTOR3& lookAt = scene->lookAt();
  const VECTOR3& up     = scene->up();
  const VECTOR3& worldCenter = scene->worldCenter();

  for (unsigned int x = 0; x < 3; x++)
  {
    Eye[x] = eye[x];
    LookAtCntr[x] = lookAt[x];
    Up[x] = up[x];
    WorldCenter[x] = worldCenter[x];
  }

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 1000.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);
  glvu.SetWorldCenter(WorldCenter);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  //scene = new STRAND_SCENE();
  // scene = new STRAND_L_SCENE();
  //scene = new STRAND_DVT_SCENE();
  // scene = new STRAIGHT_STRAND_SCENE();
  //scene = new MULTI_STRAND_SCENE();
  //scene = new STRAND_LINE_SCENE();
  //scene = new STRAND_WAVE_SCENE();
  //scene = new STRAND_COLLISION_SCENE();
  
  // scene = new STRAND_GRID_SCENE();
  // scene = new NET_SCENE();
  scene = new NET_SCENE_CORNER();
  // scene = new NET_CAN();
  // scene = new NET_SCENE_DROP();
  //scene = new HAIRBALL();
  //scene = new CURLY_HAIRBALL();
  // scene = new BEND_TWIST_VALIDATION();

  //scene = new SINGLE_STRAND_INSTABILITY();
  // scene = new SANDBOX();

  //scene = new COMPRESSION_INSTABILITY();
  //scene = new HAIRBALL_UNSTABLE();

  // make sure to print out what we should expect from this scene
  scene->printSceneDescription();

  bool success = scene->buildScene();

  //finiteDiffForces();
  //finiteDiffTwistHessian();

  // prepare scene for JSON output
  //initializeSceneJSON(*scene);

  if (!success)
  {
    cout << " Failed to build scene-> Exiting ..." << endl;
    return 1;
  }

  glutInit(&argc, argv);
  glvuWindow();

  return 0;
}
