#ifndef SHELL_SCENE_V_H
#define SHELL_SCENE_V_H

#include "SIMULATION_SCENE.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/TRIANGLE_MESH_FASTER.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/QUADRATIC_F_BENDING.h"
#include <filesystem>

namespace HOBAK {

#ifndef USING_HOMEBREW_GCC
namespace fs = std::__fs::filesystem;
#else
namespace fs = std::filesystem;
#endif

class SHELL_SCENE_V : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Shell V bending scene" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {

    // V shape parameters
    n0 = 31; n1 = 21; vTheta = M_PI / 4.0;

    //_triangleMeshFilename = string("../data/objs/curtain_10x10.obj");
    // _triangleMeshFilename = string("../data/objs/curtain_25x25.obj");
    // //_triangleMeshFilename = string("../data/objs/curtain_50x50.obj");
    // vector<VECTOR3> vertices;
    // vector<VECTOR3I> triangles;
    // bool success = TRIANGLE_MESH::readObjFile(_triangleMeshFilename, vertices, triangles);
    // if (!success)
    // {
    //   cout << " Failed to open file " << _triangleMeshFilename.c_str() << endl;
    //   return false;
    // }
    // _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

    buildVCloth1();
    _triangleMesh->setCollisionEps(spacing / 8.0);
    //_gravity = VECTOR3(0, 0, -1.0);
    _gravity = VECTOR3(0, -1.0, 0.0);

    const REAL stretchingStiffness = 8.0;
    const REAL bendingStiffness = 1e5;
    _strechingEnergy = new SHELL::ARAP(stretchingStiffness, 0.0);
    _bendingEnergy = new SHELL::QUADRATIC_F_BENDING(bendingStiffness);


    // this will determine the MOV and JSON filenames
    string method;
    // string method = "ours";
    if(USING_BRUTE_FORCE_CLAMP)
      method = "bruteCholesky";
    if(USING_RANK4_CORRECTION)
      method = "rank4Cholesky";
    if(USING_GAUSS_NEWTON_BENDING)
      method = "GNCholesky";
    if(USING_UNFILTERED)
      method = "unFIlCholesky";
    char buffer[100];
    sprintf(buffer, "_S%.1f_B%.1f", stretchingStiffness, bendingStiffness);
    _sceneName = "shell_bending_V_" + method + buffer;

    // build the time integrator
    _clothSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy);
    // _clothSolver -> edgeEdgeSelfCollisionsOn() = false;
    // _clothSolver -> vertexFaceSelfCollisionsOn() = false;
    // _clothSolver -> setRayeligh(0.001, 0.0);
    // _clothSolver->setDt(1.0 / 60.0);

    // cube on top
    // VECTOR3 center0(0.0, 0.0, 0.95);
    //_kinematicShapes.push_back(new CUBE(center0, 1.0));
    //_clothSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    //_clothSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    const REAL cubeScale = 8;
    VECTOR3 center0(-0.5 * cubeScale, -0.4 * cubeScale, 0.4 * cubeScale);
    // VECTOR3 center0(-0.5 * cubeScale, 0.4 * cubeScale, 0.4 * cubeScale);
    //VECTOR3 center1(0.0, -1.0, -0.01);
    _kinematicShapes.push_back(new CUBE(center0, cubeScale));
    // _kinematicShapes.push_back(new SPHERE(center1, 0.25));
    _clothSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _clothSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    // _clothSolver->addKinematicCollisionObject(_kinematicShapes.back());

    //_eye    = VECTOR3(1.24999, 2.07096, 0.502227);
    //_lookAt = VECTOR3(0.777846, 1.20965, 0.314523);
    //_up     = VECTOR3(-0.0859995, -0.166908, 0.982213);

    // build
    // _eye    = VECTOR3(0.5, 0.8, 6.0);
    // _lookAt = VECTOR3(0.3, 0.0, 0.0);
    // _up     = VECTOR3(0.0, 1.0, 0.0);

    // build 1
    _eye    = VECTOR3(5.0, 3.0, 1.0);
    _lookAt = VECTOR3(0.3, 0.0, 0.0);
    _up     = VECTOR3(0.0, 1.0, 0.0);

    _worldCenter = _triangleMesh->getRestTranslation();
    //_worldCenter = VECTOR3(0, 0, 0);
    _pauseFrame = 400;

    // try scaling see if F follows
    //for (unsigned int x = 0; x < _triangleMesh->vertices().size(); x++)
    //  _triangleMesh->vertex(x)[2] *= 2.0;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {

    _clothSolver->externalForces().setZero();
    _clothSolver->addGravity(_gravity);
    _clothSolver->solve(verbose);


    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.01, 0.0);
    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.02, 0.0);
    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    drawTriangleMesh(*_triangleMesh, true);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    drawKinematicConstraints(_triangleMesh, _clothSolver);
    //drawPlaneConstraints(_triangleMesh, _clothSolver);

    drawVertex(*_triangleMesh, _arrowCounter);
    /*
    SIMULATION_SCENE::drawScene();

    if (_drawFeature)
    {
      glDisable(GL_DEPTH_TEST);
      glPointSize(10.0);
      drawVertexFacePairs(*_tetMesh, _arrowCounter);
      //drawVertexFacePairs(*_tetMesh);
    }
    */
  };

  void buildVCloth(){
    // the origin is at one corner.
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool alter = true, alterRow = false;

    spacing = 1.0/(n1 - 1.0);
    // first part
    const REAL zStep0 = spacing;
    const REAL xStep0 = zStep0;
    VECTOR3 vert(0.0, 0.0, 0.0);
    int vId = 0; // next to add
    //first row
    for(int i = 0; i < n0; i++) {
      vertices.push_back(vert);
      vert(2) += zStep0;
    }
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(0) += xStep0; vert(2) = 0.0;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(2) += zStep0;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }

    // bended part
    const REAL zStep1 = zStep0;
    const REAL xStep1 = -xStep0 * cos(vTheta);
    const REAL yStep1 = xStep0 * sin(vTheta);
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(0) += xStep1; vert(2) = 0.0; vert(1) += yStep1;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(2) += zStep1;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }
    // cout<<"vertices: "<<endl;
    // for(int i = 0; i < vertices.size(); i++)
    //   cout<<vertices[i]<<endl;
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

  }

  void buildVCloth1(){
    // the origin is at one corner.
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool alter = true, alterRow = false;

    spacing = 1.0/(n1 - 1.0);
    // first part
    const REAL xStep0 = spacing;
    const REAL zStep0 = xStep0 * sin(0.5 * vTheta);
    const REAL yStep0 = - xStep0 * cos(0.5 * vTheta);
    VECTOR3 vert(0.0, 0.0, 0.0);
    int vId = 0; // next to add
    //first row
    for(int i = 0; i < n0; i++) {
      vertices.push_back(vert);
      vert(0) += xStep0;
    }
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(2) += zStep0; vert(1) += yStep0; vert(0) = 0.0;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(0) += xStep0;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }

    // bended part
    const REAL xStep1 = xStep0;
    const REAL zStep1 = zStep0;
    const REAL yStep1 = -yStep0;
    for(int j = 0; j< (n1 - 1)/2; j++){
      vert(2) += zStep1; vert(0) = 0.0; vert(1) += yStep1;
      vertices.push_back(vert);
      // add row
      alter = alterRow;
      for(int i = 0; i < n0 - 1; i++) {
        vert(0) += xStep1;
        vId = (int)vertices.size();
        vertices.push_back(vert);
        if(alter){
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0));
          triangles.push_back(VECTOR3I(vId - 1, vId - 1 - n0, vId - n0));
        }
        else{
          triangles.push_back(VECTOR3I(vId, vId - 1, vId - n0 - 1));
          triangles.push_back(VECTOR3I(vId, vId - 1 - n0, vId - n0));
        }
        alter = !alter;
      }
      alterRow = !alterRow;
    }
    // cout<<"vertices: "<<endl;
    // for(int i = 0; i < vertices.size(); i++)
    //   cout<<vertices[i]<<endl;
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

  }
#endif

protected:
  TRIANGLE_MESH* _triangleMesh;
  SHELL::TIMESTEPPER* _clothSolver;
  SHELL::STRETCHING* _strechingEnergy;
  SHELL::BENDING* _bendingEnergy;

  string _shellFilePath;

  // V shape parameters
  int n0, n1; // number of vertices along the edge. n1 should be odd
  REAL vTheta; // angle of V
  REAL spacing;

  string _triangleMeshFilename;
};

}

#endif
