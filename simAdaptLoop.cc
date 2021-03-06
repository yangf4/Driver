#include "chefPhasta.h"
#include "samSz.h"
#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <SimField.h>
#include <MeshSimAdapt.h>
#include <apfSIM.h>
#include <SimPartitionedMesh.h>

#define debug

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  apf::Field* getField(apf::Mesh* m) {
    /* if the value of the fldIdx'th index from the fldName
     * field is greater than fldLimit then multiply the current
     * isotropic mesh size at the vertex by szFactor */
    const unsigned fldIdx = 1;
    const double fldLimit = 1.0;
    const double szFactor = 0.5;
    const char* fldName = "solution";
    return sam::specifiedIso(m,fldName,fldIdx,fldLimit,szFactor);
  }
  
  apf::Field* getConstSF(apf::Mesh* m, double factor) {
    apf::Field* newSz = apf::createFieldOn(m,"constSz",apf::SCALAR);
    double h = 0.0; 
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      h = factor;
      apf::setScalar(newSz,vtx,0,h);
    }
    m->end(itr);
    return newSz;
  }
  
  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
  }

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
    }
  }
  
  bool overwriteMeshCoord(apf::Mesh2* m) { 
    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[f->countComponents()];
    assert(f->countComponents() == 3);
    apf::MeshEntity* vtx;
    apf::Vector3 points; 
    apf::MeshIterator* itr = m->begin(0);
    int debug = 0; 
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(f, vtx, 0, vals);
#ifdef debug
      m->getPoint(vtx, 0, points); 
      double err = (points[0] - vals[0])*(points[0] - vals[0])
                 + (points[1] - vals[1])*(points[1] - vals[1])
                 + (points[2] - vals[2])*(points[2] - vals[2]); 
      if ( err > 2.0 ) fprintf(stderr, "Node %d bigger than tolerance\n", debug);
//      std::cout << "node: " << debug << " ;Coordinates: " << points; 
//      std::cout << " ;Com: (" << vals[0] << ", " << vals[1] << ", " << vals[2] << ")" << '\n';
#endif
      for ( int i = 0; i < 3; i++ )  points[i] = vals[i];  
      m->setPoint(vtx, 0, points);
      debug++;
    }
    m->end(itr); 
    fprintf(stderr, "total number of vertex: %d\n", debug);
    delete [] vals;
    return true;  
  }
   
  void writeSequence (apf::Mesh2* m, int step) {
    const std::string filename = "test_";
    std::ostringstream oss; 
    oss << filename << step;
    const std::string tmp = oss.str();
    apf::writeVtkFiles(tmp.c_str(),m);
  }

#ifdef debug
  void writeFirstCoord (apf::Mesh2* m) {
    apf::Vector3 points; 
    apf::MeshIterator* itr = m->begin(0);
    apf::MeshEntity* vtx = m->iterate(itr);
    m->getPoint(vtx, 0, points); 
    std::cout << "First Node Coordinates: " << points << '\n'; 
  }
#endif

  void simAdapt (apf::Mesh2* m, apf::Field* szFld) {
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh sim_pm = apf_msim->getMesh();
    pMSAdapt adapter = MSA_new(sim_pm, 1);
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while ((v = m->iterate(it))) {
      double size = apf::getScalar(szFld, v, 0);
      MSA_setVertexSize(adapter, (pVertex) v, size);
    }
    m->end(it);
    apf::Field* res_fld;
    if (res_fld = m->findField("solution"))
//      apf::Field* res_fld = apf_m->findField("solution");
//      pField sim_res_fld = apf::getSIMField(res_fld);
      PList_append(sim_fld_lst, apf::getSIMField(res_fld));
    if (res_fld = m->findField("time derivative of solution"))
      PList_append(sim_fld_lst, apf::getSIMField(res_fld));
    if (res_fld = m->findField("mesh_vel"))
      PList_append(sim_fld_lst, apf::getSIMField(res_fld));
    MSA_setMapFields(adapter, sim_fld_lst);
    PList_delete(sim_fld_lst);
    pProgress progress = Progress_new();
    MSA_adapt(adapter, progress);
    Progress_delete(progress);
    MSA_delete(adapter);
  }  
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("samAdaptLoop.inp");
  /* setup file reading */
  ctrl.openfile_read = openfile_read;
  /* load the model and mesh */
  apf::Mesh2* m = apf::loadMdsMesh(
      ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str());
  chef::preprocess(m,ctrl,grs);
  rstream rs = makeRStream();
  /* setup stream reading */
  ctrl.openfile_read = openstream_read;
  ctrl.rs = rs;
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  int loop = 0;
  int seq  = 0;
  writeSequence(m,seq); seq++; 
  do {
    /* take the initial mesh as size field */
    apf::Field* szFld = samSz::isoSize(m);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs; 
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    overwriteMeshCoord(m);
    writeSequence(m,seq); seq++; 
    apf::synchronize(szFld);
    apf::synchronize(m->getCoordinateField());
    assert(szFld);
//    chef::adapt(m,szFld);
    simAdapt(m,szFld);
    writeSequence(m,seq); seq++; 
    apf::destroyField(szFld);
    chef::preprocess(m,ctrl,grs);
    clearRStream(rs);
    loop++; 
  } while( loop < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
