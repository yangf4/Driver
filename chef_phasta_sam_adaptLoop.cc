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

#ifndef WRITE_VTK
#define WRITE_VTK
#endif

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
    const char* fldName = "motion_coords";
    return sam::specifiedIso(m,fldName,fldIdx,fldLimit,szFactor);
  }
  
  apf::Field* getPreSF(apf::Mesh* m, int step) {
    apf::Field* newSz = apf::createFieldOn(m,"preSz",apf::SCALAR);
    apf::Field* coord = m->findField("motion_coords");
    double* vals = new double[coord->countComponents()];
    double dis[444];
    for ( int i =   0; i < 300; i++ )  dis[i] = 0.0;
    dis[300] = 1.1e-6;
    for ( int i = 301; i < 375; i++ )  dis[i] = 1.8e-6;
    for ( int i = 375; i < 424; i++ )  dis[i] = dis[i-1]+0.0142857e-6;
    dis[424] = 2.5e-6;
    for ( int i = 425; i < 444; i++ )  dis[i] = dis[i-1]-0.1052631e-6;
    double cen[] = {0.0, 0.0, 0.0};
    for ( int i = 300; i < step-1;i++ )  cen[0] = cen[0] + dis[i];
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(coord,vtx,0,vals);
      double dist = sqrt((vals[0]-cen[0])*(vals[0]-cen[0]) +
                         (vals[1]-cen[1])*(vals[1]-cen[1]) +
                         (vals[2]-cen[2])*(vals[2]-cen[2]))- 9.5e-6; 
      if ( dist < 0 )
        apf::setScalar(newSz,vtx,0, 1.0e-6);
      else if ( dist < 18e-6 )
        apf::setScalar(newSz,vtx,0, 1.0e-6+dist/2.0);
      else 
        apf::setScalar(newSz,vtx,0, 1e-5);
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
    ctrl.recursivePtn = 0;
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
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(f, vtx, 0, vals);
      for ( int i = 0; i < 3; i++ )  points[i] = vals[i];  
      m->setPoint(vtx, 0, points);
    }
    m->end(itr); 
    delete [] vals;
    return true;  
  }
   
  void writeSequence (apf::Mesh2* m, int step, const char* filename) {
    std::ostringstream oss; 
    oss << filename << step;
    const std::string tmp = oss.str();
#ifdef WRITE_VTK
    apf::writeVtkFiles(tmp.c_str(),m);
#endif
  }

  void writeFirstCoord (apf::Mesh2* m) {
    apf::Vector3 points; 
    apf::MeshIterator* itr = m->begin(0);
    apf::MeshEntity* vtx = m->iterate(itr);
    m->getPoint(vtx, 0, points); 
    std::cout << "First Node Coordinates: " << points << '\n'; 
  }

  void attachSizeField (apf::Field* sf, apf::Mesh2* m) {
    int out_size = 1;
    apf::Field* f = apf::createPackedField(m, "Size Field", out_size);
    double* data = new double[out_size]; 
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while ((e = m->iterate(it))) {
      apf::getComponents(sf,e, 0, data);
      apf::setComponents(f, e, 0, data);
    }
    m->end(it);
  }

  void writePHTfiles (int step, int nstep, int nproc) {
    std::ostringstream oss;
    oss << "solution_" << step << ".pht";
    const std::string tp = oss.str();
    const char* filename = tp.c_str();
    FILE* sFile = fopen (filename, "w");
    fprintf (sFile, "<?xml version=\"1.0\" ?>\n");
    fprintf (sFile, "<PhastaMetaFile number_of_pieces=\"%d\">\n", nproc);
    fprintf (sFile, "  <GeometryFileNamePattern pattern=\"%d-procs_case/geombc.%d.%%d\"\n",nproc,step);
    fprintf (sFile, "                           has_piece_entry=\"1\"\n");
    fprintf (sFile, "                           has_time_entry=\"0\"/>\n");
    fprintf (sFile, "  <FieldFileNamePattern pattern=\"%d-procs_case/restart.%%d.%%d\"\n",nproc);
    fprintf (sFile, "                        has_piece_entry=\"1\"\n");
    fprintf (sFile, "                        has_time_entry=\"1\"/>\n");
    fprintf (sFile, "  <TimeSteps number_of_steps=\"%d\"\n", nstep);
    fprintf (sFile, "             auto_generate_indices=\"1\"\n");
    fprintf (sFile, "             start_index=\"%d\"\n", step+1);
    fprintf (sFile, "             increment_index_by=\"1\"\n");
    fprintf (sFile, "             start_value=\"0.0\"\n");
    fprintf (sFile, "             increment_value_by=\"1.0e-6\">\n");
    fprintf (sFile, "  </TimeSteps>\n");
    fprintf (sFile, "  <Fields number_of_fields=\"5\">\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"pressure\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"1\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"velocity\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"1\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"temperature\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"4\"\n");
    fprintf (sFile, "           number_of_components=\"1\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"motion_coords\"\n");
    fprintf (sFile, "           phasta_field_tag=\"motion_coords\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"mesh_vel\"\n");
    fprintf (sFile, "           phasta_field_tag=\"mesh_vel\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "  </Fields>\n");
    fprintf (sFile, "</PhastaMetaFile>\n");
    fclose (sFile);
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
//  ctrl.openfile_read = openfile_read;
  /* load the model and mesh */
  apf::Mesh2* m = apf::loadMdsMesh(
      ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str());
  chef::preprocess(m,ctrl);
  chef::preprocess(m,ctrl,grs);
  rstream rs = makeRStream();
  /* setup stream reading */
  ctrl.openfile_read = openstream_read;
  ctrl.rs = rs;
  if (!ctrl.writeVizFiles)  ctrl.writeVizFiles = 2; 
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  int loop = 0;
  int seq  = 0;
  writeSequence(m,seq,"test_"); seq++; 
  writePHTfiles (0, 300, 8);
  do {
    /* take the initial mesh as size field */
//    apf::Field* szFld = samSz::isoSize(m);
    step = phasta(inp,grs,rs);
    if ( step >= 300 && step<425 && step%5==0 )
      writePHTfiles (step, 5, 8);
    else if ( step >= 425 )
      writePHTfiles (step, 1, 8);
    ctrl.rs = rs; 
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    overwriteMeshCoord(m);
    if ( (step%5==0) || (step>=425) )
      writeSequence(m,seq,"test_"); seq++; 
//    apf::Field* szFld = getField(m);
    apf::Field* szFld = getPreSF(m, step);
    apf::synchronize(szFld);
    apf::synchronize(m->getCoordinateField());
    assert(szFld);
    if ( (step>=425) || (step>300 && step%5==0))
      chef::adapt(m,szFld);
    if ( (step%5==0) || (step>=425) )
      writeSequence(m,seq,"test_"); seq++; 
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
