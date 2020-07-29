/*
 * class: I3Juliet.C
 *
 * Version $Id: I3Juliet.cxx,v 1.1.1.1 2006/03/11 04:31:23 syoshida Exp $
 *
 * Date 07 Nov 2003
 *
 * (c) 2003 IceCube Collaboration
 */

// class header
#include "I3Juliet.h"

#include <iostream>
#include <iomanip>
#include <string>

const double energy = 1e10;
const int flavorID  = 2; //[0,1,2:e,mu,tau]
const int doubletID = 1; //[0,1:neutrino,charged]
const int matID     = 0; //[0,1:ice,rock]
const int doCC      = 1;
const int doNC      = 1; 
const int doMuBrem  = 1;
const int doTauBrem = 1;
const int doMuKnock = 1;
const int doTauKnock= 1;
const int doMu2e    = 1;
const int doTau2e   = 1;
const int doMu2mu   = 1;
const int doTau2mu  = 1;
const int doMu2tau  = 1;
const int doTau2tau = 1;
const int doMuPN    = 1;
const int doTauPN   = 1;
const int doMuDecay = 1;
const int doTauDecay= 1;
const int startlocationID = 1; //[1=from earth entrance: 2=800m away from IceCube origin]
const long random_seed = 1337;

using namespace std;

// constructor and destructor

//===================================================================
//* constructor -----------------------------------------------------
I3Juliet::I3Juliet() 
{
  //---------------------------------------------
  // clear MethodIDs...
  //---------------------------------------------

  fJULIeT_getStartLocationAlongTheAxisID = 0;
  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID = 0;
  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID = 0;
  fJULIeT_wherePrimaryParticleEndsInIceCubeCoordinateID = 0;
  fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID = 0;
  fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID = 0;
  fJULIeT_definePropagatingParticleID = 0;
  fJULIeT_definePropagationGeometryID = 0;
  fJULIeT_definePropagationGeometry_deg_ID = 0;
  fJULIeT_configurePropagationGeometryID = 0;
  fJULIeT_runSingleEventID = 0;           // Method ID runSingleEvent()
  fJULIeT_getParticleIteratorID = 0;      // Method ID getParticleIterator()
  fJULIeT_getTrackParticleIteratorID = 0; // Method ID getParticleIterator()
  fJULIeT_getLocationIce3IteratorID = 0;      // Method ID getLocationIce3Iterator()
  fJULIeT_getTrackLocationIce3IteratorID = 0; // Method ID getLocationIce3Iterator()
  fJULIeT_propParticleID = 0;             // Field ID propParticle in JULIeT
  fJULIeT_primaryEnergyID = 0;            // Field ID primaryEnergy in JULIeT

  fIterator = 0;                    // iterator class (java/util/ListIterator)

  fParticle = 0;                    // the JULIeT Particle class
  fParticle_getEnergyID = 0;        // Method ID getEnergy()
  fParticle_getFlavorID = 0;        // Method ID getFlavor()
  fParticle_getDoubletID = 0;       // Method ID getDoublet()
  fParticle_particleNameID = 0;     // Method ID particleName()

  fJ3Vector = 0;                    // The J3Vector class
  fJ3Vector_getXID = 0;             // Method ID getX()
  fJ3Vector_getYID = 0;             // Method ID getY()
  fJ3Vector_getZID = 0;             // Method ID getZ()
  fJ3UnitVector = 0;                // The J3UnitVector class
  fJ3UnitVector_getXID = 0;         // Method ID getX()
  fJ3UnitVector_getYID = 0;         // Method ID getY()
  fJ3UnitVector_getZID = 0;         // Method ID getZ()
  fJ3Line = 0;                      // The J3Line class

  fDouble_doubleValueID = 0;        // Method ID Double.doubleValue()

}

//===================================================================
//* destructor ------------------------------------------------------
I3Juliet::~I3Juliet()
{ 

}

//===================================================================
//* Configure -------------------------------------------------------
void I3Juliet::Configure()
{
  // generate JavaVM, assign methodIDs, and generate
  // JULIeT object
  CreateJVM();
  GetMethodIDs();
  GenerateJULIeT();

}

//===================================================================
//* Reconfigure -----------------------------------------------------
void I3Juliet::Reconfigure()
{
  // disallow 'Reconfigure' because it makes things more simple.
  cerr << "Reconfigure:Cannot reconfigure this module" << endl;
  exit(1);
}

//===================================================================
//* create JVM  -----------------------------------------------------
void I3Juliet::CreateJVM()
{
  JavaVMInitArgs vm_args;
  JavaVMOption options[4];

  // IMPORTANT: specify vm_args version # if you use JDK1.4 and beyond
  vm_args.version = JNI_VERSION_1_4;
  options[0].optionString = "-Djava.class.path=../../../../classes";
  options[1].optionString = "-Xms256m";
  options[2].optionString = "-Xmx512m";
  options[3].optionString = "-verbose:jni";
  vm_args.options = options;
  //vm_args.nOptions = 4;
  vm_args.nOptions = 3;

  // Create the Java VM 
  jint res = JNI_CreateJavaVM(&fJvmService,(void **)&fJNIEnv,&vm_args);
  if(res < 0) {
    cerr << "Cannot create Java VM" << endl;
    exit(1);
  }

}

//===================================================================
//* Generate JULIeT object: run constructor -------------------------
void I3Juliet::GenerateJULIeT(bool isInteractive)
{

  // delete the former object if exists
  if(fJULIeTObj != 0)  fJNIEnv->DeleteLocalRef(fJULIeTObj);


  if (isInteractive) {

     fJULIeTObj = 
        fJNIEnv-> NewObject(fJULIeT, fJNIEnv->GetMethodID(fJULIeT, "<init>", "()V"));

  } else {

     fJULIeTObj = 
        fJNIEnv->NewObject(fJULIeT, fJNIEnv->GetMethodID(fJULIeT, "<init>", 
                       "(IIDIIIIIIIIIIIIIIIIIIL)V"), 
                       flavorID, doubletID, energy, matID,
                       doCC, doNC, doMuBrem, doTauBrem, doMuKnock, doTauKnock,
                       doMu2e, doTau2e, doMu2mu, doTau2mu,
                       doMu2tau, doTau2tau, doMuPN, doTauPN,
                       doMuDecay, doTauDecay, startlocationID,
                       random_seed);
  }

}

//====================================================================
//* run JULIeT.definePropagatingParticle() by JNI --------------------
//  int flavor      : particle flavor (defined by the JULIeT particle class)
//  int doublet     : particle doublet (defined by the JULIeT particle class)
//  double energy   : initial particle energy [GeV]
//===================================================================
void I3Juliet::RunJULIeTDefinePropagatingParticle(int flavor, int doublet, double energy)
{

  cerr << "DefinePropagationParticle called!" << endl;
  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    fJNIEnv->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagatingParticleID,
                        flavor,doublet,energy);
  }
  cerr << "DefinePropagationParticle End!" << endl;
}

//===================================================================
//* run JULIeT.definePropagationGeometry() by JNI --------------------
//  double x, y, z  : Initial location of the particle [cm]. Use them for J3Line
//  double nx,ny,nz    : Initial direction of the particle. Use them for J3Line
// 
// All these argumets above must be the ones represened by the IceCube coordinate.
//===================================================================
void I3Juliet::RunJULIeTDefinePropagationGeometry(
                                          double x_ice3, double y_ice3, double z_ice3,
                                          double nx, double ny, double nz)
{

  // Generate J3Vector to define the initial location
  jobject initialLocationObj =
    fJNIEnv->NewObject(fJ3Vector, fJNIEnv->GetMethodID(fJ3Vector, "<init>", "(DDD)V"),
                   x_ice3,y_ice3,z_ice3);
  // Generate J3UnitVector to define the direction
  jobject directionObj =
    fJNIEnv->NewObject(fJ3UnitVector, fJNIEnv->GetMethodID(fJ3UnitVector, "<init>", "(DDD)V"),
                       nx,ny,nz);
  // Generate J3Line to define the propagating particle axis
  jobject initialAxisObj =
    fJNIEnv->NewObject(fJ3Line, 
    fJNIEnv->GetMethodID(fJ3Line, "<init>", 
                     "(Lgeometry/J3Vector;Lgeometry/J3UnitVector;)V"),
                   initialLocationObj,directionObj);

  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    fJNIEnv->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagationGeometryID,
                        initialAxisObj);
    // delete references to the java objects
    fJNIEnv->DeleteLocalRef(initialLocationObj);
    fJNIEnv->DeleteLocalRef(directionObj);
    fJNIEnv->DeleteLocalRef(initialAxisObj);
  }else{
    // delete references to the java objects
    fJNIEnv->DeleteLocalRef(initialLocationObj);
    fJNIEnv->DeleteLocalRef(directionObj);
    fJNIEnv->DeleteLocalRef(initialAxisObj);
    return;
  }
}

//===================================================================
//* run JULIeT.definePropagationGeometry() by JNI --------------------
//  double x, y, z  : Initial location of the particle [cm]. Use them for J3Line
//  double nadirDeg   : nadir [deg] in IceCube coordinate
//  double azimuthDeg : azimuth [deg] in IceCube coordinate
// 
// All these argumets above must be the ones represened by the IceCube coordinate.
//===================================================================
void I3Juliet::RunJULIeTDefinePropagationGeometry(
                                          double x_ice3, double y_ice3, double z_ice3,
                                          double nadirDeg, double azimuthDeg)
{

  cerr << "DefinePropagationGeometry called!" << endl;
  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    fJNIEnv->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagationGeometry_deg_ID,
                        x_ice3,y_ice3,z_ice3,nadirDeg,azimuthDeg);
  }else{
    return;
  }
}

//==================== ===============================================
//* run JULIeT.configurePropagationGeometry() by JNI --------------------
//
// Configure the primary particle geometry. See JulietEventGenerator.java
// and its API document for details.
//===================================================================
void I3Juliet::RunJULIeTConfigurePropagationGeometry()
{

  cerr << "configurePropagationGeometry called!" << endl;
  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    fJNIEnv->CallVoidMethod(fJULIeTObj,fJULIeT_configurePropagationGeometryID);
  }else{
    return;
  }
}

//===================================================================
//* run JULIeT.runSingleEvent() by JNI ---------------------------
void I3Juliet::RunJULIeTRunSingleEvent()
{

  // Run JULIeT.RunSingleEvent()
  if(fJULIeTObj != 0){
    fJNIEnv->CallVoidMethod(fJULIeTObj,fJULIeT_runSingleEventID);
  }else{
    return;
  }
}

//===================================================================
//* Create objects from JULILeT ----
void I3Juliet::CreateLocalRefOfJULIeT()
{

  cerr << "CreateLocalRefOfJULIeT is called!" << endl;

  // obtain particleIterator in the JulietEventGenerator object by running
  // JULIeT.getParticleIterator()
  fJULIeT_particleIteratorObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj, fJULIeT_getParticleIteratorID);

  // obtain trackParticleIterator in the JulietEventGenerator object by running
  // JULIeT.getTrackParticleIterator()
  fJULIeT_trackParticleIteratorObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj, fJULIeT_getTrackParticleIteratorID);

  // obtain locationIterator in the JulietEventGenerator object by running
  // JULIeT.getLocationIce3Iterator()
  fJULIeT_locationIce3IteratorObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj, fJULIeT_getLocationIce3IteratorID);

  // obtain trackLocationIce3Iterator in the JulietEventGenerator object by running
  // JULIeT.getTrackLocationIce3Iterator()
  fJULIeT_trackLocationIce3IteratorObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj, fJULIeT_getTrackLocationIce3IteratorID);

  // obtain "propParticle" (of Particle class) object 
  // in the JULIeTEventGenerator
  fJULIeT_propParticleObj = 
    fJNIEnv->GetObjectField(fJULIeTObj, fJULIeT_propParticleID);

  // obtain startLocation_J3Vector_center in the JulietEventGenerator object 
  // by running JULIeT.wherePrimaryParticleStartsInEarthCenterCoordinate()
  fJULIeT_startLocationCenterObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj,
                  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID);

  // obtain startLocation_J3Vector_ice3 in the JulietEventGenerator object 
  // by running JULIeT.wherePrimaryParticleStartsInIceCubeCoordinate()
  fJULIeT_startLocationIce3Obj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj,
                          fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID);

  // obtain startLocation_J3Vector_ice3 in the JulietEventGenerator object 
  // by running JULIeT.wherePrimaryParticleStartsInIceCubeCoordinate()
  fJULIeT_endLocationIce3Obj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj,
                          fJULIeT_wherePrimaryParticleEndsInIceCubeCoordinateID);

  // obtain direction_J3UnitVector_center in the JulietEventGenerator object 
  // by running JULIeT.getPrimaryParticleDirectionInEarthCenterCoordinate()
  fJULIeT_directionCenterObj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj,
                  fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID);

  // obtain direction_J3UnitVector_ice3 in the JulietEventGenerator object 
  // by running JULIeT.getPrimaryParticleDirectionInIceCubeCoordinate()
  fJULIeT_directionIce3Obj = 
    fJNIEnv->CallObjectMethod(fJULIeTObj,
                  fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID);

}

//===================================================================
//* Delete LocalRef of JULILeT ----
void I3Juliet::DeleteLocalRefOfJULIeT() 
{
  // delete references to the java list iterator
  fJNIEnv->DeleteLocalRef(fJULIeT_particleIteratorObj);
  fJNIEnv->DeleteLocalRef(fJULIeT_locationIce3IteratorObj);
  fJNIEnv->DeleteLocalRef(fJULIeT_trackParticleIteratorObj);
  fJNIEnv->DeleteLocalRef(fJULIeT_trackLocationIce3IteratorObj);

  // delete references to the java "propParticle" object
  fJNIEnv->DeleteLocalRef(fJULIeT_propParticleObj);

  // delete references to all the J3Vector objects we refered
  fJNIEnv->DeleteLocalRef(fJULIeT_startLocationCenterObj);
  fJNIEnv->DeleteLocalRef(fJULIeT_startLocationIce3Obj);
  fJNIEnv->DeleteLocalRef(fJULIeT_endLocationIce3Obj);

  // delete references to all the J3UnitVector objects we refered
  fJNIEnv->DeleteLocalRef(fJULIeT_directionCenterObj);
  fJNIEnv->DeleteLocalRef(fJULIeT_directionIce3Obj);
}

//===================================================================
//* Take and Display all the particle objects stored in the list ----
void I3Juliet::TakeParticlesInList()
{

  cerr << "TakeParticlesInList called!" << endl;

  CreateLocalRefOfJULIeT();  // create LocalRef from JULIeT to data members

  //
  // obtain primary energy and energy after the propagation
  //

  jdouble endEnergy = 
    fJNIEnv->CallDoubleMethod(fJULIeT_propParticleObj, fParticle_getEnergyID);
  jdouble primaryEnergy =
    fJNIEnv->GetDoubleField(fJULIeTObj,fJULIeT_primaryEnergyID);

  // display energy
  cerr << "Primary Energy " << primaryEnergy << " : Final Energy " 
       << endEnergy << " [GeV]" << endl;

  //
  // obtain primary's start position and direction 
  //

  jdouble x_center = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationCenterObj,fJ3Vector_getXID);
  jdouble y_center = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationCenterObj,fJ3Vector_getYID);
  jdouble z_center = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationCenterObj,fJ3Vector_getZID);
  jdouble x_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationIce3Obj,fJ3Vector_getXID);
  jdouble y_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationIce3Obj,fJ3Vector_getYID);
  jdouble z_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_startLocationIce3Obj,fJ3Vector_getZID);
  jdouble locationAlongTheAxis =
      fJNIEnv->CallDoubleMethod(fJULIeTObj, fJULIeT_getStartLocationAlongTheAxisID);

  jdouble nx_center =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionCenterObj,fJ3Vector_getXID);
  jdouble ny_center =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionCenterObj,fJ3Vector_getYID);
  jdouble nz_center =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionCenterObj,fJ3Vector_getZID);
  jdouble nx_ice3 =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionIce3Obj,fJ3Vector_getXID);
  jdouble ny_ice3 =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionIce3Obj,fJ3Vector_getYID);
  jdouble nz_ice3 =
      fJNIEnv->CallDoubleMethod(fJULIeT_directionIce3Obj,fJ3Vector_getZID);

  // display start position 
  cerr << "Primary pos/dir(EarthCenter) : (" 
       << x_center << ", " << y_center << ", " << z_center << ") [cm] / ("
       << nx_center << ", " << ny_center << ", " << nz_center << ")" << endl;

  cerr << "Primary pos/dir(IceCube)     : (" 
       << x_ice3 << ", " << y_ice3 << ", " << z_ice3 << ") [cm] / ("
       << nx_ice3 << ", " << ny_ice3 << ", " << nz_ice3 << ")" << endl;

  cerr << "Distance from the Earth entrance point : " << locationAlongTheAxis
       << " [cm]" << endl;

  //
  // obtain primary's end position (in IceCube coordinte)
  //

  jdouble xend_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_endLocationIce3Obj,fJ3Vector_getXID);
  jdouble yend_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_endLocationIce3Obj,fJ3Vector_getYID);
  jdouble zend_ice3 = 
      fJNIEnv->CallDoubleMethod(fJULIeT_endLocationIce3Obj,fJ3Vector_getZID);
  double tot_length = sqrt( (xend_ice3 - x_ice3) * (xend_ice3 - x_ice3) +
                            (yend_ice3 - y_ice3) * (yend_ice3 - y_ice3) +
                            (zend_ice3 - z_ice3) * (zend_ice3 - z_ice3) );

  // display start position 
  cerr << "Propagation ends at pos(IceCube) : (" 
       << xend_ice3 << ", " << yend_ice3 << ", " << zend_ice3 << ")[cm] : Length " 
       << tot_length << "[cm]" << endl;

  //
  // Now iterate the track particle list 
  //

  int i = 0;
  cerr << "*** Track Loop Start ***" << endl;
  while(fJNIEnv->CallBooleanMethod(fJULIeT_trackParticleIteratorObj,
                               fJNIEnv->GetMethodID(fIterator,"hasNext","()Z"))){
    // Get the particle object representing primary cascade partcle
    jobject particleObj =
      fJNIEnv->CallObjectMethod(fJULIeT_trackParticleIteratorObj,
                   fJNIEnv->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble energy =
      fJNIEnv->CallDoubleMethod(particleObj, fParticle_getEnergyID);
    jint flavor =
      fJNIEnv->CallIntMethod(particleObj, fParticle_getFlavorID);
    jint doublet =
      fJNIEnv->CallIntMethod(particleObj, fParticle_getDoubletID);
    jstring particleName = (jstring)
      fJNIEnv->CallStaticObjectMethod(fParticle, fParticle_particleNameID,flavor,doublet);
    const char* name = fJNIEnv->GetStringUTFChars(particleName,NULL);

    // Get the J3Vector object accomodating particle location info
    jobject locationObj =
      fJNIEnv->CallObjectMethod(fJULIeT_trackLocationIce3IteratorObj,
                     fJNIEnv->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble x =
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getXID);
    jdouble y =
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getYID);
    jdouble z =
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getZID);

    double dist = sqrt( (x - x_ice3) * (x - x_ice3) +
                        (y - y_ice3) * (y - y_ice3) +
                        (z - z_ice3) * (z - z_ice3) );

    // display
    cerr << ++i << " Flabor " << flavor << " " << std::setw(8) << name
         << " Energy " << std::setw(12) << energy  << "[GeV] Pos ("
         << std::setw(12) << x << ", " 
         << std::setw(12) << y << ", " 
         << std::setw(12) << z << ")[cm] Dist " 
         << std::setw(12) << dist << "[cm]" << endl;

    // delete object reference
    fJNIEnv->ReleaseStringUTFChars(particleName,name);
    fJNIEnv->DeleteLocalRef(particleObj);
    fJNIEnv->DeleteLocalRef(locationObj);

  }

  //
  // Now iterate the cascade particle list 
  //

  i = 0;
  cerr << "*** Cascade Loop Start ***" << endl;
  while(fJNIEnv->CallBooleanMethod(fJULIeT_particleIteratorObj,
                               fJNIEnv->GetMethodID(fIterator,"hasNext","()Z"))){

    // Get the particle object representing primary cascade partcle
    jobject particleObj = 
      fJNIEnv->CallObjectMethod(fJULIeT_particleIteratorObj,
                   fJNIEnv->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble energy = 
      fJNIEnv->CallDoubleMethod(particleObj, fParticle_getEnergyID);
    jint flavor =
      fJNIEnv->CallIntMethod(particleObj, fParticle_getFlavorID);
    jint doublet =
      fJNIEnv->CallIntMethod(particleObj, fParticle_getDoubletID);
    jstring particleName = (jstring)
      fJNIEnv->CallStaticObjectMethod(fParticle, fParticle_particleNameID,flavor,doublet);
    const char* name = fJNIEnv->GetStringUTFChars(particleName,NULL);

    // Get the J3Vector object accomodating particle location info
    jobject locationObj = 
      fJNIEnv->CallObjectMethod(fJULIeT_locationIce3IteratorObj,
                     fJNIEnv->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble x = 
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getXID);
    jdouble y = 
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getYID);
    jdouble z = 
      fJNIEnv->CallDoubleMethod(locationObj,fJ3Vector_getZID);

    double distance = sqrt( (x - x_ice3) * (x - x_ice3) +
                          (y - y_ice3) * (y - y_ice3) +
                          (z - z_ice3) * (z - z_ice3) );

    // display
    cerr << ++i << " Flabor " << flavor << " " << std::setw(8) << name
         << " Energy " << std::setw(12) << energy  << "[GeV] Pos ("
         << std::setw(12) << x << ", " 
         << std::setw(12) << y << ", " 
         << std::setw(12) << z << ")[cm] Dist " 
         << std::setw(12) << distance << "[cm]" << endl;

    // delete object reference
    fJNIEnv->ReleaseStringUTFChars(particleName,name);
    fJNIEnv->DeleteLocalRef(particleObj);
    fJNIEnv->DeleteLocalRef(locationObj);

  }

  DeleteLocalRefOfJULIeT();

}

//===================================================================
//* get IDs of the Java methods -------------------------------------
void I3Juliet::GetMethodIDs()
{

  // getting Class "JulietEventGenerator"
  fJULIeT = fJNIEnv->FindClass("iceCube/uhe/event/JulietEventGenerator");
  if(fJULIeT == 0){
    cerr << "I3Juliet::GetMethodIDs:Cannot find the class JulietEventGenerator" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.runSingleEvent()
  fJULIeT_runSingleEventID = 
    fJNIEnv->GetMethodID(fJULIeT, "runSingleEvent", "()V");
  if(fJULIeT_runSingleEventID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method runSingleEvent() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getStartLocationAlongTheAxis()
  fJULIeT_getStartLocationAlongTheAxisID = 
    fJNIEnv->GetMethodID(fJULIeT, "getStartLocationAlongTheAxis", "()D");
  if(fJULIeT_getStartLocationAlongTheAxisID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getStartLocationAlongTheAxis() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.wherePrimaryParticleStartsInEarthCenterCoordinate()
  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID = 
    fJNIEnv->GetMethodID(fJULIeT, "wherePrimaryParticleStartsInEarthCenterCoordinate", 
                     "()Lgeometry/J3Vector;");
  if(fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method wherePrimaryParticleStartsInEarthCenterCoordinate() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.wherePrimaryParticleStartsInIceCubeCoordinate()
  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID = 
    fJNIEnv->GetMethodID(fJULIeT, "wherePrimaryParticleStartsInIceCubeCoordinate", 
                     "()Lgeometry/J3Vector;");
  if(fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method wherePrimaryParticleStartsInIceCubeCoordinate() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.wherePrimaryParticleEndsInIceCubeCoordinate()
  fJULIeT_wherePrimaryParticleEndsInIceCubeCoordinateID = 
    fJNIEnv->GetMethodID(fJULIeT, "wherePrimaryParticleEndsInIceCubeCoordinate", 
                     "()Lgeometry/J3Vector;");
  if(fJULIeT_wherePrimaryParticleEndsInIceCubeCoordinateID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method wherePrimaryParticleEndsInIceCubeCoordinate() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getPrimaryParticleDirectionInEarthCenterCoordinateID
  fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID = 
    fJNIEnv->GetMethodID(fJULIeT, "getPrimaryParticleDirectionInEarthCenterCoordinate",
                     "()Lgeometry/J3UnitVector;");
  if(fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getPrimaryParticleDirectionInEarthCenterCoordinate() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getPrimaryParticleDirectionInIceCubeCoordinateID
  fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID = 
    fJNIEnv->GetMethodID(fJULIeT, "getPrimaryParticleDirectionInIceCubeCoordinate",
                     "()Lgeometry/J3UnitVector;");
  if(fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getPrimaryParticleDirectionInIceCubeCoordinate() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.definePropagatingParticle()
  fJULIeT_definePropagatingParticleID = 
    fJNIEnv->GetMethodID(fJULIeT, "definePropagatingParticle", 
                     "(IIDLgeometry/J3Line;)V");
  if(fJULIeT_definePropagatingParticleID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method definePropagatingParticle() in JULIeT" << endl; 
    exit(1);
  }

  // ID JulietEventGenerator.definePropagationGeometry()
  fJULIeT_definePropagationGeometryID = 
    fJNIEnv->GetMethodID(fJULIeT, "definePropagationGeometry", 
                     "(Lgeometry/J3Line;)V");
  if(fJULIeT_definePropagationGeometryID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method definePropagationGeometry((Lgeometry/J3Line;)V) in JULIeT" << endl; 
    exit(1);
  }

  // ID JulietEventGenerator.definePropagationGeometry()
  fJULIeT_definePropagationGeometry_deg_ID = 
    fJNIEnv->GetMethodID(fJULIeT, "definePropagationGeometry", 
                     "(DDDDD)V");
  if(fJULIeT_definePropagationGeometry_deg_ID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method definePropagationGeometry((DDDDD)V) in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.configurePropagationGeometry()
  fJULIeT_configurePropagationGeometryID = 
    fJNIEnv->GetMethodID(fJULIeT, "configurePropagationGeometry", 
                     "()V");
  if(fJULIeT_configurePropagationGeometryID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method configurePropagationGeometry() in JULIeT" << endl; 
    exit(1);
  }

  // ID JulietEventGenerator.getParticleIterator()
  fJULIeT_getParticleIteratorID = 
    fJNIEnv->GetMethodID(fJULIeT, "getParticleIterator", 
                     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getParticleIteratorID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getParticlesIterator() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getTrackParticleIterator()
  fJULIeT_getTrackParticleIteratorID = 
    fJNIEnv->GetMethodID(fJULIeT, "getTrackParticleIterator", 
                     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getTrackParticleIteratorID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getTrackParticleIterator() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getLocationIce3Iterator()
  fJULIeT_getLocationIce3IteratorID = 
    fJNIEnv->GetMethodID(fJULIeT, "getLocationIce3Iterator", 
                     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getLocationIce3IteratorID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getLocationIce3Iterator() in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.getTrackLocationIce3Iterator()
  fJULIeT_getTrackLocationIce3IteratorID = 
    fJNIEnv->GetMethodID(fJULIeT, "getTrackLocationIce3Iterator", 
                     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getTrackLocationIce3IteratorID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getTrackLocationIce3Iterator() in JULIeT" << endl;
    exit(1);
  }


  // ID JulietEventGenerator.propParticle (of JULIeT particle class)
  fJULIeT_propParticleID = 
    fJNIEnv->GetFieldID(fJULIeT, "propParticle", 
                     "LiceCube/uhe/particles/Particle;"); 
  if(fJULIeT_propParticleID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the member variable propParticle in JULIeT" << endl;
    exit(1);
  }

  // ID JulietEventGenerator.primaryEnergy (of <double>)
  fJULIeT_primaryEnergyID = 
    fJNIEnv->GetFieldID(fJULIeT, "primaryEnergy","D"); 
  if(fJULIeT_primaryEnergyID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the member variable primaryEnergy in JULIeT" << endl;
    exit(1);
  }

  // Getting the Java util class ListIterator 
  fIterator = fJNIEnv->FindClass("java/util/ListIterator");
  if(fIterator == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find ListIterator class in java.util" << endl;
    exit(1);
  }

  // getting JULIeT Particle class
  fParticle = fJNIEnv->FindClass("iceCube/uhe/particles/Particle");
  if(fParticle == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find Particle class in the JULIeT" << endl;
    exit(1);
  }

  // ID Particle.getEnergy()
  fParticle_getEnergyID = fJNIEnv->GetMethodID(fParticle, "getEnergy", "()D");
  if(fParticle_getEnergyID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getEnergy() in Particle" << endl;
    exit(1);
  }

  // ID Particle.getFlavor()
  fParticle_getFlavorID = fJNIEnv->GetMethodID(fParticle, "getFlavor", "()I");
  if(fParticle_getFlavorID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getFlavor() in Particle" << endl;
    exit(1);
  }

  // ID Particle.getFlavor()
  fParticle_getDoubletID = fJNIEnv->GetMethodID(fParticle, "getDoublet", "()I");
  if(fParticle_getDoubletID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getDoublet() in Particle" << endl;
    exit(1);
  }

  // ID Particle.particleName()
  fParticle_particleNameID = fJNIEnv->GetStaticMethodID(fParticle, 
                                                    "particleName", 
                                                    "(II)Ljava/lang/String;");
  if(fParticle_particleNameID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method particleName() in Particle" << endl;
    exit(1);
  }

  // getting J3Vector class
  fJ3Vector = fJNIEnv->FindClass("geometry/J3Vector");
  if(fJ3Vector == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find J3Vector class in the JULIeT" << endl;
    exit(1);
  }

  // ID J3Vector.getX()
  fJ3Vector_getXID = fJNIEnv->GetMethodID(fJ3Vector,"getX","()D");
  if(fJ3Vector_getXID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getX() in J3Vector" << endl;
    exit(1);
  }

  // ID J3Vector.getY()
  fJ3Vector_getYID = fJNIEnv->GetMethodID(fJ3Vector,"getY","()D");
  if(fJ3Vector_getYID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getY() in J3Vector" << endl;
    exit(1);
  }

  // ID J3Vector.getZ()
  fJ3Vector_getZID = fJNIEnv->GetMethodID(fJ3Vector,"getZ","()D");
  if(fJ3Vector_getZID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getZ() in J3Vector" << endl;
    exit(1);
  }

  // getting J3UnitVector class
  fJ3UnitVector = fJNIEnv->FindClass("geometry/J3UnitVector");
  if(fJ3UnitVector == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find J3UnitVector class in the JULIeT" << endl;
    exit(1);
  }

  // ID J3UnitVector.getX()
  fJ3UnitVector_getXID = fJNIEnv->GetMethodID(fJ3UnitVector,"getX","()D");
  if(fJ3UnitVector_getXID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getX() in J3UnitVector" << endl;
    exit(1);
  }

  // ID J3UnitVector.getY()
  fJ3UnitVector_getYID = fJNIEnv->GetMethodID(fJ3UnitVector,"getY","()D");
  if(fJ3UnitVector_getYID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getY() in J3UnitVector" << endl;
    exit(1);
  }

  // ID J3UnitVector.getZ()
  fJ3UnitVector_getZID = fJNIEnv->GetMethodID(fJ3UnitVector,"getZ","()D");
  if(fJ3UnitVector_getZID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method getZ() in J3UnitVector" << endl;
    exit(1);
  }

  // getting J3Line class
  fJ3Line = fJNIEnv->FindClass("geometry/J3Line");
  if(fJ3Line == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find J3Line class in the JULIeT" << endl;
    exit(1);
  }

  // Getting the Java number class Double
  fDouble = fJNIEnv->FindClass("java/lang/Double");
  if(fDouble == 0x00){
    cerr << "I3Juliet::GetMethodIDs:Cannot find Double class in java.lang" << endl;
    exit(1);
  }

  // ID Double.dobleValue()
  fDouble_doubleValueID = fJNIEnv->GetMethodID(fDouble, "doubleValue", "()D");
  if(fDouble_doubleValueID == 0) {
    cerr << "I3Juliet::GetMethodIDs:Cannot find the method doubleValue in Double" << endl;
    exit(1);
  }

}

//===================================================================
//* Main ---------------------------------------------------------
int main(void)
{

  I3Juliet *juliet = new I3Juliet();

  int    flavorID  = 2;
  int    doubletID = 0;
  double energy    = 1e10; // GeV

  double x0        = 0.;   // cm 
  double y0        = 0.;   // cm
  double z0        = 0.;   // cm
  double nadir     = 0.;   // degree
  double azimuth   = 0.;   // degree

  // choose constructor 
  bool isInteractive = 0;
  cout << "Interactive Mode ? [no/yes: 0/1, default:no]->";
  cin >> isInteractive;

  if (isInteractive) {

     juliet->CreateJVM();
     juliet->GetMethodIDs();

     // generate juliet with interactive constructor
     juliet->GenerateJULIeT(true);

     cout << "x [cm] in the ice3 coordinate->";
     cin >> x0;
     cout << "y [cm] in the ice3 coordinate->";
     cin >> y0;
     cout << "z [cm] in the ice3 coordinate->";
     cin >> z0;

     cout << "nadir [deg] in the ice3 coordinate->";
     cin >> nadir;
     cout << "azimuth [deg] in the ice3 coordinate->";
     cin >> azimuth;

  } else {

     juliet->Configure();

     // define particle and propagation geometry
     juliet->RunJULIeTDefinePropagatingParticle(flavorID, doubletID, energy);

  }

  // define geometry

  juliet->RunJULIeTDefinePropagationGeometry(x0, y0, z0, nadir, azimuth);

  // configure geometry

  juliet->RunJULIeTConfigurePropagationGeometry();

  // run single event

  juliet->RunJULIeTRunSingleEvent();

  // show particles

  juliet->TakeParticlesInList();

}

