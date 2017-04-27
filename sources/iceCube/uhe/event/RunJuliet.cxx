#include <jni.h>
#include <iostream>
#include "RunJuliet.h"

using namespace std;

//===================================================================
//* constructor -----------------------------------------------------
RunJuliet::RunJuliet() 
{
  fJULIeT_getStartLocationAlongTheAxisID = 0;
  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID = 0;
  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID = 0;
  fJULIeT_definePropagatingParticleID = 0;
  fJULIeT_definePropagationGeometryID = 0;
  fJULIeT_definePropagationGeometry_deg_ID = 0;
  fJULIeT_configurePropagationGeometryID = 0;
  fJULIeT_runSingleEventID = 0;      // Method ID runSingleEvent()
  fJULIeT_getParticleIteratorID = 0; // Method ID getParticleIterator()
  fJULIeT_getLocationIce3IteratorID = 0; // Method ID getLocationIce3Iterator()
  fJULIeT_propParticleID = 0;        // Field ID propParticle in JULIeT
  fJULIeT_primaryEnergyID = 0;       // Field ID primaryEnergy in JULIeT

  fIterator = 0;                    // iterator class (java/util/ListIterator)

  fParticle = 0;                    // the JULIeT Particle class
  fParticle_getEnergyID = 0;        // Method ID getEnergy()
  fParticle_getFlavorID = 0;        // Method ID getFlavor()
  fParticle_particleNameID = 0;     // Method ID particleName()

  fJ3Vector = 0;                    // The J3Vector class
  fJ3UnitVector = 0;                // The J3UnitVector class
  fJ3Line = 0;                      // The J3Line class

  fDouble_doubleValueID = 0;        // Method ID Double.doubleValue()

  CreateJVM();
}
//===================================================================
//* create JVM  -----------------------------------------------------
void RunJuliet::CreateJVM()
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
  vm_args.nOptions = 4;

  // Create the Java VM 
  jint res = JNI_CreateJavaVM(&jvmService,(void **)&env,&vm_args);
  if(res < 0) {
    cerr << "Cannot create Java VM" << endl;
    exit(1);
  }
}

//===================================================================
//* get IDs of the Java methods -------------------------------------
void RunJuliet::GetMethodIDs()
{

  // getting Class "JulietEventGenerator"
  fJULIeT = env->FindClass("iceCube/uhe/event/JulietEventGenerator");
  if(fJULIeT == 0){
    cerr << "Cannot find the class JulietEventGenerator" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.runSingleEvent()
  fJULIeT_runSingleEventID = 
    env->GetMethodID(fJULIeT, "runSingleEvent", "()V");
  if(fJULIeT_runSingleEventID == 0) {
    cerr << "Cannot find the method runSingleEvent() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.getStartLocationAlongTheAxis()
  fJULIeT_getStartLocationAlongTheAxisID = 
    env->GetMethodID(fJULIeT, "getStartLocationAlongTheAxis", "()D");
  if(fJULIeT_getStartLocationAlongTheAxisID == 0) {
    cerr << "Cannot find the method getStartLocationAlongTheAxis() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.wherePrimaryParticleStartsInEarthCenterCoordinate()
  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID = 
    env->GetMethodID(fJULIeT, "wherePrimaryParticleStartsInEarthCenterCoordinate", 
		     "()Lgeometry/J3Vector;");
  if(fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID == 0) {
    cerr << "Cannot find the method wherePrimaryParticleStartsInEarthCenterCoordinate() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.wherePrimaryParticleStartsInIceCubeCoordinate()
  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID = 
    env->GetMethodID(fJULIeT, "wherePrimaryParticleStartsInIceCubeCoordinate", 
		     "()Lgeometry/J3Vector;");
  if(fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID == 0) {
    cerr << "Cannot find the method wherePrimaryParticleStartsInIceCubeCoordinate() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.definePropagatingParticle()
  fJULIeT_definePropagatingParticleID = 
    env->GetMethodID(fJULIeT, "definePropagatingParticle", 
		     "(IID)V");
  if(fJULIeT_definePropagatingParticleID == 0) {
    cerr << "Cannot find the method definePropagatingParticle() in JULIeT" 
	 << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.definePropagationGeometry()
  fJULIeT_definePropagationGeometryID = 
    env->GetMethodID(fJULIeT, "definePropagationGeometry", 
		     "(Lgeometry/J3Line;)V");
  if(fJULIeT_definePropagationGeometryID == 0) {
    cerr << "Cannot find the method definePropagationGeometry() in JULIeT" 
	 << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.definePropagationGeometry()
  fJULIeT_definePropagationGeometry_deg_ID = 
    env->GetMethodID(fJULIeT, "definePropagationGeometry", 
		     "(DDDDD)V");
  if(fJULIeT_definePropagationGeometry_deg_ID == 0) {
    cerr << "Cannot find the method definePropagationGeometry() in JULIeT" 
	 << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.configurePropagationGeometry()
  fJULIeT_configurePropagationGeometryID = 
    env->GetMethodID(fJULIeT, "configurePropagationGeometry", 
		     "()V");
  if(fJULIeT_configurePropagationGeometryID == 0) {
    cerr << "Cannot find the method configurePropagationGeometry() in JULIeT" 
	 << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.getParticleIterator()
  fJULIeT_getParticleIteratorID = 
    env->GetMethodID(fJULIeT, "getParticleIterator", 
		     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getParticleIteratorID == 0) {
    cerr << "Cannot find the method getParticlesIterator() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.getLocationIterator()
  fJULIeT_getLocationIce3IteratorID = 
    env->GetMethodID(fJULIeT, "getLocationIce3Iterator", 
		     "()Ljava/util/ListIterator;"); 
  if(fJULIeT_getLocationIce3IteratorID == 0) {
    cerr << "Cannot find the method getLocationIce3Iterator() in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.propParticle (of JULIeT particle class)
  fJULIeT_propParticleID = 
    env->GetFieldID(fJULIeT, "propParticle", 
		     "LiceCube/uhe/particles/Particle;"); 
  if(fJULIeT_propParticleID == 0) {
    cerr << "Cannot find the member variable propParticle in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID JulietEventGenerator.primaryEnergy (of <double>)
  fJULIeT_primaryEnergyID = 
    env->GetFieldID(fJULIeT, "primaryEnergy","D"); 
  if(fJULIeT_primaryEnergyID == 0) {
    cerr << "Cannot find the member variable primaryEnergy in JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // Getting the Java util class ListIterator 
  fIterator = env->FindClass("java/util/ListIterator");
  if(fIterator == 0x00){
    cerr << "Cannot find ListIterator class in java.util" << endl;
    DeleteJVM();
    exit(1);
  }

  // getting JULIeT Particle class
  fParticle = env->FindClass("iceCube/uhe/particles/Particle");
  if(fParticle == 0x00){
    cerr << "Cannot find Particle class in the JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID Particle.getEnergy()
  fParticle_getEnergyID = env->GetMethodID(fParticle, "getEnergy", "()D");
  if(fParticle_getEnergyID == 0) {
    cerr << "Cannot find the method getEnergy() in Particle" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID Particle.getFlavor()
  fParticle_getFlavorID = env->GetMethodID(fParticle, "getFlavor", "()I");
  if(fParticle_getFlavorID == 0) {
    cerr << "Cannot find the method getFlavor() in Particle" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID Particle.particleName()
  fParticle_particleNameID = env->GetStaticMethodID(fParticle, 
						    "particleName", 
						    "(II)Ljava/lang/String;");
  if(fParticle_particleNameID == 0) {
    cerr << "Cannot find the method particleName() in Particle" << endl;
    DeleteJVM();
    exit(1);
  }


  // getting J3Vector class
  fJ3Vector = env->FindClass("geometry/J3Vector");
  if(fJ3Vector == 0x00){
    cerr << "Cannot find J3Vector class in the JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }


  // getting J3UnitVector class
  fJ3UnitVector = env->FindClass("geometry/J3UnitVector");
  if(fJ3UnitVector == 0x00){
    cerr << "Cannot find J3UnitVector class in the JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }


  // getting J3Line class
  fJ3Line = env->FindClass("geometry/J3Line");
  if(fJ3Line == 0x00){
    cerr << "Cannot find J3Line class in the JULIeT" << endl;
    DeleteJVM();
    exit(1);
  }


  // Getting the Java number class Double
  fDouble = env->FindClass("java/lang/Double");
  if(fDouble == 0x00){
    cerr << "Cannot find Double class in java.lang" << endl;
    DeleteJVM();
    exit(1);
  }

  // ID Double.dobleValue()
  fDouble_doubleValueID = env->GetMethodID(fDouble, "doubleValue", "()D");
  if(fDouble_doubleValueID == 0) {
    cerr << "Cannot find the method doubleValue in Double" << endl;
    DeleteJVM();
    exit(1);
  }

}



//===================================================================
//* Generate JULIeT object: run constructor -------------------------
void RunJuliet::GenerateJULIeT()
{
  // delete the former object if exists
  if(fJULIeTObj != 0)  env->DeleteLocalRef(fJULIeTObj);

  // run the constructor
  fJULIeTObj = 
    env->NewObject(fJULIeT, env->GetMethodID(fJULIeT, "<init>", "()V"));
}

//==================== ===============================================
//* run JULIeT.definePropagatingParticle() by JNI --------------------
//  int flavor      : particle flavor (defined by the JULIeT particle class)
//  int doublet     : particle doublet (defined by the JULIeT particle class)
//  double energy   : initial particle energy [GeV]
//===================================================================
void RunJuliet::RunJULIeTDefinePropagatingParticle(int flavor, int doublet, double energy)
{

  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    env->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagatingParticleID,
			flavor,doublet,energy);
  }
}

//===================================================================
//* run JULIeT.definePropagationGeometry() by JNI --------------------
//  double x, y, z  : Initial location of the particle. Use them for J3Line
//  double nx,ny,nz : Initial direction of the particle. Use them for J3Line
// 
// All these argumets above must be the ones represened by the IceCube coordinate.
//===================================================================
void RunJuliet::RunJULIeTDefinePropagationGeometry(
					  double x_ice3, double y_ice3, double z_ice3,
					  double nx, double ny, double nz)
{

  // Generate J3Vector to define the initial location
  jobject initialLocationObj =
    env->NewObject(fJ3Vector, env->GetMethodID(fJ3Vector, "<init>", "(DDD)V"),
		   x_ice3,y_ice3,z_ice3);
  // Generate J3UnitVector to define the direction
  jobject directionObj =
    env->NewObject(fJ3UnitVector, env->GetMethodID(fJ3UnitVector, "<init>", "(DDD)V"),
		   nx,ny,nz);
  // Generate J3Line to define the propagating particle axis
  jobject initialAxisObj =
    env->NewObject(fJ3Line, 
    env->GetMethodID(fJ3Line, "<init>", 
		     "(Lgeometry/J3Vector;Lgeometry/J3UnitVector;)V"),
		   initialLocationObj,directionObj);

  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    env->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagationGeometryID,
			initialAxisObj);
    // delete references to the java objects
    env->DeleteLocalRef(initialLocationObj);
    env->DeleteLocalRef(directionObj);
    env->DeleteLocalRef(initialAxisObj);
  }else{
    // delete references to the java objects
    env->DeleteLocalRef(initialLocationObj);
    env->DeleteLocalRef(directionObj);
    env->DeleteLocalRef(initialAxisObj);
    return;
  }
}

//===================================================================
//* run JULIeT.definePropagationGeometry() by JNI --------------------
//  double x, y, z  : Initial location of the particle. Use them for J3Line
//  double nadirDeg   : nadir [deg] in IceCube coordinate
//  double azimuthDeg : azimuth [deg] in IceCube coordinate
// 
// All these argumets above must be the ones represened by the IceCube coordinate.
//===================================================================
void RunJuliet::RunJULIeTDefinePropagationGeometry(
					  double x_ice3, double y_ice3, double z_ice3,
					  double nadirDeg, double azimuthDeg)
{

  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    env->CallVoidMethod(fJULIeTObj,fJULIeT_definePropagationGeometry_deg_ID,
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
void RunJuliet::RunJULIeTConfigurePropagationGeometry()
{

  // Run JULIeT.defineProagatingParticle()
  if(fJULIeTObj != 0){
    env->CallVoidMethod(fJULIeTObj,fJULIeT_configurePropagationGeometryID);
  }else{
    return;
  }
}


//===================================================================
//* run JULIeT.runSingleEvent() by JNI ---------------------------
void RunJuliet::RunJULIeTRunSingleEvent()
{

  // Run JULIeT.RunSingleEvent()
  if(fJULIeTObj != 0){
    env->CallVoidMethod(fJULIeTObj,fJULIeT_runSingleEventID);
  }else{
    return;
  }
}

//===================================================================
//* Take and Display all the particle objects stored in the list ----
void RunJuliet::TakeParticlesInList()
{

  // obtain startLocation_J3Vector_center in the JulietEventGenerator object 
  // by running JULIeT.wherePrimaryParticleStartsInEarthCenterCoordinate()
  fJULIeT_startLocationCenterObj =
    env->CallObjectMethod(fJULIeTObj, 
			  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID);

  // obtain startLocation_J3Vector_ice3 in the JulietEventGenerator object 
  // by running JULIeT.wherePrimaryParticleStartsInIceCubeCoordinate()
  fJULIeT_startLocationIce3Obj =
    env->CallObjectMethod(fJULIeTObj, 
			  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID);

  // obtain particleIterator in the JulietEventGenerator object by running
  // JULIeT.getParticleIterator()
  fJULIeT_particleIteratorObj = 
    env->CallObjectMethod(fJULIeTObj, fJULIeT_getParticleIteratorID);


  // obtain locationIterator in the JulietEventGenerator object by running
  // JULIeT.getLocationIterator()
  fJULIeT_locationIce3IteratorObj = 
    env->CallObjectMethod(fJULIeTObj, fJULIeT_getLocationIce3IteratorID);


  // obtain "propParticle" (of Particle class) object 
  // in the JULIeTEventGnerator
  fJULIeT_propParticleObj = 
    env->GetObjectField(fJULIeTObj, fJULIeT_propParticleID);

  // obtain methods on J3Vector ("geometry package")
  jmethodID J3Vector_getXID = env->GetMethodID(fJ3Vector,"getX","()D");
  jmethodID J3Vector_getYID = env->GetMethodID(fJ3Vector,"getY","()D");
  jmethodID J3Vector_getZID = env->GetMethodID(fJ3Vector,"getZ","()D");


  // obtain primary energy and energy after the propagation
  jdouble endEnergy = 
    env->CallDoubleMethod(fJULIeT_propParticleObj, fParticle_getEnergyID);
  jdouble primaryEnergy =
    env->GetDoubleField(fJULIeTObj,fJULIeT_primaryEnergyID);

  // display
  cout << "Primary Energy " << primaryEnergy << " Final Energy " 
       << endEnergy << " [GeV]" << endl;

  // Obtain the info on the location where this primary particle
  // starts its propagation
  fJULIeT_startLocationCenterObj = 
    env->CallObjectMethod(fJULIeTObj,
			  fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID);
  fJULIeT_startLocationIce3Obj = 
    env->CallObjectMethod(fJULIeTObj,
			  fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID);
  jdouble locationAlongTheAxis =
    env->CallDoubleMethod(fJULIeTObj,
			  fJULIeT_getStartLocationAlongTheAxisID);

  // display
  jdouble x_center = 
      env->CallDoubleMethod(fJULIeT_startLocationCenterObj,J3Vector_getXID);
  jdouble y_center = 
      env->CallDoubleMethod(fJULIeT_startLocationCenterObj,J3Vector_getYID);
  jdouble z_center = 
      env->CallDoubleMethod(fJULIeT_startLocationCenterObj,J3Vector_getZID);
  cout << "This track starts at " 
       << x_center << " [cm] " << y_center << " [cm] " << z_center 
       << " [cm] in EarthCenterCoordinate" << endl;

  jdouble x_ice3 = 
      env->CallDoubleMethod(fJULIeT_startLocationIce3Obj,J3Vector_getXID);
  jdouble y_ice3 = 
      env->CallDoubleMethod(fJULIeT_startLocationIce3Obj,J3Vector_getYID);
  jdouble z_ice3 = 
      env->CallDoubleMethod(fJULIeT_startLocationIce3Obj,J3Vector_getZID);
  cout << "This track starts at " 
       << x_ice3 << " [cm] " << y_ice3 << " [cm] " << z_ice3 
       << " [cm] in IceCubeCoordinate" << endl;

  cout << "Distance from the Earth entrance point : " << locationAlongTheAxis
       << " [cm]" << endl;


  // Now iterate the particle list 
  int i =0;
  while(env->CallBooleanMethod(fJULIeT_particleIteratorObj,
			       env->GetMethodID(fIterator,"hasNext","()Z"))){

    // Get the particle object representing primary cascade partcle
    jobject particleObj = 
      env->CallObjectMethod(fJULIeT_particleIteratorObj,
			    env->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble energy = 
      env->CallDoubleMethod(particleObj, fParticle_getEnergyID);
    jint flavor =
      env->CallIntMethod(particleObj, fParticle_getFlavorID);
    jstring particleName = (jstring)
      env->CallStaticObjectMethod(fParticle, fParticle_particleNameID,flavor,1);
    const char* name = env->GetStringUTFChars(particleName,NULL);

    // Get the J3Vector object accomodating particle location info
    jobject locationObj = 
      env->CallObjectMethod(fJULIeT_locationIce3IteratorObj,
			    env->GetMethodID(fIterator,"next","()Ljava/lang/Object;"));
    jdouble x = 
      env->CallDoubleMethod(locationObj,J3Vector_getXID);
    jdouble y = 
      env->CallDoubleMethod(locationObj,J3Vector_getYID);
    jdouble z = 
      env->CallDoubleMethod(locationObj,J3Vector_getZID);

    // display
    cout << ++i << " Flabor " << flavor << " " << name 
	 << " Energy " << energy  << " [GeV] " 
	 << x << " [cm] " << y << " [cm] " << z << " [cm]" << endl;

    // delete object reference
    env->ReleaseStringUTFChars(particleName,name);
    env->DeleteLocalRef(particleObj);
    env->DeleteLocalRef(locationObj);

  }

  // delete references to all the J3Vector objects we refered
  env->DeleteLocalRef(fJULIeT_startLocationCenterObj);
  env->DeleteLocalRef(fJULIeT_startLocationIce3Obj);
  // delete references to the java list iterator
  env->DeleteLocalRef(fJULIeT_particleIteratorObj);
  env->DeleteLocalRef(fJULIeT_locationIce3IteratorObj);
  // delete references to the java "propParticle" object
  env->DeleteLocalRef(fJULIeT_propParticleObj);


}

//===================================================================
//* Main function to run myself -------------------------------------
//* This is an example. ---------------------------------------------


int main(void){

  RunJuliet* cplus;

  // generate myself
  cplus = new RunJuliet();

  // get the method IDs
  cplus->GetMethodIDs();

  // generate the JulietEventGenerator object (Java class)
  cplus->GenerateJULIeT();

  // Setup the geometry
  double x,y,z;
  cout << "x [cm] in the ice3 coordinate->";
  cin >> x;
  cout << "y [cm] in the ice3 coordinate->";
  cin >> y;
  cout << "z [cm] in the ice3 coordinate->";
  cin >> z;

  double nx,ny,nz;
  cout << "nx in the ice3 coordinate->";
  cin >> nx;
  cout << "ny in the ice3 coordinate->";
  cin >> ny;
  cout << "nz in the ice3 coordinate->";
  cin >> nz;

  cplus->RunJULIeTDefinePropagationGeometry(x,y,z,nx,ny,nz);
  cplus->RunJULIeTConfigurePropagationGeometry();

  // generate 5 events
  for(int i=0;i<5;i++){
  // Run JulietEventGenerator.runSingleEvent()
    cplus->RunJULIeTRunSingleEvent();

  // Take all the JULIeT particle pbjects and display them
    cplus->TakeParticlesInList();
  }

  // change primay propagating particle
  int flavor;
  int doublet;
  double energy;
  cout << "flavor? (0[e]/1[mu]/2[tau])->";
  cin >> flavor;
  cout << "doublet? (0[neutrno]/1[charged lepton])->";
  cin >> doublet;
  cout << "enegry? [GeV]->";
  cin >> energy;
  cplus->RunJULIeTDefinePropagatingParticle(flavor,doublet,energy);

  // Setup geometry : A different way
  cout << "x [cm] in the ice3 coordinate->";
  cin >> x;
  cout << "y [cm] in the ice3 coordinate->";
  cin >> y;
  cout << "z [cm] in the ice3 coordinate->";
  cin >> z;

  double theta,alpha;
  cout << "nadir [deg] in the ice3 coordinate->";
  cin >> theta;
  cout << "azimuth [deg] in the ice3 coordinate->";
  cin >> alpha;

  cplus->RunJULIeTDefinePropagationGeometry(x,y,z,theta,alpha);
  cplus->RunJULIeTConfigurePropagationGeometry();

  
  // generate 5 events
  for(int i=0;i<5;i++){
    // Run JulietEventGenerator.runSingleEvent()
    cplus->RunJULIeTRunSingleEvent();

    // Take all the JULIeT particle pbjects and display them
    cplus->TakeParticlesInList();
  }

  delete(cplus);

  return 0;
}


