#include <jni.h>

#ifndef _Included_RunJuliet
#define _Included_RunJuliet

class RunJuliet {

 public:
  RunJuliet();                         // constructor

  virtual ~RunJuliet(){                // destructor
    if(fJULIeTObj!=0) env->DeleteLocalRef(fJULIeTObj);
    DeleteJVM();
  }  


  void CreateJVM();                 // Create the JAVA virtual machine
  void GetMethodIDs();              // Acquire all the relevant java method IDs
  void GenerateJULIeT();            // Generate the java JULIeT object
  void RunJULIeTDefinePropagatingParticle(int flavor, int doublet, double energy);
                                    // Run the JULIeT.definePropagatingParticle()
  void RunJULIeTDefinePropagationGeometry(double x_ice3, double y_ice3, double z_ice3,
					  double nx, double ny, double nz);
                                    // Run the JULIeT.definePropagationGeometry()
  void RunJULIeTDefinePropagationGeometry(double x_ice3, double y_ice3, double z_ice3,
					  double nadirDeg, double azimuthDeg);
                                    // Run the JULIeT.definePropagationGeometry()
  void RunJULIeTConfigurePropagationGeometry();
                                    // Run the JULIeT.configurePropagationGeometry()
  void RunJULIeTRunSingleEvent();   // Run the java JULIeT.runSingleEvent()
  void TakeParticlesInList();       // Take and Display all the JULIeT particle objects
  inline void DeleteJVM() { jvmService->DestroyJavaVM();}
  // Delete the Java virtual machine

 private:

  JNIEnv *env;                             // JNI pointer
  JavaVM *jvmService;                      // JVM Pointer


  jclass fJULIeT;                          // JULIeTEventGenerator class
  jobject fJULIeTObj;                      // JULIeTEventGenerator objects
  jmethodID fJULIeT_getStartLocationAlongTheAxisID;
                                      // MethodID getStartLocationAlongTheAxis()
  jmethodID fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID;
                  // MethodID wherePrimaryParticleStartsInEarthCenterCoordinate()
  jmethodID fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID;
                  // MethodID wherePrimaryParticleStartsInIceCubeCoordinate()
  jmethodID fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID;
  jmethodID fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID;

  jmethodID fJULIeT_definePropagatingParticleID;      
                                           // Method ID definePropagatingParticle()
  jmethodID fJULIeT_definePropagationGeometryID;      
                                           // Method ID definePropagationGeometry()
  jmethodID fJULIeT_definePropagationGeometry_deg_ID;
                                           // Method ID definePropagationGeometry()
  jmethodID fJULIeT_configurePropagationGeometryID;      
                                           // Method ID configurePropagationGeometry()
  jmethodID fJULIeT_runSingleEventID;      // Method ID runSingleEvent()
  jmethodID fJULIeT_getParticleIteratorID; // Method ID getParticleIterator()
  jmethodID fJULIeT_getLocationIce3IteratorID; // Method ID getLocationIce3Iterator()
  jfieldID  fJULIeT_propParticleID;        // Field ID propParticle in JULIeT
  jfieldID  fJULIeT_primaryEnergyID;       // Field ID primaryEnergy in JULIeT

  jobject   fJULIeT_directionCenterObj;        // J3UnitVector object of the primary direction
  jobject   fJULIeT_directionIce3Obj;          // J3UnitVector object of the primary direction (ice3)
  jobject   fJULIeT_startLocationCenterObj;// J3Vector object of the start location
  jobject   fJULIeT_startLocationIce3Obj;  // J3Vector object of the start location (ice3)
  jobject   fJULIeT_particleIteratorObj;   // particleIterator 
  jobject   fJULIeT_locationIce3IteratorObj; // locationIterator 
  jobject   fJULIeT_propParticleObj;       // the JULIeT Particle object 
                                           // for "propParticle"

  jclass    fIterator;                    // iterator class (java/util/ListIterator)

  jclass fParticle;                       // the JULIeT Particle class
  jmethodID fParticle_getEnergyID;        // Method ID Particle.getEnergy() 
  jmethodID fParticle_getFlavorID;        // Method ID Particle.getFlavor()
  jmethodID fParticle_particleNameID;     // Static Method ID Particle.particleName()

  jclass    fJ3Vector;                    // The J3Vector class
  jclass    fJ3UnitVector;                // The J3UnitVector class
  jclass    fJ3Line;                      // The J3Line class

  jclass    fDouble;                      // Double class (java/lang/Double)
  jmethodID fDouble_doubleValueID;        // Method ID Double.doubleValue()

};

#endif
