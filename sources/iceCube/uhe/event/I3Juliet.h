#ifndef I3JULIET_H
#define I3JULIET_H
/**
 * class: I3Juliet
 *
 * Version $Id: I3Juliet.h,v 1.1.1.1 2006/03/11 04:31:23 syoshida Exp $
 *
 * date:  2 Dec 2004
 *
 * (c) 2003 IceCube Collaboration
 */


// other headers
#include "jni.h"

/**
 * @version $Id: I3Juliet.h,v 1.1.1.1 2006/03/11 04:31:23 syoshida Exp $
 *
 * @brief This module simulates EHE lepton propagation in Ice/Rock.
 *
 * uses JULIeT  
 *
 * @author 
 */
class I3Juliet {
 public:
  // constructor and destructor

  /**
   * Builds an instance of this class
   * @param ctx the context with which this module's built
   */
  I3Juliet();
  
  /**
   *Destroys an instance of this class
   */
  ~I3Juliet();

  // transitions
  /**
   * This module takes a configuration parameter and so it must be configured
   */
  void Configure();

  /**
   * This function just calls 'Fatal' since I haven't rigged it to be
   * reconfigurable
   */
  void Reconfigure();

  //
  // default, assignment, and copy constructor declared private
  //


  I3Juliet(const I3Juliet&);
  I3Juliet& operator=(const I3Juliet&);

  //
  // sub-routine functions
  //

  /**
   * Create Local references from JULIeT Objects to jobjects data members
   */ 
   void CreateLocalRefOfJULIeT();

  /**
   * Delete Local references of jobjects data members 
   */ 
   void DeleteLocalRefOfJULIeT();

  /**
   * Take and Display all the particle objects stored in the list ----
   */ 
   void TakeParticlesInList();

  //
  // functions required by JAVA VM / JNI interface
  //

  /**
   * Construct JAVA VM module
   */ 
   void CreateJVM();

  /**
   * Get method IDs of JULIeT via JNI 
   */ 
   void GetMethodIDs();

  /**
   * Generate JULIeT 
   */ 
   void GenerateJULIeT(bool isInteractive = false);

  /**
   * run JULIeT.definePropagatingParticle() by JNI --------------------
   * int flavor      : particle flavor (defined by the JULIeT particle class)
   * int doublet     : particle doublet (defined by the JULIeT particle class)
   * double energy   : initial particle energy [GeV]
   */ 
   void RunJULIeTDefinePropagatingParticle(int flavor, 
                                           int doublet, 
                                           double energy);
  /**
   * run JULIeT.definePropagationGeometry() by JNI --------------------
   * double x, y, z  : Initial location of the particle. Use them for J3Line
   * double nx,ny,nz : Initial direction of the particle. Use them for J3Line
   * All these argumets above must be the ones represened by the 
   * IceCube coordinate.
   */ 
   void RunJULIeTDefinePropagationGeometry(double x_ice3, double y_ice3, 
                                           double z_ice3,
                                           double nx, double ny, double nz);
  /**
   * run JULIeT.definePropagationGeometry() by JNI --------------------
   * double x, y, z  : Initial location of the particle. Use them for J3Line
   * double nadirDeg   : nadir [deg] in IceCube coordinate
   * double azimuthDeg : azimuth [deg] in IceCube coordinate
   * All these argumets above must be the ones represened by the 
   * IceCube coordinate.
   */ 
   void RunJULIeTDefinePropagationGeometry(double x_ice3, double y_ice3, 
                                           double z_ice3,
                                           double nadirDeg, double azimuthDeg);
  /**
   * run JULIeT.configurePropagationGeometry() by JNI --------------------
   * Configure the primary particle geometry. See JulietEventGenerator.java
   * and its API document for details.
   */ 
   void RunJULIeTConfigurePropagationGeometry();

  /**
   * run JULIeT.runSingleEvent() by JNI ---------------------------
   */ 
   void RunJULIeTRunSingleEvent();

private:

   //
   // fields
   // 

   double  fT0;    // start time offset of current track 

   JavaVM *fJvmService;
   JNIEnv *fJNIEnv;
   jclass  fJULIeT;
   jobject fJULIeTObj;

   jmethodID fJULIeT_getStartLocationAlongTheAxisID;
   jmethodID fJULIeT_wherePrimaryParticleStartsInEarthCenterCoordinateID;
   jmethodID fJULIeT_wherePrimaryParticleStartsInIceCubeCoordinateID;
   jmethodID fJULIeT_wherePrimaryParticleEndsInIceCubeCoordinateID;
   jmethodID fJULIeT_getPrimaryParticleDirectionInEarthCenterCoordinateID;
   jmethodID fJULIeT_getPrimaryParticleDirectionInIceCubeCoordinateID;
   jmethodID fJULIeT_definePropagatingParticleID;
   jmethodID fJULIeT_definePropagationGeometryID;
   jmethodID fJULIeT_definePropagationGeometry_deg_ID;
   jmethodID fJULIeT_configurePropagationGeometryID;
   jmethodID fJULIeT_runSingleEventID;          // Method ID runSingleEvent()
   jmethodID fJULIeT_getParticleIteratorID;     // Method ID getParticleIterator()
   jmethodID fJULIeT_getTrackParticleIteratorID;// Method ID getTrackParticleIterator()
   jmethodID fJULIeT_getLocationIce3IteratorID;     // Method ID getLocationIce3Iterator()
   jmethodID fJULIeT_getTrackLocationIce3IteratorID;// Method ID getTrackLocationIce3Iterator()
   jfieldID  fJULIeT_propParticleID;            // Field ID propParticle in JULIeT
   jfieldID  fJULIeT_primaryEnergyID;           // Field ID primaryEnergy in JULIeT

   jobject   fJULIeT_startLocationCenterObj;    // J3Vector object of the primary start location
   jobject   fJULIeT_startLocationIce3Obj;      // J3Vector object of the primary start location (ice3)
   jobject   fJULIeT_endLocationIce3Obj;        // J3Vector object of the primary end location (ice3)
   jobject   fJULIeT_directionCenterObj;        // J3UnitVector object of the primary direction
   jobject   fJULIeT_directionIce3Obj;          // J3UnitVector object of the primary direction (ice3)
   jobject   fJULIeT_particleIteratorObj;       // particleIterator 
   jobject   fJULIeT_trackParticleIteratorObj;  // trackParticleIterator 
   jobject   fJULIeT_locationIce3IteratorObj;       // locationIterator 
   jobject   fJULIeT_trackLocationIce3IteratorObj;  // trackLocationIce3Iterator 
   jobject   fJULIeT_propParticleObj;           // the JULIeT Particle object 

   jclass    fIterator;                     // iterator class (java/util/ListIterator)
   jclass    fParticle;                     // the JULIeT Particle class
   jmethodID fParticle_getEnergyID;         // Method ID getEnergy()
   jmethodID fParticle_getFlavorID;         // Method ID getFlavor()
   jmethodID fParticle_getDoubletID;        // Method ID getFlavor()
   jmethodID fParticle_particleNameID;      // Method ID particleName()
   jclass    fJ3Vector;                     // The J3Vector class
   jmethodID fJ3Vector_getXID;              // Method ID getXID()
   jmethodID fJ3Vector_getYID;              // Method ID getYID()
   jmethodID fJ3Vector_getZID;              // Method ID getZID()
   jmethodID fJ3UnitVector_getXID;          // Method ID getXID()
   jmethodID fJ3UnitVector_getYID;          // Method ID getYID()
   jmethodID fJ3UnitVector_getZID;          // Method ID getZID()
   jclass    fJ3UnitVector;                 // The J3UnitVector class
   jclass    fJ3Line;                       // The J3Line class

   jclass    fDouble;
   jmethodID fDouble_doubleValueID;         // Method ID Double.doubleValue()

};


#endif //I3PROPAGATORTOY_H
