package iceCube.uhe.event;

import numRecipes.*;
import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.particles.*;

import java.io.*;
import java.util.*;

import hep.aida.*;

public class MakeJuliet4Gen2JaidaTree {

    /** Define the type of run **/
    static int mediumID = 0; //ice
    static int doCC = 1;
    static int doNC = 1;
    static int doMuBrems = 1; 
    static int doTauBrems = 1;
    static int doMuKnock = 0;
    static int doTauKnock = 0;
    static int doMu2e = 1; 
    static int doTau2e =1;
    static int doMu2mu = 0;
    static int doTau2mu = 0;
    static int doMu2tau = 0;
    static int doTau2tau = 0;
    static int doMuPN = 1;
    static int doTauPN = 1;
    static int doGR = 0;
    static int doMuDecay = 1;
    static int doTauDecay =1;
    static int posID = 2;

    static double gen2Xorigin = -4.510014e+04; // gen2 center in the ice3 coordinate [cm]
    static double gen2Yorigin = 1.844949e+04; // gen2 center in the ice3 coordinate [cm]
    static double gen2Zorigin = 0.0; // gen2 center in the ice3 coordinate [cm]

    static double epsilon = 1.0e-2;
    static double cm2m = 1.0e-2;

    // Main method. In order to run JULIET 
    public static void main(String[] args) throws IOException {

	int flavorID = 1;
	int doubletID = 1;
	double energy = 1.0e9; // [GeV]
	int numberOfEvents = 10;
	String jaidaTreeFileName = null;
	String jaidaTreeHeaderFileName = null;

	if(args.length!=5){
            System.out.println("Usage: MakeJuliet4Gen2JaidaTree flavorID doubletID number-of-events energy [GeV] tree-file-name");
            System.exit(0);
        }else{
            flavorID = Integer.valueOf(args[0]).intValue();
            doubletID = Integer.valueOf(args[1]).intValue();
            numberOfEvents = Integer.valueOf(args[2]).intValue();
            energy = Double.valueOf(args[3]).doubleValue();
	    jaidaTreeHeaderFileName = args[4];
	    jaidaTreeFileName = jaidaTreeHeaderFileName.concat(".aida");
        }
	System.err.format("(flavor doublet) = (%d %d) Energy=%e [GeV]\n",
			  flavorID, doubletID, energy);
	System.err.println("Tree file name " + jaidaTreeFileName);

	RandomGenerator rand = new RandomGenerator();

	// generate RunManager object
	JulietEventGenerator4Gen2 generator = 
	    new  JulietEventGenerator4Gen2(flavorID, doubletID, energy, mediumID,
					   doCC, doNC, doMuBrems, doTauBrems,
					   doMuKnock, doTauKnock, doMu2e, doTau2e,
					   doMu2mu, doTau2mu, doMu2tau, doTau2tau,
					   doMuPN, doTauPN, doGR, doMuDecay, doTauDecay,posID);

	// Jaida FreeHep objects
        IAnalysisFactory jaidaFactory = IAnalysisFactory.create();
        ITree jaidaTree = jaidaFactory.createTreeFactory().createTree(jaidaTreeFileName,
								      "xml",ITreeFactory.RECREATE);
        IHistogramFactory jaidaHistoFactory =
            jaidaFactory.createHistogramFactory(jaidaTree);

	// Making the 2D histogram
	String cascadeHistName = "cascadePosition";
	String injectionHistName = "injectionPosition";
	String endHistName = "endingPosition";
	double deltaDistance = 1.0e4; // 10,000 cm = 100 m
	double minX = gen2Xorigin - 2.0*InjectionGeometryUtils.getDefaultCylinderRadius();
	double maxX = gen2Xorigin + 2.0*InjectionGeometryUtils.getDefaultCylinderRadius();
	int dimensionX = (int )((maxX-minX)/deltaDistance + epsilon);
	double minY = gen2Yorigin - 2.0*InjectionGeometryUtils.getDefaultCylinderRadius();
	double maxY = gen2Yorigin + 2.0*InjectionGeometryUtils.getDefaultCylinderRadius();
	int dimensionY = (int )((maxY-minY)/deltaDistance + epsilon);
	double minZ = gen2Zorigin - 2.0*InjectionGeometryUtils.getDefaultCylinderHeight();
	double maxZ = gen2Zorigin + 2.0*InjectionGeometryUtils.getDefaultCylinderHeight();
	int dimensionZ = (int )((maxZ-minZ)/deltaDistance + epsilon);
        IHistogram3D cascadeXYZ =
	    jaidaHistoFactory.createHistogram3D(cascadeHistName,dimensionX,minX*cm2m,maxX*cm2m,
						dimensionY,minY*cm2m,maxY*cm2m,dimensionZ,minZ*cm2m,maxZ*cm2m);
        IHistogram3D injectionXYZ =
	    jaidaHistoFactory.createHistogram3D(injectionHistName,dimensionX,minX*cm2m,maxX*cm2m,
						dimensionY,minY*cm2m,maxY*cm2m,dimensionZ,minZ*cm2m,maxZ*cm2m);
        IHistogram3D endXYZ =
	    jaidaHistoFactory.createHistogram3D(endHistName,dimensionX,minX*cm2m,maxX*cm2m,
						dimensionY,minY*cm2m,maxY*cm2m,dimensionZ,minZ*cm2m,maxZ*cm2m);



	//
	// Event Generatiorn loop
	//
	int trial = 0;

	while(trial<numberOfEvents){

	    double nadirAngleInDeg = 180.0*rand.GetRandomDouble(); // [Deg]
	    double azimuthAngleInDeg = 360.0*rand.GetRandomDouble(); // [Deg]
	    System.out.format("event %d nadir %f azimuth %f\n",
			      trial,nadirAngleInDeg,azimuthAngleInDeg);

	    double nadirAngle = Math.toRadians(nadirAngleInDeg);
	    double azimuthAngle =  Math.toRadians(azimuthAngleInDeg);
	    J3UnitVector polarVectorInjection = new J3UnitVector(
						    Math.sin(nadirAngle)*Math.sin(azimuthAngle),
						    Math.sin(nadirAngle)*Math.cos(azimuthAngle),
						    Math.cos(nadirAngle));
	    EarthLocalCoordinate injectionCoordinate = 
		new EarthLocalCoordinate(polarVectorInjection,0.0); // origin (0,0,0)

	    //double injectionRadius = 
	    //	InjectionGeometryUtils.getInjectionRadius(Math.toRadians(nadirAngleInDeg));
	    //double injectionPointR = injectionRadius*Math.sqrt(rand.GetRandomDouble()); // [cm]
	    //double injectionPointAzimuthInRad = 2.0*Math.PI*rand.GetRandomDouble();
	    //double x_injectionCoord = injectionPointR*Math.sin(injectionPointAzimuthInRad);
	    //double y_injectionCoord = injectionPointR*Math.cos(injectionPointAzimuthInRad);
	    double injectionX = InjectionGeometryUtils.getDefaultCylinderRadius();
	    double x_injectionCoord = injectionX*(2.0*rand.GetRandomDouble()-1.0); // [cm]
	    double injectionY = InjectionGeometryUtils.getInjectionRadius(Math.toRadians(nadirAngleInDeg));
	    double y_injectionCoord = injectionY*(2.0*rand.GetRandomDouble()-1.0); // [cm]

	    J3Vector r_injectionCoord = new J3Vector(x_injectionCoord,y_injectionCoord,0.0);
	    J3Vector r_gen2 = injectionCoordinate.transformVectorToEarthCenter(r_injectionCoord);
	

	    // run JULIeT
	    generator.definePropagationGeometry(r_gen2.getX(),
						r_gen2.getY(),
						r_gen2.getZ(),
						nadirAngleInDeg,azimuthAngleInDeg);

	    generator.configurePropagationGeometry();


	    generator.runSingleEvent();


	    //
	    // now extract the primary track and secondary particle info
	    //
	    ListIterator particleIterator = generator.getParticleIterator();
	    ListIterator particleLocationIterator = generator.getLocationIce3Iterator();
	    ListIterator trackIterator = generator.getTrackParticleIterator();
	    ListIterator trackLocationIterator = generator.getTrackLocationIce3Iterator();

	    // start position (=injection position) in ice3 coordinate
	    J3Vector startPosition_ice3 = generator.wherePrimaryParticleStartsInIceCubeCoordinate();
	    J3Vector startPosition_ice3_stored_inTrack =
		(J3Vector )(trackLocationIterator.next());
	    J3Vector startPosition_gen2 = generator.wherePrimaryParticleStartsInGen2Coordinate();
	    double distanceFromEarthSurface = generator.getStartLocationAlongTheAxis();
	    J3Vector startPosition_center = generator.wherePrimaryParticleStartsInEarthCenterCoordinate();

	    injectionXYZ.fill(startPosition_ice3.getX()*cm2m,startPosition_ice3.getY()*cm2m,
			      startPosition_ice3.getZ()*cm2m);

	    //  end position in ice3 coordinate
	    J3Vector endPosition_ice3 = generator.wherePrimaryParticleEndsInIceCubeCoordinate();
	    J3Vector propagationAxis = J3Vector.subtract(endPosition_ice3,startPosition_ice3);
	    double propagationDistance = propagationAxis.getLength();

	    endXYZ.fill(endPosition_ice3.getX()*cm2m,endPosition_ice3.getY()*cm2m,
			      endPosition_ice3.getZ()*cm2m);

	    // secondary particles
	    while(particleLocationIterator.hasNext()){
		Particle particle = (Particle )(particleIterator.next());
		String particleName = particle.particleName(particle.getFlavor(),
							    particle.getDoublet());
		double cascade_energy = particle.getEnergy();
		J3Vector r = (J3Vector )(particleLocationIterator.next());
		double x = r.getX();
		double y = r.getY();
		double z = r.getZ();
		double w = cascade_energy/energy;
		cascadeXYZ.fill(x*cm2m,y*cm2m,z*cm2m,w);
	    }

	    trial++;
	}

	jaidaTree.commit();

    }



}

