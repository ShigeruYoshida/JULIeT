package iceCube.uhe.event;

import numRecipes.*;
import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.particles.*;

import java.io.*;
import java.util.*;

public class RunJuliet4Gen2 {

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

    // Main method. In order to run JULIET 
    public static void main(String[] args) throws IOException {

	int flavorID = 1;
	int doubletID = 1;
	double energy = 1.0e9; // [GeV]
	int numberOfEvents = 10;
	boolean simulateARA = false;

	if(args.length<4){
            System.out.println("Usage: RunJulet2Gen2 flavorID doubletID number-of-events energy [GeV]");
            System.exit(0);
        }else{
            flavorID = Integer.valueOf(args[0]).intValue();
            doubletID = Integer.valueOf(args[1]).intValue();
            numberOfEvents = Integer.valueOf(args[2]).intValue();
            energy = Double.valueOf(args[3]).doubleValue();
	    if(args.length==5) simulateARA = true;
        }
	System.err.format("(flavor doublet) = (%d %d) Energy=%e [GeV]\n",
			  flavorID, doubletID, energy);

	RandomGenerator rand = new RandomGenerator();

	// generate RunManager object
	JulietEventGenerator4Gen2.neutrinoCSHERAZeus = true;
	if(simulateARA) JulietEventGenerator4Gen2.setARADimensionToDefaults();
	JulietEventGenerator4Gen2 generator = 
	    new  JulietEventGenerator4Gen2(flavorID, doubletID, energy, mediumID,
					   doCC, doNC, doMuBrems, doTauBrems,
					   doMuKnock, doTauKnock, doMu2e, doTau2e,
					   doMu2mu, doTau2mu, doMu2tau, doTau2tau,
					   doMuPN, doTauPN, doGR, doMuDecay, doTauDecay,posID);

	//
	// Event Generatorn loop
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
	    System.out.format("event %d nx (%f) ny(%f) nz(%f)\n", trial,polarVectorInjection.getX(),
			      polarVectorInjection.getY(),polarVectorInjection.getZ());
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
	    System.out.format("event %d injected radius %f\n",trial,r_injectionCoord.getLength());
	    System.out.format("event %d injected radius %f\n",trial,r_gen2.getLength());
	    System.out.format("event %d primary injected position %f %f %f\n",
			      trial,r_gen2.getX(),r_gen2.getY(),r_gen2.getZ());
	

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

	    System.out.format("event %d start %f %f %f in the earth center\n",
			      trial,startPosition_center.getX(),
			      startPosition_center.getY(),startPosition_center.getZ());
	    System.out.format("event %d start %f %f %f\n",trial,startPosition_ice3.getX(),
			      startPosition_ice3.getY(),startPosition_ice3.getZ());
	    System.out.format("event %d track %f %f %f\n",trial,
			      startPosition_ice3_stored_inTrack.getX(),
			      startPosition_ice3_stored_inTrack.getY(),
			      startPosition_ice3_stored_inTrack.getZ());
	    System.out.format("event %d Distance from the earth surface %f\n",trial,distanceFromEarthSurface);
	    System.out.format("event %d Distance from the gen2 center  %f\n",trial,startPosition_gen2.getLength());

	    //  end position in ice3 coordinate
	    J3Vector endPosition_ice3 = generator.wherePrimaryParticleEndsInIceCubeCoordinate();
	    J3Vector propagationAxis = J3Vector.subtract(endPosition_ice3,startPosition_ice3);
	    double propagationDistance = propagationAxis.getLength();
	    System.out.format("event %d end %f %f %f %f\n",trial,endPosition_ice3.getX(),
			      endPosition_ice3.getY(),endPosition_ice3.getZ(),propagationDistance);

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
		System.out.println("event " + trial + " secondary " + particleName + " " + 
				   cascade_energy + " [GeV] " +
				   x + " " +  y + " " + z);
	    }

	    trial++;
	}

    }

}

