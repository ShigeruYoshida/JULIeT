package iceCube.uhe.interactions;

import numRecipes.*;
import geometry.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.event.*;

import java.io.*;
import java.util.*;

public class RunJulietAndCalcLikelihood {

    private final static double ln10 = Math.log(10.0);
    /** Define the type of run **/
    static int mediumID = 0; //ice
    static int doCC = 0;
    static int doNC = 0;
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
    static double threshold = InteractionsBaseLikelihood.logEnergyProducedMinimum;
    static boolean inelasticityBase = false;
    static double inice_initial_energy = 1.0e9; // [GeV]

    // Main method. In order to run JULIET 
    public static void main(String[] args) throws IOException {

	int flavorID = 1;
	int doubletID = 1;
	int numberOfEvents = 10;
	ParticlePoint point = null;
	List llhList = null;
	String fileName = null;

	if(args.length!=7){
            System.out.println("Usage: RunJulietAndCalcLikelihood flavorID doubletID number-of-events energy [GeV] log(threshold_e [GeV]) or log10(yThreshls)  list-file-name inelasticityBase(yes 1 no 0)");
            System.exit(0);
        }else{
            flavorID = Integer.valueOf(args[0]).intValue();
            doubletID = Integer.valueOf(args[1]).intValue();
            numberOfEvents = Integer.valueOf(args[2]).intValue();
            inice_initial_energy = Double.valueOf(args[3]).doubleValue();
            threshold = Double.valueOf(args[4]).doubleValue();
	    fileName = args[5];
	    if(Integer.valueOf(args[6]).intValue()==1) inelasticityBase = true;
        }
	if(!inelasticityBase) System.err.format("(flavor doublet) = (%d %d) Energy=%e [GeV] log(E threshold [GeV])=%6.3f\n",
						flavorID, doubletID, inice_initial_energy,threshold);
	else System.err.format("(flavor doublet) = (%d %d) Energy=%e [GeV] log(Y threshold)=%6.3f\n",
						flavorID, doubletID, inice_initial_energy,threshold);
	System.err.println("list file name " + fileName);

	RandomGenerator rand = new RandomGenerator();

	// generate RunManager object
	JulietEventGenerator.neutrinoCSHERAZeus = true;
	JulietEventGenerator generator = 
	        new  JulietEventGenerator(flavorID, doubletID, inice_initial_energy, mediumID,
	    				   doCC, doNC, doMuBrems, doTauBrems,
	    				   doMuKnock, doTauKnock, doMu2e, doTau2e,
	    				   doMu2mu, doTau2mu, doMu2tau, doTau2tau,
	    				   doMuPN, doTauPN, doGR, doMuDecay, doTauDecay,posID);
	//new  JulietEventGenerator(flavorID, doubletID, inice_initial_energy, mediumID,
	//				   0, 0, 0, 0,
	//				   0, 0, 0, 0,
	//				   1, 0, 0, 0,
	//				   0, 0, 0, doMuDecay, doTauDecay,posID);

	//
	// configure geometroy of propagating particle trajectory
	//
	
	double nadirAngleInDeg = 180.0*rand.GetRandomDouble(); // [Deg]
	double azimuthAngleInDeg = 360.0*rand.GetRandomDouble(); // [Deg]
	System.out.format("nadir %f azimuth %f\n",
			  nadirAngleInDeg,azimuthAngleInDeg);

	double nadirAngle = Math.toRadians(nadirAngleInDeg);
	double azimuthAngle =  Math.toRadians(azimuthAngleInDeg);
	J3UnitVector polarVectorInjection = new J3UnitVector(
					     Math.sin(nadirAngle)*Math.sin(azimuthAngle),
					     Math.sin(nadirAngle)*Math.cos(azimuthAngle),
					     Math.cos(nadirAngle));
	System.out.format("nx (%f) ny(%f) nz(%f)\n", polarVectorInjection.getX(),
			  polarVectorInjection.getY(),polarVectorInjection.getZ());
	EarthLocalCoordinate injectionCoordinate = 
	    new EarthLocalCoordinate(polarVectorInjection,0.0); // origin (0,0,0)

	double injectionRadius = IceCubeCoordinate.elevation;  // 8.8e4 cm
	double radius_injectionCoord = injectionRadius*rand.GetRandomDouble(); // [cm]
	double phi_injectionCoord = 2.0*Math.PI*rand.GetRandomDouble(); // [rad]
	double x_injectionCoord = radius_injectionCoord*Math.sin(phi_injectionCoord);
	double y_injectionCoord = radius_injectionCoord*Math.cos(phi_injectionCoord);

	J3Vector r_injectionCoord = new J3Vector(x_injectionCoord,y_injectionCoord,0.0);
	J3Vector r_ice3 = injectionCoordinate.transformVectorToEarthCenter(r_injectionCoord);
	System.out.format("injected radius %f\n",r_injectionCoord.getLength());
	System.out.format("injected radius %f\n",r_ice3.getLength());
	System.out.format("primary injected position %f %f %f\n",
			      r_ice3.getX(),r_ice3.getY(),r_ice3.getZ());
	

	generator.definePropagationGeometry(r_ice3.getX(),r_ice3.getY(), r_ice3.getZ(),
					    nadirAngleInDeg,azimuthAngleInDeg);

	generator.configurePropagationGeometry();

	// start position (=injection position) in ice3 coordinate
	J3Vector startPosition_ice3 = generator.wherePrimaryParticleStartsInIceCubeCoordinate();
	double distanceFromEarthSurface = generator.getStartLocationAlongTheAxis();
	J3Vector startPosition_center = generator.wherePrimaryParticleStartsInEarthCenterCoordinate();

	System.out.format("start %f %f %f in the earth center\n",
			  startPosition_center.getX(),
			  startPosition_center.getY(),startPosition_center.getZ());
	System.out.format("start %f %f %f\n",startPosition_ice3.getX(),
			  startPosition_ice3.getY(),startPosition_ice3.getZ());
	System.out.format("Distance from the earth surface %f\n",distanceFromEarthSurface);
	System.out.format("Distance from the ice3 center  %f\n",startPosition_ice3.getLength());

	double nadirAngleAtEntrance = nadirAngle;
	if(nadirAngleAtEntrance>0.5*Math.PI) nadirAngleAtEntrance = Math.PI-nadirAngle;
	point = new ParticlePoint(0.0,nadirAngleAtEntrance,mediumID);
	double localDensity = point.getMediumDensity();

	
	//
	// Internaction Likelihood object
	//
	InteractionsLikelihood intLikelihood =  new InteractionsLikelihood(threshold,inelasticityBase);
	intLikelihood.setPropagatingParticleIDs(flavorID,doubletID);
	intLikelihood.configureInteractionsForLikelihoodCalculation(
			    doCC, doNC, doMuBrems, doTauBrems, 
                                 doMuKnock, doTauKnock, doMu2e, doTau2e,
                                 doMu2mu, doTau2mu, doMu2tau, doTau2tau, doMuPN, doTauPN);
	//
	// Event Generatorn loop
	//
	int trial = 0;
	llhList = new LinkedList();
	while(trial<numberOfEvents){


	    generator.runSingleEvent();


	    //
	    // now extract the primary track and secondary particle info
	    //
	    ListIterator particleIterator = generator.getParticleIterator();
	    ListIterator particleLocationIterator = generator.getLocationIce3Iterator();
	    ListIterator trackIterator = generator.getTrackParticleIterator();
	    ListIterator trackLocationIterator = generator.getTrackLocationIce3Iterator();

	    J3Vector startPosition_ice3_stored_inTrack =
		(J3Vector )(trackLocationIterator.next());
	    System.out.format("track %f %f %f\n",
			      startPosition_ice3_stored_inTrack.getX(),
			      startPosition_ice3_stored_inTrack.getY(),
			      startPosition_ice3_stored_inTrack.getZ());

	    //  end position in ice3 coordinate
	    J3Vector endPosition_ice3 = generator.wherePrimaryParticleEndsInIceCubeCoordinate();
	    J3Vector propagationAxis = J3Vector.subtract(endPosition_ice3,startPosition_ice3);
	    double propagationDistance = propagationAxis.getLength();
	    System.out.format("event %d end %f %f %f %f\n",trial,endPosition_ice3.getX(),
			      endPosition_ice3.getY(),endPosition_ice3.getZ(),propagationDistance);

	    double logLikelihoodValue = 
		InteractionsLikelihoodBuilder.calculateInteractionsLikelihoodFromJulietLists(intLikelihood,
					       trackIterator, particleIterator, particleLocationIterator,
 					     startPosition_ice3, localDensity);
	    if(!Double.isInfinite(logLikelihoodValue)){
		System.out.format("event %d llh=%e\n",trial,logLikelihoodValue);
		Double llhObj = Double.valueOf(logLikelihoodValue);
		llhList.add(llhObj);
	    }

	    trial++;
	}

	Collections.sort(llhList);
	// make output stream
        FileOutputStream out = new FileOutputStream(fileName);
	ObjectOutputStream objectOut = new ObjectOutputStream(out);
	objectOut.writeObject(llhList);
	objectOut.flush();
	out.close();

    }

}

