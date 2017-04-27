package iceCube.uhe.analysis;

import java.io.*;
import java.util.*;
import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.event.*;
import numRecipes.*;
import geometry.*;

public class RunJulietMillipede {

    /** Define the type of run **/
    static int     numberOfRun = 1000000;
    static int     numberOfPreRun = 1000;
    static int   neutrinoWeight = 5000;
    private static boolean calEnergyDepositInsideMillipedeTrack = true; // set false if cascade hypothesis, i,e, nu-e for example.
    private static boolean normalize = true;
    private static boolean weight = false;
    private static boolean weightGZK = false;
    private static final double ln10 = Math.log(10.0);
    private static double energyBase = Math.pow(10.0,Particle.getLogEnergyMinimum());

    /**
       Calculate the total amount of the energy loss along the juliet track [startVertex endVertex].
       Needs to generator.runSingleEvent() in priori to call this method.
    */
    public static double getEnergyLossAlongTheTrack(JulietEventGenerator generator, J3Vector startVertex, J3Vector endVertex){
	ListIterator particleIterator = generator.getParticleIterator();
	ListIterator locationIterator = generator.getLocationIce3Iterator();
	double depositEnergyInThisEvent = 0.0;
	while(locationIterator.hasNext()){
	    Particle particle = (Particle )(particleIterator.next());
	    double cascadeEnergy = particle.getEnergy();
	    J3Vector r = (J3Vector )(locationIterator.next());
	    double z = r.getZ();
	    if((z-startVertex.getZ())*(z-endVertex.getZ())<0.0 || !calEnergyDepositInsideMillipedeTrack){
		                                             // within the track [startVertex endVertex]
		depositEnergyInThisEvent += cascadeEnergy;
	    }
	}

	return depositEnergyInThisEvent;
    }

    /**
       Calculate the amount of the cacsade energy loss along the juliet track [startVertex endVertex].
       Needs to generator.runSingleEvent() in priori to call this method. This method demands a given i3partocle event
       to constitute a single (electron or pion) particle in the particle list along the track, 
       which must be a single cascade event  by the JULIeT algorithm. If the event is not a single cascade,
       ruturn 1.0.
    */
    public static double getSingleCascadeEnergyLossAlongTheTrack(JulietEventGenerator generator, J3Vector startVertex, J3Vector endVertex){
	//ListIterator trackIterator = generator.getTrackParticleIterator();
	//double primaryEnergy = ((Particle )(trackIterator.next())).getEnergy();
	ListIterator particleIterator = generator.getParticleIterator();
	ListIterator locationIterator = generator.getLocationIce3Iterator();
	double depositEnergyInThisEvent = 0.0; int numberOfParticles = 0;
	while(particleIterator.hasNext()){
	    Particle particle = (Particle )(particleIterator.next());
	    double cascadeEnergy = particle.getEnergy();
	    int flavor = particle.getFlavor();
	    int doublet = particle.getDoublet();
	    J3Vector r = (J3Vector )(locationIterator.next());
	    double z = r.getZ();
	    if(doublet == 1){ // either electron or hadron
		//System.err.println(Particle.particleName(flavor,doublet)+" with energy of " + cascadeEnergy 
		//		   + " from primary energy of " + primaryEnergy);
		if((z-startVertex.getZ())*(z-endVertex.getZ())<0.0 || !calEnergyDepositInsideMillipedeTrack){
		                                             // within the track [startVertex endVertex]
		    depositEnergyInThisEvent += cascadeEnergy;
		}
	    }
	    numberOfParticles++;
	    if(numberOfParticles > 1) break; // this is not a single cacade event. Return a tiny energy as dummy
	}
	return depositEnergyInThisEvent;
    }

    /**
       Calculate the deposited energy probability. Based upon the monopod reconstruction by Jakob
     */
    public static double getDepositEnergyProbability(I3Particle thisEvent, double mc_energy){
	double depositEnergy = thisEvent.getRecoEnergy();
	double dE = (depositEnergy-mc_energy)/mc_energy;
	double dEmean = 0.0; double dEsigma = 0.1; // generic resolution.
	int eventNumber = thisEvent.getIceCubeData().getEventNumber();
	if(eventNumber == 118545){
	    //dEmean = 0.0597725;
	    //dEsigma = 0.0144999;
	    dEmean = 0.0;
	    //dEsigma = 0.35;
	    dEsigma = 0.15;
	    //System.err.println(" 2011 Aug event");
	}else if(eventNumber == 119316){
	    //dEmean = 0.0779754;
	    //dEsigma = 0.0192319;
	    dEmean = 0.0;
	    //dEsigma = 0.35;
	    dEsigma = 0.15;
	    //System.err.println(" 2012 Jan event");
	}
	double prob = SpecialFunctions.gauss(dEmean,dEsigma,dE);
	return prob;
    }


    /** Main method. In order to run JULIET, JulietEventGenerator object 
        and DataOutPutStream are generated. */
    public static void main(String[] args) throws IOException {

	String i3dataFileName = null;
	String i3dataOutName = null;
	String i3surfacedataOutName = null;
	NeutrinoFlux gzkFlux = null;
	double baseFlux = 1.0;
	boolean setDifferentParticleSpieces = false;
	boolean singleCascadeStudy = false;
	boolean shiftRecoEnergyUp = false;  // The systematic error consideration
	boolean shiftRecoEnergyDown = false;  // The systematic error consideration
	double systematicErrorUpper = 1.175;  // The +17.5 % sysmematics
	double systematicErrorLower = 0.825;  // The -17.5 % sysmematics
	int flavorToSet = 1;
	int doubletToSet = 1;

        if(args.length<3){
            System.out.println("Usage: RunJulietMillepede input-i3particle-file-name output-i3particle-name output-i3surface-name (flavor doublet bracket-vertx?(no 0))");
            System.exit(0);
        }else{
             i3dataFileName = args[0];
             i3dataOutName = args[1];
             i3surfacedataOutName = args[2];
             if(args.length >3){
	 setDifferentParticleSpieces = true;
	 flavorToSet = Integer.valueOf(args[3]).intValue();
	 doubletToSet = Integer.valueOf(args[4]).intValue();
	 System.err.println(" will be setting flavor " + flavorToSet + " doublet " + doubletToSet);
             }
             if(args.length >5){
		 if(Integer.valueOf(args[5]).intValue()==0){
		     calEnergyDepositInsideMillipedeTrack=false;
		     System.err.println(" calculate energy loss no matter where the vertecis are");
		 }
	     }
        }

        // Interactive session
        DataInputStream input = new DataInputStream(System.in);
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
        String buffer;
	System.out.print("Normalize? [yes(1)/no(0)] ->");
	buffer   = d.readLine();
        if(Integer.valueOf(buffer).intValue()==0){
	    normalize = false;
	    System.err.println(" Won't normalize to unity");
	}
	System.out.print("Weight to E^-2? [yes(1)/no(0)] ->");
	buffer   = d.readLine();
        if(Integer.valueOf(buffer).intValue()==1){
	    weight = true;
	    System.err.println(" weight to E^-2");
	}
	if(!weight){
	    System.out.print("Weight to GZK? [yes(1)/no(0)] ->");
	    buffer   = d.readLine();
	    if(Integer.valueOf(buffer).intValue()==1){
		weightGZK = true;
		gzkFlux = new NeutrinoFlux(12); // Fermi Best
		System.err.println(" weight to GZK");
	    }
	}

	System.out.print("run the special single cascade study? [yes(1)/no(0)] ->");
	buffer   = d.readLine();
	if(Integer.valueOf(buffer).intValue()==1){
	    singleCascadeStudy = true;
	    System.err.println(" single cascade study");
	}

	System.out.print("shift the energy by + 17.5%? [yes(1)/no(0)] ->");
	buffer   = d.readLine();
	if(Integer.valueOf(buffer).intValue()==1){
	    shiftRecoEnergyUp = true;
	    System.err.println(" shift UP the reco cascade study");
	}

	if(!shiftRecoEnergyUp) System.out.print("shift the energy by - 17.5%? [yes(1)/no(0)] ->");
	buffer   = d.readLine();
	if(Integer.valueOf(buffer).intValue()==1){
	    shiftRecoEnergyDown = true;
	    System.err.println(" shift DOWN the reco cascade study");
	}
	//
	// interaction IDs 
	//
	int flavorID=1, doubletID=1, mediumID=0; // muon ice
	int doCC=0, doNC=0, doMuBrems=1, doTauBrems=1; 
	int doMuKnock=0, doTauKnock=0, doMu2e=1, doTau2e=1;
	int doMu2mu=0, doTau2mu=0, doMu2tau=0, doTau2tau=0;
	int doMuPN=1, doTauPN=1, doMuDecay=1, doTauDecay=1;
	int posID=2, doGR=0; // start position 884m away from i3 center

	double energy = 1.0e6; // 10^7 GeV (dummy value)

	// generate JulietEventGenerator object
	JulietEventGenerator generator = null;
	if(!setDifferentParticleSpieces){
                    generator = 
	    new JulietEventGenerator(flavorID, doubletID, energy, mediumID, 
				     doCC, doNC, doMuBrems, doTauBrems, 
				     doMuKnock, doTauKnock, doMu2e, doTau2e,
				     doMu2mu, doTau2mu, doMu2tau, doTau2tau,
				     doMuPN, doTauPN, 
				     doGR,
				     doMuDecay, doTauDecay, posID);
	}else{
	    generator = new JulietEventGenerator();
	}

	// make output stream
	DataOutputStream out = new DataOutputStream(new FileOutputStream(i3dataOutName));
	DataOutputStream outsurface = new DataOutputStream(new FileOutputStream(i3surfacedataOutName));


	// input stream to read I3Particle
	InputStream in = ClassLoader.getSystemResourceAsStream(i3dataFileName);

	// PropagatingMatrixFactory object
	PropagationMatrixFactory matrix = new PropagationMatrixFactory();

	//
	// Now reading i3particle and run MC
	//
	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in)) !=null){

	    if(setDifferentParticleSpieces){
		iceParticle.setFlavorAndDoublet(flavorToSet,doubletToSet);
	    }
	    int flavor = iceParticle.getFlavor();
	    int doublet = iceParticle.getDoublet();
	    if(weightGZK) baseFlux = gzkFlux.getDFDLogEwzOsci(Particle.getLogEnergyMinimum(),flavor+1);

	    iceParticle.switchToReco();

	    J3Line particleAxis_ice3 = I3ParticleWrapper.getParticleAxisInIceCubeCoordinate(iceParticle);
	    J3Vector startVertex = particleAxis_ice3.getR0();
	    J3Vector endVertex = new J3Vector(particleAxis_ice3.getX(),particleAxis_ice3.getY(),particleAxis_ice3.getZ());

	    double depositEnergy = iceParticle.getRecoEnergy();
	    if(shiftRecoEnergyUp){ // The sysmetic shift
		depositEnergy = depositEnergy*systematicErrorUpper;
		iceParticle.putRecoEnergy(depositEnergy);
	    }else if(shiftRecoEnergyDown){
		depositEnergy = depositEnergy*systematicErrorLower;
		iceParticle.putRecoEnergy(depositEnergy);
	    }
		
	    double depositEnergyMax = 1.2*depositEnergy;
	    double depositEnergyMin = 0.8*depositEnergy;

	    generator.definePropagatingParticle(flavor,doublet,energy,particleAxis_ice3);
	    generator.configurePropagationGeometry();
	    if(doublet ==0 && flavor <=2){ // this is a neurtino in-ice particle
		generator.setNeutrinoInteractionWeight(neutrinoWeight);
	    }

	    double distanceFromEarthSurface = generator.getStartLocationAlongTheAxis();
	    iceParticle.setDistanceFromEarthSurfaceToIceCube(distanceFromEarthSurface);

	    System.err.format(" flavor(%d) doublet(%d) energy(%e) distance(%e)\n",
			      flavor,doublet,depositEnergy,distanceFromEarthSurface);
	    System.err.format(" start vertex x(%e) y(%e) z(%e)\n",startVertex.getX(),startVertex.getY(),startVertex.getZ());
	    System.err.format(" end vertex x(%e) y(%e) z(%e)\n",endVertex.getX(),endVertex.getY(),endVertex.getZ());

	    // pre run to define the in-ice energy rannge
	    int iStart = 0;
	    int iEnd =0;
	    boolean inTheRange = false;
	    for(int iLogE=0; iLogE<Particle.getDimensionOfLogEnergyMatrix(); iLogE++){
		double logPrimaryEnergy = 
		    Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
		double primaryEnergy = Math.pow(10.0,logPrimaryEnergy);

		int bingo = 0;
		for(int run=0;run<numberOfPreRun;run++){
		    // Run the Juliet
		    generator.definePropagatingParticle(flavor,doublet,primaryEnergy);
		    generator.runSingleEvent();
		    //generator.getListedEvents(new DataOutputStream(System.out));

		    // Extract the energy loss
		    double depositEnergyInThisEvent;
		    if(!singleCascadeStudy){
			depositEnergyInThisEvent = getEnergyLossAlongTheTrack(generator,startVertex,endVertex);
		    }else{
			depositEnergyInThisEvent = getSingleCascadeEnergyLossAlongTheTrack(generator,startVertex,endVertex);
		    }
		    //System.err.format(" run(%d) deposited energy(%e)\n",run, depositEnergyInThisEvent);

		    if(depositEnergyInThisEvent>=depositEnergyMin && depositEnergyMax>=depositEnergyInThisEvent){
			bingo++;
		    }
		}
		//System.err.format(" log(in-ice energy)=%e bingo times(%d)\n",logPrimaryEnergy,bingo);
		if((!inTheRange)&&(bingo>0)){
		    iStart=iLogE; inTheRange = true;
		}
		if((inTheRange)&&(bingo>0)){
		    iEnd=iLogE;
		}


	    }

	    iStart -= 50; if(iStart<0) iStart = 0;
	    iEnd += 150; if(iEnd>=Particle.getDimensionOfLogEnergyMatrix()) iEnd = Particle.getDimensionOfLogEnergyMatrix();

	    // Main run
	    int dimension = iEnd-iStart;
	    double[] bingo = new double[dimension]; int i=0;
	    double[] numOfNotSingleCascadeEvents = new double[dimension];
	    for(int iLogE=iStart; iLogE<iEnd; iLogE++){
		double logPrimaryEnergy = 
		    Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
		double primaryEnergy = Math.pow(10.0,logPrimaryEnergy);

		bingo[i]=0.0;
		if(singleCascadeStudy) numOfNotSingleCascadeEvents[i] = 0.0;
		for(int run=0;run<numberOfRun;run++){
		    // Run the Juliet
		    generator.definePropagatingParticle(flavor,doublet,primaryEnergy);
		    generator.runSingleEvent();
		    //generator.getListedEvents(new DataOutputStream(System.out));

		    // Extract the energy loss
		    double depositEnergyInThisEvent;
		    if(!singleCascadeStudy){
			depositEnergyInThisEvent = getEnergyLossAlongTheTrack(generator,startVertex,endVertex);
			bingo[i] += getDepositEnergyProbability(iceParticle, depositEnergyInThisEvent)*Particle.getDeltaLogEnergy();
		    }else{
			depositEnergyInThisEvent = getSingleCascadeEnergyLossAlongTheTrack(generator,startVertex,endVertex);
			if(depositEnergyInThisEvent>1.0){ // ok, this is a single cacade event
			    bingo[i] += getDepositEnergyProbability(iceParticle, depositEnergyInThisEvent)*Particle.getDeltaLogEnergy();
			}else if (depositEnergyInThisEvent == 1.0){
			    numOfNotSingleCascadeEvents[i] += 1.0;
			}
		    }
		    //System.err.format(" run(%d) deposited energy(%e)\n",run, depositEnergyInThisEvent);


		}

		System.err.format(" log(in-ice energy)=%5.3f bingo times(%8.4e)\n",logPrimaryEnergy,bingo[i]);
		if(singleCascadeStudy) System.err.format(" number of events with more than one cascade(%5.1f)\n",numOfNotSingleCascadeEvents[i]);
		if(weight){
		    bingo[i] = bingo[i]*(energyBase/primaryEnergy); // rewrighting to E^-2
		    if(singleCascadeStudy) numOfNotSingleCascadeEvents[i] = 
					       numOfNotSingleCascadeEvents[i]*(energyBase/primaryEnergy);
		}
		else if(weightGZK){
		    double  flux = gzkFlux.getDFDLogEwzOsci(logPrimaryEnergy,flavor+1);
		    bingo[i] = bingo[i]*flux/baseFlux; // rewrighting to GZK
		    if(singleCascadeStudy) numOfNotSingleCascadeEvents[i] = numOfNotSingleCascadeEvents[i]*flux/baseFlux;
		}
		i++;
	    }
	    // normalization 
	    double norm = 0.0;
	    for(i=0;i<dimension;i++){
		norm += bingo[i];
		if(singleCascadeStudy) norm += numOfNotSingleCascadeEvents[i];
	    }

	    // store the probability distribution in the array in the particle class
	    i=0;
	    iceParticle.generateLogEnergyMatrix();
	    for(int iLogE=iStart; iLogE<iEnd; iLogE++){
		double prob = bingo[i]/norm;
		if(!normalize) prob = bingo[i];
		iceParticle.putLogEnergyMatrix(iLogE,prob);
		i++;
	    }

	    // confirm the written probability
	    //for(int iLogE=0; iLogE<Particle.getDimensionOfLogEnergyMatrix(); iLogE++){
	    //	double logPrimaryEnergy = 
	    //	    Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
	    //	double prob = iceParticle.getLogEnergyMatrix(iLogE);
	    //	System.err.format(" log(in-ice energy)=%e probability(%7.5f)\n",logPrimaryEnergy,prob);
	    //}


	    //
	    //
	    // Read Propagation Matrix and fill the surface neutrino energy distribution
	    //
	    //
	    iceParticle.switchToReco();
	    J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	    double distance = 
		iceParticle.getDistanceFromEarthSurfaceToIceCube();
	    boolean isDownWard = true;
	    double cosZenith = -n_ice.getZ(); // Reversed vector
	    if(cosZenith<0.0) isDownWard = false; // this is an upward event
	    // read matrix
	    String matrixFileName = 
		I3ParticleWeightFiller.getMatrixFileName(n_ice,distance,isDownWard);
	    System.err.println(" --matrix filename " + matrixFileName);
	    DataInputStream inMtx = 
		new DataInputStream(new FileInputStream(matrixFileName));
	    matrix.readMatrix(inMtx);
	    inMtx.close( );

	    // Generate I3Particle to store the energy distribution of surface neutrinos
	    I3Particle surfaceParticle = new I3Particle(iceParticle.getFlavor(),iceParticle.getDoublet());
	    surfaceParticle.setDistanceFromEarthSurfaceToIceCube(distance);
            double iniceEnergy = iceParticle.getEnergy();
 	    surfaceParticle.putEnergy(iniceEnergy);
            if(iniceEnergy>0.0){
                double logEnergy = Math.log(iniceEnergy)/ln10;
                surfaceParticle.putLogEnergy(logEnergy);
	    }
	    surfaceParticle.putRecoEnergy(iceParticle.getRecoEnergy());
	    surfaceParticle.switchToReco();
	    surfaceParticle.setParticleAxisInIceCubeCoordinate(particleAxis_ice3);
	    surfaceParticle.transformParticleAxisToEarthCenterCoordinate();
	    I3Particle.I3Data i3Data = iceParticle.getIceCubeData();
	    // Fills the IceCube data
	    surfaceParticle.getIceCubeData().setEventNumber(i3Data.getEventNumber());
	    surfaceParticle.getIceCubeData().setNpeATWD(i3Data.getNpeATWD());
	    surfaceParticle.getIceCubeData().setNpeFADC(i3Data.getNpeFADC());
	    surfaceParticle.getIceCubeData().setBestNpe(i3Data.getBestNpe());
	    surfaceParticle.getIceCubeData().setNDOMsATWD(i3Data.getNDOMsATWD());
	    surfaceParticle.getIceCubeData().setNDOMsFADC(i3Data.getNDOMsFADC());
	    surfaceParticle.getIceCubeData().setNDOMsLaunch(i3Data.getNDOMsLaunch());

	    surfaceParticle.generateLogEnergyMatrix();

	    // now fill the surface energy distribution
	    double epsilon = 1.0e-4; // round-off margin for binning
	    double logSurfaceEnergy = Particle.getLogEnergyMinimum() + epsilon;
	    double logSurfaceMaxEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	    int iLogE = 0;
	    while(logSurfaceEnergy < logSurfaceMaxEnergy){ // E < 10^12 GeV
		double logInIceEnergy = Particle.getLogEnergyMinimum() + epsilon;
		int jLogE = 0;
		double prob = 0.0;
		while(logInIceEnergy <= logSurfaceEnergy){
		    double flux = matrix.getDF(0,0,logSurfaceEnergy, flavor,doublet,logInIceEnergy);// nu-e
		    flux += matrix.getDF(1,0,logSurfaceEnergy, flavor,doublet,logInIceEnergy); // nu-mu
		    flux += matrix.getDF(2,0,logSurfaceEnergy, flavor,doublet,logInIceEnergy); // nu-tau
		    prob += flux*iceParticle.getLogEnergyMatrix(jLogE);
		    jLogE++;
		    logInIceEnergy += Particle.getDeltaLogEnergy();
		}

		surfaceParticle.putLogEnergyMatrix(iLogE,prob);

		iLogE++;
		logSurfaceEnergy += Particle.getDeltaLogEnergy();
	    }


	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	    I3ParticleOutputStream.outputI3Particle(surfaceParticle, outsurface);
	}
    	out.close();
    	outsurface.close();
    }
}
