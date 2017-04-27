package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.event.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.points.*;
import numRecipes.*;

import java.util.*;
import java.io.*;

/**
   <pre>

   This is a run manager for the Event class. It initializes all the obejcts
   and parameters necessary for events and then run the Event objects.

   </pre>
*/

public class RunManager {
    
    /** Common logarithm */
    static final double ln10 = Math.log(10.0);

    /** Round off error */
    static final double epsilon = 10e-5;

    /** Neutrino factor */
    static final int neutrinoFactor = InteractionsBase.neutrinoFactor;

    /** Matrices for Cascade Energy */
    double[][] emgCascadeMtx;
    double[][] hadronCascadeMtx; 
    double[][] totalCascadeMtx;

    double[][] oneEmgCascadeMtx;
    double[][] oneHadronCascadeMtx;
    double[][] oneTotalCascadeMtx;

    /** Number of cascades */
    int   numberOfEmgCascade;
    int   numberOfHadronCascade;
    int   numberOfTotalCascade;
    int[] numberOfCascade;

    /** Array for storing cascade energy of each interaction */
    double[] energyOfCascade;

    /** ParticlePoint object to define the trajectory and the medium
	in the propagation */
    ParticlePoint point;
    int           materialNumber;

    /** Propagating Particle object */
    Particle   propParticle;

    /** Initial value of propagating particle */
    double     primaryEnergy;
    static int primaryFlavor;
    static int primaryDoublet;

    /** Total Distance of the particle path. Default value: 1km = 1.0e5 cm */
    static double trajectoryLength = 1.0e5;

    /** Started location of incident particle */
    double startLocation;

    /** Array of the MonteCarloBase objects */
    MonteCarloBase[] mcBases;

    /** Array of the InteractionsMatrix objects **/
    InteractionsMatrix[] intMtx;

    /** Palameter for making MonteCarloBase[] **/
    int     numOfDecay;
    boolean mudecay;
    boolean taudecay;
    int     tauDecayFlag;
    int     electronBaseFlag;
    int     mTauDecay;

    /** DecayMatrix Objects **/
    MuDecayBase  muDecayBase;
    TauDecayBase tauDecayBase;

    /** Event Object **/
    Event event;

    /** Index of bin of the log primary energy **/
    int primaryiLogE;

    /** Directory path for the dumped InteractionsMatrix objects */
    String interactionsMatrixDirectoryInIce  = "iceCube/uhe/interactions/ice/";
    String interactionsMatrixDirectoryInRock = "iceCube/uhe/interactions/rock/";
    String interactionsMatrixDirectory;

    /** Random Generator */
    RandomGenerator rand;

    /** Define the type of run **/
    static int     numberOfEvent;
    static boolean integral        = false;
    static boolean differential    = false;
    static boolean eachInteraction = false;
    static String  fileName; 
    static int     typeOfEvent;

    /** Dimension of LogEnergyMatrix */
    int dim = Particle.getDimensionOfLogEnergyMatrix();

    /**
      Constructor. Here all the obejcts including individual InteractionsBase objects
      (which implies it reads out the InteractionsMatrix objects dumped in the directory
      classes/iceCube/uhe/interactions/) and  Mu/TauDecayBase objects are generated. Each
      object is generated in an interactive way to an user.
    */
    public RunManager() throws IOException{
	
	DataInputStream input = new DataInputStream(System.in);
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer;

	System.out.print("Number of Event ->");	
	buffer        = d.readLine();
	numberOfEvent = Integer.valueOf(buffer).intValue();

	System.out.print(
	"Type of Event (0:single energy (full data) 1:single energy 2:multi energy (Matrix ) ) ->");	
	buffer      = d.readLine();
	typeOfEvent = Integer.valueOf(buffer).intValue();

	if(typeOfEvent==1) {
	    System.out.print("Each Interaction Data yes(1)/no(0)?->");
	    buffer = d.readLine();
	    if(Integer.valueOf(buffer).intValue()==1) eachInteraction = true;
	    System.out.print("Type Of Matrix (0:EnrgyDeposit/1km 1:EnergyDeposit/cascade) ->");
	    buffer = d.readLine();
	    if(Integer.valueOf(buffer).intValue()==0) integral = true;
	    else if(Integer.valueOf(buffer).intValue()==1) differential = true;
	}

	else if(typeOfEvent==2) {
	    System.out.print("Type Of Matrix (0:EnrgyDeposit/1km 1:EnergyDeposit/cascade) ->");
	    buffer = d.readLine();
	    if(Integer.valueOf(buffer).intValue()==0) integral = true;
	    else if(Integer.valueOf(buffer).intValue()==1) differential = true;
	}

	System.out.print("Name of File ->");	
	buffer   = d.readLine();
	fileName = String.valueOf(buffer).trim();

	System.out.print("Particle Flavor ->");
	buffer        = d.readLine();
	primaryFlavor = Integer.valueOf(buffer).intValue();

	System.out.print("Particle Doublet ->");
	buffer         = d.readLine();
	primaryDoublet = Integer.valueOf(buffer).intValue();

	if(typeOfEvent==0 || typeOfEvent==1) {
	    System.out.print("Particle Energy [GeV] ->");
	    buffer        = d.readLine();
	    primaryEnergy = Double.valueOf(buffer).doubleValue();
	}
	else{
	    primaryEnergy = 1.0e6;
	}

	// Kind Of Medium
	System.out.print("Medium in Propagation ice(0) rock(1)->");
	buffer         = d.readLine();
	materialNumber = Integer.valueOf(buffer).intValue();

	
	// Generate Random Generator
 	rand = new RandomGenerator();
	System.out.println("Random Generator has been generated");


	// Register Interactions and read the InteractionMatrix objects
	String[] fileName = new String[16];
	String matrixName = null;
	if(materialNumber==0)  interactionsMatrixDirectory = interactionsMatrixDirectoryInIce;
	else interactionsMatrixDirectory = interactionsMatrixDirectoryInRock;
	int n = 0; int ans = 0;

	System.out.print("Charged Current Interactions yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	electronBaseFlag = 0;
	if(ans == 1){
	    if(primaryFlavor==0) {matrixName = "ENeutrinoChargeMtx"; electronBaseFlag = 1;} 
	    else if(primaryFlavor==1) matrixName = "MuNeutrinoChargeMtx";
	    else if(primaryFlavor==2) matrixName = "TauNeutrinoChargeMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Neutral Current Interactions yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    if(primaryFlavor==0) matrixName      = "ENeutrinoNeutralMtx";
	    else if(primaryFlavor==1) matrixName = "MuNeutrinoNeutralMtx";
	    else if(primaryFlavor==2) matrixName = "TauNeutrinoNeutralMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon Bremsstrahlung yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "muBremsstrahlungMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau Bremsstrahlung yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauBremsstrahlungMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon Knock-on Electrons yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "muKnockOnElectronsMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau Knock-on Electrons yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauKnockOnElectronsMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon e+e- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "muToePairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau e+e- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauToePairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon mu+mu- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "muTomuPairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau mu+mu- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauTomuPairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon tau+tau- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "muTotauPairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau tau+tau- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauTotauPairCreationFitMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Muon Photo-nuclear interactions yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName = "muPhotoNuclearMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	System.out.print("Tau Photo-nuclear interactions yes(1)/no(0)?->");
	buffer = d.readLine();
	ans    = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    matrixName   = "tauPhotoNuclearMtx";
	    fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	}

	numOfDecay = 0;
	mudecay    = false;
	taudecay   = false;

	System.out.print("Mu decay yes(1)/no(0)?->");
	buffer = d.readLine();
	ans = Integer.valueOf(buffer).intValue();
	if(ans == 1){
	    Particle       mu         = new Particle(1,1,primaryEnergy);
	    MuDecayYMatrix muDecayMtx = new MuDecayYMatrix(mu); 
	    muDecayBase = new MuDecayBase(muDecayMtx);
	    n++; numOfDecay++; mudecay=true;
	}

	System.out.print("Tau Decay yes(1)/no(0)?->");
	buffer = d.readLine();
	ans = Integer.valueOf(buffer).intValue();
	tauDecayFlag = 0;
	if(ans == 1){
	    Particle        tau         = new Particle(2, 1, primaryEnergy);
	    TauDecayYMatrix tauDecayMtx = new TauDecayYMatrix(tau);
	    tauDecayBase = new TauDecayBase(tauDecayMtx);
	    n++; numOfDecay++; taudecay=true; tauDecayFlag=1;
	}

	// make an array for MonteCarloBase
	//mcBases = new MonteCarloBase[n - numOfDecay + electronBaseFlag + 2];
	mcBases = new MonteCarloBase[n  + electronBaseFlag];
	// regist interactionsMtx
	intMtx = new InteractionsMatrix[n - numOfDecay];
	//intMtx = new InteractionsMatrix[n];
	System.out.println("Registered Interaction Matrix....");
	for(int i=0; i<(n-numOfDecay) ;i++){       
	    System.out.println(fileName[i]);
	    //FileInputStream in = new FileInputStream(fileName[i]);
	    InputStream in = ClassLoader.getSystemResourceAsStream(fileName[i]);
	    intMtx[i]  = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    mcBases[i] = new InteractionsBase(intMtx[i]); 
	    in.close();
	}

	// register ElectronBase when primary particle is nu-e and charged current is choozen
	if(electronBaseFlag == 1) {
	    Particle  e = new Particle(0,1,primaryEnergy);
	    mcBases[n-numOfDecay] = new ElectronBase(e);
	    System.out.println("ElectronBase is registered");
	}

	// register decay matrix
	if(mudecay) {
	    mcBases[n-numOfDecay+electronBaseFlag] = muDecayBase;
	    System.out.println("MuDecayBase is registered");
	    if(taudecay) {
		mcBases[n-numOfDecay+1+electronBaseFlag] = tauDecayBase;
		mTauDecay = n - numOfDecay + 1+ electronBaseFlag;
		System.out.println("TauDecayBase is registered");
	    }
	}
	else if(taudecay) {
	    mcBases[n-numOfDecay+electronBaseFlag] = tauDecayBase;
	    mTauDecay = n - numOfDecay + electronBaseFlag;
	    System.out.println("TauDecayBase is registered");
	}

	d.close();

    }

    /** Set the location of primary particle */
    public void setStartLocation(double start) {
	startLocation = start;
    }

    /** Get the location of primary particle */
    public double getStartLocation() {
	return startLocation;
    }

    /** Method to run event for making textual full-data-file. This is mainly for
     debugging. All the profile of particle propagation is dumped in the file.
     **/
    public void runEventDump(DataOutputStream out) throws IOException {

	double sumRange      = 0.0; ////////
	double sum2Range     = 0.0; ///////
	int sumEvent      = 0; ///////
	double mean = 0.0;/////
	double sig = 0.0;/////


	for(int n=0; n<numberOfEvent; n++){
	    System.out.println(n);
	    propParticle = new Particle(primaryFlavor,primaryDoublet,primaryEnergy);
	    System.err.println(
            Particle.particleName(propParticle.getFlavor(),propParticle.getDoublet()) +
			   " has been generated with energy of " + propParticle.getEnergy() +
			   " [GeV]");    

	    point = new ParticlePoint(0.0, 0.0, materialNumber);
	    point.setIceRockBoundaryRadius(1.01*point.REarth);
	    System.out.println("Medium Density " +  point.getMediumDensity() + " [g/cm^2]");
	    
	    event = new Event(mcBases, propParticle, point);
	    System.out.println("Event object has been generated.");

	    double startlocation = point.getParticleLocation();
	    setStartLocation(startlocation);

	    System.out.println("Current location " + getStartLocation());
	    out.writeBytes("Current location " + getStartLocation());
	    out.write('\n');
	    System.out.println("Mass number in the medium " + event.getMassNumber());
	    out.writeBytes("Mass number in the medium " + event.getMassNumber());
	    out.write('\n');

	    double currentEnergy = 0.0; //////////
	    double currentX      = 0.0; //////////
	    double range         = 0.0; //////////
	    double primaryLogEnergy = Math.log(primaryEnergy)/ln10;/////////
	    

	    while(true){
		
		double transferedEnergy = 0.0;
		currentEnergy = event.propParticle.getEnergy(); //////////
		

		    double pathLength = event.getPhysicalPathLength(rand);
		    event.setStepDx(0.1*pathLength);
		    event.traceParticle(pathLength);
		    currentX = point.getParticleLocation()-getStartLocation(); //////////

		    if((point.getParticleLocation()-getStartLocation())>trajectoryLength){
			out.writeBytes("through");
		    break;   // particles passed through your DETECTOR
		    }else{
			transferedEnergy = event.collideNow(rand);

			

		    System.out.println("Colliding via " + 
				       event.interactionsNameInPlay() + " producing " + 
				       Particle.particleName(
					   event.getFlavorByInteractionsInPlay(),1));  
		    out.writeBytes("Colliding via " + 
				   event.interactionsNameInPlay() + " producing " + 
				   Particle.particleName(
				       event.getFlavorByInteractionsInPlay(),1));
		    out.write('\t');
		    
		    System.out.println("Particle Energy " + 
				       event.propParticle.getEnergy() + " [GeV]");
		    
		    out.writeBytes("Particle Energy " + 
				   event.propParticle.getEnergy() + " [GeV]");
		    
		    out.write('\t');
		    System.out.println("Transfered Energy " + transferedEnergy + " [GeV]");
		    out.writeBytes("Transfered Energy " + transferedEnergy + " [GeV]");
		    out.write('\t');
		    
		    System.out.println("Current distance the particle gained " + 
				       (point.getParticleLocation()-getStartLocation()) + " [cm]");
		    out.writeBytes("Current distance the particle gained " + 
				   (point.getParticleLocation()-getStartLocation()) + " [cm]");
		    out.write('\n');  
/*
		    //currentEnergy = event.propParticle.getEnergy(); //////////
		    //if(currentEnergy <= 1.0/Math.E*primaryEnergy+epsilon){ //////////
			    out.write('\n'); 
			    //range = currentX/1.0e2; // [m]*point.getMediumDensity();
			    System.out.println("currentX="+currentX+" range="+range);
			    //out.writeBytes("currentX="+currentX+" range="+range+" intRange="+(int)range);
			    //out.writeBytes((int)range + " ");
			    //out.write('\n');
			    //out.writeBytes("primaryLogEnergy =" + primaryLogEnergy);
			    out.write('\n');
			    sumEvent++;
			    sumRange += range;
			    sum2Range += range*range;
/*			    if(sumEvent>=12900){
				mean = sumRange/sumEvent;
				sig = Math.sqrt(sum2Range/sumEvent - mean*mean);
				out.writeBytes("mean=" + mean);
				out.writeBytes("  sigma=" + sig);
				out.write('\n');
			    }
			    break;
			}////////////////
*/
			out.write('\n');

			double logEnergy = event.propParticle.getLogEnergy(); // logEnergy after interaction
			if(logEnergy <= InteractionsBase.getLogEnergyProducedMinimum()) break;
			// Journey ends because all the primary energy has been lost
			
			if(event.mcBaseInPlay.getTypeOfInteraction() == 1 &&
			   event.mcBaseInPlay.getProducedFlavor()!= 1){ //Tau to Hadron/Electron decay
			    break;                                      // or Tau to Mu decay
			}
			
		    } //else trajectory
	    
			System.out.println();
	    }
	}
    }

    
    /** Method to run a single event with a monocromatic primary energy given by
	the constructor.**/

    public void runSingleEvent(){

	    propParticle = new Particle(primaryFlavor,primaryDoublet,primaryEnergy);
	    point        = new ParticlePoint(0.0, 0.0, materialNumber);
	    point.setIceRockBoundaryRadius(1.01*point.REarth);
	    event        = new Event(mcBases, propParticle, point);

	    double startlocation = point.getParticleLocation();
	    setStartLocation(startlocation);

	    while(true){

		double pathLength = event.getPhysicalPathLength(rand);
		event.setStepDx(0.1*pathLength); 
		event.traceParticle(pathLength);
		if((point.getParticleLocation()-getStartLocation())>trajectoryLength){
		    break;   // particles passed through your DETECTOR
		}else{
		    double transferedEnergy = event.collideNow(rand);

		    if(event.getFlavorByInteractionsInPlay()==0) {  // 0 is electron flavor
			int jLogE = (int)((Math.log(transferedEnergy)/ln10 - 
					   InteractionsBase.getLogEnergyProducedMinimum())
					  /Particle.getDeltaLogEnergy() - epsilon);
			if(jLogE<0)jLogE=0;
			setOneEmgCascadeMtx(primaryiLogE, jLogE);
			setOneTotalCascadeMtx(primaryiLogE, jLogE);
			
		    }else if(event.getFlavorByInteractionsInPlay()==3) { // 3 is hadron flavor
			int jLogE = (int)((Math.log(transferedEnergy)/ln10 - 
					   InteractionsBase.getLogEnergyProducedMinimum())
					  /Particle.getDeltaLogEnergy() - epsilon);
			if(jLogE<0)jLogE=0;
			setOneHadronCascadeMtx(primaryiLogE, jLogE);
			setOneTotalCascadeMtx(primaryiLogE, jLogE);
		    }
		}

		double logEnergy = event.propParticle.getLogEnergy(); // logEnergy after interaction
		if(logEnergy <= InteractionsBase.getLogEnergyProducedMinimum()) break;
                       // Journey ends because all the primary energy has been lost
		
		if(event.mcBaseInPlay.getTypeOfInteraction() == 1 && //Tau to Hadron/Electron decay
		   event.mcBaseInPlay.getProducedFlavor()!= 1) break; //or Tau to Mu decay
	    }
    }
    

    /** Method to run multiple events with a monocromatic primary energy given by
     the constructor. This method can make data of each interaction/decay respectively. */
    public void runMultipleEvent_Each(DataOutputStream out) throws IOException {
	
	double logEnergyMinimum = Particle.getLogEnergyMinimum();
	int    expandedDim      = dim + 
	    (int )((logEnergyMinimum-InteractionsBase.getLogEnergyProducedMinimum())/
		   Particle.getDeltaLogEnergy());

	int length = mcBases.length + tauDecayFlag;  // length of the mcBases.
                                                     // Tau decay can both of emg and hadron cascades.
	oneTotalCascadeMtx = new double[length][expandedDim];
	totalCascadeMtx    = new double[length][expandedDim];
	energyOfCascade    = new double[length];
	numberOfCascade    = new int[length];

	for(int n=0; n<numberOfEvent; n++){
	    primaryiLogE = 0;  // These --Mtx are not a matrix but an array.
	    System.out.println(n);
	    
	    // Initialization
	    for(int l=0; l<length; l++){
		energyOfCascade[l] = 0.0;
	    }
	    
	    propParticle = new Particle(primaryFlavor,primaryDoublet,primaryEnergy);
	    point        = new ParticlePoint(0.0, 0.0, materialNumber);
	    point.setIceRockBoundaryRadius(1.01*point.REarth);
	    event        = new Event(mcBases, propParticle, point);

	    double startlocation = point.getParticleLocation();
	    setStartLocation(startlocation);

	    while(true){

		double pathLength = event.getPhysicalPathLength(rand);
		event.setStepDx(0.1*pathLength);
		event.traceParticle(pathLength);
		if((point.getParticleLocation()-getStartLocation())>trajectoryLength){
		    break;   // particles passed through your DETECTOR
		}else{
		    double transferedEnergy = event.collideNow(rand);
		    int m = 10000;
		    for(int mm=0; mm<mcBases.length; mm++) {
			if(event.mcBaseInPlay.getInteractionName().equals(mcBases[mm].getInteractionName())) m=mm;  
		    }
		    
		    if(tauDecayFlag==1 && m==mTauDecay) {   // when the mcBase = tau decay
			if(event.getFlavorByInteractionsInPlay()==0) m=length-2; //tau2ElectronDecay
			else if(event.getFlavorByInteractionsInPlay()==3) m=length-1; //tau2HadronDecay
		    }

		    energyOfCascade[m] += transferedEnergy;
		    int jLogE = (int)((Math.log(transferedEnergy)/ln10 - 
				       InteractionsBase.getLogEnergyProducedMinimum())
				      /Particle.getDeltaLogEnergy() - epsilon);
		    if(jLogE>=expandedDim) jLogE = expandedDim-1;
		    if(jLogE<0) jLogE=0;
		    setOneCascadeMtx(m,jLogE);
		    
		    double logEnergy = event.propParticle.getLogEnergy(); // logEnergy after interaction
		    if(logEnergy <= InteractionsBase.getLogEnergyProducedMinimum()) break;
			// Journey ends because all the primary energy has been lost
			
		    if(event.mcBaseInPlay.getTypeOfInteraction() == 1 &&
		       event.mcBaseInPlay.getProducedFlavor()!= 1) break;		    
		}
	    }

	    // make Array for Energy deposit/1km
	    if(integral) {
		for(int nn=0; nn<length; nn++) {
		    
		    int jLogE = (int)((Math.log(energyOfCascade[nn])/ln10
				       - InteractionsBase.getLogEnergyProducedMinimum())
				      /Particle.getDeltaLogEnergy() - epsilon);
		    if(jLogE>=expandedDim) jLogE = expandedDim-1;
		    if(jLogE<0) jLogE=0;
		    if(jLogE>=0) {
			totalCascadeMtx[nn][jLogE] += 1.0;
		    }
		}
	    }
	    
	} 

	// write out Matrix 
	if(primaryDoublet==0) numberOfEvent = numberOfEvent * neutrinoFactor;
	if(integral) { // EnergyDeposit/1km
	    for(int k=0; k<mcBases.length; k++){ 
		System.out.println(mcBases[k].getInteractionName());
	    }
	    for(int nn=0; nn<length; nn++) {
		for(int jLogE=0; jLogE<expandedDim; jLogE++) {
		    totalCascadeMtx[nn][jLogE] = totalCascadeMtx[nn][jLogE]/numberOfEvent;
		    out.writeDouble(totalCascadeMtx[nn][jLogE]);
		}
	    }
	}
	if(primaryDoublet==0) numberOfEvent = numberOfEvent / neutrinoFactor;

	if(differential) { // EnergyDeposit/oneCascade
	    for(int nn=0; nn<length; nn++) {
		if(primaryDoublet==0) numberOfCascade[nn] = numberOfCascade[nn] * neutrinoFactor;
		for(int jLogE=0; jLogE<expandedDim; jLogE++) {
		    oneTotalCascadeMtx[nn][jLogE] = oneTotalCascadeMtx[nn][jLogE]/numberOfCascade[nn];
		    out.writeDouble(oneTotalCascadeMtx[nn][jLogE]);
		}
		if(primaryDoublet==0) numberOfCascade[nn] = numberOfCascade[nn] / neutrinoFactor;
	    }
	}

    }
    
    /** make histgram for EnergyDeposit/cascade */
    public void setOneCascadeMtx(int n, int jLogE) {
	oneTotalCascadeMtx[n][jLogE] += 1.0;
	numberOfCascade[n]++;
    }

    /** Method to run multiple events with a monocromatic primary energy given by
     the constructor. This method can make data of electromagnetic and hadronic cascade
     energy respectively. */
    public void runMultipleEvent(DataOutputStream out) throws IOException {

	double logEnergyMinimum = Particle.getLogEnergyMinimum();
	int    expandedDim      = dim + 
                       	          (int )((logEnergyMinimum-InteractionsBase.getLogEnergyProducedMinimum())/
				   Particle.getDeltaLogEnergy());

	emgCascadeMtx    = new double[1][expandedDim];
	hadronCascadeMtx = new double[1][expandedDim];
	totalCascadeMtx  = new double[1][expandedDim];
	
	oneEmgCascadeMtx    = new double[1][expandedDim];
	oneHadronCascadeMtx = new double[1][expandedDim];
	oneTotalCascadeMtx  = new double[1][expandedDim];
	    
	numberOfEmgCascade    = 0;
	numberOfHadronCascade = 0;
	numberOfTotalCascade  = 0;

	for(int n=0; n<numberOfEvent; n++){
	    System.out.println(n);
	    
	    primaryiLogE = 0;  // These --Mtx have only one row.
	    runSingleEvent();

	    if(integral) {  	    // make matrix for energy deposit/1km

		if(event.getCascadeEmgEnergy() != 0){
		    int jLogE = 
			(int)((Math.log(event.getCascadeEmgEnergy())/ln10 - 
			       InteractionsBase.getLogEnergyProducedMinimum())
			      /Particle.getDeltaLogEnergy() - epsilon);
		    if(jLogE>=expandedDim) jLogE = expandedDim-1;
		    if(jLogE<0) jLogE=0;
		    if(jLogE>=0) {
			emgCascadeMtx[0][jLogE] += 1.0;
		    }
		}

		if(event.getCascadeHadronEnergy()!=0){
		    int jLogE = 
			(int)((Math.log(event.getCascadeHadronEnergy())/ln10 -
                               InteractionsBase.getLogEnergyProducedMinimum())
			      /Particle.getDeltaLogEnergy() - epsilon);
		    if(jLogE>=expandedDim) jLogE = expandedDim-1;
		    if(jLogE<0) jLogE=0;
		    if(jLogE>=0) {
			hadronCascadeMtx[0][jLogE] += 1.0;
		    }
		}
		if(event.getCascadeTotalEnergy()>0.0){
		    int jLogE = 
			(int)((Math.log(event.getCascadeTotalEnergy())/ln10 - 
			       InteractionsBase.getLogEnergyProducedMinimum())
			      /Particle.getDeltaLogEnergy()  - epsilon);
		    if(jLogE>=expandedDim) jLogE = expandedDim-1;
		    if(jLogE<0) jLogE=0;
		    if(jLogE>=0) {
			totalCascadeMtx[0][jLogE] += 1.0;
		    }
		}
	    }
	}  			

	  
	if(primaryDoublet==0) numberOfEvent = numberOfEvent * neutrinoFactor;
	for(int jLogE=0; jLogE<expandedDim; jLogE++){  // normalization
	    emgCascadeMtx[0][jLogE]    = emgCascadeMtx[0][jLogE]/numberOfEvent; 
	    hadronCascadeMtx[0][jLogE] = hadronCascadeMtx[0][jLogE]/numberOfEvent; 
	    totalCascadeMtx[0][jLogE]  = totalCascadeMtx[0][jLogE]/numberOfEvent; 
	}
    	if(primaryDoublet==0) numberOfEvent = numberOfEvent / neutrinoFactor;
	// make Matrix for oneCascade
	if(differential) {
	    
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){	// Normalization
		oneEmgCascadeMtx[0][jLogE]    = oneEmgCascadeMtx[0][jLogE]/numberOfEmgCascade;
		oneHadronCascadeMtx[0][jLogE] = oneHadronCascadeMtx[0][jLogE]/numberOfHadronCascade;
		oneTotalCascadeMtx[0][jLogE]  = oneTotalCascadeMtx[0][jLogE]/numberOfTotalCascade;
	    }
	}
    
	// write out Matrix 
	if(integral) { // EnergyDeposit/1km
	    out.writeInt(primaryFlavor);
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(emgCascadeMtx[0][jLogE]);
	    }
	
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(hadronCascadeMtx[0][jLogE]);
	    }
	    
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(totalCascadeMtx[0][jLogE]);
	    }
	    
	}
	
	if(differential) { // EnergyDep/cascade
	    
	    out.writeInt(primaryFlavor);
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(oneEmgCascadeMtx[0][jLogE]);
	    }
	    
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(oneHadronCascadeMtx[0][jLogE]);
	    }
	    
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){
		out.writeDouble(oneTotalCascadeMtx[0][jLogE]);
	    }
	    
	}
    }

    /** Method to run multiple events with various primary energies from
	logE = Particle.getLogEnergyMinimum() all the way up to 10^12 GeV.
        The results are stored in the matrix form and will be handled
	by EventMatrix.class in the event package.
    */
    public void runEventOnMatrix(DataOutputStream out) throws IOException {
	double logEnergyMinimum = Particle.getLogEnergyMinimum();
	int expandedDim = dim + 
                       (int )((logEnergyMinimum-InteractionsBase.getLogEnergyProducedMinimum())/
			      Particle.getDeltaLogEnergy());
	emgCascadeMtx    = new double[dim][expandedDim];
	hadronCascadeMtx = new double[dim][expandedDim];
	totalCascadeMtx  = new double[dim][expandedDim];

	oneEmgCascadeMtx    = new double[dim][expandedDim];
	oneHadronCascadeMtx = new double[dim][expandedDim];
	oneTotalCascadeMtx  = new double[dim][expandedDim];

	numberOfEmgCascade    = 0;
	numberOfHadronCascade = 0;
	numberOfTotalCascade  = 0;
	
	for(int iLogE=0; iLogE<dim; iLogE++){
	    primaryiLogE = iLogE;
	    System.out.println(iLogE);
	    double logPrimaryEnergy = 
		Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
	    primaryEnergy = Math.pow(10.0,logPrimaryEnergy);

	    for(int n=0; n<numberOfEvent; n++){

		runSingleEvent();

		// make matrix for energy deposit/1km
		if(integral) {
		    
		    if(event.getCascadeEmgEnergy()!=0){
			int jLogE = 
			    (int)((Math.log(event.getCascadeEmgEnergy())/ln10 - 
				   InteractionsBase.getLogEnergyProducedMinimum())
				  /Particle.getDeltaLogEnergy() - epsilon);
			if(jLogE>=expandedDim) jLogE = expandedDim-1;
			if(jLogE<0) jLogE=0;
			if(jLogE>=0) {
			    emgCascadeMtx[iLogE][jLogE] += 1.0;
			}
		    }
		    
		    if(event.getCascadeHadronEnergy()!=0){
			int jLogE = 
			    (int)((Math.log(event.getCascadeHadronEnergy())
				   /ln10 - InteractionsBase.getLogEnergyProducedMinimum())
				  /Particle.getDeltaLogEnergy() - epsilon);
			if(jLogE>=expandedDim) jLogE = expandedDim-1;
			if(jLogE<0) jLogE=0;
			if(jLogE>=0) {
			    hadronCascadeMtx[iLogE][jLogE] += 1.0;
			}
		    }
		    
		    if(event.getCascadeTotalEnergy()>0.0){
			int jLogE = 
			    (int)((Math.log(event.getCascadeTotalEnergy())/ln10 - 
				   InteractionsBase.getLogEnergyProducedMinimum())
				  /Particle.getDeltaLogEnergy() - epsilon);
			if(jLogE>=expandedDim) jLogE = expandedDim-1;
			if(jLogE<0) jLogE=0;
			if(jLogE>=0) {
			    totalCascadeMtx[iLogE][jLogE] += 1.0;
			}
		    }
		}
	    }
	    if(primaryDoublet==0) numberOfEvent = numberOfEvent * neutrinoFactor;
	    for(int jLogE=0; jLogE<expandedDim; jLogE++){  // normalization
		emgCascadeMtx[iLogE][jLogE]    = emgCascadeMtx[iLogE][jLogE]/numberOfEvent; 
		hadronCascadeMtx[iLogE][jLogE] = hadronCascadeMtx[iLogE][jLogE]/numberOfEvent; 
		totalCascadeMtx[iLogE][jLogE]  = totalCascadeMtx[iLogE][jLogE]/numberOfEvent; 
	    }
	    if(primaryDoublet==0) numberOfEvent = numberOfEvent / neutrinoFactor;
	}
	
	// make Matrix for oneCascade
	if(differential) {
	    
	    for(int iLogE=0; iLogE<dim; iLogE++){  // Normalization
		for(int jLogE=0; jLogE<expandedDim; jLogE++){	
		    oneEmgCascadeMtx[iLogE][jLogE]    = oneEmgCascadeMtx[iLogE][jLogE]/numberOfEmgCascade;
		    oneHadronCascadeMtx[iLogE][jLogE] = oneHadronCascadeMtx[iLogE][jLogE]/numberOfHadronCascade;
		    oneTotalCascadeMtx[iLogE][jLogE]  = oneTotalCascadeMtx[iLogE][jLogE]/numberOfTotalCascade;
		}
	    }
	}
	// write out Matrix 
	if(integral) { // EnergyDeposit/1km
	    out.writeInt(primaryFlavor);
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(emgCascadeMtx[iLogE][jLogE]);
		}
	    }
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(hadronCascadeMtx[iLogE][jLogE]);
		}
	    }
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(totalCascadeMtx[iLogE][jLogE]);
		}
	    }
	}
	
	if(differential) { // EnergyDep/cascade
	    
	    out.writeInt(primaryFlavor);
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(oneEmgCascadeMtx[iLogE][jLogE]);
		}
	    }
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(oneHadronCascadeMtx[iLogE][jLogE]);
		}
	    }
	    for(int iLogE=0; iLogE<dim; iLogE++){
		for(int jLogE=0; jLogE<expandedDim; jLogE++){
		    out.writeDouble(oneTotalCascadeMtx[iLogE][jLogE]);
		}
	    }
	}
    }
	
    public void setOneEmgCascadeMtx(int iLogE, int jLogE) {
	oneEmgCascadeMtx[iLogE][jLogE] += 1.0;
	numberOfEmgCascade++;
    }

    public void setOneHadronCascadeMtx(int iLogE, int jLogE) {
	oneHadronCascadeMtx[iLogE][jLogE] += 1.0;
	numberOfHadronCascade++;
    }

    public void setOneTotalCascadeMtx(int iLogE, int jLogE) {
	oneTotalCascadeMtx[iLogE][jLogE] += 1.0;
	numberOfTotalCascade++;
    }


    /** Main method. In order to run JULIET, it is need that RunManager object 
        and DataOutPutStream are generated. */
    public static void main(String[] args) throws IOException {

	// generate RunManager object
	RunManager run = new RunManager(); 

	// make output stream
	DataOutputStream out = new DataOutputStream(new FileOutputStream(fileName));

	if(typeOfEvent==0) {
	    run.runEventDump(out);
	}

	else if(typeOfEvent==1) {
	    if(eachInteraction) run.runMultipleEvent_Each(out);
	    else  run.runMultipleEvent(out);
	}

	else if(typeOfEvent==2) {
	    run.runEventOnMatrix(out);
	}
    	out.close();
    }

}



