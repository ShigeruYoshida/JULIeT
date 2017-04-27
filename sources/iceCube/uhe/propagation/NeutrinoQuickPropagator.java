package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.propagation.*;
import java.io.*;

/**

This class provides the methods to cauculate energy distribution
of particles resulted from the neutrino propagation
in the handy approximated way. You can let PropagationMatrix
do the same job, but it would take more CPU time.
Instead this class uses the analytical expressons for
neutrino propagation and the precalculated PropagationMatrix
for the secondary generated muon and tau propagation,
assuming that the regeneration of neutrinos and their
successive interactions can be neglected as a first order
approximation.

<pre>

     |-----------neutrino -----|----------- charged lepton --------------|
     |-------------------------XX----------------------------------------|
     <----    analytical ------>
                                <---  PropagationMatrix mutiplication --->

</pre>

As illustrated above, the charged lepton part is taken care of by
multiplying the matrix with step size Delta X. 
The PropagationMatrixFactory object, one of the protected member
in this class,  is responsible for interface with the matrix data
with step size of Delta X.

The vertex point XX
is integrated from the begining though the end inside the method
propagateNeutrino(). 

*/

public class NeutrinoQuickPropagator {

    private static final double ln10 = Math.log(10.0);

    protected static String[] pathname =
    {"../data/neutrino_earth/ice/","../data/neutrino_earth/rock/"};

    protected static String[] intMtxPathname =
    {"iceCube/uhe/interactions/ice/","iceCube/uhe/interactions/rock/"};

    protected static final String[] matrixFile = 
    {"1e6_00cm.data","1e5_00cm.data", "1e4_00cm.data"};

    protected static double[] stepSizeBase = {1.0e6,1.0e5,1.0e4}; //[cm]

    protected static double distanceMaximum = 5.0e7; // [cm]

    protected static double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level


    /** Propagation Matrix to handle the charged lepton propgation 
	after the neutrino interaction */
    protected PropagationMatrixFactory matrix = null;
    private boolean matrixHasBeenRead = false;
    /** Propagated medium class */
    protected ParticlePoint s = null;
    /** The Propagation Matrix elements from mu/tau to muon - 
	the neutrino interaction vertex to the end */
    double[][] FmuToMu,FtauToMu;
    /** The Propagation Matrix elements from mu/tau to tau - 
	the neutrino interaction vertex to the end */
    double[][] FmuToTau,FtauToTau;
    /** The Propagation Matrix elements from neutrino to the charged lepton
	the start point to the end */
    double[][] FnuEToMu,FnuEToTau,FnuMuToMu,FnuTauToMu,FnuMuToTau,FnuTauToTau;
    /** The Propgation Matrix elements from neutrino to neutrino
	the start point to the end */
    double[] FnuToNu,FnuEToNuE;
    private boolean hasPropagated = false;


    /** The default neutrino matrix charged current interaction 
	to generate the charged lepton */
    protected InteractionsMatrix nuCCMtx = null; 
    private String nuCCMtxObjectFile = "ENeutrinoChargeMtx";
    /** The Leptonic Glashow Resonance interaction to generate the charged lepton */
    protected InteractionsMatrix grLeptonMtx = null;
    private String grLeptonMtxObjectFile = "eGlashowResonanceLeptonicMtx";
    /** The Hadronic Glashow Resonance interaction to generate the charged lepton */
    private InteractionsMatrix grHadronMtx = null;
    private String grHadronMtxObjectFile = "glashowResonanceHadronicMtx";

    public double nuCCEnhancementFactor = 1.0;

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();

    /** 
	Constructor.
	<pre>
	(1). Register the particle point class
	(2). Read the neutrino interaction matrix 
	(3). Generate the PropagationMatrixFactory object to calculate
             the muon/tau propagation
	</pre>
    */

    public NeutrinoQuickPropagator(ParticlePoint s) throws IOException{
	this.s = s;
	generateNeutrinoInteractionMatrix();
	generatePropagationMatrixArray();
	matrix = new PropagationMatrixFactory();
	initLeptonMatrix();
	initNeutrinoMatrix();
    }

    /** Allocate memory for the propagation matrix array */
    private void generatePropagationMatrixArray(){
        FmuToMu= new double[dimension][dimension];
        FmuToTau= new double[dimension][dimension];
        FtauToMu= new double[dimension][dimension];
        FtauToTau= new double[dimension][dimension];
        FnuEToMu= new double[dimension][dimension];
        FnuEToTau= new double[dimension][dimension];
        FnuMuToMu= new double[dimension][dimension];
        FnuMuToTau= new double[dimension][dimension];
        FnuTauToMu= new double[dimension][dimension];
        FnuTauToTau= new double[dimension][dimension];

	FnuToNu = new double[dimension];
	FnuEToNuE = new double[dimension];
    }

    /** generate the neutrino interaction matrices */
    private void generateNeutrinoInteractionMatrix() throws IOException {
	// The Leptonic Glashow
	String objectFile = 
	    intMtxPathname[s.getMaterialNumber()].concat(grLeptonMtxObjectFile);
	InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	grLeptonMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	in.close( );
	// The Hadronic Glashow
	objectFile = intMtxPathname[s.getMaterialNumber()].concat(grHadronMtxObjectFile);
	in = ClassLoader.getSystemResourceAsStream(objectFile);
	grHadronMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	// The Charged Current neutrino-nucleon
	objectFile = intMtxPathname[s.getMaterialNumber()].concat(nuCCMtxObjectFile);
	in = ClassLoader.getSystemResourceAsStream(objectFile);
	nuCCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
    }

    /** Read the propagation matrix data for handing the charged lepton 
	propagation. The PropagationMatrixFactory object is responsible
	for this task.
    */
    public void readMatrix(DataInputStream in) throws IOException {
	matrix.readMatrix(in);
	matrixHasBeenRead = true;
    }

    /** initialized the propagation matrix elements */
    protected void initLeptonMatrix(){
	// Initialization
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<dimension;jLogE++){
		FmuToMu[iLogE][jLogE] = 0.0;
		FmuToTau[iLogE][jLogE] = 0.0;
		FtauToMu[iLogE][jLogE] = 0.0;
		FtauToTau[iLogE][jLogE] = 0.0;
	    }
	}

	for(iLogE=0;iLogE<dimension;iLogE++){
	    FmuToMu[iLogE][iLogE] = 1.0;
	    FtauToTau[iLogE][iLogE] = 1.0;
	}
    }

    /** initialized the propagation matrix elements */
    protected void initNeutrinoMatrix(){
	// Initialization
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    FnuToNu[iLogE] = 0.0;
	    FnuEToNuE[iLogE] = 0.0;
	    for(jLogE=0;jLogE<dimension;jLogE++){
		FnuEToMu[iLogE][jLogE] = 0.0; 
		FnuEToTau[iLogE][jLogE] = 0.0;
		FnuMuToMu[iLogE][jLogE] = 0.0;
		FnuMuToTau[iLogE][jLogE] = 0.0;
		FnuTauToMu[iLogE][jLogE] = 0.0;
		FnuTauToTau[iLogE][jLogE] = 0.0;
	    }
	}
	hasPropagated = false;
    }

    /** If you need to read the matrix data without the Glashow Resonance 
	call this method first with flag=false before readMatrix(in).
     */
    public void whetherPropagationMatrixWithGlashowResonance(boolean flag){
	matrix.whetherPropagationMatrixWithGlashowResonance(flag);
    }

    /** Set the particle point class. You may need to call this method, for instance,
        when you switch the material from ice to rock, or run the method
        of propagateNeutrinoToIceCubeDepth(nadir, nuCCenhancement) with various nadir angles.
    */
    public void setParticlePoint(ParticlePoint s){
	this.s = s;
    }


    /** 
	Calculate the neutrino propagation generating leptons 
	in the appriximated handy method.
        The propagation is divided into the two part and calculated separately.
        <pre>
     double nuCCEnhancementFactor : multiplication factor to the standard CC cross section
     double slantDepth            :depth [g/cm^2] of the neutrino interaction vertex.
     double deltaX                : step size [g/cm^2] of the particle propagation

               |-----------neutrino -----|----------- charged lepton ---------|
               |-------------------------XX-----------------------------------|
<-----    slantDepth ------------------->
               <-------------------  totalPropagationLength ------------------>
<-offSetLength->
	</pre>
    */
    public void propagateNeutrinoToLepton(double totalPropagationLength, 
				  double neutrinoOffSetLength,
				  double deltaX,
				  double nuCCEnhancementFactor){

	if(!matrixHasBeenRead){ // Propagation Matrix Factory has not read
	                        // any matrix data yet!
	    System.err.println("No matrix data has been read!");
	    System.err.println("You must call readPropagationMatrix(in).");
	    System.exit(0);
	}

	// initialization 
	this.nuCCEnhancementFactor = nuCCEnhancementFactor;

	//
	// Propagate Neutrino
	//

	// neutrino to charged lepton
	double chargedLeptonPropagationLength = deltaX;
	double slantDepth = totalPropagationLength + neutrinoOffSetLength - deltaX;
	int numberOfSteps = (int )((totalPropagationLength)/deltaX + 0.01);
	System.err.println("Following the neutrino propagation by " + numberOfSteps +
			   " times steps");
	for(int times =1;times<=numberOfSteps;times++){
	    System.err.println(" " + times + " th step");
	    propagateChargedLepton(); // Propagate the charged lepton
	    System.err.println(" Charged Lepton Propagated");
	    calculateNeutrinoToLeptonTransfer(slantDepth,deltaX);
	    System.err.println(" Propagated Neutrino at depth " + slantDepth +
			       " [g/cm^2]");
	    slantDepth -= deltaX;
	}
	System.err.println("done");

	hasPropagated = true;
    }

    /** 
	Calculate the charged leptron generation matrix 
	which plays a role of propagation matrix from neutrino to charged lepton.
        <pre>
     double slantDepth            :depth [g/cm^2] of the neutrino interaction vertex.
     double deltaX                : step size [g/cm^2] of the particle propagation

     |-----------neutrino -----|----------- charged lepton --------------|
     |-------------------------XX----------------------------------------|
     <----    slantDepth ------>
	</pre>
    */
    protected void calculateNeutrinoToLeptonTransfer(double slantDepth, 
						     double deltaX){

	// Mutiply the matrix by the neutrino interaction probablity
	for(int iLogE=0;iLogE<dimension;iLogE++){
	    double sigmaCC = nuCCMtx.getSigmaMatrix(iLogE)*nuCCEnhancementFactor;
	    double sigmaGRLepton = grLeptonMtx.getSigmaMatrix(iLogE);
	    double sigmaGRHadron = grHadronMtx.getSigmaMatrix(iLogE);
	    double dumpingNuFactor = Math.exp(-ParticlePoint.NA*sigmaCC*slantDepth);
	    double dumpingNuEFactor = 
		Math.exp(-ParticlePoint.NA*
			 (sigmaCC+0.5*(3.0*sigmaGRLepton+sigmaGRHadron))*slantDepth);

	    // nu-e to muon
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuEToMu = dumpingNuEFactor*ParticlePoint.NA*
			 0.5*grLeptonMtx.getLeptonTransferMatrix(iLogE,kLogE);
		    double nuEToTau = nuEToMu;
		    element += 
			nuEToMu*FmuToMu[kLogE][jLogE]+nuEToTau*FtauToMu[kLogE][jLogE];
		}
		FnuEToMu[iLogE][jLogE] += element*deltaX;
	    }


	    // nu-mu to muon
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuMuToMu = dumpingNuFactor*ParticlePoint.NA*
			nuCCMtx.getLeptonTransferMatrix(iLogE,kLogE)*
			nuCCEnhancementFactor;
		    element += 
			nuMuToMu*FmuToMu[kLogE][jLogE];
		}
		FnuMuToMu[iLogE][jLogE] += element*deltaX;
	    }

	    // nu-tau to muon
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuTauToTau = dumpingNuFactor*ParticlePoint.NA*
			nuCCMtx.getLeptonTransferMatrix(iLogE,kLogE)*
			nuCCEnhancementFactor;
		    element += 
			nuTauToTau*FtauToMu[kLogE][jLogE];
		}
		FnuTauToMu[iLogE][jLogE] += element*deltaX;
	    }

	    // nu-e to tau
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuEToTau = dumpingNuEFactor*ParticlePoint.NA*
			 0.5*grLeptonMtx.getLeptonTransferMatrix(iLogE,kLogE);
		    double nuEToMu = nuEToTau;
		    element += 
			nuEToTau*FtauToTau[kLogE][jLogE]+nuEToMu*FmuToTau[kLogE][jLogE];
		}
		FnuEToTau[iLogE][jLogE] += element*deltaX;
	    }


	    // nu-mu to tau
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuMuToMu = dumpingNuFactor*ParticlePoint.NA*
			nuCCMtx.getLeptonTransferMatrix(iLogE,kLogE)*
			nuCCEnhancementFactor;
		    element += 
			nuMuToMu*FmuToTau[kLogE][jLogE];
		}
		FnuMuToTau[iLogE][jLogE] += element*deltaX;
	    }

	    // nu-tau to tau
	    for(int jLogE=0;jLogE<=iLogE;jLogE++){
		double element = 0.0;
		for(int kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    double nuTauToTau = dumpingNuFactor*ParticlePoint.NA*
			nuCCMtx.getLeptonTransferMatrix(iLogE,kLogE)*
			nuCCEnhancementFactor;
		    element += 
			nuTauToTau*FtauToTau[kLogE][jLogE];
		}
		FnuTauToTau[iLogE][jLogE] += element*deltaX;
	    }
	}
    }

    /** 
	propagate the charged leptons by the propagation matrix step size deltaX.
	by Mutiplication of the  lepton propagation matrix 
        mu->mu, tau->mu, mu->tau and tau->tau 
    */
    protected void propagateChargedLepton(){

	int iLogE,jLogE,kLogE;

	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){

		// Mu to Mu
		double elementMuToMu = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    elementMuToMu += 
			// neglect
			//FmuToNuE[iLogE][kLogE]*matrix.FnuEToMu[kLogE][jLogE]+
			//FmuToNuMu[iLogE][kLogE]*matrix.FnuMuToMu[kLogE][jLogE]+
			//FmuToNuTau[iLogE][kLogE]*matrix.FnuTauToMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*matrix.FmuToMu[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*matrix.FtauToMu[kLogE][jLogE];
		}

		// Tau to Mu
		double elementTauToMu = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    elementTauToMu += 
			// neglect
			//FtauToNuE[iLogE][kLogE]*matrix.FnuEToMu[kLogE][jLogE]+
			//FtauToNuMu[iLogE][kLogE]*matrix.FnuMuToMu[kLogE][jLogE]+
			//FtauToNuTau[iLogE][kLogE]*matrix.FnuTauToMu[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*matrix.FmuToMu[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*matrix.FtauToMu[kLogE][jLogE];
		}

		// Mu to Tau
		double elementMuToTau = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    elementMuToTau += 
			// neglect
			//FmuToNuE[iLogE][kLogE]*matrix.FnuEToTau[kLogE][jLogE]+
			//FmuToNuMu[iLogE][kLogE]*matrix.FnuMuToTau[kLogE][jLogE]+
			//FmuToNuTau[iLogE][kLogE]*matrix.FnuTauToTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*matrix.FmuToTau[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*matrix.FtauToTau[kLogE][jLogE];
		}

		// Tau to Tau
		double elementTauToTau = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    elementTauToTau += 
			// neglect
			//FtauToNuE[iLogE][kLogE]*matrix.FnuEToTau[kLogE][jLogE]+
			//FtauToNuMu[iLogE][kLogE]*matrix.FnuMuToTau[kLogE][jLogE]+
			//FtauToNuTau[iLogE][kLogE]*matrix.FnuTauToTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*matrix.FmuToTau[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*matrix.FtauToTau[kLogE][jLogE];
		}

		FmuToMu[iLogE][jLogE] = elementMuToMu;
		FtauToMu[iLogE][jLogE] = elementTauToMu;
		FmuToTau[iLogE][jLogE] = elementMuToTau;
		FtauToTau[iLogE][jLogE] = elementTauToTau;

	    }
	}
    }

    /** 
	Calculate the neutrino propagation in the appriximated handy method.
        The propagation is divided into the two part and calculated separately.
        <pre>
     double nuCCEnhancementFactor : multiplication factor to the standard CC cross section
     double slantDepth            :depth [g/cm^2] of the neutrino interaction vertex.

          |-----------neutrino -------------------------|
          <-----------    slantDepth ------------------->
	</pre>
    */
    public void propagateNeutrinoToNeutrino(double slantDepth,
				  double nuCCEnhancementFactor){

	this.nuCCEnhancementFactor = nuCCEnhancementFactor;

	// neutrino to neutrino
	for(int iLogE=0;iLogE<dimension;iLogE++){ 
	    double sigmaCC = nuCCMtx.getSigmaMatrix(iLogE)*nuCCEnhancementFactor;
	    double sigmaGRLepton = grLeptonMtx.getSigmaMatrix(iLogE);
	    double sigmaGRHadron = grHadronMtx.getSigmaMatrix(iLogE);
	    double dumpingNuFactor = Math.exp(-ParticlePoint.NA*sigmaCC*slantDepth);
	    double dumpingNuEFactor = 
		Math.exp(-ParticlePoint.NA*
			 (sigmaCC+0.5*(3.0*sigmaGRLepton+sigmaGRHadron))*slantDepth);
	    FnuToNu[iLogE] = dumpingNuFactor;
	    FnuEToNuE[iLogE] = dumpingNuEFactor;
	}
    }


    /**
       Let the propagationMatrixFactory read the propagation matrix data
       corresponds to the "fileName" 
    */
    protected void readPropagationMatrix(String filename)
    throws IOException {
	String matrixFileName = pathname[s.getMaterialNumber()].concat(filename);
	System.err.println(" --matrix filename " + matrixFileName);
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(matrixFileName));
	matrix.readMatrix(in);
	matrixHasBeenRead = true;
    }

    /** Calculation of the neutrino propagation with the methods
        of propagateNeutrinoToNeutrino() and propagateNeutrinoToLepton() 
      <pre>
     double propagationDepth   : the neutrino-lepton total propagation distance [g/cm^2]
     </pre>
    */
    public void propagateNeutrino(double propagationDepth,
				  double nuCCEnhancementFactor)	throws IOException {


	// initialization
	generateNeutrinoInteractionMatrix();
	initLeptonMatrix();
	initNeutrinoMatrix();

       	int stepIndex = 0;
	for(int i=0;i<stepSizeBase.length;i++){
	    stepIndex = i;
	    if(propagationDepth > stepSizeBase[i]*s.getMediumDensity()){
		break;
	    }
	    
	}

	// Transportation from Neutrino to Neutrino
	propagateNeutrinoToNeutrino(propagationDepth,nuCCEnhancementFactor);

	//
	// Transportation from Neutrino to muon and taus
	double distanceToPropagate = propagationDepth;
	double maximumDepth = distanceMaximum*s.getMediumDensity();

	if(propagationDepth>maximumDepth){     // Bayond the maximum limit
	    distanceToPropagate = maximumDepth;// of the mu/tau propagation length
	}

	// calculate the charged lepton propagation length in each of stepSize
	double[]  distanceInThisStep = new double[stepSizeBase.length];
	for(int i=stepIndex;i<stepSizeBase.length;i++){
	    double stepSize = stepSizeBase[i]*s.getMediumDensity();
	    distanceInThisStep[i] = stepSize*(double )
		((int )(distanceToPropagate/stepSize+0.01));
	    if(i != (stepSizeBase.length-1)) distanceInThisStep[i] -= stepSize;
	    distanceToPropagate -= distanceInThisStep[i];
	}

	distanceToPropagate = propagationDepth;
	double distancePropagated = 0.0;
	if(propagationDepth>maximumDepth){   // Bayond the maximum limit
	    distancePropagated = propagationDepth-maximumDepth; 
	}

	// Integral of the paths from neutrino to leptons
	for(int i=(stepSizeBase.length-1);i>=stepIndex;i--){
	    // read the matrix
	    readPropagationMatrix(matrixFile[i]);
	    System.err.println(" reading the matrix data done");

	    double stepSize = stepSizeBase[i]*s.getMediumDensity();
	    propagateNeutrinoToLepton(distanceInThisStep[i],
				      distanceToPropagate-distanceInThisStep[i],
				      stepSize, nuCCEnhancementFactor);

	    distanceToPropagate -= distanceInThisStep[i];
	    distancePropagated += distanceInThisStep[i];
	}


	System.err.println(" Neutrino propagation completed");
	System.err.println(" Distance propagated = " + distancePropagated + 
			   " [g/cm^2] " + distancePropagated/s.getMediumDensity() +
			   " [cm]");

    }


    /** Calculation of the neutrino propagation with the methods
        of propagateNeutrinoToNeutrino() and propagateNeutrinoToLepton() 
      <pre>
     double nadirAngle   : nadire (upgoing) zenith (down) of the neutrino-lepton track [deg]

     </pre>
     Note : nadir_angle must be registered in the ParticlePoint object member s
    */
    public void propagateNeutrinoToIceCubeDepth(double nadirAngle,
				  double nuCCEnhancementFactor)	throws IOException {

	// calculate the distance [cm]
	double trajectoryLength;
	if(s.getMaterialNumber()==0){ // ice i.e., downgoing
	    double zenithAngle =  nadirAngle;  // Zenith Angle  0deg for vertical
	    double cos_zenith = Math.cos(zenithAngle*Math.PI/180.0);
	    double sq_term = Math.sqrt((ParticlePoint.REarth-detectorDepth)*
				       (ParticlePoint.REarth-detectorDepth)
				       *cos_zenith*cos_zenith + 
				       2.0*ParticlePoint.REarth*detectorDepth-detectorDepth*detectorDepth);
	    double cos_nadir = sq_term/ParticlePoint.REarth;
	    trajectoryLength = sq_term - (ParticlePoint.REarth-detectorDepth)*cos_zenith;
	}else{ // rock i.e., upgoing
	    trajectoryLength = 2.0*s.REarth*Math.cos(nadirAngle*Math.PI/180.0);
	}

	// calculate the slant depth [g/cm^2];
	double slantDepth = 0.0;          // [g/cm^2]
	double distancePropagated = 0.0; // [cm]
	double dX = 10.0; // 10 g
	s.setParticleLocation(0.0);
	while(distancePropagated< trajectoryLength){
	    double deltaL = dX/s.getMediumDensity( ); 
	    distancePropagated += deltaL;
	    slantDepth += dX;
	    s.setParticleLocation(distancePropagated);
	}

	System.err.println(" angle(" + nadirAngle + ") [deg] distance=" +
			   trajectoryLength + " [cm] slant depth=" +
			   slantDepth);

	propagateNeutrino(slantDepth,nuCCEnhancementFactor);

    }

    /** 
	Returns dF/dLogE * deltaLogE (inputNeutrino ---> outputParticle)
	calculated by th method propagateNeutrino().

	<pre>
	int   neutrinoFlavor   : flavor of the NEUTRINO entering into the earth.
	double logEneutrino    : logE [GeV] of the NEUTRINO entering into the earth.
	int   outputFlavor   : flavor of the Particle object after the propagation.
	int   outputDoublet  : doublet of the Particle object after the propagation.
	double logEoutput    : logE [GeV] of the particle after the propagation.
	</pre>

	delta LogE is the bin size of the propagation matrix. It is defined by
	Particle.getDimensionOfLogEnergyMatrix()
    */
    public double getDF(int neutrinoFlavor, double logEneutrino,
			int outputFlavor, int outputDoublet, double logEoutput){

	if(!hasPropagated){ // the propagation calculation has not done yet!
	    System.err.println("You must call the method propagateNeutrino()");
	    System.exit(0);
	}

	int iLogE = (int)((logEneutrino + 0.1*Particle.getDeltaLogEnergy()
		   - Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy());
	int jLogE = (int)((logEoutput + 0.1*Particle.getDeltaLogEnergy()
		   - Particle.getLogEnergyMinimum())/Particle.getDeltaLogEnergy());

	double count = 0.0;
	if(((0 <= iLogE) && (iLogE<dimension)) && 
	   ((jLogE<= iLogE) && (0 <= jLogE))){ // In the valid energy range

	    if(outputDoublet == 1){ // charged lepton

		if(outputFlavor == 1) { // muon
		    if(neutrinoFlavor==0) count = FnuEToMu[iLogE][jLogE];// nu_e to mu
		    else if(neutrinoFlavor==1) count = FnuMuToMu[iLogE][jLogE];// nu_mu to mu
		    else if(neutrinoFlavor==2) count = FnuTauToMu[iLogE][jLogE];// nu_tau to mu
		}else if(outputFlavor == 2) { // tau
		    if(neutrinoFlavor==0) count = FnuEToTau[iLogE][jLogE];// nu_e to tau
		    else if(neutrinoFlavor==1) count = FnuMuToTau[iLogE][jLogE];// nu_mu to tau
		    else if(neutrinoFlavor==2) count = FnuTauToTau[iLogE][jLogE];// nu_tau to tau
		}
	    }else if(outputDoublet == 0){ // neutrino as output particle

		if(neutrinoFlavor == outputFlavor){ // with no flavor flip
		    if(neutrinoFlavor!=0){   // nu_mu/tau to nu_mu/tau 
			if(iLogE == jLogE) count = FnuToNu[iLogE];
		    }else{ // nu_e to nu_e
			if(iLogE == jLogE) count = FnuEToNuE[iLogE];
		    }
		}
	    }

	}

	return count;

    }

}


