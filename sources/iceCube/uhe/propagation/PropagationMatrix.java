package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.points.*;
import java.io.*;

/**
<pre>
   Matrix of the Energy and particle spiecies Transfer 
    by the Interactions and Particle Decays.
    The matrix elements are a priori calculated by the methods supplied
    by the InteractionMatrix in iceCube/uhe/interactions/ and
    by the Decay classes in iceCube/uhe/decay/ .


The infinitesimal propagation over dX is calculated as

nuToNu        : by the Neutral Current Intereactions.
nuEToLepton   : by the Charged Current Intereactions.
nuEToHadron   : by the Neutral/Charged currents.
nuEToNuE      : by the Electronic Glashow Resonance.

For Glashow Resonance
nuEToNuMu     : by the Muonic Glashow Resonance.
nuEToNuTau    : by the Tauonic Glashow Resonance.
nuEToHadron   : by the Hadrinic Glashow Resonance.

muToNuE       : by                                     Decay
muToNuMu      : by the Charged Current Intereactions   Decay
muToE         : by Pair Creation and Bremsstrahlung    Decay
muToMu        : by Pair Creation and survival.         Decay
muToTau       : by Pair Creation.
muToHadron    : by PhotoNuclear Interactions.
tauToNuE      : by                                     Decay
tauToNuMu     : by                                     Decay
tauToNuTau    : by the Charged Current Intereactions   Decay
tauToE        : by Pair Creation and Bremsstrahlung    Decay
tauToMu       : by Pair Creation.                      Decay
tauToTau      : by Pair Creation and survival
tauToHadron   : by PhotoNuclear Interactions.          Decay

Allocated memory for these array can also be used for
storing the finite propagation length.

The Main finite Propagation Matrix are formed such as

<channels to nuE>
FnuEToNuE(logEin,logEout),FnuMuToNuE(logEin,logEout),FnuTauToNuE(logEin,logEout),
FmuToNuE(logEin,logEout),FtauToNuE(logEin,logEout)

<channels to Taus>
FnuMuToTau(logEin,logEout),FnuTauToTau(logEin,logEout),
FmuToTau(logEin,logEout),FtauToTau(logEin,logEout)

etc.



Switching on/off an individual interactions channels
are done by the bits patemeter interactionsSwitch and decaySwitch.

interactionsSwitch
 bit8    bit7     bit6      bit5     bit4     bit3     bit2     bit1     bit0
GlaRes  LepWeak  PhotoNucl Bremss  KnockOn  PairCHeavy PairC  Neutral  Charged

decaySwitch
bit7-2   bit 1     bit0
Reserv.  TauDecay  MuDecay


These switches are read when generate this object.

</pre>
*/

public class PropagationMatrix {

    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    double massNumber = 0.0; // total target mass number in the medium
    /** For Glashow Resonance 
	dimension at resonance **/
    final int dimResonance = 180;
    /** Step size of the propagation [g/cm^2].
        Initialized in the constructor.*/
    public double dX = 1.0; 
    /** Step size of the propagation [g/cm^2] determined by 
	the shortest decay length.
        Initialized in the constructor.*/
    public double dXDecay = 1.0;
    /** Speed of light [cm/sec].*/
    public final static double c = 2.99792e10;
    protected double[][][] temp;
    protected double[]   intProbNeutrino,intProbMu,intProbTau;
    /** For Glashow Resonance **/
    protected double[]   intProbNuE;
    //protected double[][] nuToNu,nuToLepton,nuToHadron;
    protected double[][] nuEToNuE,nuMuToNuE,nuTauToNuE,muToNuE,tauToNuE;
    protected double[][] nuEToNuMu,nuMuToNuMu,nuTauToNuMu,muToNuMu,tauToNuMu;
    protected double[][] nuEToNuTau,nuMuToNuTau,nuTauToNuTau,muToNuTau,tauToNuTau;
    protected double[][] nuEToE,nuMuToE,nuTauToE,muToE,tauToE;
    protected double[][] nuEToMu,nuMuToMu,nuTauToMu,muToMu,tauToMu;
    protected double[][] nuEToTau,nuMuToTau,nuTauToTau,muToTau,tauToTau;
    protected double[][] nuEToHadron,nuMuToHadron,nuTauToHadron,muToHadron,tauToHadron;
    protected double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    protected double[][] FnuEToNuMu,FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    protected double[][] FnuEToNuTau,FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    protected double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    protected double[][] FnuEToMu,FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    protected double[][] FnuEToTau,FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    protected double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;
    protected double[][] SnuEToNuE,SnuMuToNuE,SnuTauToNuE,SmuToNuE,StauToNuE;
    protected double[][] SnuEToNuMu,SnuMuToNuMu,SnuTauToNuMu,SmuToNuMu,StauToNuMu;
    protected double[][] SnuEToNuTau,SnuMuToNuTau,SnuTauToNuTau,SmuToNuTau,StauToNuTau;
    protected double[][] SnuEToE,SnuMuToE,SnuTauToE,SmuToE,StauToE;
    protected double[][] SnuEToMu,SnuMuToMu,SnuTauToMu,SmuToMu,StauToMu;
    protected double[][] SnuEToTau,SnuMuToTau,SnuTauToTau,SmuToTau,StauToTau;
    protected double[][] SnuEToHadron,SnuMuToHadron,SnuTauToHadron,SmuToHadron,StauToHadron;

    protected final static int CHARGED_FLAG = 1; // Allows the Charged Current.
    protected final static int NEUTRAL_FLAG = 2; // Allows the neutral Current.
    protected final static int PAIRC_FLAG = 4;// Allows the e+e- pair Creation.
    protected final static int PAIRCH_FLAG = 8;// Allows the mu+mu- pair Creation.
    protected final static int KNOCK_FLAG = 16; // Allows the Knock on electrons.
    protected final static int BREMSS_FLAG = 32;// Allows the Bremsstrahlung.
    protected final static int PHOTO_FLAG = 64; // Allows the Photo-nuclear.
    protected final static int LEPTW_FLAG = 128; // Allows the charged leptons
                                       //   involved with the charged current.

    /** For Glashow Resonance **/
    protected final static int GR_FLAG = 256; // Allows the Glashow Resonance

    protected final static int MUDECAY_FLAG = 1; // Allows the mu decay.
    protected final static int TAUDECAY_FLAG = 2;// Allows the tau decay.

    /** For Glashow Resonance **/
    //protected final static int ALL_FLAG = 255; 
    protected final static int ALL_FLAG = 511; 
                    // involving all the intereaction and decay channels;

    // Directory path for dumped InteractionMatrix objects.
    protected String[] pathName = {
	"iceCube/uhe/interactions/ice/","iceCube/uhe/interactions/rock/"
    };


    // Week IntereactionMatrix generated by MakeNeutrinoChargeMtx.java
    // and MakeNeutrinoNeutralMtx.java. The Interactions Object is
    // NeutrinoCharge.java and NeutrinoNeutral.java 
    private InteractionsMatrix nuCCMtx = null; 
    private static String nuCCMtxObjectCCH5File = "ENeutrinoChargeMtx";
    private InteractionsMatrix nuNCMtx = null; 
    private static String nuNCMtxObjectCCH5File = "ENeutrinoNeutralMtx";
    private double neutrinoFactor = 1.0; // Neutrino CC enhancement Factor

    // Pair Creation InteractionMatrix generated by MakePairCreationFitMtx.java.
    // The Interaction Object is PairCreationFit.java
    private InteractionsMatrix muToEPairCMtx = null; 
    private String muToEPairCMtxObjectFile = "muToePairCreationFitMtx";  // e+e-
    private InteractionsMatrix tauToEPairCMtx = null; 
    private String tauToEPairCMtxObjectFile = "tauToePairCreationFitMtx"; // e+e-
    private InteractionsMatrix muToMuPairCMtx = null; 
    private String muToMuPairCMtxObjectFile = "muTomuPairCreationFitMtx"; // mu+mu-
    private InteractionsMatrix tauToMuPairCMtx = null; 
    private String tauToMuPairCMtxObjectFile = "tauTomuPairCreationFitMtx";// mu+mu-
    private InteractionsMatrix muToTauPairCMtx = null; 
    private String muToTauPairCMtxObjectFile = "muTotauPairCreationFitMtx";// tau+tau-
    private InteractionsMatrix tauToTauPairCMtx = null; 
    private String tauToTauPairCMtxObjectFile = "tauTotauPairCreationFitMtx";// tau+tau-

    // Bremsstrahlung InteractionMatrix generated by MakeBremsstrahlungMtx.java.
    // The Interaction Object is Bremsstrahlung.java
    private InteractionsMatrix muBremssMtx = null; 
    private String muBremssMtxObjectFile = "muBremsstrahlungMtx";
    private InteractionsMatrix tauBremssMtx = null; 
    private String tauBremssMtxObjectFile = "tauBremsstrahlungMtx";

    // Knock-on Electron InteractionMatrix generated by MakeKnockOnElectronsMtx.java.
    // The Interaction Object is KnockOnElectrons.java
    private InteractionsMatrix muKnockOnMtx = null; 
    private String muKnockOnMtxObjectFile = "muKnockOnElectronsMtx";
    private InteractionsMatrix tauKnockOnMtx = null; 
    private String tauKnockOnMtxObjectFile = "tauKnockOnElectronsMtx";

    // Photo-nuclear InteractionMatrix generated by MakePhotoNuclearFitMtx.java
    // The Interaction Object is PhotoNuclearFit.java
    private InteractionsMatrix muPhotoNuclMtx = null; 
    private String muPhotoNuclMtxObjectFile = "muPhotoNuclearMtx";
    private InteractionsMatrix tauPhotoNuclMtx = null;
    private String tauPhotoNuclMtxObjectFile = "tauPhotoNuclearMtx";


    /** For Glashow Resonance -begin **/

    // Leptonic Glashow Resonance InteractionMatrix generated by MakeGlashowResonanceLeptonicMtx.java
    // The Interaction Object is GlashowResonanceLeptonic.java
    private InteractionsMatrix grLeptonMtx = null;
    private String grLeptonMtxObjectFile = "eGlashowResonanceLeptonicMtx";

    // Hadronic Glashow Resonance InteractionMatrix generated by MakeGlashowResonanceLeptonicMtx.java
    // The Interaction Object is GlashowResonanceHadronic.java
    private InteractionsMatrix grHadronMtx = null;
    private String grHadronMtxObjectFile = "glashowResonanceHadronicMtx";

    /** For Glashow Resonance -end **/


    // DecayMatrix of mu's and tau's
    private MuDecayMatrix muDecayMtx = null;
    private TauDecayMatrix tauDecayMtx = null;



    // Particles
    Particle nuE = null;   // Electron-Neutrino
    Particle nuMu = null;  // Muon-Neutrino
    Particle nuTau = null; // Tau-Neutrino
    Particle e = null;     // Electrons
    Particle mu = null;    // Muons
    Particle tau = null;   // Taus
    Particle pi = null;    // pions (hadrons)

    // Particle Point
    ParticlePoint s;// Point of the propagation


    // Switches for the interactions and decays involved
    int interactionsSwitch;
    int decaySwitch;



    /** Constructor. Reading all the InteractiosMatrix objects
	and generating the DecayMatrix objects. The infinitesimal propagation 
	distance dX [g/cm^2] is also determined here.
	String nuCCMtxObjectFile, String nuNCMtxObjectFile in the arguments*/
    public PropagationMatrix(Particle nuE, Particle nuMu, Particle nuTau,
			     Particle e,   Particle mu,   Particle tau,
			     Particle pi,  ParticlePoint s,
			     int interactionsSwitch, int decaySwitch,
			     double neutrinoFactor,
			     String nuCCMtxObjectFile, String nuNCMtxObjectFile) 
    throws IOException{

	this.interactionsSwitch = interactionsSwitch;
	this.decaySwitch = decaySwitch;
	this.s = s;
	this.neutrinoFactor = neutrinoFactor;
	System.err.println("Will construct PropagationMatrix. ");

	// Checking the Particle flavors
	if(nuE.getFlavor( )== 0 &&  nuE.getDoublet( ) == 0){// Check nuE
	    System.err.println("Checked nuE.");
	    this.nuE = nuE;
	}else{
	    System.err.println("This is NOT e-nu !!");
	    System.exit(0);
	}
	if(nuMu.getFlavor( )== 1 &&  nuMu.getDoublet( ) == 0){// Check nuMu
	    this.nuMu = nuMu;
	}else{
	    System.err.println("This is NOT mu-nu !!");
	    System.exit(0);
	}
	if(nuTau.getFlavor( )== 2 &&  nuTau.getDoublet( ) == 0){// Check nuTau
	    this.nuTau = nuTau;
	}else{
	    System.err.println("This is NOT tau-nu !!");
	    System.exit(0);
	}

	if(e.getFlavor( )== 0 &&  e.getDoublet( ) == 1){// Check e
	    this.e = e;
	}else{
	    System.err.println("This is NOT e !!");
	    System.exit(0);
	}
	if(mu.getFlavor( )== 1 &&  mu.getDoublet( ) == 1){// Check mu
	    this.mu = mu;
	}else{
	    System.err.println("This is NOT mu-nu !!");
	    System.exit(0);
	}
	if(tau.getFlavor( )== 2 &&  tau.getDoublet( ) == 1){// Check Tau
	    this.tau = tau;
	}else{
	    System.err.println("This is NOT tau !!");
	    System.exit(0);
	}

	if(pi.getFlavor( )== 3 &&  pi.getDoublet( ) == 1){// Check pi
	    this.pi = pi;
	}else{
	    System.err.println("This is NOT pi !!");
	    System.exit(0);
	}

	    


	// Reading the InteractionMatrix objects

	// the Charged Current
	if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(nuCCMtxObjectFile);
	    System.err.println("Set object file name." + objectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    System.err.println("Set input stream.");
	    nuCCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    System.err.println("Set interaction Matrix.");
	    in.close( );
	    System.err.println("Closed input file.");
	}else{
	    System.err.println("Charged Current Interactions switched off.");
	}

	// the Neutral Current
	if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(nuNCMtxObjectFile);
	    System.err.println("Set object file name." + objectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    nuNCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Neutral Current Interactions switched off.");
	}


	// the e+e- Pair Creation
	if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(muToEPairCMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muToEPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(tauToEPairCMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauToEPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("e+e- Pair Creations switched off.");
	}


	// the mu+mu- and tau+tau- Pair Creation
	if((interactionsSwitch & PAIRCH_FLAG) == PAIRCH_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(muToMuPairCMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muToMuPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( ); 
	    objectFile = pathName[s.getMaterialNumber()].concat(tauToMuPairCMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauToMuPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(muToTauPairCMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muToTauPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(tauToTauPairCMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauToTauPairCMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Heavier Leptons Pair Creations switched off.");
	}


	// Bremsstrahlung
	if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(muBremssMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muBremssMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(tauBremssMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauBremssMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Bremsstrahlung switched off.");
	}


	// KnockOn Electrons
	if((interactionsSwitch & KNOCK_FLAG) == KNOCK_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(muKnockOnMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muKnockOnMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(tauKnockOnMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauKnockOnMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Knock-on Electrons switched off.");
	}


	// Photo-nucelar
	if((interactionsSwitch & PHOTO_FLAG) == PHOTO_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(muPhotoNuclMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    muPhotoNuclMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(tauPhotoNuclMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    tauPhotoNuclMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Photo-Nuclear Interactions switched off.");
	}

	/** For Glashow Resonance -begin **/
	// Glashow Resonance
	if((interactionsSwitch & GR_FLAG) == GR_FLAG){
	    String objectFile = pathName[s.getMaterialNumber()].concat(grLeptonMtxObjectFile);
	    InputStream in = ClassLoader.getSystemResourceAsStream(objectFile);
	    grLeptonMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    objectFile = pathName[s.getMaterialNumber()].concat(grHadronMtxObjectFile);
	    in = ClassLoader.getSystemResourceAsStream(objectFile);
	    grHadronMtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}else{
	    System.err.println("Glashow Resonance Interactions switched off.");
	}
	/** For Glashow Resonance -end **/


	System.err.println("Reading the Interaction Matrix objects done.");

	// Mu-Decay
	if((decaySwitch & MUDECAY_FLAG) == MUDECAY_FLAG){
	    muDecayMtx = new MuDecayMatrix(mu);
	    int iLogE;
	    for(iLogE=0;iLogE<mu.getDimensionOfLogEnergyMatrix();iLogE++){
		// LifeTime [sec]
		muDecayMtx.setLifeTimeMatrix(iLogE);
		// Decay Matrix
		int jLogE;
		for(jLogE=0;jLogE<mu.getDimensionOfLogEnergyMatrix();jLogE++){
		    muDecayMtx.setMuDecayMatrix(iLogE,jLogE);
		}
	    }
	}else{
	    System.err.println("Mu decay switched off.");
	}


	// Tau-Decay
	if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
	    tauDecayMtx = new TauDecayMatrix(tau);
	    int iLogE;
	    for(iLogE=0;iLogE<tau.getDimensionOfLogEnergyMatrix();iLogE++){
		// LifeTime [sec]
		tauDecayMtx.setLifeTimeMatrix(iLogE);
		// Decay Matrix
		int jLogE;
		for(jLogE=0;jLogE<tau.getDimensionOfLogEnergyMatrix();jLogE++){
		    tauDecayMtx.setTauDecayMatrix(iLogE,jLogE);
		}

	    }
	}else{
	    System.err.println("Tau decay switched off.");
	}

	System.err.println("Decay Matrix calculation done.");

	// Generate the propagation matrix

	/** For Glashow Resonance **/
	System.err.println("Will generate probMtx for nue done.");
	intProbNuE = new double[dimension];
	System.err.println("Generate probMtx for nue done.");

	intProbNeutrino = new double[dimension];
	intProbMu = new double[dimension];
	intProbTau = new double[dimension];
	/** For Glashow Resonance **/
	//nuToNu = new double[dimension][dimension];
        //nuToLepton = new double[dimension][dimension];
	//nuToHadron = new double[dimension][dimension];
	System.err.println("Generate interaction probavirity Mtx done.");

	/** For Glashow Resonance **/
	//temp = new double[31][dimension][dimension];
	temp = new double[35][dimension][dimension];

        nuEToNuE= new double[dimension][dimension];
        nuEToNuMu= new double[dimension][dimension];
        nuEToNuTau= new double[dimension][dimension];
        nuEToE= new double[dimension][dimension];
        nuEToMu= new double[dimension][dimension];
        nuEToTau= new double[dimension][dimension];
        nuEToHadron= new double[dimension][dimension];
	System.err.println("Generate propMtx for nuE done.");
        nuMuToNuE= new double[dimension][dimension];
        nuMuToNuMu= new double[dimension][dimension];
        nuMuToNuTau= new double[dimension][dimension];
        nuMuToE= new double[dimension][dimension];
        nuMuToMu= new double[dimension][dimension];
        nuMuToTau= new double[dimension][dimension];
        nuMuToHadron= new double[dimension][dimension];
	System.err.println("Generate propMtx for nuMu done.");
        nuTauToNuE= new double[dimension][dimension];
        nuTauToNuMu= new double[dimension][dimension];
        nuTauToNuTau= new double[dimension][dimension];
        nuTauToE= new double[dimension][dimension];
        nuTauToMu= new double[dimension][dimension];
        nuTauToTau= new double[dimension][dimension];
        nuTauToHadron= new double[dimension][dimension];
	System.err.println("Generate propMtx for nuTau done.");
        muToNuE= new double[dimension][dimension];
        muToNuMu= new double[dimension][dimension];
        muToNuTau= new double[dimension][dimension];
        muToE= new double[dimension][dimension];
        muToMu= new double[dimension][dimension];
        muToTau= new double[dimension][dimension];
        muToHadron= new double[dimension][dimension];
	System.err.println("Generate propMtx for muon done.");
        tauToNuE= new double[dimension][dimension];
        tauToNuMu= new double[dimension][dimension];
        tauToNuTau= new double[dimension][dimension];
        tauToE= new double[dimension][dimension];
        tauToMu= new double[dimension][dimension];
        tauToTau= new double[dimension][dimension];
        tauToHadron= new double[dimension][dimension];
	System.err.println("Generate propMtx for tau done.");
        FnuEToNuE= new double[dimension][dimension];
        FnuEToNuMu= new double[dimension][dimension];
        FnuEToNuTau= new double[dimension][dimension];
        FnuEToE= new double[dimension][dimension];
        FnuEToMu= new double[dimension][dimension];
        FnuEToTau= new double[dimension][dimension];
        FnuEToHadron= new double[dimension][dimension];
        FnuMuToNuE= new double[dimension][dimension];
        FnuMuToNuMu= new double[dimension][dimension];
        FnuMuToNuTau= new double[dimension][dimension];
        FnuMuToE= new double[dimension][dimension];
        FnuMuToMu= new double[dimension][dimension];
        FnuMuToTau= new double[dimension][dimension];
        FnuMuToHadron= new double[dimension][dimension];
        FnuTauToNuE= new double[dimension][dimension];
        FnuTauToNuMu= new double[dimension][dimension];
        FnuTauToNuTau= new double[dimension][dimension];
        FnuTauToE= new double[dimension][dimension];
        FnuTauToMu= new double[dimension][dimension];
        FnuTauToTau= new double[dimension][dimension];
        FnuTauToHadron= new double[dimension][dimension];
        FmuToNuE= new double[dimension][dimension];
        FmuToNuMu= new double[dimension][dimension];
        FmuToNuTau= new double[dimension][dimension];
        FmuToE= new double[dimension][dimension];
        FmuToMu= new double[dimension][dimension];
        FmuToTau= new double[dimension][dimension];
        FmuToHadron= new double[dimension][dimension];
        FtauToNuE= new double[dimension][dimension];
        FtauToNuMu= new double[dimension][dimension];
        FtauToNuTau= new double[dimension][dimension];
        FtauToE= new double[dimension][dimension];
        FtauToMu= new double[dimension][dimension];
        FtauToTau= new double[dimension][dimension];
        FtauToHadron= new double[dimension][dimension];
        SnuEToNuE= new double[dimension][dimension];
        SnuEToNuMu= new double[dimension][dimension];
        SnuEToNuTau= new double[dimension][dimension];
        SnuEToE= new double[dimension][dimension];
        SnuEToMu= new double[dimension][dimension];
        SnuEToTau= new double[dimension][dimension];
        SnuEToHadron= new double[dimension][dimension];
        SnuMuToNuE= new double[dimension][dimension];
        SnuMuToNuMu= new double[dimension][dimension];
        SnuMuToNuTau= new double[dimension][dimension];
        SnuMuToE= new double[dimension][dimension];
        SnuMuToMu= new double[dimension][dimension];
        SnuMuToTau= new double[dimension][dimension];
        SnuMuToHadron= new double[dimension][dimension];
        SnuTauToNuE= new double[dimension][dimension];
        SnuTauToNuMu= new double[dimension][dimension];
        SnuTauToNuTau= new double[dimension][dimension];
        SnuTauToE= new double[dimension][dimension];
        SnuTauToMu= new double[dimension][dimension];
        SnuTauToTau= new double[dimension][dimension];
        SnuTauToHadron= new double[dimension][dimension];
        SmuToNuE= new double[dimension][dimension];
        SmuToNuMu= new double[dimension][dimension];
        SmuToNuTau= new double[dimension][dimension];
        SmuToE= new double[dimension][dimension];
        SmuToMu= new double[dimension][dimension];
        SmuToTau= new double[dimension][dimension];
        SmuToHadron= new double[dimension][dimension];
        StauToNuE= new double[dimension][dimension];
        StauToNuMu= new double[dimension][dimension];
        StauToNuTau= new double[dimension][dimension];
        StauToE= new double[dimension][dimension];
        StauToMu= new double[dimension][dimension];
        StauToTau= new double[dimension][dimension];
        StauToHadron= new double[dimension][dimension];


	// Matrix Initialization.
	initALL( );
	System.err.println("Matrix Initialization done.");

	// total atomic number in the initial medium
        for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
            massNumber += s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
        }

	// Propagation Step size
	if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
	    dX = 1.0e-1/(s.NA/massNumber*muToEPairCMtx.getSigmaMatrix(dimension-1));
	    // Pair Creation  mean free path
	}else if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
	    dX = 1.0e-1/(s.NA/massNumber*muBremssMtx.getSigmaMatrix(dimension-1));
	    // Bremsstrahlung  mean free path
	}else if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
	    dX = 1.0e-2/(s.NA*nuCCMtx.getSigmaMatrix(dimension-1)*neutrinoFactor);
	    // Charged current interaction
	}else if((interactionsSwitch & GR_FLAG) == GR_FLAG){
	    /** For Glashow Resonance 
		dX is tempolarily changed to smaller value [g/cm^2] **/
	    dX = 1.0e-6/(s.NA*grLeptonMtx.getSigmaMatrix(dimResonance));
	    //dX = 1.0;
	    // Glashow Resonance 
	}else{
	    dX = 1.0e6;
	}

	if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
	    s.setIceRockBoundaryRadius(1.01*s.REarth);//Rock
	    dXDecay = 0.1*c*tauDecayMtx.getLifeTimeMatrix(0)*s.getMediumDensity( );
	}else{
	    dXDecay = 1.0e6;
	}
	System.err.println("dX = " + dX + " dXDecay = " + dXDecay);
	if(dXDecay<dX) dX = dXDecay;
	
	calculateTransferMatrix( );
	System.err.println("Elementary transfer matrix calculation done");

    }

    /** Constructor. Reading all the InteractiosMatrix objects
	and generating the DecayMatrix objects. The infinitesimal propagation 
	distance dX [g/cm^2] is also determined here. The neutrino-nuclen interaction files are 
	given by the default: CTECH5 */
    public PropagationMatrix(Particle nuE, Particle nuMu, Particle nuTau,
			     Particle e,   Particle mu,   Particle tau,
			     Particle pi,  ParticlePoint s,
			     int interactionsSwitch, int decaySwitch,
			     double neutrinoFactor)
    throws IOException{
	this(nuE, nuMu, nuTau, e, mu, tau, pi, s,interactionsSwitch, decaySwitch, neutrinoFactor,
	     nuCCMtxObjectCCH5File,nuNCMtxObjectCCH5File);
    }


    /** Constructor. Reading all the InteractiosMatrix objects
	and generating the DecayMatrix objects. The infinitesimal propagation 
	distance dX [g/cm^2] is also determined here. The neutrino-nuclen interaction files are given by
	String nuCCMtxObjectFile, String nuNCMtxObjectFile in the arguments. The default neutrino factor*/
    public PropagationMatrix(Particle nuE, Particle nuMu, Particle nuTau,
			     Particle e,   Particle mu,   Particle tau,
			     Particle pi,  ParticlePoint s,
			     int interactionsSwitch, int decaySwitch,
			     String nuCCMtxObjectFile, String nuNCMtxObjectFile) throws IOException{
	this(nuE, nuMu, nuTau, e, mu, tau, pi, s,interactionsSwitch, decaySwitch, 1.0,
	     nuCCMtxObjectFile,nuNCMtxObjectFile);
    }

    /** Constructor. Reading all the InteractiosMatrix objects
	and generating the DecayMatrix objects. The infinitesimal propagation 
	distance dX [g/cm^2] is also determined here.
	Use the default value of the neutrino CC/NC enhancement foctor and the default
	neutrino-nucleon intertaction matrix, CTECH5
    */
    public PropagationMatrix(Particle nuE, Particle nuMu, Particle nuTau,
			     Particle e,   Particle mu,   Particle tau,
			     Particle pi,  ParticlePoint s,
			     int interactionsSwitch, int decaySwitch) throws IOException{
	this(nuE, nuMu, nuTau, e, mu, tau, pi, s,interactionsSwitch, decaySwitch, 1.0,
	     nuCCMtxObjectCCH5File,nuNCMtxObjectCCH5File);
    }

    /** Initialize the propagation matrices. */
    public void init( ){
	// Initialization
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<dimension;jLogE++){
		FnuEToNuE[iLogE][jLogE]= 0.0;
		FnuEToNuMu[iLogE][jLogE]= 0.0;
		FnuEToNuTau[iLogE][jLogE]= 0.0;
		FnuEToE[iLogE][jLogE]= 0.0;
		FnuEToMu[iLogE][jLogE]= 0.0;
		FnuEToTau[iLogE][jLogE]= 0.0;
		FnuEToHadron[iLogE][jLogE]= 0.0;
		FnuMuToNuE[iLogE][jLogE]= 0.0;
		FnuMuToNuMu[iLogE][jLogE]= 0.0;
		FnuMuToNuTau[iLogE][jLogE]= 0.0;
		FnuMuToE[iLogE][jLogE]= 0.0;
		FnuMuToMu[iLogE][jLogE]= 0.0;
		FnuMuToTau[iLogE][jLogE]= 0.0;
		FnuMuToHadron[iLogE][jLogE]= 0.0;
		FnuTauToNuE[iLogE][jLogE]= 0.0;
		FnuTauToNuMu[iLogE][jLogE]= 0.0;
		FnuTauToNuTau[iLogE][jLogE]= 0.0;
		FnuTauToE[iLogE][jLogE]= 0.0;
		FnuTauToMu[iLogE][jLogE]= 0.0;
		FnuTauToTau[iLogE][jLogE]= 0.0;
		FnuTauToHadron[iLogE][jLogE]= 0.0;
		FmuToNuE[iLogE][jLogE]= 0.0;
		FmuToNuMu[iLogE][jLogE]= 0.0;
		FmuToNuTau[iLogE][jLogE]= 0.0;
		FmuToE[iLogE][jLogE]= 0.0;
		FmuToMu[iLogE][jLogE]= 0.0;
		FmuToTau[iLogE][jLogE]= 0.0;
		FmuToHadron[iLogE][jLogE]= 0.0;
		FtauToNuE[iLogE][jLogE]= 0.0;
		FtauToNuMu[iLogE][jLogE]= 0.0;
		FtauToNuTau[iLogE][jLogE]= 0.0;
		FtauToE[iLogE][jLogE]= 0.0;
		FtauToMu[iLogE][jLogE]= 0.0;
		FtauToTau[iLogE][jLogE]= 0.0;
		FtauToHadron[iLogE][jLogE]= 0.0;
		
	    }
	}
	for(iLogE=0;iLogE<dimension;iLogE++){
	    FnuEToNuE[iLogE][iLogE] = 1.0;
	    FnuMuToNuMu[iLogE][iLogE] = 1.0;
	    FnuTauToNuTau[iLogE][iLogE] = 1.0;
	    FmuToMu[iLogE][iLogE] = 1.0;
	    FtauToTau[iLogE][iLogE] = 1.0;
	}
    }


    /** Initialize ALL the propagation matrices including the store matrix. */
    public void initALL( ){
	// Initialization
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<dimension;jLogE++){
		FnuEToNuE[iLogE][jLogE]= 0.0;
		FnuEToNuMu[iLogE][jLogE]= 0.0;
		FnuEToNuTau[iLogE][jLogE]= 0.0;
		FnuEToE[iLogE][jLogE]= 0.0;
		FnuEToMu[iLogE][jLogE]= 0.0;
		FnuEToTau[iLogE][jLogE]= 0.0;
		FnuEToHadron[iLogE][jLogE]= 0.0;
		FnuMuToNuE[iLogE][jLogE]= 0.0;
		FnuMuToNuMu[iLogE][jLogE]= 0.0;
		FnuMuToNuTau[iLogE][jLogE]= 0.0;
		FnuMuToE[iLogE][jLogE]= 0.0;
		FnuMuToMu[iLogE][jLogE]= 0.0;
		FnuMuToTau[iLogE][jLogE]= 0.0;
		FnuMuToHadron[iLogE][jLogE]= 0.0;
		FnuTauToNuE[iLogE][jLogE]= 0.0;
		FnuTauToNuMu[iLogE][jLogE]= 0.0;
		FnuTauToNuTau[iLogE][jLogE]= 0.0;
		FnuTauToE[iLogE][jLogE]= 0.0;
		FnuTauToMu[iLogE][jLogE]= 0.0;
		FnuTauToTau[iLogE][jLogE]= 0.0;
		FnuTauToHadron[iLogE][jLogE]= 0.0;
		FmuToNuE[iLogE][jLogE]= 0.0;
		FmuToNuMu[iLogE][jLogE]= 0.0;
		FmuToNuTau[iLogE][jLogE]= 0.0;
		FmuToE[iLogE][jLogE]= 0.0;
		FmuToMu[iLogE][jLogE]= 0.0;
		FmuToTau[iLogE][jLogE]= 0.0;
		FmuToHadron[iLogE][jLogE]= 0.0;
		FtauToNuE[iLogE][jLogE]= 0.0;
		FtauToNuMu[iLogE][jLogE]= 0.0;
		FtauToNuTau[iLogE][jLogE]= 0.0;
		FtauToE[iLogE][jLogE]= 0.0;
		FtauToMu[iLogE][jLogE]= 0.0;
		FtauToTau[iLogE][jLogE]= 0.0;
		FtauToHadron[iLogE][jLogE]= 0.0;
		
		SnuEToNuE[iLogE][jLogE]= 0.0;
		SnuEToNuMu[iLogE][jLogE]= 0.0;
		SnuEToNuTau[iLogE][jLogE]= 0.0;
		SnuEToE[iLogE][jLogE]= 0.0;
		SnuEToMu[iLogE][jLogE]= 0.0;
		SnuEToTau[iLogE][jLogE]= 0.0;
		SnuEToHadron[iLogE][jLogE]= 0.0;
		SnuMuToNuE[iLogE][jLogE]= 0.0;
		SnuMuToNuMu[iLogE][jLogE]= 0.0;
		SnuMuToNuTau[iLogE][jLogE]= 0.0;
		SnuMuToE[iLogE][jLogE]= 0.0;
		SnuMuToMu[iLogE][jLogE]= 0.0;
		SnuMuToTau[iLogE][jLogE]= 0.0;
		SnuMuToHadron[iLogE][jLogE]= 0.0;
		SnuTauToNuE[iLogE][jLogE]= 0.0;
		SnuTauToNuMu[iLogE][jLogE]= 0.0;
		SnuTauToNuTau[iLogE][jLogE]= 0.0;
		SnuTauToE[iLogE][jLogE]= 0.0;
		SnuTauToMu[iLogE][jLogE]= 0.0;
		SnuTauToTau[iLogE][jLogE]= 0.0;
		SnuTauToHadron[iLogE][jLogE]= 0.0;
		SmuToNuE[iLogE][jLogE]= 0.0;
		SmuToNuMu[iLogE][jLogE]= 0.0;
		SmuToNuTau[iLogE][jLogE]= 0.0;
		SmuToE[iLogE][jLogE]= 0.0;
		SmuToMu[iLogE][jLogE]= 0.0;
		SmuToTau[iLogE][jLogE]= 0.0;
		SmuToHadron[iLogE][jLogE]= 0.0;
		StauToNuE[iLogE][jLogE]= 0.0;
		StauToNuMu[iLogE][jLogE]= 0.0;
		StauToNuTau[iLogE][jLogE]= 0.0;
		StauToE[iLogE][jLogE]= 0.0;
		StauToMu[iLogE][jLogE]= 0.0;
		StauToTau[iLogE][jLogE]= 0.0;
		StauToHadron[iLogE][jLogE]= 0.0;
		
	    }
	}
	for(iLogE=0;iLogE<dimension;iLogE++){
	    FnuEToNuE[iLogE][iLogE] = 1.0;
	    FnuMuToNuMu[iLogE][iLogE] = 1.0;
	    FnuTauToNuTau[iLogE][iLogE] = 1.0;
	    FmuToMu[iLogE][iLogE] = 1.0;
	    FtauToTau[iLogE][iLogE] = 1.0;

	    SnuEToNuE[iLogE][iLogE] = 1.0;
	    SnuMuToNuMu[iLogE][iLogE] = 1.0;
	    SnuTauToNuTau[iLogE][iLogE] = 1.0;
	    SmuToMu[iLogE][iLogE] = 1.0;
	    StauToTau[iLogE][iLogE] = 1.0;
	}
    }

    /**** Calculate the elementary interaction/decay transfer matrix. 
	  The mass density of the proagation medium is approxed to be
	  constant and equal to the cuurent particle location
	  for saving CPU time. 
	  The infinitesimal propagation over dX is a subject to consider
	  here. All the calculated matrix elements are stored
	  in nuToNu[kLogE][jLogE], tauToNuE[kLogE][jLogE] etc.
	  The bin width and the relation between energy index kLogE (integer)
	  and energy are defined in Particle.class as usual since it relies
	  on InteractionMatrix object which is pre-calculated and stored
	  in the file. InteractionMatrix object for each interaction
	  is read by the constructor PropagationMatrix( ).
    */
    public void calculateTransferMatrix( ){
	int kLogE;
	for(kLogE=0;kLogE<dimension;kLogE++){
	    // total intereaction probability including decay.
	    /** For Glashow Resonance **/
	    intProbNuE[kLogE] = 0.0;
	    intProbNeutrino[kLogE] = 0.0;
	    intProbMu[kLogE] = 0.0;
	    intProbTau[kLogE] = 0.0;

	    if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
		intProbNeutrino[kLogE] += dX*s.NA*nuCCMtx.getSigmaMatrix(kLogE)*
		    neutrinoFactor;
	    }
	    if((interactionsSwitch & LEPTW_FLAG) == LEPTW_FLAG){
		intProbMu[kLogE] += dX*s.NA*nuCCMtx.getSigmaMatrix(kLogE)*neutrinoFactor;
		intProbTau[kLogE] += dX*s.NA*nuCCMtx.getSigmaMatrix(kLogE)*neutrinoFactor;
	    }
	    if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
		intProbNeutrino[kLogE] += dX*s.NA*nuNCMtx.getSigmaMatrix(kLogE)*
		    neutrinoFactor;
	    }
	    if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
		intProbMu[kLogE] += dX*s.NA/massNumber*muToEPairCMtx.getSigmaMatrix(kLogE);
		intProbTau[kLogE] += dX*s.NA/massNumber*tauToEPairCMtx.getSigmaMatrix(kLogE);
	    }
	    if((interactionsSwitch & PAIRCH_FLAG) == PAIRCH_FLAG){
		intProbMu[kLogE] += 
		    dX*s.NA/massNumber*(
					muToMuPairCMtx.getSigmaMatrix(kLogE)+
					muToTauPairCMtx.getSigmaMatrix(kLogE));
		intProbTau[kLogE] += 
		    dX*s.NA/massNumber*(
					tauToMuPairCMtx.getSigmaMatrix(kLogE)+
					tauToTauPairCMtx.getSigmaMatrix(kLogE));
	    }
	    if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
		intProbMu[kLogE] += dX*s.NA/massNumber*muBremssMtx.getSigmaMatrix(kLogE);
		intProbTau[kLogE] += dX*s.NA/massNumber*tauBremssMtx.getSigmaMatrix(kLogE);
	    }
	    if((interactionsSwitch & KNOCK_FLAG) == KNOCK_FLAG){
		intProbMu[kLogE] += dX*s.NA/massNumber*muKnockOnMtx.getSigmaMatrix(kLogE);
		intProbTau[kLogE] += dX*s.NA/massNumber*tauKnockOnMtx.getSigmaMatrix(kLogE);
	    }
	    if((interactionsSwitch & PHOTO_FLAG) == PHOTO_FLAG){
		intProbMu[kLogE] += dX*s.NA/massNumber*muPhotoNuclMtx.getSigmaMatrix(kLogE);
		intProbTau[kLogE] += dX*s.NA/massNumber*tauPhotoNuclMtx.getSigmaMatrix(kLogE);
	    }
	    /** For Glashow Resonance -begin **/
	    if((interactionsSwitch & GR_FLAG) == GR_FLAG){
		intProbNuE[kLogE] += 3.0*0.5*dX*s.NA*grLeptonMtx.getSigmaMatrix(kLogE);
		intProbNuE[kLogE] += 0.5*dX*s.NA*grHadronMtx.getSigmaMatrix(kLogE);
	    }
	    /** For Glashow Resonance -end **/
	    if((decaySwitch & MUDECAY_FLAG) == MUDECAY_FLAG){
		intProbMu[kLogE] += 
		    dX/(c*muDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
	    }
	    if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		intProbTau[kLogE] += 
		    dX/(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
	    }


	    int jLogE;
	    for(jLogE=0;jLogE<=kLogE;jLogE++){


		// To NuE
		nuEToNuE[kLogE][jLogE] = 0.0;
                if((interactionsSwitch & GR_FLAG) == GR_FLAG){
                    nuEToNuE[kLogE][jLogE] += 0.5*
                        dX*s.NA*grLeptonMtx.getTransferMatrix(kLogE,jLogE);
                }
		if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
		    nuEToNuE[kLogE][jLogE] += 
			dX*s.NA*nuNCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		muToNuE[kLogE][jLogE] = 0.0;
		if((decaySwitch & MUDECAY_FLAG) == MUDECAY_FLAG){
		    muToNuE[kLogE][jLogE] +=
  		        dX*muDecayMtx.getMuToNuEDecayMatrix(kLogE,jLogE)
			/(c*muDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		tauToNuE[kLogE][jLogE] = 0.0;
		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToNuE[kLogE][jLogE] += 
		        dX*tauDecayMtx.getTauToNuDecayMatrix(kLogE,jLogE)
		        /(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		// To NuMu
		nuEToNuMu[kLogE][jLogE] = 0.0;
                if((interactionsSwitch & GR_FLAG) == GR_FLAG){
                    nuEToNuMu[kLogE][jLogE] += 0.5*
                        dX*s.NA*grLeptonMtx.getTransferMatrix(kLogE,jLogE);
                }

		nuMuToNuMu[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
		    nuMuToNuMu[kLogE][jLogE] += 
			dX*s.NA*nuNCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		muToNuMu[kLogE][jLogE] = 0.0;
		if((decaySwitch & MUDECAY_FLAG) == MUDECAY_FLAG){
		    muToNuMu[kLogE][jLogE] +=
  		        dX*muDecayMtx.getMuToNuMuDecayMatrix(kLogE,jLogE)
			/(c*muDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}
		if((interactionsSwitch & LEPTW_FLAG) == LEPTW_FLAG){
		    muToNuMu[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		tauToNuMu[kLogE][jLogE] = 0.0;
		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToNuMu[kLogE][jLogE] += 
		        dX*tauDecayMtx.getTauToNuDecayMatrix(kLogE,jLogE)
		        /(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		// To NuTau
		nuEToNuTau[kLogE][jLogE] = 0.0;
		nuEToNuTau[kLogE][jLogE] += nuEToNuMu[kLogE][jLogE];

		nuTauToNuTau[kLogE][jLogE] = 0.0;
		nuTauToNuTau[kLogE][jLogE] += nuMuToNuMu[kLogE][jLogE];

		tauToNuTau[kLogE][jLogE] = 0.0;
		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToNuTau[kLogE][jLogE] +=
  		        dX*tauDecayMtx.getTauToNuTauDecayMatrix(kLogE,jLogE)
			/(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}
		if((interactionsSwitch & LEPTW_FLAG) == LEPTW_FLAG){
		    tauToNuTau[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		// To E
		nuEToE[kLogE][jLogE] = 0.0;
                if((interactionsSwitch & GR_FLAG) == GR_FLAG){
                    nuEToE[kLogE][jLogE] += 0.5*
                        dX*s.NA*grLeptonMtx.getLeptonTransferMatrix(kLogE,jLogE);
                }
		if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
		    nuEToE[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		muToE[kLogE][jLogE] = 0.0; tauToE[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
		    muToE[kLogE][jLogE] += 
			dX*s.NA/massNumber*muBremssMtx.getTransferMatrix(kLogE,jLogE);
		    tauToE[kLogE][jLogE] += 
			dX*s.NA/massNumber*tauBremssMtx.getTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & KNOCK_FLAG) == KNOCK_FLAG){
		    muToE[kLogE][jLogE] += 
			dX*s.NA/massNumber*muKnockOnMtx.getTransferMatrix(kLogE,jLogE);
		    tauToE[kLogE][jLogE] += 
			dX*s.NA/massNumber*tauKnockOnMtx.getTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
		    muToE[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*muToEPairCMtx.getTransferMatrix(kLogE,jLogE);
		    tauToE[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*tauToEPairCMtx.getTransferMatrix(kLogE,jLogE);
		}

		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToE[kLogE][jLogE] += 
		        dX*tauDecayMtx.getTauToChargedLeptonDecayMatrix(kLogE,jLogE)
		        /(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		if((decaySwitch & MUDECAY_FLAG) == MUDECAY_FLAG){
		    muToE[kLogE][jLogE] +=
  		        dX*muDecayMtx.getMuToEDecayMatrix(kLogE,jLogE)
			/(c*muDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		// To Mu
		nuEToMu[kLogE][jLogE] = 0.0;
                if((interactionsSwitch & GR_FLAG) == GR_FLAG){
                    nuEToMu[kLogE][jLogE] += 0.5*
                        dX*s.NA*grLeptonMtx.getLeptonTransferMatrix(kLogE,jLogE);
                }

		nuMuToMu[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
		    nuMuToMu[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getLeptonTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}
		nuTauToMu[kLogE][jLogE] = 0.0; // reserved for the new physics

		muToMu[kLogE][jLogE] = 0.0; tauToMu[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & PAIRCH_FLAG) == PAIRCH_FLAG){
		    muToMu[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*muToMuPairCMtx.getTransferMatrix(kLogE,jLogE);
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muToMuPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muToTauPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		    tauToMu[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*tauToMuPairCMtx.getTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muToEPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muBremssMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & KNOCK_FLAG) == KNOCK_FLAG){
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muKnockOnMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & PHOTO_FLAG) == PHOTO_FLAG){
		    muToMu[kLogE][jLogE] += dX*s.NA/massNumber*
			muPhotoNuclMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToMu[kLogE][jLogE] += 
			dX*tauDecayMtx.getTauToChargedLeptonDecayMatrix(kLogE,jLogE)/
			(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}

		// To Tau
		nuEToTau[kLogE][jLogE] = 0.0;
		nuEToTau[kLogE][jLogE] += nuEToMu[kLogE][jLogE];

		nuMuToTau[kLogE][jLogE] = 0.0; // reserved for the new physics

		nuTauToTau[kLogE][jLogE] = 0.0;
		nuTauToTau[kLogE][jLogE] += nuMuToMu[kLogE][jLogE];

		muToTau[kLogE][jLogE] = 0.0; tauToTau[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & PAIRCH_FLAG) == PAIRCH_FLAG){
		    muToTau[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*muToTauPairCMtx.getTransferMatrix(kLogE,jLogE);
		    tauToTau[kLogE][jLogE] += 2.0*
			dX*s.NA/massNumber*tauToTauPairCMtx.getTransferMatrix(kLogE,jLogE);
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauToTauPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauToMuPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & PAIRC_FLAG) == PAIRC_FLAG){
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauToEPairCMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & BREMSS_FLAG) == BREMSS_FLAG){
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauBremssMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & KNOCK_FLAG) == KNOCK_FLAG){
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauKnockOnMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		if((interactionsSwitch & PHOTO_FLAG) == PHOTO_FLAG){
		    tauToTau[kLogE][jLogE] += dX*s.NA/massNumber*
			tauPhotoNuclMtx.getLeptonTransferMatrix(kLogE,jLogE);
		}

		// To Hadron
		nuEToHadron[kLogE][jLogE]= 0.0;
                if((interactionsSwitch & GR_FLAG) == GR_FLAG){
                    nuEToHadron[kLogE][jLogE] += 0.5*
			dX*s.NA*grHadronMtx.getTransferMatrix(kLogE,jLogE);
                }
		if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
		    nuEToHadron[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}
		if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
		    nuEToHadron[kLogE][jLogE] += 
			dX*s.NA*nuNCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		nuMuToHadron[kLogE][jLogE]= 0.0;
		if((interactionsSwitch & CHARGED_FLAG) == CHARGED_FLAG){
		    nuMuToHadron[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}
		if((interactionsSwitch & NEUTRAL_FLAG) == NEUTRAL_FLAG){
		    /** Modified ! **/
		    nuMuToHadron[kLogE][jLogE] += 
			dX*s.NA*nuNCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}

		nuTauToHadron[kLogE][jLogE] = 0.0;
		nuTauToHadron[kLogE][jLogE] += nuMuToHadron[kLogE][jLogE];

		muToHadron[kLogE][jLogE] = 0.0;tauToHadron[kLogE][jLogE] = 0.0;
		if((interactionsSwitch & PHOTO_FLAG) == PHOTO_FLAG){
		    muToHadron[kLogE][jLogE] += dX*s.NA/massNumber*
			muPhotoNuclMtx.getTransferMatrix(kLogE,jLogE);
		    tauToHadron[kLogE][jLogE] += dX*s.NA/massNumber*
			tauPhotoNuclMtx.getTransferMatrix(kLogE,jLogE);
		}
		if((decaySwitch & TAUDECAY_FLAG) == TAUDECAY_FLAG){
		    tauToHadron[kLogE][jLogE] += 
		        dX*tauDecayMtx.getTauToHadronDecayMatrix(kLogE,jLogE)
		        /(c*tauDecayMtx.getLifeTimeMatrix(kLogE)*s.getMediumDensity( ));
		}
		if((interactionsSwitch & LEPTW_FLAG) == LEPTW_FLAG){
		    tauToHadron[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		    muToHadron[kLogE][jLogE] += 
			dX*s.NA*nuCCMtx.getTransferMatrix(kLogE,jLogE)*neutrinoFactor;
		}
	    }
	}
    }




    /**** Propagate the particles involved over a dX [g/cm^2].
     The calculation is performed by multiplication of the elemental
     transfer matrix calculated by calculateTransferMatrix( ) called in
     the constructor PropagationMatrix( ). The particle energy distribution
     after propagation over n x dX [g] can be, for example, obtained
     by calling this method n times. The resultant energy distributions
     are stored in FmuToNuE[iLogE][kLogE] etc.
    */
    public void propagateDX( ){

	double dNFromNuE,dNFromNuMu,dNFromNuTau,dNFromMu,dNFromTau;

	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		int kLogE;

		//Interactions and Decays to NuE
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    // NuE to NuE
		    /** For Glashow Resonance **/
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		    //NuMu to NuE
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE] +
			FnuMuToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] +
			FnuMuToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		    //NuTau to NuE
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] + 
			FnuTauToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		    //Mu to NuE
		    dNFromMu += 
			FmuToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE] +
			FmuToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] +
			FmuToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		    //tau to nuE
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE] +
			FtauToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] + 
			FtauToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];

		}

		/** For Glashow Resonance **/
		//dNFromNuE += FnuEToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		dNFromNuE += FnuEToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]-intProbNuE[jLogE]);
		FnuEToNuE[iLogE][jLogE] = dNFromNuE;

	        dNFromNuMu += FnuMuToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]-intProbNuE[jLogE]);
		FnuMuToNuE[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]-intProbNuE[jLogE]);
		FnuTauToNuE[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]-intProbNuE[jLogE]);
		FmuToNuE[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToNuE[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]-intProbNuE[jLogE]);
		FtauToNuE[iLogE][jLogE] = dNFromTau;

		

		//Interactions and Decays to NuMu
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    /** For Glashow Resonance **/
		    //NuE to NuMu
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE] + 
			FnuEToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		    //NuMu to NuMu
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE] +
			FnuMuToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE] +
			FnuMuToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE] +
			FnuMuToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		    //NuTau to NuMu
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE] +
			FnuTauToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE] +
			FnuTauToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE] +
			FnuTauToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		    //Mu to NuMu
		    dNFromMu += 
			FmuToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE] +
			FmuToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE] +
			FmuToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE] +
			FmuToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		    //tau to nuMu
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE] +
			FtauToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE] +
			FtauToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE] +
			FtauToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];

		}

		/** For Glashow Resonance **/
		dNFromNuE += FnuEToNuMu[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuEToNuMu[iLogE][jLogE] = dNFromNuE;

		dNFromNuMu += FnuMuToNuMu[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuMuToNuMu[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToNuMu[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuTauToNuMu[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToNuMu[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FmuToNuMu[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToNuMu[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FtauToNuMu[iLogE][jLogE] = dNFromTau;



		//Interactions and Decays to NuTau
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    /** For Glashow Resonance **/
		    //NuE to NuTau
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];

		    //NuMu to NuTau
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE] +
			FnuMuToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE] +
			FnuMuToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		    //NuTau to NuTau
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE] +
			FnuTauToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE] +
			FnuTauToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		    //Mu to NuTau
		    dNFromMu += 
			FmuToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE] +
			FmuToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE] +
			FmuToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		    //tau to nuTau
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE] +
			FtauToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE] +
			FtauToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];

		}

		/** For Glashow Resonance **/
		dNFromNuE += FnuEToNuTau[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuEToNuTau[iLogE][jLogE] = dNFromNuE;

		dNFromNuMu += FnuMuToNuTau[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuMuToNuTau[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToNuTau[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FnuTauToNuTau[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToNuTau[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FmuToNuTau[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToNuTau[iLogE][jLogE]*(1.0-intProbNeutrino[jLogE]);
		FtauToNuTau[iLogE][jLogE] = dNFromTau;



		//Interactions and Decays to Electrons
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    
		    /** For Glashow Resonance **/
		    // NuE to Electrons
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		    //NuMu to Electrons
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToE[kLogE][jLogE] + 
			FnuMuToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		    //NuTau to Electrons
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToE[kLogE][jLogE] + 
			FnuTauToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		    //Mu to Electrons
		    dNFromMu +=
			FmuToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+ 
			FmuToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		    //tau to Electrons
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+ 
			FtauToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];

		}

		/** For Glashow Resonance **/
	        dNFromNuE += FnuEToE[iLogE][jLogE];
		FnuEToE[iLogE][jLogE] = dNFromNuE;

	        dNFromNuMu += FnuMuToE[iLogE][jLogE];
		FnuMuToE[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToE[iLogE][jLogE];
		FnuTauToE[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToE[iLogE][jLogE];
		FmuToE[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToE[iLogE][jLogE];
		FtauToE[iLogE][jLogE] = dNFromTau;



		//Interactions and Decays to Muons
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    /** For Glashow Resonance **/
		    //NuE to Muons
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToMu[kLogE][jLogE] + 
			FnuEToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];

		    // NuMu to Muons
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE] +
			FnuMuToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToMu[kLogE][jLogE] + 
			FnuMuToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];

		    //NuTau to Muons
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE] +
			FnuTauToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+ // reserved for new physics
			FnuTauToMu[iLogE][kLogE]*muToMu[kLogE][jLogE] + 
			FnuTauToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		    //Mu to Muons
		    dNFromMu += 
			FmuToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE] +
			FmuToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToMu[kLogE][jLogE] +
			FmuToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		    //tau to Muons
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE] +
			FtauToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+ 
			FtauToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+ 
			FtauToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];

		}


		/** For Glashow Resonance **/
		dNFromNuE += FnuEToMu[iLogE][jLogE]*(1.0-intProbMu[jLogE]);
		FnuEToMu[iLogE][jLogE] = dNFromNuE;

	        dNFromNuMu += FnuMuToMu[iLogE][jLogE]*(1.0-intProbMu[jLogE]);
		FnuMuToMu[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToMu[iLogE][jLogE]*(1.0-intProbMu[jLogE]);
		FnuTauToMu[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToMu[iLogE][jLogE]*(1.0-intProbMu[jLogE]);
		FmuToMu[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToMu[iLogE][jLogE]*(1.0-intProbMu[jLogE]);
		FtauToMu[iLogE][jLogE] = dNFromTau;



		//Interactions and Decays to Taus
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    /** For Glashow Resonance **/
		    //NuE to Tauons
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+ 
			FnuEToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];

		    // NuTau to Tauons
		    dNFromNuTau += 
			FnuTauToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE] +
			FnuTauToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+ 
			FnuTauToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];

		    //NuMu to Tauons
		    dNFromNuMu += 
			FnuMuToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE] +
			FnuMuToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+ // reserved for new physics
			FnuMuToMu[iLogE][kLogE]*muToTau[kLogE][jLogE] + 
			FnuMuToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		    //Mu to Tauons
		    dNFromMu += 
			FmuToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE] +
			FmuToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+ 
			FmuToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		    //tau to Tauons
		    dNFromTau += 
			FtauToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE] +
			FtauToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+ 
			FtauToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];

		}

		/** For Glashow Resonance **/
		dNFromNuE += FnuEToTau[iLogE][jLogE]*(1.0-intProbTau[jLogE]);
		FnuEToTau[iLogE][jLogE] = dNFromNuE;

	        dNFromNuMu += FnuMuToTau[iLogE][jLogE]*(1.0-intProbTau[jLogE]);
		FnuMuToTau[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToTau[iLogE][jLogE]*(1.0-intProbTau[jLogE]);
		FnuTauToTau[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToTau[iLogE][jLogE]*(1.0-intProbTau[jLogE]);
		FmuToTau[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToTau[iLogE][jLogE]*(1.0-intProbTau[jLogE]);
		FtauToTau[iLogE][jLogE] = dNFromTau;



		//Interactions and Decays to Hadrons
		dNFromNuE = 0.0;  //From NuE
		dNFromNuMu = 0.0; //From NuMu
		dNFromNuTau = 0.0;//From NuTau
		dNFromMu = 0.0;   //From Muons
		dNFromTau = 0.0;  //From Tauons
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){

		    /** For Glashow Resonance **/
		    // Neutrinos to Hadrons
		    dNFromNuE += 
			FnuEToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] + 
			FnuEToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] + 
			FnuEToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] + 
			FnuEToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];

		    dNFromNuMu +=
			FnuMuToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FnuMuToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] + 
			FnuMuToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] + 
			FnuMuToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] + 
			FnuMuToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];

		    dNFromNuTau +=
			FnuTauToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FnuTauToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] + 
			FnuTauToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] + 
			FnuTauToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] + 
			FnuTauToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];

		    //Mu to Hadrons
		    dNFromMu +=
			FmuToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FmuToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] + 
			FmuToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] + 
			FmuToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] + 
			FmuToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];

		    //tau to Hadrons
		    dNFromTau +=
			FtauToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FtauToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] + 
			FtauToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] + 
			FtauToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] + 
			FtauToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}

	        dNFromNuE += FnuEToHadron[iLogE][jLogE];
		FnuEToHadron[iLogE][jLogE] = dNFromNuE;

	        dNFromNuMu += FnuMuToHadron[iLogE][jLogE];
		FnuMuToHadron[iLogE][jLogE] = dNFromNuMu;

	        dNFromNuTau += FnuTauToHadron[iLogE][jLogE];
		FnuTauToHadron[iLogE][jLogE] = dNFromNuTau;

	        dNFromMu += FmuToHadron[iLogE][jLogE];
		FmuToHadron[iLogE][jLogE] = dNFromMu;

	        dNFromTau += FtauToHadron[iLogE][jLogE];
		FtauToHadron[iLogE][jLogE] = dNFromTau;



	    }
	}

    }




    /**** Propagate the particles involved over a 2^n Delta x [g/cm^2] 
	  where Delta x is the propagation distance 
	  for the finite propagation matrix.
	  "n" is the times to call this method.
	  It doubles the finit propagation matrix (FnuEtoNuE[iLogE][jLog] for example)
	  calculated by propagateDX( ). The resultant matrix element describes
	  the energy distribution of particles.
    */
    public void propagateDXpowered( ){
	int iLogE,jLogE,kLogE;

	/** For Glashow Resonance -begin **/
	// NuE to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[0][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[0][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[1][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[1][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[2][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[2][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[3][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[3][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[4][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[4][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}


	/** For Glashow Resonance -begin **/
	// NuE to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[5][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[5][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[6][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[6][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[7][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[7][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[8][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[8][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[9][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[9][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[10][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[10][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -begin **/

	// NuMu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[11][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[11][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[12][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[12][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[13][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[13][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[14][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[14][iLogE][jLogE] += 
                        /** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}





	/** For Glashow Resonance -begin **/
	// NuE to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[15][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[15][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[15][iLogE][jLogE] += FnuEToE[iLogE][jLogE];
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[16][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[16][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[16][iLogE][jLogE] += FnuMuToE[iLogE][jLogE];
	    }
	}

	// NuTau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[17][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[17][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[17][iLogE][jLogE] += FnuTauToE[iLogE][jLogE];
	    }
	}

	// Mu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[18][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[18][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[18][iLogE][jLogE] += FmuToE[iLogE][jLogE];
	    }
	}

	// Tau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[19][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[19][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[19][iLogE][jLogE] += FtauToE[iLogE][jLogE];
	    }
	}





	/** For Glashow Resonance -begin **/
	// NuE to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[20][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[20][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[21][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[21][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[22][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[22][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// Mu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[23][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[23][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[24][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[24][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[25][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[25][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[26][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[26][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[27][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[27][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[28][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[28][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[29][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[29][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[30][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[30][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[30][iLogE][jLogE] += FnuEToHadron[iLogE][jLogE];
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[31][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[31][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[31][iLogE][jLogE] += FnuMuToHadron[iLogE][jLogE];
	    }
	}

	// NuTau to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[32][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[32][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[32][iLogE][jLogE] += FnuTauToHadron[iLogE][jLogE];
	    }
	}

	// Mu to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[33][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[33][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[33][iLogE][jLogE] += FmuToHadron[iLogE][jLogE];
	    }
	}

	// Tau to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[34][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[34][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[34][iLogE][jLogE] += FtauToHadron[iLogE][jLogE];
	    }
	}


	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		// To NuE
		FnuEToNuE[iLogE][jLogE] = temp[0][iLogE][jLogE];
		FnuMuToNuE[iLogE][jLogE] = temp[1][iLogE][jLogE];
		FnuTauToNuE[iLogE][jLogE] = temp[2][iLogE][jLogE];
		FmuToNuE[iLogE][jLogE] = temp[3][iLogE][jLogE];
		FtauToNuE[iLogE][jLogE] = temp[4][iLogE][jLogE];
		// To NuMu
		/** For Glashow Resonance **/
		FnuEToNuMu[iLogE][jLogE] = temp[5][iLogE][jLogE];
		FnuMuToNuMu[iLogE][jLogE] = temp[6][iLogE][jLogE];
		FnuTauToNuMu[iLogE][jLogE] = temp[7][iLogE][jLogE];
		FmuToNuMu[iLogE][jLogE] = temp[8][iLogE][jLogE];
		FtauToNuMu[iLogE][jLogE] = temp[9][iLogE][jLogE];
		// To NuTau
		/** For Glashow Resonance **/
		FnuEToNuTau[iLogE][jLogE] = temp[10][iLogE][jLogE];
		FnuMuToNuTau[iLogE][jLogE] = temp[11][iLogE][jLogE];
		FnuTauToNuTau[iLogE][jLogE] = temp[12][iLogE][jLogE];
		FmuToNuTau[iLogE][jLogE] = temp[13][iLogE][jLogE];
		FtauToNuTau[iLogE][jLogE] = temp[14][iLogE][jLogE];
		// To E
		FnuEToE[iLogE][jLogE] = temp[15][iLogE][jLogE];
		FnuMuToE[iLogE][jLogE] = temp[16][iLogE][jLogE];
		FnuTauToE[iLogE][jLogE] = temp[17][iLogE][jLogE];
		FmuToE[iLogE][jLogE] = temp[18][iLogE][jLogE];
		FtauToE[iLogE][jLogE] = temp[19][iLogE][jLogE];
		// To Mu
		/** For Glashow Resonance **/
		FnuEToMu[iLogE][jLogE] = temp[20][iLogE][jLogE];
		FnuMuToMu[iLogE][jLogE] = temp[21][iLogE][jLogE];
		FnuTauToMu[iLogE][jLogE] = temp[22][iLogE][jLogE];
		FmuToMu[iLogE][jLogE] = temp[23][iLogE][jLogE];
		FtauToMu[iLogE][jLogE] = temp[24][iLogE][jLogE];
		// To Tau
		/** For Glashow Resonance **/
		FnuEToTau[iLogE][jLogE] = temp[25][iLogE][jLogE];
		FnuMuToTau[iLogE][jLogE] = temp[26][iLogE][jLogE];
		FnuTauToTau[iLogE][jLogE] = temp[27][iLogE][jLogE];
		FmuToTau[iLogE][jLogE] = temp[28][iLogE][jLogE];
		FtauToTau[iLogE][jLogE] = temp[29][iLogE][jLogE];
		// To Hadron
		FnuEToHadron[iLogE][jLogE] = temp[30][iLogE][jLogE];
		FnuMuToHadron[iLogE][jLogE] = temp[31][iLogE][jLogE];
		FnuTauToHadron[iLogE][jLogE] = temp[32][iLogE][jLogE];
		FmuToHadron[iLogE][jLogE] = temp[33][iLogE][jLogE];
		FtauToHadron[iLogE][jLogE] = temp[34][iLogE][jLogE];
	    }
	}

    }


    /**** Propagate the particles involved over X [g/cm^2] 
	  where Delta x is the propagation distance 
	  for the finite propagation matrix.
	  It multiples the finit propagation matrix calculated by
	  propagateDX( ) and proagateDXpowerd( ) and storerd
	  in nuETonuE[][] etc by copyTransferMatrix( ).
    */
    public void propagateX( ){
	int iLogE,jLogE,kLogE;


	/** For Glashow Resonance -begin **/
	// NuE to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[0][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[0][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToNuE[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*nuTauToNuE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[1][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[1][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToNuE[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToNuE[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[2][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[2][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToNuE[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToNuE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[3][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[3][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToNuE[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToNuE[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[4][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[4][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*nuEToNuE[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToNuE[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToNuE[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToNuE[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToNuE[kLogE][jLogE];
		}
	    }
	}


	/** For Glashow Resonance -begin **/
	// NuE to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[5][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[5][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*nuTauToNuMu[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[6][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[6][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToNuMu[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[7][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[7][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToNuMu[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[8][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[8][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToNuMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[9][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[9][iLogE][jLogE] += 
			/** Modified ! **/
			FtauToNuE[iLogE][kLogE]*nuEToNuMu[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToNuMu[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToNuMu[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToNuMu[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToNuMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[10][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[10][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*nuMuToNuTau[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*muToNuTau[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[11][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[11][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToNuTau[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToNuTau[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[12][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[12][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToNuTau[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToNuTau[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[13][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[13][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToNuTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToNuTau[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[14][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[14][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*nuEToNuTau[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToNuTau[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToNuTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToNuTau[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToNuTau[kLogE][jLogE];
		}
	    }
	}





	/** For Glashow Resonance -begin **/
	// NuE to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[15][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[15][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToE[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*nuTauToE[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToE[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		}
		temp[15][iLogE][jLogE] += FnuEToE[iLogE][jLogE];
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[16][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[16][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToE[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToE[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		}
		temp[16][iLogE][jLogE] += FnuMuToE[iLogE][jLogE];
	    }
	}

	// NuTau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[17][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[17][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToE[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToE[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		}
		temp[17][iLogE][jLogE] += FnuTauToE[iLogE][jLogE];
	    }
	}

	// Mu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[18][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[18][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToE[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToE[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		}
		temp[18][iLogE][jLogE] += FmuToE[iLogE][jLogE];
	    }
	}

	// Tau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[19][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[19][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*nuEToE[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToE[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToE[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToE[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToE[kLogE][jLogE];
		}
		temp[19][iLogE][jLogE] += FtauToE[iLogE][jLogE];
	    }
	}





	/** For Glashow Resonance -begin **/
	// NuE to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[20][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[20][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[21][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[21][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		}
	    }
	}

	// NuMu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[22][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[22][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[23][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[23][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[24][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[24][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*nuEToMu[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToMu[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToMu[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToMu[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[25][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[25][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE]+
			FnuEToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+
			FnuEToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuEToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+
			FnuEToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[26][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[26][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuMuToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[27][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[27][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FnuTauToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[28][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[28][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FmuToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[29][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[29][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			FtauToNuE[iLogE][kLogE]*nuEToTau[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToTau[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToTau[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToTau[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToTau[kLogE][jLogE];
		}
	    }
	}




	// NuE to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[30][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[30][iLogE][jLogE] += 
			FnuEToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE] +
			FnuEToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE] +
			FnuEToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE] +
			FnuEToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE] +
			FnuEToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}
		temp[30][iLogE][jLogE] += FnuEToHadron[iLogE][jLogE];
	    }
	}

	// NuMu to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[31][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[31][iLogE][jLogE] += 
			FnuMuToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE]+
			FnuMuToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE]+
			FnuMuToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE]+
			FnuMuToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE]+
			FnuMuToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}
		temp[31][iLogE][jLogE] += FnuMuToHadron[iLogE][jLogE];
	    }
	}

	// NuTau to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[32][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[32][iLogE][jLogE] += 
			FnuTauToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE]+
			FnuTauToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE]+
			FnuTauToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE]+
			FnuTauToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE]+
			FnuTauToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}
		temp[32][iLogE][jLogE] += FnuTauToHadron[iLogE][jLogE];
	    }
	}

	// Mu to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[33][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[33][iLogE][jLogE] += 
			FmuToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE]+
			FmuToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE]+
			FmuToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE]+
			FmuToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE]+
			FmuToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}
		temp[33][iLogE][jLogE] += FmuToHadron[iLogE][jLogE];
	    }
	}

	// Tau to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[34][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[34][iLogE][jLogE] += 
			FtauToNuE[iLogE][kLogE]*nuEToHadron[kLogE][jLogE]+
			FtauToNuMu[iLogE][kLogE]*nuMuToHadron[kLogE][jLogE]+
			FtauToNuTau[iLogE][kLogE]*nuTauToHadron[kLogE][jLogE]+
			FtauToMu[iLogE][kLogE]*muToHadron[kLogE][jLogE]+
			FtauToTau[iLogE][kLogE]*tauToHadron[kLogE][jLogE];
		}
		temp[34][iLogE][jLogE] += FtauToHadron[iLogE][jLogE];
	    }
	}


	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		// To NuE
		FnuEToNuE[iLogE][jLogE] = temp[0][iLogE][jLogE];
		FnuMuToNuE[iLogE][jLogE] = temp[1][iLogE][jLogE];
		FnuTauToNuE[iLogE][jLogE] = temp[2][iLogE][jLogE];
		FmuToNuE[iLogE][jLogE] = temp[3][iLogE][jLogE];
		FtauToNuE[iLogE][jLogE] = temp[4][iLogE][jLogE];
		// To NuMu
		/** For Glashow Resonence **/
		FnuEToNuMu[iLogE][jLogE] = temp[5][iLogE][jLogE];
		FnuMuToNuMu[iLogE][jLogE] = temp[6][iLogE][jLogE];
		FnuTauToNuMu[iLogE][jLogE] = temp[7][iLogE][jLogE];
		FmuToNuMu[iLogE][jLogE] = temp[8][iLogE][jLogE];
		FtauToNuMu[iLogE][jLogE] = temp[9][iLogE][jLogE];
		// To NuTau
		/** For Glashow Resonence **/
		FnuEToNuTau[iLogE][jLogE] = temp[10][iLogE][jLogE];
		FnuMuToNuTau[iLogE][jLogE] = temp[11][iLogE][jLogE];
		FnuTauToNuTau[iLogE][jLogE] = temp[12][iLogE][jLogE];
		FmuToNuTau[iLogE][jLogE] = temp[13][iLogE][jLogE];
		FtauToNuTau[iLogE][jLogE] = temp[14][iLogE][jLogE];
		// To E
		FnuEToE[iLogE][jLogE] = temp[15][iLogE][jLogE];
		FnuMuToE[iLogE][jLogE] = temp[16][iLogE][jLogE];
		FnuTauToE[iLogE][jLogE] = temp[17][iLogE][jLogE];
		FmuToE[iLogE][jLogE] = temp[18][iLogE][jLogE];
		FtauToE[iLogE][jLogE] = temp[19][iLogE][jLogE];
		// To Mu
		/** For Glashow Resonence **/
		FnuEToMu[iLogE][jLogE] = temp[20][iLogE][jLogE];
		FnuMuToMu[iLogE][jLogE] = temp[21][iLogE][jLogE];
		FnuTauToMu[iLogE][jLogE] = temp[22][iLogE][jLogE];
		FmuToMu[iLogE][jLogE] = temp[23][iLogE][jLogE];
		FtauToMu[iLogE][jLogE] = temp[24][iLogE][jLogE];
		// To Tau
		/** For Glashow Resonence **/
		FnuEToTau[iLogE][jLogE] = temp[25][iLogE][jLogE];
		FnuMuToTau[iLogE][jLogE] = temp[26][iLogE][jLogE];
		FnuTauToTau[iLogE][jLogE] = temp[27][iLogE][jLogE];
		FmuToTau[iLogE][jLogE] = temp[28][iLogE][jLogE];
		FtauToTau[iLogE][jLogE] = temp[29][iLogE][jLogE];
		// To Hadron
		FnuEToHadron[iLogE][jLogE] = temp[30][iLogE][jLogE];
		FnuMuToHadron[iLogE][jLogE] = temp[31][iLogE][jLogE];
		FnuTauToHadron[iLogE][jLogE] = temp[32][iLogE][jLogE];
		FmuToHadron[iLogE][jLogE] = temp[33][iLogE][jLogE];
		FtauToHadron[iLogE][jLogE] = temp[34][iLogE][jLogE];
	    }
	}

    }


    /** Copy the transfer matrix */
    public void copyTransferMatrix( ){
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		// To NuE
		nuEToNuE[iLogE][jLogE] = FnuEToNuE[iLogE][jLogE];
		nuMuToNuE[iLogE][jLogE] = FnuMuToNuE[iLogE][jLogE];
		nuTauToNuE[iLogE][jLogE] =FnuTauToNuE[iLogE][jLogE];
		muToNuE[iLogE][jLogE] = FmuToNuE[iLogE][jLogE];
		tauToNuE[iLogE][jLogE] = FtauToNuE[iLogE][jLogE]; 
		// To NuMu
		/** For Glashow Resonance **/
		nuEToNuMu[iLogE][jLogE] = FnuEToNuMu[iLogE][jLogE];
		nuMuToNuMu[iLogE][jLogE] = FnuMuToNuMu[iLogE][jLogE];
		nuTauToNuMu[iLogE][jLogE] = FnuTauToNuMu[iLogE][jLogE];
		muToNuMu[iLogE][jLogE] = FmuToNuMu[iLogE][jLogE];
		tauToNuMu[iLogE][jLogE] = FtauToNuMu[iLogE][jLogE];
		// To NuTau
		/** For Glashow Resonance **/
		nuEToNuTau[iLogE][jLogE] = FnuEToNuTau[iLogE][jLogE];
		nuMuToNuTau[iLogE][jLogE] = FnuMuToNuTau[iLogE][jLogE];
		nuTauToNuTau[iLogE][jLogE] = FnuTauToNuTau[iLogE][jLogE];
		muToNuTau[iLogE][jLogE] = FmuToNuTau[iLogE][jLogE]; 
		tauToNuTau[iLogE][jLogE] = FtauToNuTau[iLogE][jLogE];
		// To E
		nuEToE[iLogE][jLogE] = FnuEToE[iLogE][jLogE];
		nuMuToE[iLogE][jLogE] = FnuMuToE[iLogE][jLogE];
		nuTauToE[iLogE][jLogE] = FnuTauToE[iLogE][jLogE];
		muToE[iLogE][jLogE] = FmuToE[iLogE][jLogE];
		tauToE[iLogE][jLogE] = FtauToE[iLogE][jLogE];
		// To Mu
		/** For Glashow Resonance **/
		nuEToMu[iLogE][jLogE] = FnuEToMu[iLogE][jLogE];
		nuMuToMu[iLogE][jLogE] = FnuMuToMu[iLogE][jLogE];
		nuTauToMu[iLogE][jLogE] = FnuTauToMu[iLogE][jLogE];
		muToMu[iLogE][jLogE] = FmuToMu[iLogE][jLogE];
		tauToMu[iLogE][jLogE] = FtauToMu[iLogE][jLogE];
		// To Tau
		/** For Glashow Resonance **/
		nuEToTau[iLogE][jLogE] = FnuEToTau[iLogE][jLogE];
		nuMuToTau[iLogE][jLogE] = FnuMuToTau[iLogE][jLogE];
		nuTauToTau[iLogE][jLogE] = FnuTauToTau[iLogE][jLogE];
		muToTau[iLogE][jLogE] = FmuToTau[iLogE][jLogE];
		tauToTau[iLogE][jLogE] = FtauToTau[iLogE][jLogE];
		// To Hadron
		nuEToHadron[iLogE][jLogE] = FnuEToHadron[iLogE][jLogE];
		nuMuToHadron[iLogE][jLogE] = FnuMuToHadron[iLogE][jLogE];
		nuTauToHadron[iLogE][jLogE] = FnuTauToHadron[iLogE][jLogE];
		muToHadron[iLogE][jLogE] =  FmuToHadron[iLogE][jLogE];
		tauToHadron[iLogE][jLogE] = FtauToHadron[iLogE][jLogE];
	    }
	}

    }





    /****
	 Store the propagation matrix calculated so far
	 to the store matrix which save energy distribution
	 of neutrinos and leptons propagating to the current
	 location.
	 You can then initialize the propagation matrix
	 by init( ) and start another propagation calculation
	 in a different section of the trajectory.
	 You can get the results stored here back
	 to the propagation matrix
	 by calling the method copyTransferMatrixFromStore( ).
    */
    public void storePropagateMatrix( ){
	int iLogE,jLogE,kLogE;


	/** For Glahow Resonance -begin **/
	// NuE to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[0][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[0][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE] +
			SnuEToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE] +
			SnuEToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE] +
			SnuEToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE] +
			SnuEToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}
	/** For Glahow Resonance -end **/

	// NuMu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[1][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[1][iLogE][jLogE] += 
			SnuMuToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[2][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[2][iLogE][jLogE] += 
			SnuTauToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[3][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[3][iLogE][jLogE] += 
			SmuToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[4][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[4][iLogE][jLogE] += 
			StauToNuE[iLogE][kLogE]*FnuEToNuE[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToNuE[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToNuE[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToNuE[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToNuE[kLogE][jLogE];
		}
	    }
	}


	/** For Glashow Resonance -begin **/
	// NuE to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[5][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[5][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			SnuEToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			SnuEToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			SnuEToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			SnuEToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glashow Resonance -begin **/

	// NuMu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[6][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[6][iLogE][jLogE] += 
			/** For Glshow Resonance **/
			SnuMuToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[7][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[7][iLogE][jLogE] += 
			/** For Glshow Resonance **/
			SnuTauToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[8][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[8][iLogE][jLogE] += 
			/** For Glshow Resonance **/
			SmuToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[9][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[9][iLogE][jLogE] += 
			/** For Glshow Resonance **/
			StauToNuE[iLogE][kLogE]*FnuEToNuMu[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToNuMu[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToNuMu[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToNuMu[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToNuMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glshow Resonance -begin **/
	// NuE to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[10][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[10][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			SnuEToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			SnuEToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			SnuEToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			SnuEToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glshow Resonance -end **/

	// NuMu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[11][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[11][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			SnuMuToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[12][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[12][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			SnuTauToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[13][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[13][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			SmuToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[14][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[14][iLogE][jLogE] += 
			/** For Glashow Resonance **/
			StauToNuE[iLogE][kLogE]*FnuEToNuTau[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToNuTau[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToNuTau[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToNuTau[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToNuTau[kLogE][jLogE];
		}
	    }
	}




	/** For Glashow Resonance -begin **/
	// NuE to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[15][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[15][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE] +
			SnuEToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE] +
			SnuEToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE] +
			/** Modified ! **/
			SnuEToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE] +
			SnuEToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[15][iLogE][jLogE] += FnuEToE[iLogE][jLogE];
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[16][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[16][iLogE][jLogE] += 
			SnuMuToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[16][iLogE][jLogE] += FnuMuToE[iLogE][jLogE];
	    }
	}

	// NuTau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[17][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[17][iLogE][jLogE] += 
			SnuTauToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[17][iLogE][jLogE] += FnuTauToE[iLogE][jLogE];
	    }
	}

	// Mu to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[18][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[18][iLogE][jLogE] += 
			SmuToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[18][iLogE][jLogE] += FmuToE[iLogE][jLogE];
	    }
	}

	// Tau to E
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[19][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[19][iLogE][jLogE] += 
			StauToNuE[iLogE][kLogE]*FnuEToE[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToE[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToE[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToE[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToE[kLogE][jLogE];
		}
		temp[19][iLogE][jLogE] += FtauToE[iLogE][jLogE];
	    }
	}





	/** For Glahow Resonance -begin **/
	// NuE to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[20][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[20][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			SnuEToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			SnuEToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			SnuEToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			SnuEToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}
	/** For Glahow Resonance -end **/

	// NuMu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[21][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[21][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SnuMuToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[22][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[22][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SnuTauToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// Mu to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[23][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[23][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SmuToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}

	// Tau to Mu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[24][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[24][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			StauToNuE[iLogE][kLogE]*FnuEToMu[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToMu[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToMu[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToMu[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToMu[kLogE][jLogE];
		}
	    }
	}




	/** For Glahow Resonance -begin **/
	// NuE to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[25][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[25][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			SnuEToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			SnuEToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			SnuEToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			SnuEToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}
	/** For Glahow Resonance -end **/

	// NuMu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[26][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[26][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SnuMuToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// NuTau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[27][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[27][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SnuTauToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// Mu to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[28][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[28][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			SmuToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}

	// Tau to Tau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[29][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[29][iLogE][jLogE] += 
			/** For Glahow Resonance **/
			StauToNuE[iLogE][kLogE]*FnuEToTau[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToTau[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToTau[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToTau[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToTau[kLogE][jLogE];
		}
	    }
	}





	/** For Glashow Resonance -begin **/
	// NuE to Hadron
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[30][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[30][iLogE][jLogE] += 
			SnuEToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE] +
			SnuEToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE] +
			SnuEToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE] +
			SnuEToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE] +
			SnuEToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[30][iLogE][jLogE] += FnuEToHadron[iLogE][jLogE];
	    }
	}
	/** For Glashow Resonance -end **/

	// NuMu to Hadron 
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[31][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[31][iLogE][jLogE] += 
			SnuMuToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			SnuMuToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			SnuMuToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			SnuMuToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			SnuMuToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[31][iLogE][jLogE] += FnuMuToHadron[iLogE][jLogE];
	    }
	}

	// NuTau to Hadron 
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[32][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[32][iLogE][jLogE] += 
			SnuTauToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			SnuTauToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			SnuTauToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			SnuTauToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			SnuTauToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[32][iLogE][jLogE] += FnuTauToHadron[iLogE][jLogE];
	    }
	}

	// Mu to Hadron 
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[33][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[33][iLogE][jLogE] += 
			SmuToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			SmuToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			SmuToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			SmuToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			SmuToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[33][iLogE][jLogE] += FmuToHadron[iLogE][jLogE];
	    }
	}

	// Tau to Hadron 
	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		temp[34][iLogE][jLogE] = 0.0;
		for(kLogE=iLogE;kLogE>=jLogE;kLogE--){
		    temp[34][iLogE][jLogE] += 
			StauToNuE[iLogE][kLogE]*FnuEToHadron[kLogE][jLogE]+
			StauToNuMu[iLogE][kLogE]*FnuMuToHadron[kLogE][jLogE]+
			StauToNuTau[iLogE][kLogE]*FnuTauToHadron[kLogE][jLogE]+
			StauToMu[iLogE][kLogE]*FmuToHadron[kLogE][jLogE]+
			StauToTau[iLogE][kLogE]*FtauToHadron[kLogE][jLogE];
		}
		temp[34][iLogE][jLogE] += FtauToHadron[iLogE][jLogE];
	    }
	}


	for(iLogE=0;iLogE<dimension;iLogE++){
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		// To NuE
		SnuEToNuE[iLogE][jLogE] = temp[0][iLogE][jLogE];
		SnuMuToNuE[iLogE][jLogE] = temp[1][iLogE][jLogE];
		SnuTauToNuE[iLogE][jLogE] = temp[2][iLogE][jLogE];
		SmuToNuE[iLogE][jLogE] = temp[3][iLogE][jLogE];
		StauToNuE[iLogE][jLogE] = temp[4][iLogE][jLogE];
		// To NuMu
		/** For Glashow Resonance **/
		SnuEToNuMu[iLogE][jLogE] = temp[5][iLogE][jLogE];
		SnuMuToNuMu[iLogE][jLogE] = temp[6][iLogE][jLogE];
		SnuTauToNuMu[iLogE][jLogE] = temp[7][iLogE][jLogE];
		SmuToNuMu[iLogE][jLogE] = temp[8][iLogE][jLogE];
		StauToNuMu[iLogE][jLogE] = temp[9][iLogE][jLogE];
		// To NuTau
		/** For Glashow Resonance **/
		SnuEToNuTau[iLogE][jLogE] = temp[10][iLogE][jLogE];
		SnuMuToNuTau[iLogE][jLogE] = temp[11][iLogE][jLogE];
		SnuTauToNuTau[iLogE][jLogE] = temp[12][iLogE][jLogE];
		SmuToNuTau[iLogE][jLogE] = temp[13][iLogE][jLogE];
		StauToNuTau[iLogE][jLogE] = temp[14][iLogE][jLogE];
		// To E
		SnuEToE[iLogE][jLogE] = temp[15][iLogE][jLogE];
		SnuMuToE[iLogE][jLogE] = temp[16][iLogE][jLogE];
		SnuTauToE[iLogE][jLogE] = temp[17][iLogE][jLogE];
		SmuToE[iLogE][jLogE] = temp[18][iLogE][jLogE];
		StauToE[iLogE][jLogE] = temp[19][iLogE][jLogE];
		// To Mu
                /** For Glashow Resonance **/
		SnuEToMu[iLogE][jLogE] = temp[20][iLogE][jLogE];
		SnuMuToMu[iLogE][jLogE] = temp[21][iLogE][jLogE];
		SnuTauToMu[iLogE][jLogE] = temp[22][iLogE][jLogE];
		SmuToMu[iLogE][jLogE] = temp[23][iLogE][jLogE];
		StauToMu[iLogE][jLogE] = temp[24][iLogE][jLogE];
		// To Tau
                /** For Glashow Resonance **/
		SnuEToTau[iLogE][jLogE] = temp[25][iLogE][jLogE];
		SnuMuToTau[iLogE][jLogE] = temp[26][iLogE][jLogE];
		SnuTauToTau[iLogE][jLogE] = temp[27][iLogE][jLogE];
		SmuToTau[iLogE][jLogE] = temp[28][iLogE][jLogE];
		StauToTau[iLogE][jLogE] = temp[29][iLogE][jLogE];
		// To Hadron
		SnuEToHadron[iLogE][jLogE] = temp[30][iLogE][jLogE];
		SnuMuToHadron[iLogE][jLogE] = temp[31][iLogE][jLogE];
		SnuTauToHadron[iLogE][jLogE] = temp[32][iLogE][jLogE];
		SmuToHadron[iLogE][jLogE] = temp[33][iLogE][jLogE];
		StauToHadron[iLogE][jLogE] = temp[34][iLogE][jLogE];
	    }
	}

    }


    /** Copy the transfer matrix from Store matrix*/
    public void copyTransferMatrixFromStore( ){
	int iLogE;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    int jLogE;
	    for(jLogE=0;jLogE<=iLogE;jLogE++){
		// To NuE
		FnuEToNuE[iLogE][jLogE] = SnuEToNuE[iLogE][jLogE];
		FnuMuToNuE[iLogE][jLogE] = SnuMuToNuE[iLogE][jLogE];
		FnuTauToNuE[iLogE][jLogE] =SnuTauToNuE[iLogE][jLogE];
		FmuToNuE[iLogE][jLogE] = SmuToNuE[iLogE][jLogE];
		FtauToNuE[iLogE][jLogE] = StauToNuE[iLogE][jLogE]; 
		// To NuMu
                /** For Glashow Resonance **/
		FnuEToNuMu[iLogE][jLogE] = SnuEToNuMu[iLogE][jLogE];
		FnuMuToNuMu[iLogE][jLogE] = SnuMuToNuMu[iLogE][jLogE];
		FnuTauToNuMu[iLogE][jLogE] = SnuTauToNuMu[iLogE][jLogE];
		FmuToNuMu[iLogE][jLogE] = SmuToNuMu[iLogE][jLogE];
		FtauToNuMu[iLogE][jLogE] = StauToNuMu[iLogE][jLogE];
		// To NuTau
                /** For Glashow Resonance **/
		/** Modified ! **/
		FnuEToNuTau[iLogE][jLogE] = SnuEToNuTau[iLogE][jLogE];
		FnuMuToNuTau[iLogE][jLogE] = SnuMuToNuTau[iLogE][jLogE];
		FnuTauToNuTau[iLogE][jLogE] = SnuTauToNuTau[iLogE][jLogE];
		FmuToNuTau[iLogE][jLogE] = SmuToNuTau[iLogE][jLogE]; 
		FtauToNuTau[iLogE][jLogE] = StauToNuTau[iLogE][jLogE];
		// To E
		FnuEToE[iLogE][jLogE] = SnuEToE[iLogE][jLogE];
		FnuMuToE[iLogE][jLogE] = SnuMuToE[iLogE][jLogE];
		FnuTauToE[iLogE][jLogE] = SnuTauToE[iLogE][jLogE];
		FmuToE[iLogE][jLogE] = SmuToE[iLogE][jLogE];
		FtauToE[iLogE][jLogE] = StauToE[iLogE][jLogE];
		// To Mu
                /** For Glashow Resonance **/
		FnuEToMu[iLogE][jLogE] = SnuEToMu[iLogE][jLogE];
		FnuMuToMu[iLogE][jLogE] = SnuMuToMu[iLogE][jLogE];
		FnuTauToMu[iLogE][jLogE] = SnuTauToMu[iLogE][jLogE];
		FmuToMu[iLogE][jLogE] = SmuToMu[iLogE][jLogE];
		FtauToMu[iLogE][jLogE] = StauToMu[iLogE][jLogE];
		// To Tau
                /** For Glashow Resonance **/
		FnuEToTau[iLogE][jLogE] = SnuEToTau[iLogE][jLogE];
		FnuMuToTau[iLogE][jLogE] = SnuMuToTau[iLogE][jLogE];
		FnuTauToTau[iLogE][jLogE] = SnuTauToTau[iLogE][jLogE];
		FmuToTau[iLogE][jLogE] = SmuToTau[iLogE][jLogE];
		FtauToTau[iLogE][jLogE] = StauToTau[iLogE][jLogE];
		// To Hadron
		FnuEToHadron[iLogE][jLogE] = SnuEToHadron[iLogE][jLogE];
		FnuMuToHadron[iLogE][jLogE] = SnuMuToHadron[iLogE][jLogE];
		FnuTauToHadron[iLogE][jLogE] = SnuTauToHadron[iLogE][jLogE];
		FmuToHadron[iLogE][jLogE] =  SmuToHadron[iLogE][jLogE];
		FtauToHadron[iLogE][jLogE] = StauToHadron[iLogE][jLogE];
	    }
	}

    }


    /****** Change the infinitesimal propagation length */
    public void setDx(double dX){
	this.dX = dX;
    }



    /***** Get the propagation matrix element from nuE with energy iLogE to
     nuE with energy jLogE ***/
    public double getFnuEToNuE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToNuE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /** For Glashow Resonance **/
    /***** Get the propagation matrix element from nuE with energy iLogE to
     nuMu with energy jLogE ***/
    public double getFnuEToNuMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToNuMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /** For Glashow Resonance **/
    /***** Get the propagation matrix element from nuE with energy iLogE to
     nuTau with energy jLogE ***/
    public double getFnuEToNuTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToNuTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuE with energy iLogE to
     Electron with energy jLogE ***/
    public double getFnuEToE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /** For Glashow Resonance **/
    /***** Get the propagation matrix element from nuE with energy iLogE to
     Muon with energy jLogE ***/
    public double getFnuEToMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /** For Glashow Resonance **/
    /***** Get the propagation matrix element from nuE with energy iLogE to
     Tau with energy jLogE ***/
    public double getFnuEToTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuE with energy iLogE to
     Hadron with energy jLogE ***/
    public double getFnuEToHadron(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuEToHadron[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     nuE with energy jLogE ***/
    public double getFnuMuToNuE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToNuE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     nuMu with energy jLogE ***/
    public double getFnuMuToNuMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToNuMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     nuTau with energy jLogE ***/
    public double getFnuMuToNuTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToNuTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     Electron with energy jLogE ***/
    public double getFnuMuToE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     Muon with energy jLogE ***/
    public double getFnuMuToMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
    Tau with energy jLogE ***/
    public double getFnuMuToTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuMu with energy iLogE to
     Hadron with energy jLogE ***/
    public double getFnuMuToHadron(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuMuToHadron[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     nuE with energy jLogE ***/
    public double getFnuTauToNuE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToNuE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     nuMu with energy jLogE ***/
    public double getFnuTauToNuMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToNuMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     nuTau with energy jLogE ***/
    public double getFnuTauToNuTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToNuTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     Electron with energy jLogE ***/
    public double getFnuTauToE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     Muon with energy jLogE ***/
    public double getFnuTauToMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     Tau with energy jLogE ***/
    public double getFnuTauToTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from nuTau with energy iLogE to
     Hadron with energy jLogE ***/
    public double getFnuTauToHadron(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FnuTauToHadron[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    } 
    /***** Get the propagation matrix element from Muon with energy iLogE to
     nuE with energy jLogE ***/
    public double getFmuToNuE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToNuE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     nuMu with energy jLogE ***/
    public double getFmuToNuMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToNuMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     nuTau with energy jLogE ***/
    public double getFmuToNuTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToNuTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     Electron with energy jLogE ***/
    public double getFmuToE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     Muon with energy jLogE ***/
    public double getFmuToMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     Tau with energy jLogE ***/
    public double getFmuToTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Muon with energy iLogE to
     Hadron with energy jLogE ***/
    public double getFmuToHadron(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FmuToHadron[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     nuE with energy jLogE ***/
    public double getFtauToNuE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToNuE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     nuMu with energy jLogE ***/
    public double getFtauToNuMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToNuMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     nuTau with energy jLogE ***/
    public double getFtauToNuTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToNuTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     Electron with energy jLogE ***/
    public double getFtauToE(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToE[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     Muon with energy jLogE ***/
    public double getFtauToMu(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToMu[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     Tau with energy jLogE ***/
    public double getFtauToTau(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToTau[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }
    /***** Get the propagation matrix element from Tauon with energy iLogE to
     Hadron with energy jLogE ***/
    public double getFtauToHadron(int iLogE,int jLogE){
	if((0<=iLogE && iLogE<dimension)&&(0<=jLogE && jLogE<dimension)){
	    return FtauToHadron[iLogE][jLogE];
	}
	else{
	    return 0.0;
	}
    }


}

