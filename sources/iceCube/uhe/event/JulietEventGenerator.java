package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.points.*;
import iceCube.uhe.event.*;
import iceCube.uhe.geometry.*;
import numRecipes.*;
import geometry.*;

import java.util.*;
import java.io.*;

/**
   <pre>

   This is a run manager for the Event class. It initializes all the obejcts
   and parameters necessary for events and then run the Event objects.

   Initial version written by Shigeru Yoshida November 28 2004.
   Ice3 geometry calculation is implemented by S.Yoshida and T.Noda December 29 2004.
   Modified to adopt IceTray framework by K.Hoshina January 27 2005.
      added functions:
        - J3UnitVector getPrimaryParticleDirectionInEarthCenterCoordinate() 
        - J3UnitVector getPrimaryParticleDirectionInIceCubeCoordinate() 
        - J3Vector     getEndLocationInIceCubeCoordinate() 
        - J3Vector     wherePrimaryParticleEndsInIceCubeCoordinate()
      added field:
        - J3Vector endLocation_J3Vector_ice3
   Modified for neutrinos induced simulations by M.Ono Aplel 1 2007.
      modified a function:
        - public void runSingleEvent()
      added :
        - double neutrinoMinimumEnergyInTravel
   Modified for Glashow Resonance by M.Ono December 2 2007.
   Modified for function to shift Geometry by S.Yoshida December 15 2007.
   </pre>
*/

public class JulietEventGenerator {
    

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

    /** Geometry of propagating particle represended by EarthCenterCoordinate*/
    J3Line particleAxis_J3Line_center = null;
    /** Geometry of propagating particle represended by IceCubeCoordinate */
    J3Line particleAxis_J3Line_ice3 = null;
    /** Progation starting point location in IceCube coordinate */
    J3Vector startLocation_J3Vector_ice3 = null;
    /** Progation starting point location in Earth center coordinate */
    J3Vector startLocation_J3Vector_center = null;
    /** Progation end point location in IceCube coordinate */
    J3Vector endLocation_J3Vector_ice3 = null;
    /** Geometry Shift (if exists) of the particle trajectory 
	in IceCube cooridinate */
    J3Vector shift_J3Vector_ice3 = null;

    /** IceCube coordinate */
    IceCubeCoordinate ice3Coordinate = null;
    EarthCenterCoordinate earthCoordinate = null;

    /** IceCube Volume and IceCube outer volume */
    IceCubeVolume ice3Vol = null;
    Volume ice3OuterVol = null;

    /** Starting location of incident particle.
	For particlePoint/J3Line axislength. 
	Distance along the shower Axis from earth [cm] */
    double startLocation;

    /** Nadir angle [rad] when the particle enters the earth.*/
    double nadirAngleAtEntrance;

    /** start flag 1.. start propagation from the particle enterance point to the earth.
	flag 2.. start the particle run from 2km away from the ice3 coordinate origin.
    */
    int propagationFlag;

    /** Array of the MonteCarloBase objects */
    MonteCarloBase[] mcBases;

    /** Array of the InteractionsMatrix objects **/
    InteractionsMatrix[] intMtx;

    /** Palameter for making MonteCarloBase[] **/
    int     numOfDecay;
    /** For Glashow Resonance */
    int     numOfGR;
    boolean mudecay;
    boolean taudecay;
    int     tauDecayFlag;
    int     electronBaseFlag;

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

    /** Neutrino interaction matrix name */ // for IceCube-Gen2
    private static String eNuCCMtxObjectCCH5File = "ENeutrinoChargeMtx";
    private static String eNuCCMtxObjectZeusFile = "ENeutrinoChargeZeusNewMtx";
    private static String muNuCCMtxObjectCCH5File = "MuNeutrinoChargeMtx";
    private static String muNuCCMtxObjectZeusFile = "MuNeutrinoChargeZeusNewMtx";
    private static String tauNuCCMtxObjectCCH5File = "TauNeutrinoChargeMtx";
    private static String tauNuCCMtxObjectZeusFile = "TauNeutrinoChargeZeusNewMtx";

    private static String eNuNCMtxObjectCCH5File = "ENeutrinoNeutralMtx";
    private static String eNuNCMtxObjectZeusFile = "ENeutrinoNeutralZeusNewMtx";
    private static String muNuNCMtxObjectCCH5File = "MuNeutrinoNeutralMtx";
    private static String muNuNCMtxObjectZeusFile = "MuNeutrinoNeutralZeusNewMtx";
    private static String tauNuNCMtxObjectCCH5File = "TauNeutrinoNeutralMtx";
    private static String tauNuNCMtxObjectZeusFile = "TauNeutrinoNeutralZeusNewMtx";
    public static boolean neutrinoCSHERAZeus = true; // neutrino HERA-Zeus-PDF based cross section

    /** Random Generator */
    RandomGenerator rand;
    long seed;
    long[] random_state;

    /** List of the cascade particles, energy deposite and 
	the interaction points along the track */
    List particleList = null;
    ListIterator particleIterator = null;
    List locationIce3List = null;
    ListIterator locationIce3Iterator = null;
    List particleInteractionsList = null;

    /** List of the track particles and 
        the interaction points */
    List         trackParticleList     = null;
    ListIterator trackParticleIterator = null;
    List         trackLocationIce3List     = null;
    ListIterator trackLocationIce3Iterator = null;

    /** Dimension of LogEnergyMatrix */
    int dim = Particle.getDimensionOfLogEnergyMatrix();

    /**
      Constructor. Here all the obejcts including individual InteractionsBase objects
      (which implies it reads out the InteractionsMatrix objects dumped in the directory
      classes/iceCube/uhe/interactions/) and  Mu/TauDecayBase objects are generated. Each
      object is generated in an interactive way to an user.
    */
    public JulietEventGenerator(int flavorID, int doubletID, double energy, int mediumID, 
                                int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                int doMuPN, int doTauPN, int doGR, int doMuDecay,
                                int doTauDecay, int posID) throws IOException{

        long[] random_state = null;
        long seed = -1;

        configureJULIeT(flavorID, doubletID, energy, mediumID, 
                        doCC, doNC, doMuBrems, doTauBrems, 
                        doMuKnock, doTauKnock, doMu2e, doTau2e,
                        doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                        doMuPN, doTauPN, doGR, doMuDecay,
                        doTauDecay, posID, seed, random_state);

    }


    /** The same constructor, but without the glashow resoanance. 
	This constructor exists to maintain the backward compatibility
    */
    public JulietEventGenerator(int flavorID, int doubletID, double energy, int mediumID, 
                                int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                int doMuPN, int doTauPN, int doMuDecay,
                                int doTauDecay, int posID) throws IOException{

	    int doGR = 0;
        long[] random_state = null;
        long seed = -1;

        configureJULIeT(flavorID, doubletID, energy, mediumID, 
                        doCC, doNC, doMuBrems, doTauBrems, 
                        doMuKnock, doTauKnock, doMu2e, doTau2e,
                        doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                        doMuPN, doTauPN, doGR, doMuDecay,
                        doTauDecay, posID, seed, random_state);

    }

    /**
      Constructor using the seed for the random generator.
    */
    public JulietEventGenerator(int flavorID, int doubletID, double energy, int mediumID, 
                                int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                int doMuPN, int doTauPN, int doGR, int doMuDecay,
                                int doTauDecay, int posID, long seed) throws IOException{

        long[] random_state = null;

        configureJULIeT(flavorID, doubletID, energy, mediumID, 
                        doCC, doNC, doMuBrems, doTauBrems, 
                        doMuKnock, doTauKnock, doMu2e, doTau2e,
                        doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                        doMuPN, doTauPN, doGR, doMuDecay,
                        doTauDecay, posID, seed, random_state);

    }

    /** The same constructor, but without the glashow resoanance. 
    This constructor exists to maintain the backward compatibility
    */
    public JulietEventGenerator(int flavorID, int doubletID, double energy, int mediumID, 
                                int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                int doMuPN, int doTauPN, int doMuDecay,
                                int doTauDecay, int posID, long seed) throws IOException{

        int doGR = 0;
        long[] random_state = null;

        configureJULIeT(flavorID, doubletID, energy, mediumID, 
                        doCC, doNC, doMuBrems, doTauBrems, 
                        doMuKnock, doTauKnock, doMu2e, doTau2e,
                        doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                        doMuPN, doTauPN, doGR, doMuDecay,
                        doTauDecay, posID, seed, random_state);

    }

    /**
      Constructor using the random_state for the generator.
    */
    public JulietEventGenerator(int flavorID, int doubletID, double energy, int mediumID, 
                                int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                int doMuPN, int doTauPN, int doGR, int doMuDecay,
                                int doTauDecay, int posID, long[] random_state) throws IOException{

        long seed = -1;

        configureJULIeT(flavorID, doubletID, energy, mediumID, 
                        doCC, doNC, doMuBrems, doTauBrems, 
                        doMuKnock, doTauKnock, doMu2e, doTau2e,
                        doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                        doMuPN, doTauPN, doGR, doMuDecay,
                        doTauDecay, posID, seed, random_state);

    }



    public JulietEventGenerator() throws IOException{

        int flavorID, doubletID, mediumID; 
        int doCC, doNC, doMuBrems, doTauBrems; 
        int doMuKnock, doTauKnock, doMu2e, doTau2e;
        int doMu2mu, doTau2mu, doMu2tau, doTau2tau;
        int doMuPN, doTauPN, doMuDecay, doTauDecay, posID;
        long seed;
        long[] random_state = null;

    	/** For Glashow Resonance */
    	int doGR;

        double energy;
	
	DataInputStream input = new DataInputStream(System.in);
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input));
	String buffer;

	System.out.print("Particle Flavor [e(0)/mu(1)/tau(2)] ->");
	buffer   = d.readLine();
	flavorID = Integer.valueOf(buffer).intValue();

	System.out.print("Particle Doublet [neutrino(0)/charged(1)] ->");
	buffer    = d.readLine();
	doubletID = Integer.valueOf(buffer).intValue();

	System.out.print("Particle Energy [GeV] ->");
	buffer = d.readLine();
	energy = Double.valueOf(buffer).doubleValue();

	System.out.print("Medium in Propagation ice(0)/rock(1)?->");
	buffer   = d.readLine();
	mediumID = Integer.valueOf(buffer).intValue();

	System.out.print("Charged Current Interactions yes(1)/no(0)/allFlavor(2)?->");
	buffer = d.readLine();
	doCC   = Integer.valueOf(buffer).intValue();

	System.out.print("Neutral Current Interactions yes(1)/no(0)/allFlavor(2)?->");
	buffer = d.readLine();
	doNC   = Integer.valueOf(buffer).intValue();

	System.out.print("Muon Bremsstrahlung yes(1)/no(0)?->");
	buffer    = d.readLine();
	doMuBrems = Integer.valueOf(buffer).intValue();

	System.out.print("Tau Bremsstrahlung yes(1)/no(0)?->");
	buffer     = d.readLine();
	doTauBrems = Integer.valueOf(buffer).intValue();

	System.out.print("Muon Knock-on Electrons yes(1)/no(0)?->");
	buffer    = d.readLine();
	doMuKnock = Integer.valueOf(buffer).intValue();

	System.out.print("Tau Knock-on Electrons yes(1)/no(0)?->");
	buffer     = d.readLine();
	doTauKnock = Integer.valueOf(buffer).intValue();

	System.out.print("Muon e+e- Pair Creation yes(1)/no(0)?->");
	buffer = d.readLine();
	doMu2e = Integer.valueOf(buffer).intValue();

	System.out.print("Tau e+e- Pair Creation yes(1)/no(0)?->");
	buffer  = d.readLine();
	doTau2e = Integer.valueOf(buffer).intValue();

	System.out.print("Muon mu+mu- Pair Creation yes(1)/no(0)?->");
	buffer  = d.readLine();
	doMu2mu = Integer.valueOf(buffer).intValue();

	System.out.print("Tau mu+mu- Pair Creation yes(1)/no(0)?->");
	buffer   = d.readLine();
	doTau2mu = Integer.valueOf(buffer).intValue();

	System.out.print("Muon tau+tau- Pair Creation yes(1)/no(0)?->");
	buffer   = d.readLine();
	doMu2tau = Integer.valueOf(buffer).intValue();

	System.out.print("Tau tau+tau- Pair Creation yes(1)/no(0)?->");
	buffer    = d.readLine();
	doTau2tau = Integer.valueOf(buffer).intValue();

	System.out.print("Muon Photo-nuclear interactions yes(1)/no(0)?->");
	buffer = d.readLine();
	doMuPN = Integer.valueOf(buffer).intValue();

	System.out.print("Tau Photo-nuclear interactions yes(1)/no(0)?->");
	buffer  = d.readLine();
	doTauPN = Integer.valueOf(buffer).intValue();

	/** For Glashow Resonance -begin*/
	System.out.print("Glashow Resonance yes(1)/no(0)?->");
	buffer     = d.readLine();
	doGR = Integer.valueOf(buffer).intValue();
	/** For Glashow Resonance -end */

	System.out.print("Mu decay yes(1)/no(0)?->");
	buffer    = d.readLine();
	doMuDecay = Integer.valueOf(buffer).intValue();

	System.out.print("Tau Decay yes(1)/no(0)?->");
	buffer     = d.readLine();
	doTauDecay = Integer.valueOf(buffer).intValue();

	System.out.print("Propagation starts at the earth entrance point yes(1)/from 880m away(2)?->");
	buffer = d.readLine();
	posID = Integer.valueOf(buffer).intValue();

    System.out.print("Seed for the RandomGenerator");
    buffer = d.readLine();
    seed = Long.valueOf(buffer).longValue();

	/** For Glashow Resonance */
    configureJULIeT(flavorID, doubletID, energy, mediumID, 
                    doCC, doNC, doMuBrems, doTauBrems, 
                    doMuKnock, doTauKnock, doMu2e, doTau2e,
                    doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                    doMuPN, doTauPN, doGR, doMuDecay,
                    doTauDecay, posID, seed, random_state);
    }

    /**
      configureJULIeT
      Initialize primary particle, generate interactionsMatrix and generate 
      IceCube geometry with coordinate systems. Primary particle information 
      can be updated after the configuration, while others are fixed here.
      Called by constructors.
    */
    private void configureJULIeT(int flavorID, int doubletID, double energy, int mediumID, 
                                 int doCC, int doNC, int doMuBrems, int doTauBrems, 
                                 int doMuKnock, int doTauKnock, int doMu2e, int doTau2e,
                                 int doMu2mu, int doTau2mu, int doMu2tau, int doTau2tau,
                                 int doMuPN, int doTauPN, int doGR, int doMuDecay,
                                 int doTauDecay, int posID, long seed,
                                 long[] random_state) throws IOException{


        primaryFlavor  = flavorID;
        primaryDoublet = doubletID;
        primaryEnergy  = energy;
        materialNumber = mediumID;

        // Generate Random Generator
        if(random_state != null){
          System.out.println("Using random_state to setup Random Generator");
          rand = new RandomGenerator(random_state);
        }
        else if(seed != -1){
          System.out.println("Using seed to setup Random Generator");
          rand = new RandomGenerator(seed);
        }
        else{
          System.out.println("Using system time to setup Random Generator");
          rand = new RandomGenerator();
        }

        // Register Interactions and read the InteractionMatrix objects
	    /** For Glashow Resonance 16->20*/
        String[] fileName = new String[20];
        String matrixName = null;
        if(materialNumber==0)  interactionsMatrixDirectory = interactionsMatrixDirectoryInIce;
        else interactionsMatrixDirectory = interactionsMatrixDirectoryInRock;
        int n = 0; 
        electronBaseFlag = 0;

        // Charged Current Interactions [yes(1)/no(0)/allFlavor(2)]
        if(doCC == 1){
	    if(!neutrinoCSHERAZeus){
		if(primaryFlavor==0) {matrixName = eNuCCMtxObjectCCH5File; electronBaseFlag = 1;} 
		else if(primaryFlavor==1) matrixName = muNuCCMtxObjectCCH5File;
		else if(primaryFlavor==2) matrixName = tauNuCCMtxObjectCCH5File;
	    }else{
		if(primaryFlavor==0) {matrixName = eNuCCMtxObjectZeusFile; electronBaseFlag = 1;} 
		else if(primaryFlavor==1) matrixName = muNuCCMtxObjectZeusFile;
		else if(primaryFlavor==2) matrixName = tauNuCCMtxObjectZeusFile;
	    }
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }else if(doCC == 2){
	    if(!neutrinoCSHERAZeus){
		matrixName = eNuCCMtxObjectCCH5File; electronBaseFlag = 1;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = muNuCCMtxObjectCCH5File;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = tauNuCCMtxObjectCCH5File;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    }else{
		matrixName = eNuCCMtxObjectZeusFile; electronBaseFlag = 1;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = muNuCCMtxObjectZeusFile;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = tauNuCCMtxObjectZeusFile;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    }
        }

        // Neutral Current Interactions [yes(1)/no(0)/allFlavor(2)]
        if(doNC == 1){
	    if(!neutrinoCSHERAZeus){
		if(primaryFlavor==0) matrixName      = eNuNCMtxObjectCCH5File;
		else if(primaryFlavor==1) matrixName = muNuNCMtxObjectCCH5File;
		else if(primaryFlavor==2) matrixName = tauNuNCMtxObjectCCH5File;
	    }else{
		if(primaryFlavor==0) matrixName      = eNuNCMtxObjectZeusFile;
		else if(primaryFlavor==1) matrixName = muNuNCMtxObjectZeusFile;
		else if(primaryFlavor==2) matrixName = tauNuNCMtxObjectZeusFile;
	    }
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }else if(doNC == 2){
	    if(!neutrinoCSHERAZeus){
		matrixName = eNuNCMtxObjectCCH5File;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = muNuNCMtxObjectCCH5File;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = tauNuNCMtxObjectCCH5File;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    }else{
		matrixName = eNuNCMtxObjectZeusFile;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = muNuNCMtxObjectZeusFile;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
		matrixName = tauNuNCMtxObjectZeusFile;
		fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
	    }
        }

        // Muon Bremsstrahlung [yes(1)/no(0)]
        if(doMuBrems == 1){
            matrixName   = "muBremsstrahlungMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Bremsstrahlung [yes(1)/no(0)]
        if(doTauBrems == 1){
            matrixName   = "tauBremsstrahlungMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon Knock-on Electrons [yes(1)/no(0)]
        if(doMuKnock == 1){
            matrixName   = "muKnockOnElectronsMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Knock-on Electrons [yes(1)/no(0)]
        if(doTauKnock == 1){
            matrixName   = "tauKnockOnElectronsMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon e+e- Pair Creation [yes(1)/no(0)]
        if(doMu2e == 1){
            matrixName   = "muToePairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau e+e- Pair Creation [yes(1)/no(0)]
        if(doTau2e == 1){
            matrixName   = "tauToePairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon mu+mu- Pair Creation [yes(1)/no(0)]
        if(doMu2mu == 1){
            matrixName   = "muTomuPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau mu+mu- Pair Creation [yes(1)/no(0)]
        if(doTau2mu == 1){
            matrixName   = "tauTomuPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon tau+tau- Pair Creation [yes(1)/no(0)]
        if(doMu2tau == 1){
            matrixName   = "muTotauPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau tau+tau- Pair Creation [yes(1)/no(0)]
        if(doTau2tau == 1){
            matrixName   = "tauTotauPairCreationFitMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Muon Photo-nuclear interactions [yes(1)/no(0)]
        if(doMuPN == 1){
            matrixName = "muPhotoNuclearMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

        // Tau Photo-nuclear interactions [yes(1)/no(0)]
        if(doTauPN == 1){
            matrixName   = "tauPhotoNuclearMtx";
            fileName[n++]= interactionsMatrixDirectory.concat(matrixName);
        }

	/** For Glashow Resonance */
        // Glashow Resonance Electronic [yes(1)/no(0)]
	if(doGR == 1) electronBaseFlag = 1;

        numOfDecay = 0;
        mudecay    = false;
        taudecay   = false;

        // Mu decay [yes(1)/no(0)]
        if(doMuDecay == 1){
            Particle       mu         = new Particle(1,1,primaryEnergy);
            MuDecayYMatrix muDecayMtx = new MuDecayYMatrix(mu); 
            muDecayBase = new MuDecayBase(muDecayMtx);
            n++; numOfDecay++; mudecay=true;
        }

        // Tau Decay [yes(1)/no(0)]
        tauDecayFlag = 0;
        if(doTauDecay == 1){
            Particle        tau         = new Particle(2, 1, primaryEnergy);
            TauDecayYMatrix tauDecayMtx = new TauDecayYMatrix(tau);
            tauDecayBase = new TauDecayBase(tauDecayMtx);
            n++; numOfDecay++; taudecay=true; tauDecayFlag=1;
        }

        // create an array for MonteCarloBase
        //mcBases = new MonteCarloBase[n - numOfDecay + electronBaseFlag + 2];
        mcBases = new MonteCarloBase[n  + electronBaseFlag + Particle.NumberOfFlavor*doGR];

        // create an array for interactionsMatrix
        intMtx = new InteractionsMatrix[n - numOfDecay];
        //intMtx = new InteractionsMatrix[n];

        // regist interactionsMtx
        System.out.println("Registered Interaction Matrix....");

	for(int i=0; i<(n-numOfDecay) ;i++){       
            System.out.println(fileName[i]);
            //FileInputStream in = new FileInputStream(fileName[i]);
            InputStream in = ClassLoader.getSystemResourceAsStream(fileName[i]);
            intMtx[i]  = InteractionsMatrixInput.inputInteractionsMatrix(in);
            mcBases[i] = new InteractionsBase(intMtx[i]); 
            in.close();
        }

	// register interactionMtx of Glashow Resonance
	for(int i=0; i<Particle.NumberOfFlavor*doGR ;i++){       
            mcBases[n-numOfDecay+i] = new GlashowResonanceBase(i, materialNumber); 
	    }

        // register ElectronBase when primary particle is nu-e and charged current is choozen
        if(electronBaseFlag == 1) {
            Particle  e = new Particle(0,1,primaryEnergy);
            mcBases[n-numOfDecay+Particle.NumberOfFlavor*doGR] = new ElectronBase(e);
            System.out.println("ElectronBase is registered");
       }

        // register decay matrix
        int m = n-numOfDecay+Particle.NumberOfFlavor*doGR+electronBaseFlag;
        if(mudecay) {
            mcBases[m++] = muDecayBase;
            System.out.println("MuDecayBase is registered");
        }if(taudecay) {
            mcBases[m++] = tauDecayBase;
            System.out.println("TauDecayBase is registered");
        }

        // Generate the coordinate system 
        ice3Coordinate = new IceCubeCoordinate();
        earthCoordinate = new EarthCenterCoordinate();

        // Generate the icecube volume
        ice3Vol = new IceCubeVolume();
        ice3OuterVol = new Volume(2.0*ice3Coordinate.elevation); // 2x884m x 2x884m x 2x884m
        System.out.println("IceCube Coordinates and ice3volume generated");

        // Propagation starts at the earth entrance point [1:from EarthEntrancePoint/2:from 884m away]
        if(posID == 1) propagationFlag=1;
        else if(posID == 2) propagationFlag=2;
	else propagationFlag=3;

    }

    /** Set the start location of primary particle along its trajectory axis
     (particle point/J3Line)*/
    public void setStartLocationAlongTheAxis(double start) {
	if(start>0.0) startLocation = start;
    }

    /** Get the start location of primary particle. Return the distance [cm]
     along the particle axis from the Earth(+glacier) surface to
    the start location. */
    public double getStartLocationAlongTheAxis() {
	return startLocation;
    }

    /** Get the point location where the particle run starts.
	Represented by Earth center coordinate */
    public J3Vector wherePrimaryParticleStartsInEarthCenterCoordinate() {
	if(startLocation_J3Vector_center!=null){
	    return startLocation_J3Vector_center;
	}else {
	    return ice3Coordinate.getOrigin();
	}
    }


    /** Get the point location where the particle run starts.
	Represented by the IceCube coordinate */
    public J3Vector wherePrimaryParticleStartsInIceCubeCoordinate() {
	if(startLocation_J3Vector_ice3!=null){
	    return startLocation_J3Vector_ice3;
	}else {
	    return ice3Coordinate.getOrigin();
	}
    }

    /** Get the point location where the particle run starts.
	Represented by the IceCube coordinate */
    public J3Vector wherePrimaryParticleEndsInIceCubeCoordinate() {
	if(endLocation_J3Vector_ice3!=null){
	    return endLocation_J3Vector_ice3;
	}else {
	    System.err.println("Warning: endLocation_J3Vector_ice3 is null! return 0 vector");
	    return ice3Coordinate.getOrigin();
	}
    }

    /** Get of propagating particle direction represented by Earth center coordinate */
    public J3UnitVector getPrimaryParticleDirectionInEarthCenterCoordinate() {
	if(particleAxis_J3Line_center != null && particleAxis_J3Line_center.getDirection() != null){
	    return particleAxis_J3Line_center.getDirection();
	}else {
	    System.err.println("Warning: direction_J3UnitLine_center.getDirection() failed! return 0 vector");
	    return new J3UnitVector(0., 0., 0.);
	}
    }

    /** Get of propagating particle direction represented by IceCube coordinate */
    public J3UnitVector getPrimaryParticleDirectionInIceCubeCoordinate() {
	if(particleAxis_J3Line_ice3 != null && particleAxis_J3Line_ice3.getDirection() != null){
	    return particleAxis_J3Line_ice3.getDirection();
	}else {
	    System.err.println("Warning: direction_J3UnitLine_ice3.getDirection() failed! return 0 vector");
	    return new J3UnitVector(0., 0., 0.);
	}
    }

    /** Set the geometrical shift in IceCube coordinate */
    public void setGeometryShift(J3Vector shift){
	shift_J3Vector_ice3 = shift;
    }

    /** Get the geometrical shift in IceCube coordinate */
    public J3Vector getGeometryShift(){
	if(shift_J3Vector_ice3 != null){
	    return shift_J3Vector_ice3;
	}else{
	    return new J3Vector(0.0,0.0,0.0);
	}
    }


    /** Method to set primary energy, flavor, doublet for defining
	the input propagating particle, and its geomerty(for future use).
	If you do not call this method, the generator use the default
	values set by the constructor JulietEventGenerator().
	<pre>

	int flavor       :  flavor of the particle (c.f. Particle.java)
	int doublet      :  doublet of the particle (c.f. Particle.java)
	double energy    :  primary energy of the particle [GeV]
	J3Line axis      :  geometry of the particle trajectory. Must be represented
                            by the IceCubeCoordinate.
    */
    public void definePropagatingParticle(int flavor, int doublet,
					  double energy, J3Line axis){
	if(primaryFlavor!=flavor){
	    System.err.println("Warning: defferent flavor is set");
	    System.err.println("You have to load the corresponging interaction matrix");
	}
	primaryFlavor = flavor;
	primaryDoublet = doublet;
	primaryEnergy = energy;

	particleAxis_J3Line_ice3 = axis;
    }
    public void definePropagatingParticle(int flavor, int doublet,
					  double energy){
	if(primaryFlavor!=flavor){
	    System.err.println("Warning: defferent flavor is set");
	    System.err.println("You have to load the corresponging interaction matrix");
	}
	primaryFlavor = flavor;
	primaryDoublet = doublet;
	primaryEnergy = energy;
    }

    /** define the geometry of propagating particles 
	<pre>
	double x_ice3  : x [cm] of a particle trajectory point in the IceCube coordinate
	double y_ice3  : y [cm] of a particle trajectory point in the IceCube coordinate
	double z_ice3  : z [cm] of a particle trajectory point in the IceCube coordinate
	double nadirAngle_ice3   : nadir angle [deg] define in the IceCube coodrinate
	double azimuthAngle_ice3 : nadir angle [deg] define in the IceCube coodrinate
       </pre>
    */
    public void definePropagationGeometry(double x_ice3, double y_ice3, double z_ice3,
					  double nadirAngle_ice3, double azimuthAngle_ice3){
	J3Vector passingPoint_J3Vector_ice3 = new J3Vector(x_ice3,y_ice3,z_ice3);
	J3Vector n_ice3 = ice3Coordinate.getPointVectorFromPolarCoordinate
	    (1.0, Math.toRadians(nadirAngle_ice3), Math.toRadians(azimuthAngle_ice3));
	J3UnitVector direction_J3UnitVector_ice3 = 
	    new J3UnitVector(n_ice3.getX(),n_ice3.getY(),n_ice3.getZ());
	particleAxis_J3Line_ice3 = new J3Line(passingPoint_J3Vector_ice3,
					      direction_J3UnitVector_ice3);
    }


    /** define the geometry of propagating particles 
	<pre>
	double x_ice3  : x [cm] of a particle trajectory point in the IceCube coordinate
	double y_ice3  : y [cm] of a particle trajectory point in the IceCube coordinate
	double z_ice3  : z [cm] of a particle trajectory point in the IceCube coordinate
	double nx,ny,nz : unit direction vector represented by the IceCube coordinate
       </pre>
    */
    public void definePropagationGeometry(double x_ice3, double y_ice3, double z_ice3,
					  double nx, double ny, double nz){
	J3Vector passingPoint_J3Vector_ice3 = new J3Vector(x_ice3,y_ice3,z_ice3);
	J3UnitVector direction_J3UnitVector_ice3 = 
	    new J3UnitVector(nx,ny,nz);
	particleAxis_J3Line_ice3 = new J3Line(passingPoint_J3Vector_ice3,
					      direction_J3UnitVector_ice3);
    }


    /** define the geometry of propagating particles 
	<pre>
	J3Line particleAxisIce3 : The primary particle trajectory. Should be
                                  defined by the IceCubeCoordinate.
       </pre>
    */
    public void definePropagationGeometry(J3Line particleAxisIce3){
	particleAxis_J3Line_ice3 = particleAxisIce3;
    }



    /** Configure the particle progation geometry. 
	<pre>
	1. Setup the trajectory axis (J3Line) in the Earth center coordinate
	2. Calculate the point where the primary particle enters the earth
	3. Calculate the start location where the MC propgation run begines.
           It depnds on the propagationFlag setup by the constructor.
	   when 1, the location is equal to the earth entrance point.
	   when 2, it is IceCubeCoordinate.elevation away from the ice3 coordinate origin.
	   when 3, it is equal to the original R0 of J3Line.
	   You can change this location by calling the method 
	   setStartLocationAlongTheAxis() if necessary
       </pre>
    */
    public void configurePropagationGeometry(){

	// 1. locate the start point in Earth center coordinate
	J3Vector passingPoint_J3Vector_center = 
	    ice3Coordinate.transformVectorToEarthCenter(particleAxis_J3Line_ice3.getR0());
	J3UnitVector direction_J3UnitVector_center =
            ice3Coordinate.transformUnitVectorToEarthCenter(particleAxis_J3Line_ice3.getDirection());
	particleAxis_J3Line_center = new J3Line(passingPoint_J3Vector_center,
						direction_J3UnitVector_center);

        // 2. calculate initial entrance point of primary particle
	ParticleTracker.setInitialPoint(particleAxis_J3Line_center,ice3Coordinate);
	J3Vector particleEntrance = new J3Vector(
						 particleAxis_J3Line_center.getX(),
						 particleAxis_J3Line_center.getY(),
						 particleAxis_J3Line_center.getZ());
	particleAxis_J3Line_center.setR0(particleEntrance);
	particleAxis_J3Line_center.setAxisLength(0.0); // defines axis length =0
	                                               // when enters the earth

	nadirAngleAtEntrance = 
	    Math.PI - J3Vector.getAngleInRadian(particleAxis_J3Line_center,
						direction_J3UnitVector_center);

        // 3. calculate the start location
	if(propagationFlag==2){ // starting location is 880 m away 
                                // from the ice3 coordinate origin
	    if(shift_J3Vector_ice3 == null){
		J3Utility.setJ3LineNegativeAxisLengthForGivenLength(
	        particleAxis_J3Line_center,ice3Coordinate.getOrigin(),ice3Coordinate.elevation);
	    }else{
		J3Vector shift_J3Vector_center // absolute shifted point in earth center coordinate
		    = ice3Coordinate.transformVectorToEarthCenter(shift_J3Vector_ice3);
		J3Utility.setJ3LineNegativeAxisLengthForGivenLength(
                particleAxis_J3Line_center,shift_J3Vector_center,ice3Coordinate.elevation);
	    }
	    if(particleAxis_J3Line_center.getLength()>
	       (point.REarth+ice3Coordinate.getGlacierDepth())){// outside earth+glacier!!
		particleAxis_J3Line_center.setAxisLength(0.0); // back to the earth entrance
	    }

	}else if(propagationFlag==3){ // starting from the original vertex location
	    J3Vector startVertex = J3Vector.subtract(passingPoint_J3Vector_center,particleEntrance);
	    startLocation = startVertex.getLength();
	    particleAxis_J3Line_center.setAxisLength(startLocation);
	}
	    


	startLocation = particleAxis_J3Line_center.getAxisLength();

	startLocation_J3Vector_center = new J3Vector(
						     particleAxis_J3Line_center.getX(),
						     particleAxis_J3Line_center.getY(),
						     particleAxis_J3Line_center.getZ());

	startLocation_J3Vector_ice3 = 
	    ice3Coordinate.transformVectorToThisCoordinate(
		   particleAxis_J3Line_center,earthCoordinate);

	J3Vector particleEntrance_ice3 =
	    ice3Coordinate.transformVectorToThisCoordinate(
		   particleEntrance,earthCoordinate);
	particleAxis_J3Line_ice3.setR0(particleEntrance_ice3);
	particleAxis_J3Line_ice3.setAxisLength(startLocation);

    }


    /** Method to run a single event with a monocromatic primary energy given by
	the constructor.**/

    public void runSingleEvent( ){

	// generate the data list
	particleList          =    new LinkedList();
    particleInteractionsList = new LinkedList();
	locationIce3List      =    new LinkedList();
	trackParticleList     =    new LinkedList();
	trackLocationIce3List =    new LinkedList();

	// generate a propagating particle
	propParticle = new Particle(primaryFlavor,primaryDoublet,primaryEnergy);
	//System.out.println(
	//	   Particle.particleName(propParticle.getFlavor(),propParticle.getDoublet()) +
	//	   " has been generated with energy of " + propParticle.getEnergy() +
	//	   " [GeV]");    

	// Put the primary particle to the start location
	particleAxis_J3Line_ice3.setAxisLength(startLocation);
	particleAxis_J3Line_center.setAxisLength(startLocation);

    // add propParticle to trackParticleList and trackLocationIce3List
    // *** It's an initial primary track ***
    trackParticleList.add(new Particle(primaryFlavor, primaryDoublet, primaryEnergy)); 
    trackLocationIce3List.add(startLocation_J3Vector_ice3);

	// generate the particle point
	point = new ParticlePoint(0.0, nadirAngleAtEntrance, materialNumber);
	point.setIceRockBoundaryRadius(1.01*point.REarth);
	point.setParticleLocation(startLocation);

	//System.out.println("Propagation starts at x = " +
	//	   startLocation_J3Vector_center.getX() + " y = " +
	//	   startLocation_J3Vector_center.getY() + " z = " +
	//	   startLocation_J3Vector_center.getZ());
	//System.out.println("Distance from Earth center " + 
	//	   startLocation_J3Vector_center.getLength() + " [cm]");
	//System.out.println("Distance from IceCube coordinate origin " + 
	//	   startLocation_J3Vector_ice3.getLength() + " [cm]");

	event        = new Event(mcBases, propParticle, point);
	int haveBeenInsideIce3OuterVolume = 0;
	// outer volume flag. 

	// flag to see if the collision is the week interaction/decay process.
	boolean wasWeakInt = false;
    boolean wasNC = false;
	boolean wasDecay = false;
	/** For Glashow Resonance */
	boolean wasGR = false;
	boolean wasGRLepton = false;

	// 100GeV is the minimum energy of neutrino cross-section 
	// which you need in filling the interaction weight
	double neutrinoMinimumEnergyInTravel = 100.0;

	while(true){
	    double axisLengthBefore = particleAxis_J3Line_center.getAxisLength();

	    // running the particle to the interaction point
	    double pathLength = event.getPhysicalPathLength(rand);
	    moveParticleAxis(pathLength);
	    J3Vector particleLocation_J3Vector_ice3 =
		ice3Coordinate.transformVectorToThisCoordinate(
			       particleAxis_J3Line_center,earthCoordinate);
	    double axisLengthNow = particleAxis_J3Line_center.getAxisLength();


	    if(haveBeenInsideIce3OuterVolume==0){ // Not through the outer volume yet
		if((shift_J3Vector_ice3==null &&
		    ice3OuterVol.isJ3LineInsideVolume(particleAxis_J3Line_ice3,
						      axisLengthBefore,
						      axisLengthNow))||
		   (shift_J3Vector_ice3!=null &&
		    ice3OuterVol.isJ3LineInsideVolume(particleAxis_J3Line_ice3,
						      shift_J3Vector_ice3,
						      axisLengthBefore,
						      axisLengthNow))){
			haveBeenInsideIce3OuterVolume = 1;
		}
	    }

            if((shift_J3Vector_ice3==null && 
		!ParticleTracker.isInsideEarth(particleAxis_J3Line_center,
					       ice3Coordinate,
					       earthCoordinate,
					       ice3OuterVol,
					       haveBeenInsideIce3OuterVolume)) ||
	       (shift_J3Vector_ice3!=null && 
		!ParticleTracker.isInsideEarth(particleAxis_J3Line_center,
					       shift_J3Vector_ice3,
					       ice3Coordinate,
					       earthCoordinate,
					       ice3OuterVol,
					       haveBeenInsideIce3OuterVolume))){
		//System.out.println("Out of the Ice3Outer Volume");
		endLocation_J3Vector_ice3 = particleLocation_J3Vector_ice3;
		break;   // particles passed through your DETECTOR

	    }else{ // Collision occurs and the paticle's energy is lost

		double transferedEnergy  = event.collideNow(rand);

        // get current interaction's name
        String   curInteractionsName  = event.interactionsNameInPlay();

        // Collision occured by NC, WeakInteraction(CC, NC), Decay or Glashow Resonance?
        wasWeakInt = curInteractionsName.startsWith("Neutrino-Nuclen ");
        wasNC      = curInteractionsName.startsWith("Neutrino-Nuclen N");
		wasDecay   = (event.mcBaseInPlay.getTypeOfInteraction() == 1);
		/** For Glashow Resonance */
        wasGR = curInteractionsName.startsWith("Glashow ");
        wasGRLepton = curInteractionsName.startsWith("Glashow Resonance Leptonic ");

		// get produced particle's name
		//String   producedParticleName = Particle.particleName(
		//				event.getFlavorByInteractionsInPlay(),1);  
		//if(wasGRLepton){
		//   producedParticleName = Particle.particleName(
		//		           event.getFlavorByInteractionsInPlay(),0);
		//}

        // get current propagation particle
        Particle curPropParticle        = event.propParticle;
        int      curPropParticleFlavor  = curPropParticle.getFlavor();
        int      curPropParticleDoublet = curPropParticle.getDoublet();
        String   curPropParticleName    = Particle.particleName(curPropParticleFlavor,
                                                                curPropParticleDoublet);

                
		//System.out.println("Colliding via " + curInteractionsName +
		//            " and producing " + producedParticleName);
		//System.out.println("  -Current propagation particle : " + curPropParticleName);

        // add prpagation track to the lists
        boolean addTrack = false;
        if (curPropParticleFlavor == 1) { // muon flavor

		    if (wasWeakInt || wasDecay || wasGR) addTrack = true;

               } else if (curPropParticleFlavor == 2) { // tau flavor

                   if (wasWeakInt || wasGR) addTrack = true;

                } else if (curPropParticleFlavor == 0) { // e flavor
		    if(wasNC){
			addTrack = true;
		    }
		}

        if (addTrack) {
           trackParticleList.add(new Particle(curPropParticle.getFlavor(), 
                                              curPropParticle.getDoublet(),
                                              curPropParticle.getEnergy())); 
           trackLocationIce3List.add(particleLocation_J3Vector_ice3);
        }

                // add secondary cascade to the lists
		if ((event.getFlavorByInteractionsInPlay()==0 && !wasGR) || 
                    event.getFlavorByInteractionsInPlay()==3) { // 0 is electron flavor
		                                                // 3 is hadron flavor

		    particleList.add(new Particle(event.getFlavorByInteractionsInPlay(),
						  1,transferedEnergy));
		    locationIce3List.add(particleLocation_J3Vector_ice3);
            particleInteractionsList.add(curInteractionsName);
		}else if(wasGRLepton){ // The GR produced neurtino as a propagating particle
		    particleList.add(new Particle(event.getFlavorByInteractionsInPlay(),
						  0,transferedEnergy)); // neutrino
		    locationIce3List.add(particleLocation_J3Vector_ice3);
            particleInteractionsList.add(curInteractionsName);
		}else{  // producing secondary track. current version of JULIeT
                        // ignore secondary tracks.

		    //System.out.println("  -" + producedParticleName +
		    //                 " secondary track is now ignored");
		}
	    }

	    double afterInteractionLogEnergy = event.propParticle.getLogEnergy(); 
                                           // logEnergy after interaction
        double afterInteractionEnergy = event.propParticle.getEnergy();
                                           // Energy after interaction
        
        //System.out.println("logEnergy after interaction is  "+ afterInteractionLogEnergy);
        //System.out.println("Energy after interaction is  "+ afterInteractionEnergy);

	    if ((afterInteractionLogEnergy <= InteractionsBase.getLogEnergyProducedMinimum()) || 
                (wasNC && afterInteractionEnergy < neutrinoMinimumEnergyInTravel)){
		// Journey ends because all the primary energy has been lost
		endLocation_J3Vector_ice3 = particleLocation_J3Vector_ice3;
		//System.out.println("Energy has been lost...");
		break;
	    }

	    if (event.mcBaseInPlay.getTypeOfInteraction() == 1 && 
	       event.mcBaseInPlay.getProducedFlavor()!= 1){
               // Tau to Hadron/Electron decay
		endLocation_J3Vector_ice3 = particleLocation_J3Vector_ice3;
		break;
	    }

	} // while loop end


	//System.out.println("Journey ends with energy of " + 
	//		   event.propParticle.getEnergy() + " [GeV]");

    }


    /** move the propagating axis by given slant depth [g/cm^2] */
    public void moveParticleAxis(double depth){
	double l = point.getParticleLocation(); 
	           // The current particle location along the trajectory.
	double deltaL;
	double xSum = 0.0;
	double stepDx = 1.0e-1*depth;
	while(xSum<depth){
	    deltaL = stepDx/point.getMediumDensity(); 
	    l     += deltaL; 
	    xSum  += stepDx;
	    point.setParticleLocation(l);
	    particleAxis_J3Line_center.setAxisLength(l);
	}
    }


    /** Method to acquire all the listed variables: the cascade particles,
	their energies, and their locations which occured along the track.
	Returns number of the interactions. */
    public int getListedEvents(DataOutputStream out) throws IOException {

	out.writeBytes(primaryEnergy + " [GeV] " + 
		       event.propParticle.getEnergy() + " [GeV]");
	out.writeBytes(" Distance from the earth surface " +
		       startLocation + " [cm]");
	out.write('\n');
	particleIterator = particleList.listIterator();
	locationIce3Iterator = locationIce3List.listIterator();

	while(locationIce3Iterator.hasNext()){
	    Particle particle = (Particle )(particleIterator.next());
	    String particleName = particle.particleName(particle.getFlavor(),
							particle.getDoublet());
	    double energy = particle.getEnergy();
	    J3Vector r = (J3Vector )(locationIce3Iterator.next());
	    double x = r.getX();
	    double y = r.getY();
	    double z = r.getZ();
	    out.writeBytes(particleName + " " + energy + " [GeV] " + 
			   x + " [cm] " +  y + " [cm] " + z + " [cm]");
	    out.write('\n');
	}

	return locationIce3List.size();
    }


    /** Return particle hit iterator which allows
	an external object to access each Particle 
	stored in this object */
    public ListIterator getParticleIterator(){
	return particleList.listIterator();
    }
    public ListIterator getTrackParticleIterator(){
	return trackParticleList.listIterator();
    }
    public ListIterator getParticleInteractionsIterator(){
    return particleInteractionsList.listIterator();
    }


    /** Return location hit iterator which allows
	an external object to access each location
	stored in this object */
    public ListIterator getLocationIce3Iterator(){
	return locationIce3List.listIterator();
    }
    public ListIterator getTrackLocationIce3Iterator(){
	return trackLocationIce3List.listIterator();
    }

    public long[] getRandomState(){
        long[] state = rand.GetState();
        return state;
    }

    public void setRandomState(long[] state){
        rand = new RandomGenerator(state);
    }

    /** Method to run multiple events (numberOfEvent) 
        with various primary energies from
	logE = Particle.getLogEnergyMinimum() all the way up to 10^12 GeV.
        The results are stored in the matrix form and will be handled
	by EventMatrix.class in the event package.
    */
    public void runEventOnMatrix(int numberOfEvent, DataOutputStream out) throws IOException {
	
	for(int iLogE=0; iLogE<dim; iLogE++){
	    primaryiLogE = iLogE;
        //System.out.println(iLogE);
	    double logPrimaryEnergy = 
		Particle.getDeltaLogEnergy()*(double )iLogE + Particle.getLogEnergyMinimum();
	    primaryEnergy = Math.pow(10.0,logPrimaryEnergy);

	    for(int n=0; n<numberOfEvent; n++){
		runSingleEvent();
		int numberOfcollisions = getListedEvents(out);
	    }

	}
    }

    /** Change the Neutrino interaction weight in the InteractionBase */
    public static void setNeutrinoInteractionWeight(int weight){

	/** Neutrino factor */
	if(weight>=1){
	    InteractionsBase.neutrinoFactor = weight;
	    GlashowResonanceBase.neutrinoFactor = weight;
	}
    }


    /** Get the present neutrino interaction weight in the InteractionBase */
    public static int getNeutrinoInteractionWeight(){

	return  InteractionsBase.neutrinoFactor;

    }

	
}

