package iceCube.uhe.particles;

import java.io.*;


/** 
    <pre>
    The Particle class to define their mass, names, flavors, lifetime
    It also provides the instance variables such as energy. 

             flavor  0      1      2        3
doublet 
    0              e-nu   mu-nu   tau-nu    hadron(pi0)
    1              e      muon    tauon     hadron(pi+)

    In adition this class can yield a maxtrix to accomodate dN/dlogE, 
    the energy distribution of particles. The matrix "logEnergyMatrix"
    corresponds to (dN/dlogE)_i (0&lt i = &lt DimensionOflogEnergyMatrix).
    The lowest logE and bin width in this table are defined by 
    "deltaLogE and alogEnergyMinimu", respectively.
    </pre>
*/

public class Particle implements Serializable {


    private double energy;        // Particle Energy [GeV].
    private double logEnergy;     // log(Particle Energy [GeV]).
    private double mass = -1.0;   // Particle mass [GeV].
    private double lifetime;      // Partilcle Life Time [sec].
    private int flavor;           // Particle Flavor.
    private int doublet;          // Lepton Doublet index.

    private static final long serialVersionUID = -8403086164562645730L;
    /** Number of the "Flavor" valuables to define the particle. */
    public final static int NumberOfFlavor = 4;

    /** Number of the "Doublet" valuables to define the particle.
	<pre>
     * 0 .. Neutrinos (or pi0 if Flabor = 3)
     * 1 .. Charged lepton (or pi+) 
        </pre>
    */
    public final static int NumberOfDoublet = 2;

    /** Particle Mass [GeV]*/
    public final static double[][] particleMasses = 
    {
	{0.0, 510.99906e-6},
	{0.0, 105.658389e-3},
	{0.0, 1.7841},
	{134.9743e-3, 139.5679e-3}
    };

    private final static double[][] particleLifeTimes =
    {
	{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY},
	{Double.POSITIVE_INFINITY, 2.19703e-6},
	{Double.POSITIVE_INFINITY, 3.05e-13},
	{8.4e-17,2.603e-8},
    };

    private final double ln10 = Math.log(10.0);


    // Variables to handle the energy distribution of the particles
    //  dN/dlogE. The matrix is used in numerical calculation of
    // the transport equation for particles propgatin in rock/ice.

    private static int DimensionOflogEnergyMatrix = 700;
                                              // Dimension of logEnergyMatrix 
    private static double deltaLogE = 0.01;   // Bin Width of logEnergyMatrix
    private static double logEnergyMinimum = 5.0;    // E >= 100TeV
    private double[] logEnergyMatrix;



    /** Default initial energy 1EeV = 1.0e9 GeV */
    public Particle(int initialFlavor, int initialDoublet) {
	this(initialFlavor, initialDoublet, 1.0e9);
    }


    /** Constructor for Checking the given flavor and doublet index.
	   <pre>
	   initialFlabor ... flavor valuable
	   initialDoblet ... doublet valuable
	   initialEnergy ... initial Energy [GeV]
	   </pre>
    */
    public Particle(int initialFlavor, int initialDoublet, 
		    double initialEnergy) {
	if(isValidFlavor(initialFlavor) && isValidDoublet(initialDoublet)){
	    this.flavor = initialFlavor;
	    this.doublet = initialDoublet;
	    this.mass = particleMasses[initialFlavor][initialDoublet];
	    this.lifetime = particleLifeTimes[initialFlavor][initialDoublet];
	    if(isValidEnergy(initialEnergy)){
		this.energy = initialEnergy;
		this.logEnergy = Math.log(initialEnergy)/ln10;
	    } else{
		System.err.println("Illegal Energy!");
		System.exit(0);
	    }
	} else{
	    System.err.println("Illegal Flavor/Doublet!");
	    System.exit(0);
	}

    }


    public static boolean isValidFlavor(int initialFlavor) {
	if(0 <= initialFlavor && initialFlavor < NumberOfFlavor) {
	    return true;
	} else{
	    return false;
	}
    }

    public static boolean isValidDoublet(int initialDoublet) {
	if(0 <= initialDoublet && initialDoublet < NumberOfDoublet) {
	    return true;
	} else{
	    return false;
	}
    }


    public boolean isValidEnergy(double initialEnergy) {
	if(mass>=0.0 && mass <= initialEnergy){ // Energy must be greater then
	    return true;                        // mass.
	} else{
	    return false;
	}
    }

    /** Sets the flavor and doublet. This is a new method to allow you to
     re-define your particle speice.*/
    public void setFlavorAndDoublet(int flavor, int doublet){
	if(isValidFlavor(flavor) && isValidDoublet(doublet)){
	    this.flavor = flavor;
	    this.doublet = doublet;
	    this.mass = particleMasses[flavor][doublet];
	    this.lifetime = particleLifeTimes[flavor][doublet];
	}
    }


    public int getFlavor( ){
	return flavor;
    }

    public int getDoublet( ){
	return doublet;
    }

    public double getMass( ){
	return mass;
    }

    public double getLifeTime( ){
	return lifetime;
    }

    public double getEnergy( ){
	return energy;
    }

    public double getLogEnergy( ){
	return logEnergy;
    }

    public void putEnergy(double newEnergy){
	if(isValidEnergy(newEnergy)){
	    energy = newEnergy;
	}
    }

    public void putLogEnergy(double newLogEnergy){
	if(newLogEnergy>=logEnergyMinimum){
	    logEnergy = newLogEnergy;
	}
    }

    /** Display the particle name. 
	Returns the String where the particle name is written. */
    public static String particleName(int flavor, int doublet){
	String name = "";
	String finalName;
	switch (flavor) {
	case 0: name = "Electron"; break;
	case 1: name = "Muon";break;
	case 2: name = "Tau";break;
	case 3: name = "Pion";break;
	default: System.err.println("Illegal Flavor Number!"); break;
	}

	if(flavor <3){
	    if(doublet == 0){
		finalName = name.concat(" Neutrino");
	    } else{
		finalName = new String(name);
	    }
	}else{
	    if(doublet == 1){
		finalName = name.concat("+");
	    } else{
		finalName = name.concat("0");
	    }
	}

	return finalName;
    }



    public static int getDimensionOfLogEnergyMatrix( ){
	return DimensionOflogEnergyMatrix;
    }

    public void putDimensionOfLogEnergyMatrix(int dimension){
	DimensionOflogEnergyMatrix = dimension;
    }

    public static double getLogEnergyMinimum( ){
	return logEnergyMinimum;
    }

    public void putLogEnergyMinimum(double newLogEnergy){
	logEnergyMinimum = newLogEnergy;
    }


    public static double getDeltaLogEnergy( ){
	return deltaLogE;
    }

    public void putDeltaLogEnergy(double newDeltaLogE){
	deltaLogE = newDeltaLogE;
    }

    public void generateLogEnergyMatrix( ){
	logEnergyMatrix = new double[DimensionOflogEnergyMatrix];
    }

    public double getLogEnergyMatrix(int ilogE){
	return logEnergyMatrix[ilogE];
    }


    public void putLogEnergyMatrix(double logE, double matrixElement){
	int ilogE;

	ilogE = (int )((logE-logEnergyMinimum)/deltaLogE+0.01);
	if(0<=ilogE && ilogE<DimensionOflogEnergyMatrix){
	    logEnergyMatrix[ilogE] = matrixElement;
	}else{
	    System.err.println("Illegal Energy Range "+ logE + "(" +ilogE+")");
	}

    }


    public void putLogEnergyMatrix(int ilogE, double matrixElement){
	logEnergyMatrix[ilogE] = matrixElement;
    }



}
