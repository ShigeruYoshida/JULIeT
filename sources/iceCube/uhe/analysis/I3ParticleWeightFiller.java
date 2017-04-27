package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 
    This class provides the methods to fill I3Particles with GZK and 
    Atmospheric muon weights. Every methods are static so that you can
    call them directly from your object.

    Written by S. Yoshida 2007 February 1st
*/

public class I3ParticleWeightFiller {

    /** Data file names of the calculated propagation matricis. */
    protected static String[][] matrixFileName = 
    {
	{
    "1e5_00cm.data", "1e5_05cm.data", "1e5_10cm.data", "1e5_15cm.data", "1e5_20cm.data", 
    "1e5_25cm.data", "1e5_30cm.data", "1e5_35cm.data", "1e5_40cm.data", "1e5_45cm.data", 
    "1e5_50cm.data", "1e5_55cm.data", "1e5_60cm.data", "1e5_65cm.data", "1e5_70cm.data", 
    "1e5_75cm.data", "1e5_80cm.data", "1e5_85cm.data", "1e5_90cm.data", "1e5_95cm.data", 
    "1e6_00cm.data", "1e6_10cm.data", "1e6_20cm.data", "1e6_30cm.data", "1e6_40cm.data", 
    "1e6_50cm.data", "1e6_60cm.data", "1e6_70cm.data", "1e6_80cm.data", "1e6_90cm.data", 
    "1e7_00cm.data", "1e7_10cm.data", "1e7_20cm.data", "1e7_30cm.data", "1e7_40cm.data", "1e7_50cm.data"
	},
 
	{ 
      "50Deg.data",  "60Deg.data",    "65Deg.data",   "70Deg.data",  "72_5Deg.data",   
      "75Deg.data",  "77_5Deg.data",  "80Deg.data",   "81Deg.data",  "82Deg.data",   
      "83Deg.data",  "84Deg.data",    "85Deg.data",   "86Deg.data",  "87Deg.data",  
      "88Deg.data",  "89Deg.data",    "89_5Deg.data",  "89_8Deg.data"
	},

	{
	"0Deg.data", "18_19Deg.data", "25_84Deg.data", "31_79Deg.data", "36_87Deg.data",  
	"41_41Deg.data", "45_57Deg.data", "49_46Deg.data", "53_13Deg.data", "56_63Deg.data",  
	"60_00Deg.data", "63_26Deg.data", "66_42Deg.data", "69_51Deg.data", "72_54Deg.data",  
	"75_52Deg.data", "77_00Deg.data", "78_46Deg.data", "79_92Deg.data", "81_37Deg.data", 
	"82_82Deg.data", "84_26Deg.data", "85_70Deg.data", "87_13Deg.data", "88_56Deg.data"  
	}

    };

    protected static String[] pathname =
    {"../data/neutrino_earth/ice/","../data/neutrino_earth/rock/"};


    protected static double[] rockrange ={//nadir
      0.0, 55.0, 62.5, 67.5, 71.25, 
      73.75, 76.25, 78.75, 80.5, 81.5, 
      82.5, 83.5, 84.5, 85.5, 86.5, 
      87.5, 88.5, 89.25, 89.7, 90.0};

    protected static double[] icerange ={//zenith
    0.0, 18.19, 25.84, 31.79, 36.87, 41.41, 45.57, 49.46, 53.13, 56.63, 
    60.0, 63.26, 66.42, 69.51, 72.54, 75.52, 77.0, 78.46, 79.92, 81.37, 
    82.82, 84.26, 85.70, 87.13, 88.56, 90.0};

    
    protected static double[] icedist ={//distance in log10(distance[cm])
      5.00,  5.025, 5.075, 5.125, 5.175, 5.225, 5.275, 5.325, 5.375, 5.425, 
      5.475, 5.525, 5.575, 5.625, 5.675, 5.725, 5.775, 5.825, 5.875, 5.925, 
      5.975, 6.05,  6.15,  6.25,  6.35,  6.45,  6.55,  6.65,  6.75,  6.85, 
      6.95,  7.05,  7.15,  7.25,  7.35,  7.45, 7.50};

    /** A minimum threshold to fill weights. For the faster processing */
    protected static int minNDOMsToFill = 60;

    private final static double ln10 = Math.log(10.0);

    /** Fill the weight besed on PropagatingAtmMuonFlux in the muonModel package */
    public static void fillPropagatingAtmMuonFluxWeight(I3Particle iceParticle, 
						    PropagatingAtmMuonFlux muonFlux,
						    String fluxName) throws IOException{


	    J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	    double distance = 
		iceParticle.getDistanceFromEarthSurfaceToIceCube();
	    boolean isDownWard = true;
	    double cosZenith = -n_ice.getZ(); // Reversed vector
	    if(cosZenith<0.0) isDownWard = false; // this is an upward event

	    double flux = 0.0;
	    if(iceParticle.getIceCubeData().getNDOMsLaunch()>= minNDOMsToFill){
		// OK, let's get the flux from the matrix data
		readPropagationMatrix(n_ice,distance,isDownWard,muonFlux);
		double logE = iceParticle.getLogEnergy();
		flux =  muonFlux.getDFMuDLogE(logE,cosZenith);
	    }
	    iceParticle.setAtmosphericMuonFlux(flux,fluxName);
    }

    /** Fill the CR energy distribution  besed on PropagatingAtmMuonFlux in the muonModel package */
    public static void fillCRFluxWeight(I3Particle iceParticle, 
					PropagatingAtmMuonFlux muonFlux) 
	throws IOException{

	J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	double distance = 
	    iceParticle.getDistanceFromEarthSurfaceToIceCube();
	boolean isDownWard = true;
	double cosZenith = -n_ice.getZ(); // Reversed vector
	if(cosZenith<0.0) isDownWard = false; // this is an upward event

	if(iceParticle.getIceCubeData().getNDOMsLaunch()>= minNDOMsToFill){
	    // OK, let's get the CR energy distribution
	    readPropagationMatrix(n_ice,distance,isDownWard,muonFlux);
	    iceParticle.generateLogEnergyMatrix();

	    double logInIceMuonEnergy = iceParticle.getLogEnergy();

	    double epsilon = 1.0e-4; // round-off margin for binning
	    double logCREnergy = Particle.getLogEnergyMinimum() + epsilon;
	    double logCRMaxEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	    int iLogE = 0;
	    while(logCREnergy < logCRMaxEnergy){ // E < 10^12 GeV
		double flux =  muonFlux.getDFMuDLogCREDLogE(logCREnergy, Particle.getDeltaLogEnergy()*0.5,
							    logInIceMuonEnergy, cosZenith)*
		    Particle.getDeltaLogEnergy();

		iceParticle.putLogEnergyMatrix(iLogE,flux);

		iLogE++;
		logCREnergy += Particle.getDeltaLogEnergy();
	    }
	}
    }

    /** Fill the weight besed on CosmicRayFlux in the muonModel package.
	This method should be used for Corsika-mmc I3Particles in which
	the primary cosmic ray energy is available and stored by
	I3Particle.putRecoEnergy(energy).
     */
    public static void fillPropagatingAtmMuonFluxWeight(I3Particle iceParticle, 
						    CosmicRayFlux crFlux,
						    String fluxName) {


	double logCosmicRayEnergy = iceParticle.getLogRecoEnergy();
	double flux =  crFlux.getDFDLogE(logCosmicRayEnergy);
	iceParticle.setAtmosphericMuonFlux(flux,fluxName);
    }


    /** Fill the weight besed on PropagatingNeutrinoFlux in the neutrinoModel package */
    public static void fillPropagatingNeutrinoFluxWeight(I3Particle iceParticle, 
						    PropagatingNeutrinoFlux leptonFlux,
						    String fluxName) throws IOException{

	    J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	    double distance = 
		iceParticle.getDistanceFromEarthSurfaceToIceCube();
	    boolean isDownWard = true;
	    double cosZenith = -n_ice.getZ(); // Reversed vector
	    if(cosZenith<0.0) isDownWard = false; // this is an upward event

	    double flux = 0.0;
	    if(iceParticle.getIceCubeData().getNDOMsLaunch()>= minNDOMsToFill){
		// OK, let's get the flux from the matrix data
		readPropagationMatrix(n_ice,distance,isDownWard,leptonFlux);
		double logE = iceParticle.getLogEnergy();
		if(iceParticle.getFlavor()== 1 && iceParticle.getDoublet() == 1){
		    flux =  leptonFlux.getDFMuDLogE(logE);
		}else if(iceParticle.getFlavor()== 1 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuMuDLogE(logE);
		}else if(iceParticle.getFlavor()== 2 && iceParticle.getDoublet() == 1){
		    flux =  leptonFlux.getDFTauDLogE(logE);
		}else if(iceParticle.getFlavor()== 2 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuTauDLogE(logE);
		}else if(iceParticle.getFlavor()== 0 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuEDLogE(logE);
		}
	    }
	    iceParticle.setGZKNeutrinoFlux(flux,fluxName);
    }

    /** Fill the weight besed on PropagatingNeutrinoFlux in the neutrinoModel package. 
        The flux specified by the model number given in the argument is calculated
        and put in the weight map. Refer NeutrinoFlux.java for the model number definition.
        You have to call readPropagationMatrix() first before usng this method.
     */
    public static void fillPropagatingNeutrinoFluxWeight(I3Particle iceParticle, 
						 PropagatingNeutrinoFlux leptonFlux,
						 int modelNumber,
						 String fluxName) throws IOException{

	    double flux = 0.0;
	    if(iceParticle.getIceCubeData().getNDOMsLaunch()>= minNDOMsToFill){

		double logE = iceParticle.getLogEnergy();
		if(iceParticle.getFlavor()== 1 && iceParticle.getDoublet() == 1){
		    flux =  leptonFlux.getDFMuDLogE(modelNumber,logE);
		}else if(iceParticle.getFlavor()== 1 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuMuDLogE(modelNumber,logE);
		}else if(iceParticle.getFlavor()== 2 && iceParticle.getDoublet() == 1){
		    flux =  leptonFlux.getDFTauDLogE(modelNumber,logE);
		}else if(iceParticle.getFlavor()== 2 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuTauDLogE(modelNumber,logE);
		}else if(iceParticle.getFlavor()== 0 && iceParticle.getDoublet() == 0){
		    flux =  leptonFlux.getDFNuEDLogE(modelNumber,logE);
		}
	    }
	    iceParticle.setGZKNeutrinoFlux(flux,fluxName);
    }

    /** 
	Reading the propagation matrix stored in the file corresponsing
	to the track geometry. 
	<pre>
	J3UnitVector n_ice : Direction of I3Particle trajectory in the IceCube Coordinate
	double distance    : Propagation distance of I3Particle from the earth surface [cm]
	boolean isDownWard :   whether this event is downward-going or not
	</pre>
    */
    private static void readPropagationMatrix(J3UnitVector n_ice, double distance,
					      boolean isDownWard,
					      PropagatingAtmMuonFlux muonFlux) 
	throws IOException {

	String matrixFileName = getMatrixFileName(n_ice,distance,isDownWard);
	System.err.println(" --matrix filename " + matrixFileName);
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(matrixFileName));
	muonFlux.readMatrix(in);
	in.close( );

    }

    /** 
	Reading the propagation matrix stored in the file corresponsing
	to the track geometry. 
	<pre>
	J3UnitVector n_ice    : Direction of I3Particle trajectory in the IceCube Coordinate
	double distance    : Propagation distance of I3Particle from the earth surface [cm]
	boolean isDownWard :   whether this event is downward-going or not
	</pre>
    */
    private static void readPropagationMatrix(J3UnitVector n_ice, double distance,
					      boolean isDownWard,
					      PropagatingNeutrinoFlux neutrinoFlux) 
	throws IOException {

	String matrixFileName = getMatrixFileName(n_ice,distance,isDownWard);
	System.err.println(" --matrix filename " + matrixFileName);
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(matrixFileName));
	neutrinoFlux.readMatrix(in);
	in.close( );

    }


    protected static String getMatrixFileName(J3UnitVector n_ice, double distanceFromEarth,
				     boolean isDownWard){

	String fileName = null;
	if(isDownWard){
	    double logDistance = Math.log(distanceFromEarth)/ln10;
	    System.err.println(" Distance = " + distanceFromEarth + 
			       " logDistance = " + logDistance);

      if((logDistance>=icedist[0]      )&&(logDistance<icedist[1] )) fileName = pathname[0] + matrixFileName[0][0];
      else if((logDistance>=icedist[1] )&&(logDistance<icedist[2] )) fileName = pathname[0] + matrixFileName[0][1];
      else if((logDistance>=icedist[2] )&&(logDistance<icedist[3] )) fileName = pathname[0] + matrixFileName[0][2];
      else if((logDistance>=icedist[3] )&&(logDistance<icedist[4] )) fileName = pathname[0] + matrixFileName[0][3];
      else if((logDistance>=icedist[4] )&&(logDistance<icedist[5] )) fileName = pathname[0] + matrixFileName[0][4];
      else if((logDistance>=icedist[5] )&&(logDistance<icedist[6] )) fileName = pathname[0] + matrixFileName[0][5];
      else if((logDistance>=icedist[6] )&&(logDistance<icedist[7] )) fileName = pathname[0] + matrixFileName[0][6];
      else if((logDistance>=icedist[7] )&&(logDistance<icedist[8] )) fileName = pathname[0] + matrixFileName[0][7];
      else if((logDistance>=icedist[8] )&&(logDistance<icedist[9] )) fileName = pathname[0] + matrixFileName[0][8];
      else if((logDistance>=icedist[9] )&&(logDistance<icedist[10])) fileName = pathname[0] + matrixFileName[0][9];
      else if((logDistance>=icedist[10])&&(logDistance<icedist[11])) fileName = pathname[0] + matrixFileName[0][10];
      else if((logDistance>=icedist[11])&&(logDistance<icedist[12])) fileName = pathname[0] + matrixFileName[0][11];
      else if((logDistance>=icedist[12])&&(logDistance<icedist[13])) fileName = pathname[0] + matrixFileName[0][12];
      else if((logDistance>=icedist[13])&&(logDistance<icedist[14])) fileName = pathname[0] + matrixFileName[0][13];
      else if((logDistance>=icedist[14])&&(logDistance<icedist[15])) fileName = pathname[0] + matrixFileName[0][14];
      else if((logDistance>=icedist[15])&&(logDistance<icedist[16])) fileName = pathname[0] + matrixFileName[0][15];
      else if((logDistance>=icedist[16])&&(logDistance<icedist[17])) fileName = pathname[0] + matrixFileName[0][16];
      else if((logDistance>=icedist[17])&&(logDistance<icedist[18])) fileName = pathname[0] + matrixFileName[0][17];
      else if((logDistance>=icedist[18])&&(logDistance<icedist[19])) fileName = pathname[0] + matrixFileName[0][18];
      else if((logDistance>=icedist[19])&&(logDistance<icedist[20])) fileName = pathname[0] + matrixFileName[0][19];
      else if((logDistance>=icedist[20])&&(logDistance<icedist[21])) fileName = pathname[0] + matrixFileName[0][20];
      else if((logDistance>=icedist[21])&&(logDistance<icedist[22])) fileName = pathname[0] + matrixFileName[0][21];
      else if((logDistance>=icedist[22])&&(logDistance<icedist[23])) fileName = pathname[0] + matrixFileName[0][22];
      else if((logDistance>=icedist[23])&&(logDistance<icedist[24])) fileName = pathname[0] + matrixFileName[0][23];
      else if((logDistance>=icedist[24])&&(logDistance<icedist[25])) fileName = pathname[0] + matrixFileName[0][24];
      else if((logDistance>=icedist[25])&&(logDistance<icedist[26])) fileName = pathname[0] + matrixFileName[0][25];
      else if((logDistance>=icedist[26])&&(logDistance<icedist[27])) fileName = pathname[0] + matrixFileName[0][26];
      else if((logDistance>=icedist[27])&&(logDistance<icedist[28])) fileName = pathname[0] + matrixFileName[0][27];
      else if((logDistance>=icedist[28])&&(logDistance<icedist[29])) fileName = pathname[0] + matrixFileName[0][28];
      else if((logDistance>=icedist[29])&&(logDistance<icedist[30])) fileName = pathname[0] + matrixFileName[0][29];
      else if((logDistance>=icedist[30])&&(logDistance<icedist[31])) fileName = pathname[0] + matrixFileName[0][30];
      else if((logDistance>=icedist[31])&&(logDistance<icedist[32])) fileName = pathname[0] + matrixFileName[0][31];
      else if((logDistance>=icedist[32])&&(logDistance<icedist[33])) fileName = pathname[0] + matrixFileName[0][32];
      else if((logDistance>=icedist[33])&&(logDistance<icedist[34])) fileName = pathname[0] + matrixFileName[0][33];
      else if((logDistance>=icedist[34])&&(logDistance<icedist[35])) fileName = pathname[0] + matrixFileName[0][34];
      else if((logDistance>=icedist[35])&&(logDistance<icedist[36])) fileName = pathname[0] + matrixFileName[0][35];
      else{
	  double cosZenith = -n_ice.getZ(); // Reversed vector
	  double zenithAngle = Math.acos(cosZenith);
	  double zenithAngleDeg = Math.toDegrees(zenithAngle);
	  if((zenithAngleDeg>= icerange[10]) && (zenithAngleDeg <icerange[11] )) 
	      fileName = pathname[0] + matrixFileName[2][10];
	  else if((zenithAngleDeg>= icerange[11]) && (zenithAngleDeg <icerange[12] )) 
	      fileName = pathname[0] + matrixFileName[2][11];
	  else if((zenithAngleDeg>= icerange[12]) && (zenithAngleDeg <icerange[13] )) 
	      fileName = pathname[0] + matrixFileName[2][12];
	  else if((zenithAngleDeg>= icerange[13]) && (zenithAngleDeg <icerange[14] )) 
	      fileName = pathname[0] + matrixFileName[2][13];
	  else if((zenithAngleDeg>= icerange[14]) && (zenithAngleDeg <icerange[15] )) 
	      fileName = pathname[0] + matrixFileName[2][14];
	  else if((zenithAngleDeg>= icerange[15]) && (zenithAngleDeg <icerange[16] )) 
	      fileName = pathname[0] + matrixFileName[2][15];
	  else if((zenithAngleDeg>= icerange[16]) && (zenithAngleDeg <icerange[17] )) 
	      fileName = pathname[0] + matrixFileName[2][16];
	  else if((zenithAngleDeg>= icerange[17]) && (zenithAngleDeg <icerange[18] )) 
	      fileName = pathname[0] + matrixFileName[2][17];	
	  else if((zenithAngleDeg>= icerange[18]) && (zenithAngleDeg <icerange[19] )) 
	      fileName = pathname[0] + matrixFileName[2][18];	
	  else if((zenithAngleDeg>= icerange[19]) && (zenithAngleDeg <icerange[20] )) 
	      fileName = pathname[0] + matrixFileName[2][19];	
	  else if((zenithAngleDeg>= icerange[20]) && (zenithAngleDeg <icerange[21] )) 
	      fileName = pathname[0] + matrixFileName[2][20];	
	  else if((zenithAngleDeg>= icerange[21]) && (zenithAngleDeg <icerange[22] )) 
	      fileName = pathname[0] + matrixFileName[2][21];	
	  else if((zenithAngleDeg>= icerange[22]) && (zenithAngleDeg <icerange[23] )) 
	      fileName = pathname[0] + matrixFileName[2][22];	
	  else if((zenithAngleDeg>= icerange[23]) && (zenithAngleDeg <icerange[24] )) 
	      fileName = pathname[0] + matrixFileName[2][23];	
	  else if((zenithAngleDeg>= icerange[24]) && (zenithAngleDeg <icerange[25] )) 
	      fileName = pathname[0] + matrixFileName[2][24];
	  else{
	      System.err.println("error getting filename!");
	      System.err.println("NadirAngle(" + zenithAngleDeg + ")[deg] distance="
				 + distanceFromEarth);
	      System.exit(0);
	  }
      }

	}else{

	    double cosZenith = -n_ice.getZ(); // Reversed vector
	    double zenithAngle = Math.acos(cosZenith);
	    double nadirAngle = Math.PI - zenithAngle;
	    double nadirAngleDeg = Math.toDegrees(nadirAngle);
	    System.err.println(" NadirAngle = " + nadirAngleDeg);

      if((nadirAngleDeg>= rockrange[0]     ) && (nadirAngleDeg <rockrange[1] )) fileName = pathname[1] + matrixFileName[1][0];
      else if((nadirAngleDeg>=rockrange[1] ) && (nadirAngleDeg<rockrange[2] )) fileName = pathname[1] + matrixFileName[1][1];
      else if((nadirAngleDeg>=rockrange[2] ) && (nadirAngleDeg<rockrange[3] )) fileName = pathname[1] + matrixFileName[1][2];
      else if((nadirAngleDeg>=rockrange[3] ) && (nadirAngleDeg<rockrange[4] )) fileName = pathname[1] + matrixFileName[1][3];
      else if((nadirAngleDeg>=rockrange[4] ) && (nadirAngleDeg<rockrange[5] )) fileName = pathname[1] + matrixFileName[1][4];
      else if((nadirAngleDeg>=rockrange[5] ) && (nadirAngleDeg<rockrange[6] )) fileName = pathname[1] + matrixFileName[1][5];
      else if((nadirAngleDeg>=rockrange[6] ) && (nadirAngleDeg<rockrange[7] )) fileName = pathname[1] + matrixFileName[1][6];
      else if((nadirAngleDeg>=rockrange[7] ) && (nadirAngleDeg<rockrange[8] )) fileName = pathname[1] + matrixFileName[1][7];
      else if((nadirAngleDeg>=rockrange[8] ) && (nadirAngleDeg<rockrange[9] )) fileName = pathname[1] + matrixFileName[1][8];
      else if((nadirAngleDeg>=rockrange[9] ) && (nadirAngleDeg<rockrange[10])) fileName = pathname[1] + matrixFileName[1][9];
      else if((nadirAngleDeg>=rockrange[10]) && (nadirAngleDeg<rockrange[11])) fileName = pathname[1] + matrixFileName[1][10];
      else if((nadirAngleDeg>=rockrange[11]) && (nadirAngleDeg<rockrange[12])) fileName = pathname[1] + matrixFileName[1][11];
      else if((nadirAngleDeg>=rockrange[12]) && (nadirAngleDeg<rockrange[13])) fileName = pathname[1] + matrixFileName[1][12];
      else if((nadirAngleDeg>=rockrange[13]) && (nadirAngleDeg<rockrange[14])) fileName = pathname[1] + matrixFileName[1][13];
      else if((nadirAngleDeg>=rockrange[14]) && (nadirAngleDeg<rockrange[15])) fileName = pathname[1] + matrixFileName[1][14];
      else if((nadirAngleDeg>=rockrange[15]) && (nadirAngleDeg<rockrange[16])) fileName = pathname[1] + matrixFileName[1][15];
      else if((nadirAngleDeg>=rockrange[16]) && (nadirAngleDeg<rockrange[17])) fileName = pathname[1] + matrixFileName[1][16];
      else if((nadirAngleDeg>=rockrange[17]) && (nadirAngleDeg<rockrange[18])) fileName = pathname[1] + matrixFileName[1][17];
      else if((nadirAngleDeg>=rockrange[18]) && (nadirAngleDeg<rockrange[19])) fileName = pathname[1] + matrixFileName[1][18];
      else {
	  System.err.println("error getting filename!");
	  System.err.println("NadirAngle(" + nadirAngleDeg + ")[deg] cosZenith="
			      + cosZenith);
	  System.exit(0);
      }
	}

      return fileName;

    }


    /** Main method -- Reading out the stored I3Particles and fills the weight */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	int neutrinoModel = 4; // Model parameter for NeutrinoFlux class
	String neutrinoModelName = null;   // model name of Neutrino Flux
	String atmMuonModelName = null;    // model name of Atmospheic Muon Flux
	String bundleCELModelName = null;  // model name of CEL approximated bundle model name
	double alpha = 1.99;              // multiplicity parameter of the bundle
	double muEthreshold = 1175.0;     // 1175 GeV - muon threhold in a bundle
	boolean cutOffExists = false;
        ParticlePoint s = null;


	PropagatingNeutrinoFlux leptonFlux = null;
	PropagatingAtmMuonFlux muonFlux = null;
	AtmMuonBundleFlux bundleFlux = null;

	if(args.length<2){
	    System.out.println("Usage: I3ParticleWeightFiller input-file-name output-file-name");
	    System.exit(0);
        }else { 
            inputFileName = args[0];
            outputFileName = args[1];
        }


	// Interactive session to set the flux(es) to be filled
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 

	boolean fillGZK = false;
        System.out.print("GZK fill [yes(1)/no(0)] ->"); 
        buffer   = d.readLine(); 
        if(Integer.valueOf(buffer).intValue()==1) fillGZK = true; 
	if(fillGZK){
	    System.out.print(" GZK model name ->"); 
	    neutrinoModelName   = d.readLine(); 
	    System.out.print(" model number (default:4 GZK m=4 Zmax=4) ->"); 
	    buffer   = d.readLine(); 
	    neutrinoModel = Integer.valueOf(buffer).intValue();
	}

	boolean fillAtmMuon = false;
	boolean fillCELBundle = false;
        System.out.print("Atm Muon fill [yes(1)/no(0)] ->"); 
        buffer   = d.readLine(); 
        if(Integer.valueOf(buffer).intValue()==1) fillAtmMuon = true; 
	if(fillAtmMuon){
	    System.out.print(" Atm Muon model name ->"); 
	    atmMuonModelName   = d.readLine();
	    System.out.print(" alpha  ->"); 
	    buffer   = d.readLine(); 
	    alpha = Double.valueOf(buffer).doubleValue();
	    System.out.print(" Muon Energy Threshold  [GeV]->"); 
	    buffer   = d.readLine(); 
	    muEthreshold = Double.valueOf(buffer).doubleValue();
	    System.out.print(" with GZK cutoff feature? [yes(1)/no(0)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) cutOffExists = true; 
	    System.out.print(" CEL approximation? [yes(1)/no(0)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) fillCELBundle = true; 

	}


	// Dislay the model setting you entered for double-checking
	if(fillGZK){
	    System.out.println("Neutrino model Name " + neutrinoModelName +
			       " model number(" + neutrinoModel + ")");
	}
	if(fillAtmMuon && ! fillCELBundle){
	    System.out.println("PropagatingNeutrinoFlux name " + atmMuonModelName +
			       "alpha(" + alpha + ") Eth=" + muEthreshold);
	}
	if(fillAtmMuon &&  fillCELBundle){
	    System.out.println("AtmMuonBunxleFlux name " + atmMuonModelName +
			       "alpha(" + alpha + ") Eth=" + muEthreshold);
	}


	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);


	// PropagatingNeutrinoFlux object
	if(fillGZK)
	leptonFlux = new PropagatingNeutrinoFlux(neutrinoModel);
	// PropagatingAtmMuonFlux object
	if(fillAtmMuon && ! fillCELBundle)
	muonFlux = new PropagatingAtmMuonFlux(alpha,muEthreshold, cutOffExists);
	// AtmMuonBundleFlux with the CEL approximation object
	if(fillAtmMuon && fillCELBundle){
	    bundleFlux = new AtmMuonBundleFlux();
	    bundleFlux.setAlpha(alpha); 
	    bundleFlux.setMuonThresholdEnergy(muEthreshold); 
	    bundleFlux.setCutOffFeature(cutOffExists); 
	    s = new ParticlePoint(0.0,0.0,0); // ice 
	}

	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    // Remove the previousely stored weight by the same model name
	    // if exists
	    if(fillGZK) iceParticle.removeGZKNeutrinoFlux(neutrinoModelName);
	    if(fillAtmMuon) iceParticle.removeAtmosphericMuonFlux(atmMuonModelName);

	    // Fills the flux weights
	    if(fillGZK){
		I3ParticleWeightFiller.fillPropagatingNeutrinoFluxWeight(iceParticle,
									 leptonFlux,
									 neutrinoModelName);
	    }
	    if(fillAtmMuon && ! fillCELBundle){
		I3ParticleWeightFiller.fillPropagatingAtmMuonFluxWeight(iceParticle,
									muonFlux,
									atmMuonModelName);
	    }
	    if(fillAtmMuon && fillCELBundle){
		AtmMuonBundleFitter.fluxWeightName=atmMuonModelName;
		AtmMuonBundleFitter.fillAtmMuonBundleFluxWeight(iceParticle,s,bundleFlux);
	    }

	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }
}
