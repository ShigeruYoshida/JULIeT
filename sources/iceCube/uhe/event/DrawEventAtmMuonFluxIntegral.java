package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

public class DrawEventAtmMuonFluxIntegral {

    private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix() + 
   (int )(( Particle.getLogEnergyMinimum()-InteractionsBase.getLogEnergyProducedMinimum())/Particle.getDeltaLogEnergy());
    private static double[][] integralFlux;

    /** Data file names of the calculated propagation matricis. */

    private static String pathname = "../data/neutrino_earth/ice/";

    private static String[] matrixFileName = 
    { 
      "0Deg.data", "18_19Deg.data", "25_84Deg.data", "31_79Deg.data", "36_87Deg.data", 
      "41_41Deg.data", "45_57Deg.data", "49_46Deg.data", "53_13Deg.data", "56_63Deg.data", 
      "60_00Deg.data", "63_26Deg.data", "66_42Deg.data", "69_51Deg.data", "72_54Deg.data", 
      "75_52Deg.data", "77_00Deg.data", "78_46Deg.data", "79_92Deg.data", "81_37Deg.data",
      "82_82Deg.data", "84_26Deg.data", "85_70Deg.data", "87_13Deg.data", "88_56Deg.data" 
    };

/*
    private static String pathname = "../data/neutrino_earth/rock/";

    private static String[] matrixFileName = 
    { 
      "89_8Deg.data", "89_5Deg.data", "89Deg.data", "88Deg.data", "87Deg.data", "86Deg.data",
      "85Deg.data", "84Deg.data", "83Deg.data", "82Deg.data", "81Deg.data", "80Deg.data",
      "77_5Deg.data", "75Deg.data", "72_5Deg.data", "70Deg.data", "65Deg.data", "60Deg.data",
      "50Deg.data" 
    };

*/



    /** Zenith angle bound [deg] for solid angle estimation. */

    private static double[] zenithBound = 
    {
	0.0, 18.19, 25.84, 31.79, 36.87, 41.41, 45.57, 49.46, 53.13, 56.63,
	60.0, 63.26, 66.42, 69.51, 72.54, 75.52, 77.0, 78.46, 79.92, 81.37,
	82.82, 84.26, 85.70, 87.13, 88.56, 90.0
    };

/*
    private static double[] zenithBound = 
    {
	90.0, 89.7, 89.3, 88.5, 87.5, 86.5, 85.5, 84.5, 83.5, 82.5, 81.5, 80.5,
	78.75, 76.25, 73.75, 71.25, 67.5, 62.5, 55.0, 40.0
    };
*/




    public static void main(String[] args) throws IOException {

	int iLogE,ip;
	double logE;

	integralFlux = new double[3][dimension];
	for(ip=0;ip<3;ip++){
	    for(iLogE=0;iLogE<dimension;iLogE++){
		integralFlux[ip][iLogE] = 0.0;
	    }
	}

       int model = 1;
       String eventMatrixfileName = null;
       if(args.length!=1){
            System.out.println("Usage: DrawEventAtmMuonFluxIntegral eventMatrixName");
            System.exit(0);
        }else{
            eventMatrixfileName = args[0];
	}

        // Read the Event Matrix file
       DataInputStream in = new DataInputStream(new FileInputStream(eventMatrixfileName));
       EventAtmMuonFlux cascadeFlux = new EventAtmMuonFlux(in);
       in.close( );




       // Integration Loop
       int itheta;
       for(itheta=0;itheta<matrixFileName.length;itheta++){

	   // Read the serialized object of the Neutrino Charged Interaction Matrix
	   String fileName = pathname.concat(matrixFileName[itheta]);
	   in = new DataInputStream(new FileInputStream(fileName));
	   cascadeFlux.propMuonFlux.readMatrix(in);
	   in.close( );
	   System.err.println("Reading the matrix from " + fileName + " done.");

	   // Solid angle calculation
	   double radiansUp = Math.toRadians(zenithBound[itheta]);
	   double radiansDown = Math.toRadians(zenithBound[itheta+1]);
	   double solidAngle = 2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	   System.err.println("Zenith " + zenithBound[itheta] + " Solid angle " + solidAngle);
	   double cosTheta = Math.cos(radiansUp);

	   // EMG cascades
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		   Particle.getDeltaLogEnergy()*(double)iLogE;
	       integralFlux[0][iLogE] +=  
		   cascadeFlux.getDFEmgCascadeDLogE(logE,cosTheta)*solidAngle;
	   }

	   // Hadron cascades
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[1][iLogE] +=  
		   cascadeFlux.getDFHadronCascadeDLogE(logE,cosTheta)*solidAngle;
	   }

	   // Total cascades
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[2][iLogE] +=  
		   cascadeFlux.getDFTotalCascadeDLogE(logE,cosTheta)*solidAngle;
	   }


       }



       // Drawing
        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1])");
        System.out.println("scal 6.0 12.0 -14.0 -6.0");


	// EMG cascades
        System.out.println("lncl 2");
	double integral_emg = 0.0;
	double integral_emg1PeV = 0.0;
	double integral_emg10PeV = 0.0;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[0][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		if(iLogE>=400) integral_emg1PeV += Flux*Particle.getDeltaLogEnergy( );
		if(iLogE>=500) integral_emg10PeV += Flux*Particle.getDeltaLogEnergy( );
		integral_emg += Flux*Particle.getDeltaLogEnergy( );
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Hadron
        System.out.println("lncl 4");
	double integral_hadron = 0.0;
	double integral_hadron1PeV = 0.0;
	double integral_hadron10PeV = 0.0;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[1][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		if(iLogE>=400) integral_hadron1PeV += Flux*Particle.getDeltaLogEnergy( );
		if(iLogE>=500) integral_hadron10PeV += Flux*Particle.getDeltaLogEnergy( );
		integral_hadron += Flux*Particle.getDeltaLogEnergy( );
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Total
        System.out.println("lncl 1");
	double integral = 0.0;
	double integral_1PeV = 0.0;
	double integral_10PeV = 0.0;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[2][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		if(iLogE>=400) integral_1PeV += Flux*Particle.getDeltaLogEnergy( );
		if(iLogE>=500) integral_10PeV += Flux*Particle.getDeltaLogEnergy( );
		integral += Flux*Particle.getDeltaLogEnergy( );
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("mssg I_emg(1PeV) = " + integral_emg1PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_emg(10PeV) = " + integral_emg10PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_had(1PeV) = " + integral_hadron1PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_had(10PeV) = " + integral_hadron10PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I(1PeV) = " + integral_1PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I(10PeV) = " + integral_10PeV + " [cm^-2! sec^-1!]");

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
    }
}
