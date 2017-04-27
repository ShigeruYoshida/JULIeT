package iceCube.uhe.neutrinoModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.neutrinoModel.*;

import java.io.*;

public class DrawNeutrinoFluxIntegral {

    private static final double ln10 = Math.log(10.0);
    private static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    private static double[][] integralFlux;

    /** Data file names of the calculated propagation matricis. */

    private static String[] matrixFileName = 
    { 
      "0Deg.data", "18_19Deg.data", "25_84Deg.data", "31_79Deg.data", "36_87Deg.data", 
      "41_41Deg.data", "45_57Deg.data", "49_46Deg.data", "53_13Deg.data", "56_63Deg.data", 
      "60_00Deg.data", "63_26Deg.data", "66_42Deg.data", "69_51Deg.data", "72_54Deg.data", 
      "75_52Deg.data", "77_00Deg.data", "78_46Deg.data", "79_92Deg.data", "81_37Deg.data",
      "82_82Deg.data", "84_26Deg.data", "85_70Deg.data", "87_13Deg.data", "88_56Deg.data" 
    };
/*

    private static String[] matrixFileName = 
    { 
      "89_8Deg.data", "89_5Deg.data", "89Deg.data", "88Deg.data", "87Deg.data", "86Deg.data",
      "85Deg.data", "84Deg.data", "83Deg.data", "82Deg.data", "81Deg.data", "80Deg.data",
      "77_5Deg.data", "75Deg.data", "72_5Deg.data", "70Deg.data", "65Deg.data", "60Deg.data",
      "50Deg.data" 
    };


    private static String[] matrixFileName = 
    { 
      "89_8DegNoPhotoNucl.data", "89_5DegNoPhotoNucl.data", "89DegNoPhotoNucl.data", 
      "88DegNoPhotoNucl.data", "87DegNoPhotoNucl.data", "86DegNoPhotoNucl.data",
      "85DegNoPhotoNucl.data", "84DegNoPhotoNucl.data", "83DegNoPhotoNucl.data", 
      "82DegNoPhotoNucl.data", "81DegNoPhotoNucl.data", "80DegNoPhotoNucl.data",
      "77_5DegNoPhotoNucl.data", "75DegNoPhotoNucl.data", "72_5DegNoPhotoNucl.data", 
      "70DegNoPhotoNucl.data", "65DegNoPhotoNucl.data", "60DegNoPhotoNucl.data",
      "50DegNoPhotoNucl.data" 
    };
*/


    //private static String pathname = "../data/neutrino_earth/rock/";

    private static String pathname = "../data/neutrino_earth/ice/";


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

	integralFlux = new double[5][dimension];
	for(ip=0;ip<5;ip++){
	    for(iLogE=0;iLogE<dimension;iLogE++){
		integralFlux[ip][iLogE] = 0.0;
	    }
	}

       int model = 1;
       if(args.length!=1){
            System.out.println("Usage: DrawNeutrinoFluxIntegral model-parameter");
            System.exit(0);
        }else{
	    model = Integer.valueOf(args[0]).intValue();
	}

       PropagatingNeutrinoFlux leptonFlux = new PropagatingNeutrinoFlux(model);




       // Integration Loop
       int itheta;
       for(itheta=0;itheta<matrixFileName.length;itheta++){

	   // Read the serialized object of the Neutrino Charged Interaction Matrix
	   String fileName = pathname.concat(matrixFileName[itheta]);
	   DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	   leptonFlux.readMatrix(in);
	   in.close( );
	   System.err.println("Reading the matrix from " + fileName + " done.");

	   // Solid angle calculation
	   double radiansUp = Math.toRadians(zenithBound[itheta]);
	   double radiansDown = Math.toRadians(zenithBound[itheta+1]);
	   double solidAngle = 2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	   System.err.println("Zenith " + zenithBound[itheta] + " Solid angle " + solidAngle);

	   // NuE
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = Particle.getLogEnergyMinimum( ) + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[0][iLogE] += leptonFlux.getDFNuEDLogE(iLogE)*solidAngle;
	   }

	   // NuMu
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = Particle.getLogEnergyMinimum( ) + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[1][iLogE] += leptonFlux.getDFNuMuDLogE(iLogE)*solidAngle;
	   }

	   // NuTau
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = Particle.getLogEnergyMinimum( ) + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[2][iLogE] += leptonFlux.getDFNuTauDLogE(iLogE)*solidAngle;
	   }

	   // Mu
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = Particle.getLogEnergyMinimum( ) + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[3][iLogE] += leptonFlux.getDFMuDLogE(iLogE)*solidAngle;
	   }

	   // Tau
	   for(iLogE=0;iLogE<dimension;iLogE++){
	       logE = Particle.getLogEnergyMinimum( ) + 
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       integralFlux[4][iLogE] += leptonFlux.getDFTauDLogE(iLogE)*solidAngle;
	   }

       }



       // Drawing
        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1])");
        System.out.println("scal 6.0 12.0 -15.0 -7.0");


	// NuE
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[0][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// NuMu
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[1][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// NuTau
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[2][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Mu
	double integral_mu = 0.0;
	double integral_mu1PeV = 0.0;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[3][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		if(iLogE>=100) integral_mu += Flux*Particle.getDeltaLogEnergy( );
		integral_mu1PeV += Flux*Particle.getDeltaLogEnergy( );
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Tau
	double integral_tau = 0.0;
	double integral_tau1PeV = 0.0;
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = Particle.getLogEnergyMinimum( ) + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[4][iLogE];
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		if(iLogE>=100) integral_tau += Flux*Particle.getDeltaLogEnergy( );
		integral_tau1PeV += Flux*Particle.getDeltaLogEnergy( );
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("mssg I_mu(10PeV) = " + integral_mu + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_mu(1PeV) = " + integral_mu1PeV + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_tau(10PeV) = " + integral_tau + " [cm^-2! sec^-1!]");
	System.out.println("mssg I_tau(1PeV) = " + integral_tau1PeV + " [cm^-2! sec^-1!]");

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
    }
}
