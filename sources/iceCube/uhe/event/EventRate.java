package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

public class EventRate {

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
    private static String pathname = "../data/neut_earth/rock/";

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
	double logE,logEcascadeThreshold=6.0;

	integralFlux = new double[3][Particle.getDimensionOfLogEnergyMatrix()];
	for(ip=0;ip<3;ip++){
	    for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
		integralFlux[ip][iLogE] = 0.0;
	    }
	}

       int model = 1;
       String eventMuMatrixfileName = null;
       String eventTauMatrixfileName = null;
       if(args.length!=3){
            System.out.println(
            "Usage: EventRate logEcascadeThreshold MuEventMatrixName TauEventMatrixName");
            System.exit(0);
        }else{
	    logEcascadeThreshold = Double.valueOf(args[0]).doubleValue();
            eventMuMatrixfileName = args[1];
            eventTauMatrixfileName = args[2];
	}

        // Read the Event Matrix file
       
       DataInputStream in = new DataInputStream(new FileInputStream(eventMuMatrixfileName));
       EventMonoEnergyFlux cascadeMuFlux = new EventMonoEnergyFlux(in);
       in.close( );
       in = new DataInputStream(new FileInputStream(eventTauMatrixfileName));
       EventMonoEnergyFlux cascadeTauFlux = new EventMonoEnergyFlux(in);
       in.close( );




       // Integration Loop
       int itheta;
       for(itheta=0;itheta<matrixFileName.length;itheta++){

	   // Read the serialized object of the Neutrino Charged Interaction Matrix
	   String fileName = pathname.concat(matrixFileName[itheta]);
	   in = new DataInputStream(new FileInputStream(fileName));
	   cascadeMuFlux.readMatrix(in);
	   in.close( );
	   in = new DataInputStream(new FileInputStream(fileName));
	   cascadeTauFlux.readMatrix(in);
	   in.close( );
	   System.err.println("Reading the matrix from " + fileName + " done.");

	   // Solid angle calculation
	   double radiansUp = Math.toRadians(zenithBound[itheta]);
	   double radiansDown = Math.toRadians(zenithBound[itheta+1]);
	   double solidAngle = 2.0*Math.PI*Math.abs(Math.cos(radiansDown)-Math.cos(radiansUp));
	   System.err.println("Zenith " + zenithBound[itheta] + " Solid angle " + solidAngle);

	   // Total cascades
	   for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
	       logE = Particle.getLogEnergyMinimum() +
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	       double logEcascade = logEcascadeThreshold;
	       while(logEcascade<=logE){
		   integralFlux[1][iLogE] +=  //nu-mu
		       cascadeMuFlux.getDFTotalCascadeDLogE(1,logE,logEcascade)*solidAngle
		       + cascadeTauFlux.getDFTotalCascadeDLogE(1,logE,logEcascade)*solidAngle;
		   integralFlux[2][iLogE] +=  //nu-tau
		       cascadeMuFlux.getDFTotalCascadeDLogE(2,logE,logEcascade)*solidAngle
		       + cascadeTauFlux.getDFTotalCascadeDLogE(2,logE,logEcascade)*solidAngle;
		   logEcascade += Particle.getDeltaLogEnergy( );
	       }
	   }


       }



	// Total
	for(iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
	    logE = Particle.getLogEnergyMinimum() +
		   Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = integralFlux[2][iLogE];

	    System.out.println(logE + " " +  integralFlux[1][iLogE] + " " + integralFlux[2][iLogE]);
	}

    }
}
