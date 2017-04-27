package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

import numRecipes.*;


public class WriteEventAtmFluxOverAllAngle {

    private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix() + 
   (int )(( Particle.getLogEnergyMinimum()-InteractionsBase.getLogEnergyProducedMinimum())/Particle.getDeltaLogEnergy());

    private static double[][] flux;
    private static double[] cosTable;

    /** Data file names of the calculated propagation matricis. */
    private static String[][] matrixFileName = 
    {
	{ 
      "0Deg.data", "18_19Deg.data", "25_84Deg.data", "31_79Deg.data", "36_87Deg.data", 
      "41_41Deg.data", "45_57Deg.data", "49_46Deg.data", "53_13Deg.data", "56_63Deg.data", 
      "60_00Deg.data", "63_26Deg.data", "66_42Deg.data", "69_51Deg.data", "72_54Deg.data", 
      "75_52Deg.data", "77_00Deg.data", "78_46Deg.data", "79_92Deg.data", "81_37Deg.data",
      "82_82Deg.data", "84_26Deg.data", "85_70Deg.data", "87_13Deg.data", "88_56Deg.data" 
	},

	{ 
      "89_8Deg.data", "89_5Deg.data", "89Deg.data", "88Deg.data", "87Deg.data", "86Deg.data",
      "85Deg.data", "84Deg.data", "83Deg.data", "82Deg.data", "81Deg.data", "80Deg.data",
      "77_5Deg.data", "75Deg.data", "72_5Deg.data", "70Deg.data", "65Deg.data", "60Deg.data",
      "50Deg.data" 
	}
    };

    private static String[] pathname =
    {"../data/neutrino_earth/ice/","../data/neutrino_earth/rock/"};


    /** Zenith angle [deg] */
    private static double[][] zenithAngle =
    {
	{
	0.0, 18.19, 25.84, 31.79, 36.87, 41.41, 45.57, 49.46, 53.13, 56.63,
	60.0, 63.26, 66.42, 69.51, 72.54, 75.52, 77.0, 78.46, 79.92, 81.37,
	82.82, 84.26, 85.70, 87.13, 88.56
	},

	{
       89.8, 89.5, 89.0, 88.0, 87.0, 86.0, 85.0, 84.0, 83.0, 82.0, 81.0, 80.0,
       77.5, 75.0, 72.5, 70.0, 65.0, 60.0, 50.0
	}
    };






    public static void main(String[] args) throws IOException {

	int iLogE,iLogEth,ip,upDown,index;
	double logE;
	double logEth = 6.0;

	int model = 1;
       String eventMatrixfileName = null;

	if(args.length!=2){
            System.out.println("Usage: WriteEventAtmFluxOverAllAngle Log(threshold E[GeV]) eventMatrixName");
            System.exit(0);
        }else{
	    logEth = Double.valueOf(args[0]).doubleValue();
            eventMatrixfileName = args[1];
	}

	iLogEth = (int )((logEth-InteractionsBase.getLogEnergyProducedMinimum())/
			 Particle.getDeltaLogEnergy( ));
	if(iLogEth<0){
	    System.err.println("Log(Threshold Energy[GeV]) " + logEth +
			       "must be grater than " + 
			       InteractionsBase.getLogEnergyProducedMinimum());
	    System.exit(0);
	}

	int numberOfCos = zenithAngle[0].length;

	flux = new double[dimension-iLogEth][numberOfCos];
	cosTable = new double[numberOfCos];

        // Read the Event Matrix file
       DataInputStream in = new DataInputStream(new FileInputStream(eventMatrixfileName));
       EventAtmMuonFlux cascadeFlux = new EventAtmMuonFlux(in);
       in.close( );




       // Angular Loop
       int itheta;
       for(itheta=0;itheta<zenithAngle[1].length;itheta++) 
	   zenithAngle[1][itheta] = 180.0-zenithAngle[1][itheta];
                                    // Nadir angle to Zenith angle
       index = 0;
       for(upDown=0;upDown<1;upDown++){
	   for(itheta=0;itheta<matrixFileName[upDown].length;itheta++){

	       // Read the serialized object of the Neutrino Charged Interaction Matrix
	       String fileName = pathname[upDown].concat(matrixFileName[upDown][itheta]);
	       in = new DataInputStream(new FileInputStream(fileName));
	       cascadeFlux.propMuonFlux.readMatrix(in);
	       in.close( );
	       System.err.println("Reading the matrix from " + fileName + " done.");

	       // Solid angle calculation
	       double cosTheta = Math.cos(Math.toRadians(zenithAngle[upDown][itheta]));
	       cosTable[index] = cosTheta;
	       System.err.println("Zenith " + zenithAngle[upDown][itheta] + " Costh " + cosTheta);

	       // obtain dF/dLogE_totalCascade
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = InteractionsBase.getLogEnergyProducedMinimum() + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   flux[iLogE-iLogEth][index] =  
		       cascadeFlux.getDFTotalCascadeDLogE(logE,cosTheta);
	       }

	       index++;


	   }
       }


       // Write the dF/dLogE over all phase space bin of solid angle
       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
	   logE = InteractionsBase.getLogEnergyProducedMinimum() + 
	       Particle.getDeltaLogEnergy( )*(double )iLogE;
	   for(int icos = 0; icos< 100; icos++){
	       double cosTheta = 1.0-0.01*(double )icos;
	       double dFdLogE = 
		   Interpolation.mThPolynominalInterpolate(cosTable,
							   flux[iLogE-iLogEth],
							   cosTable.length,
							   cosTheta,6);
	       System.out.println(logE + " " + cosTheta + " " + dFdLogE);
	   }
       }



    }
}
