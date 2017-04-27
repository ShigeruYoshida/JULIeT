package iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.muonModel.*;

import java.io.*;

public class DrawAngularAtmMuonFlux {

    private static final double ln10 = Math.log(10.0);
    private static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    private static double[] integralFlux;

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
    {"../data/neutrino_earth/ice/","../data/neutrino_earth/rock"};


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
	double alpha = 1.757;
	double muonEthreshold = 100.0; // 100 GeV

	integralFlux = new double[50];

	int model = 1;

	if(args.length!=3){
            System.out.println("Usage: DrawAtmMuonFluxIntegral alpha muonEth Log(threshold E[GeV])");
            System.exit(0);
        }else{
            alpha = Double.valueOf(args[0]).doubleValue();  
            muonEthreshold = Double.valueOf(args[1]).doubleValue(); // with the given parameters 
	    logEth = Double.valueOf(args[2]).doubleValue();
	}

	iLogEth = (int )((logEth-Particle.getLogEnergyMinimum( ))/Particle.getDeltaLogEnergy( ));
	if(iLogEth<0){
	    System.err.println("Log(Threshold Energy[GeV]) " + logEth +
			       "must be grater than " + Particle.getLogEnergyMinimum( ));
	    System.exit(0);
	}

	PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux(alpha,muonEthreshold,true);




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
	       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	       muonFlux.readMatrix(in);
	       in.close( );
	       System.err.println("Reading the matrix from " + fileName + " done.");

	       // Solid angle calculation
	       double cosTheta = Math.cos(Math.toRadians(zenithAngle[upDown][itheta]));
	       System.err.println("Zenith " + zenithAngle[upDown][itheta] + " Costh " + cosTheta);

	       //Initialization
	       integralFlux[index] = 0.0;

	       // Mu
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[index] += 
		       muonFlux.getDFMuDLogE(logE,cosTheta)*Particle.getDeltaLogEnergy( );
	       }

	       index++;
	   }


       }



       System.out.println("titx Cos(theta)");
       System.out.println("tity Flux [cm^-2! sec^-1! sr^-1!])");
       System.out.println("scal -1.0 1.0 0.0 3.0E-19");


       index = 0;
       for(upDown=0;upDown<2;upDown++){
	   for(itheta=0;itheta<matrixFileName[upDown].length;itheta++){
	       double cosTheta = Math.cos(Math.toRadians(zenithAngle[upDown][itheta]));
	       // Drawing
	       double Flux = integralFlux[index];
	       System.out.println("data " + cosTheta + " 0.0 " + Flux + " 0.0");
	       index++;
	   }
       }

       System.out.println("join");
       System.out.println("disp");
       System.out.println("endg");

	    
    }
}
