package iceCube.uhe.neutrinoModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.neutrinoModel.*;

import java.io.*;

public class DrawAngularQuickNeutrinoFlux {

    private static final double ln10 = Math.log(10.0);
    private static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    private static double[][] integralFlux;

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
	double nuCCEnhancementFactor = 1.0;

	integralFlux = new double[5][50];

	int model = 1;

	if(args.length!=3){
            System.out.println(
	       "Usage: DrawFlux model-parameter Log(threshold E[GeV]) nuCCenhanceFactor");
            System.exit(0);
        }else{
	    model = Integer.valueOf(args[0]).intValue();
	    logEth = Double.valueOf(args[1]).doubleValue();
	    nuCCEnhancementFactor = Double.valueOf(args[2]).doubleValue();
	}

	iLogEth = (int )((logEth-Particle.getLogEnergyMinimum( ))/Particle.getDeltaLogEnergy( ));
	if(iLogEth<0){
	    System.err.println("Log(Threshold Energy[GeV]) " + logEth +
			       "must be grater than " + Particle.getLogEnergyMinimum( ));
	    System.exit(0);
	}

	QuickPropagatingNeutrinoFlux leptonFlux = new QuickPropagatingNeutrinoFlux(model);
	ParticlePoint s = new ParticlePoint(0.0,0.0,0);
	leptonFlux.propagator = new NeutrinoQuickPropagator(s);

       // Angular Loop
       int itheta;
       index = 0;
       for(upDown=0;upDown<2;upDown++){
	   for(itheta=0;itheta<zenithAngle[upDown].length;itheta++){
	       double zenith = zenithAngle[upDown][itheta];

	       // propagate Neutrino
               s= new ParticlePoint(0.0,Math.toRadians(zenith),upDown);
	       leptonFlux.propagator.setParticlePoint(s);
	       leptonFlux.propagator.propagateNeutrinoToIceCubeDepth(zenith,
								     nuCCEnhancementFactor);


	       // Solid angle calculation
	       if(upDown == 1) zenith = 180.0-zenith;
	       double cosTheta = Math.cos(Math.toRadians(zenith));
	       System.err.println("Zenith " + zenith + " Costh " + cosTheta);


	       //Initialization
	       for(ip=0;ip<5;ip++){
		   integralFlux[ip][index] = 0.0;
	       }

	       // NuE
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[0][index] += 
		       leptonFlux.getDFNuEDLogE(iLogE)*Particle.getDeltaLogEnergy( );
	       }

	       // NuMu
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[1][index] += 
		       leptonFlux.getDFNuMuDLogE(iLogE)*Particle.getDeltaLogEnergy( );
	       }

	       // NuTau
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[2][index] += 
		       leptonFlux.getDFNuTauDLogE(iLogE)*Particle.getDeltaLogEnergy( );
	       }

	       // Mu
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[3][index] += 
		       leptonFlux.getDFMuDLogE(iLogE)*Particle.getDeltaLogEnergy( );
	       }

	       // Tau
	       for(iLogE=iLogEth;iLogE<dimension;iLogE++){
		   logE = Particle.getLogEnergyMinimum( ) + 
		       Particle.getDeltaLogEnergy( )*(double )iLogE;
		   integralFlux[4][index] += 
		       leptonFlux.getDFTauDLogE(iLogE)*Particle.getDeltaLogEnergy( );
	       }

	       index++;
	   }


       }



       System.out.println("titx Cos(theta)");
       System.out.println("tity Flux [cm^-2! sec^-1! sr^-1!])");
       System.out.println("scal -1.0 1.0 0.0 3.0E-19");


       for(ip=0;ip<5;ip++){
	   index = 0;
	   for(upDown=0;upDown<2;upDown++){
	       for(itheta=0;itheta<zenithAngle[upDown].length;itheta++){
		   double zenith = zenithAngle[upDown][itheta];
		   if(upDown == 1) zenith = 180.0-zenith;
		   double cosTheta = Math.cos(Math.toRadians(zenith));
		   // Drawing
		   double Flux = integralFlux[ip][index];
		   System.out.println("data " + cosTheta + " 0.0 " + Flux + " 0.0");
		   index++;
	       }
	   }


	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("cont");
       }
       System.out.println("endg");

	    
    }
}
