package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import java.io.*;

public class DrawEventFlux {

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix() + 
   (int )(( Particle.getLogEnergyMinimum()-InteractionsBase.getLogEnergyProducedMinimum())/Particle.getDeltaLogEnergy());


   public static void main(String[] args) throws IOException {

       String propMatrixfileName = null;
       String eventMatrixfileName = null;
       int model = 1;
 
       if(args.length!=3){
            System.out.println("Usage: DrawEventFlux PropagationMatrixFileName model-parameter EventMatrixName");
            System.exit(0);
        }else{
            propMatrixfileName = args[0];
	    model = Integer.valueOf(args[1]).intValue();
            eventMatrixfileName = args[2];
	}

        // Read the Event Matrix file
       DataInputStream in = new DataInputStream(new FileInputStream(eventMatrixfileName));
       EventFlux cascadeFlux = new EventFlux(model,in);
       in.close( );

        // Read the Propagation Matrix file
       in = new DataInputStream(new FileInputStream(propMatrixfileName));
       cascadeFlux.propLeptonFlux.readMatrix(in);
       in.close( );

        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
        System.out.println("scal 2.0 12.0 -15.0 -7.0");


	int iLogE;
	double logE;
	double[][] fluxTable = new double[2][dimension];

	// EMG cascades
        System.out.println("lncl 2");
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = cascadeFlux.getDFEmgCascadeDLogE(logE);
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		fluxTable[0][iLogE] = Flux;
	    }else{
		logEFlux = -15.0;
		fluxTable[0][iLogE] = 0.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");

	// Hadron cascades
        System.out.println("lncl 4");
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = cascadeFlux.getDFHadronCascadeDLogE(logE);
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
		fluxTable[1][iLogE] = Flux;
	    }else{
		logEFlux = -15.0;
		fluxTable[1][iLogE] = 0.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");


	// Total cascades
        System.out.println("lncl 1");
	for(iLogE=0;iLogE<dimension;iLogE++){
	    logE = InteractionsBase.getLogEnergyProducedMinimum() + 
                    Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double Flux = cascadeFlux.getDFTotalCascadeDLogE(logE);
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux/ln10)/ln10 + logE;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	}
        System.out.println("txcl 2");
        System.out.println("txtl 2.5 -8.0 Emg Cascades");
        System.out.println("txcl 4");
        System.out.println("txtl 2.5 -7.5 Hadron Cascades");

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
   }
}
