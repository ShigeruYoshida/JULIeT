package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import java.io.*;

public class DrawEventMatrix_single{

   private static final double ln10 = Math.log(10.0);
   private static int dimension = Particle.getDimensionOfLogEnergyMatrix();

   public static void main(String[] args) throws IOException {

       String fileName = null;
       double sum;
       double logE;

       double logEprimary = 0.0;
 
       if(args.length!=2){
            System.out.println("Usage: DrawNeutrinoFlux MatrixFileName PrimaryEnergy");
            System.exit(0);
        }else{
            fileName = args[0];
	    logEprimary = Double.valueOf(args[1]).doubleValue();
	}

       // Read the data file
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       EventMatrix_single eventMtx = new EventMatrix_single(in);
       in.close( );

        System.out.println("titx Log E-cascade [GeV]");
        System.out.println("tity Log dF/dLogE");
        System.out.println("scal 2.0 12.0 -4.0 -1.0");


	// Emg Cascade
	logE = iceCube.uhe.interactions.InteractionsBase.getLogEnergyProducedMinimum();
	sum = 0.0;
        System.out.println("lncl 2");
	while(logE<logEprimary){
	    double Flux = eventMtx.getEmgCascadeFlux(logEprimary, logE);
	    sum += Flux;
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux)/ln10;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	    logE += Particle.getDeltaLogEnergy();
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");
	System.err.println("Sum of the Emg cascade matrix " + sum);

	// Hadron Cascade
	logE = iceCube.uhe.interactions.InteractionsBase.getLogEnergyProducedMinimum();
	sum = 0.0;
        System.out.println("lncl 4");
	while(logE<logEprimary){
	    double Flux = eventMtx.getHadronCascadeFlux(logEprimary, logE);
	    sum += Flux;
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux)/ln10;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	    logE += Particle.getDeltaLogEnergy();
	}
	System.err.println("Sum of the Hadron cascade matrix " + sum);

        System.out.println("txcl 1");
        System.out.println("tfon 18");
        System.out.println("txtl 2.5 -1.3 "  + Particle.particleName(eventMtx.getPropFlavor(),1));
        System.out.println("txcl 2");
        System.out.println("txtl 2.5 -1.6 Emg Cascades");
        System.out.println("txcl 4");
        System.out.println("txtl 2.5 -1.9 Hadron Cascades");

	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");


	// Total Cascade
	logE = iceCube.uhe.interactions.InteractionsBase.getLogEnergyProducedMinimum();
	sum = 0.0;
        System.out.println("lncl 1");
	double E = 0.0, sumE = 0.0;
	double E2 = 0.0, sumE2 = 0.0;
	double variance = 0.0;
	double energy = 0.0;
	double sumFlux = 0.0;
	while(logE<logEprimary){
	    double Flux = eventMtx.getTotalCascadeFlux(logEprimary, logE);
	    sum += Flux;
	    double logEFlux;
	    if(Flux > 0.0){
		logEFlux = Math.log(Flux)/ln10;
		sumFlux += Flux;
		E = Flux*logE;       E2 = Flux*logE*logE;
		sumE += E;           sumE2 += E2;
	    }else{
		logEFlux = -15.0;
	    }
	    System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	    logE += Particle.getDeltaLogEnergy();
	}

	System.err.println("Sum of the total cascade matrix " + sum);
	double mean = sumE/sum;  double mean2 = sumE2/sum;
	variance = (mean2) - (mean*mean);
	double sig = Math.sqrt(variance);
	System.err.println("Primary logE: " + logEprimary);
	System.err.println("data "+ logEprimary + " 0.0 " + (mean-sig) + " 0.0 "); 
	System.err.println("data "+ logEprimary + " 0.0 " + (mean+sig) + " 0.0 "); 
	System.err.println("data "+ logEprimary + " 0.0 " + (mean) + " 0.0 "); 
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");
	System.out.println("lncl 5");
	System.out.println("data "+ mean + " " + sig + " -15.0 " + " 0.0 "); 
	System.out.println("data "+ mean + " " + sig + " -1.0 " + " 0.0 "); 
	System.out.println("join");
	System.out.println("disp");
	System.out.println("cont");
	System.out.println("mksz 1");
	System.out.println("mkcl 5");
	System.out.println("mkpt 2");
	System.out.println("data "+ (mean-sig) + " 0.0 " + " -3.0 " + " 0.0 ");
	System.out.println("data "+ (mean+sig) + " 0.0 " + " -3.0 " + " 0.0 ");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

	    
   }
}

