package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import java.io.*;

public class DrawEventArrayForEachInteraction{

    private static final int mc  = EventArrayForEachInteraction.mc; // length of mcBases[]
    private static final double ln10 = Math.log(10.0);
    private static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    
   public static void main(String[] args) throws IOException {

       String fileName = null;
       int model = 1;
       double sum;
       double logE;
       double logEprimary=7;


       if(args.length!=2){
            System.out.println("Usage: DrawEventArrayForEachInteraction MatrixFileName primary-logE");
            System.exit(0);
        }else{
            fileName = args[0];
	    logEprimary = Double.valueOf(args[1]).doubleValue();
	}

        // Read the serialized object of the Neutrino Charged Interaction Matrix
       DataInputStream in = new DataInputStream(new FileInputStream(fileName));
       EventArrayForEachInteraction eventArray = new EventArrayForEachInteraction(in);
       in.close( );

       System.out.println("titx Log E-cascade [GeV]");
       System.out.println("tity Log dF/dLogE");
       System.out.println("scal 2.0 12.0 -10.0 -1.0");

       for(int n=0; n<mc; n++) {
	   logE = iceCube.uhe.interactions.InteractionsBase.getLogEnergyProducedMinimum();
	   sum = 0.0;
	   System.out.println("lncl ");
	   System.out.println("lnpt ");
	   while(logE<logEprimary){
	       double Flux = eventArray.getCascadeFlux(n, logE);
	       sum += Flux;
	       double logEFlux;
	       if(Flux > 0.0){
		   logEFlux = Math.log(Flux)/ln10;
	       }else{
		   logEFlux = -10.0;
	       }
	       System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	       logE += Particle.getDeltaLogEnergy();
	   }
	   
	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("cont");
       }


        System.out.println("txcl 1");
        System.out.println("tfon 18");
	System.out.println("endg");

	    
   }
}




