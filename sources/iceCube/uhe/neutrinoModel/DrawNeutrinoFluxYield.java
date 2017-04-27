package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawNeutrinoFluxYield {

   public static void main(String[] args) throws IOException {

       double zSource = 2.0;
       double alpha = 2.5;
       double neutrinoEnergy = 1.0e9;
       boolean plotDifferential = false;
 
       if(args.length<3){
            System.out.println("Usage: DrawNeutrinoFluxFunction zSource powerIndex  neutrinoEnergy");
            System.exit(0);
       }else if (args.length>=3) {
	   zSource = Double.valueOf(args[0]).doubleValue();
	   alpha = Double.valueOf(args[1]).doubleValue();
	   neutrinoEnergy = Double.valueOf(args[2]).doubleValue();
	   System.err.println(" neutrinoEnergy" + neutrinoEnergy);
       }

       NeutrinoFluxFunction neutFlux = new NeutrinoFluxFunction();

       System.out.println("titx redshift to yield neutrinos");
       System.out.println("tity dF/dLogE [cm^-2! sec^-1! sr^-1!])");
       //System.out.println("scal 0.0 2.0 1.0e-25 1.0e-16");
       System.out.println("scal 1.0e-3 3.0 1.0e-25 1.0e-16");
       System.out.println("gwin 0.2 0.9 0.2 0.9");
       System.out.println("logx");
       System.out.println("logy");

       double z = zSource;
       double zMin = zSource-0.1;
       if(zMin<0.0) zMin = 0.0;
       while(z>=zMin){
	   double flux = neutFlux.getNeutrinoYieldFromAsource(zSource,z,alpha,neutrinoEnergy)*neutrinoEnergy*Math.log(10.0);
	   if(flux>0.0){
	       System.out.println("data " + z + " 0.0 " + flux + " 0.0");
	   }
	   z -= 0.001;
       }

       System.out.println("join");
       System.out.println("disp");
       System.out.println("endg");
   }
}
