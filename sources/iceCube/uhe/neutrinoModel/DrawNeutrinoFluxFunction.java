package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawNeutrinoFluxFunction {

   public static void main(String[] args) throws IOException {

       double zMax = 2.0;
       double zConst = 2.0;
       double alpha = 2.5;
       double m = 2.0;
       double I_10PeV = 0.0;
       boolean plotDifferential = false;
       boolean assumeConstEvolution = false;
       boolean addFromSource = false;
 
       if(args.length<3){
            System.out.println("Usage: DrawNeutrinoFluxFunction zMax powerIndex m (log differential? yes 1 no 0) (zconst)");
            System.exit(0);
       }else if (args.length>=3) {
	   zMax = Double.valueOf(args[0]).doubleValue();
	   alpha = Double.valueOf(args[1]).doubleValue();
	   m = Double.valueOf(args[2]).doubleValue();
       }if(args.length >= 4){
	   if(Integer.valueOf(args[3]).intValue() == 1){
	       plotDifferential = true;
	       System.err.println(" plot differential spectrum");
	   }
       }if(args.length == 5){
	   zConst = Double.valueOf(args[4]).doubleValue();
	   System.err.println(" evolution costant above" + zConst);
	   assumeConstEvolution = true;
       }

       NeutrinoFluxFunction neutFlux = new NeutrinoFluxFunction();

       // set up the parameters
       DataInputStream input = new DataInputStream(System.in); 
       BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
       String buffer; 
       System.err.print("add from source ? (yes 0) ->"); 
       buffer   = d.readLine(); 
       if(Integer.valueOf(buffer).intValue()==0) addFromSource = true;

       NeutrinoFluxFromSource neutFluxSource = null;
       if(addFromSource){
	   neutFluxSource = new NeutrinoFluxFromSource();
	   System.err.print("maximum neutrino energy from source [GeV] ->"); 
	   buffer   = d.readLine(); 
	   double neutrinoEnergyMax = Double.valueOf(buffer).doubleValue();
	   neutFluxSource.setMaximumNeutrinoEnergy(neutrinoEnergyMax);
	   System.err.println("minimum target photon energy [GeV] " + neutFluxSource.getTargetPhotonEnergyMinimum()); 
	   System.err.print("power law index of the target photon field spectrum ->"); 
	   buffer   = d.readLine(); 
	   double gamma = Double.valueOf(buffer).doubleValue();
	   neutFluxSource.setPowerLawIndexOfRadiation(gamma);
       }

       if(!plotDifferential){
	   System.out.println("titx Log E[GeV]");
	   System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
	   if(!addFromSource) System.out.println("scal 6.0 12.0 -14.0 -4.0");
	   else System.out.println("scal 4.0 12.0 -14.0 -4.0");

	   double logE = 6.0;
	   if(addFromSource) logE = 4.0;
	   while(logE<12.0){
	       double neutrinoEnergy = Math.pow(10.0,logE);
	       double EFlux = 0.0;
	       if(!assumeConstEvolution)
		   EFlux = neutFlux.getDFDE(zMax,m,alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	       else
		   EFlux = neutFlux.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;

	       if(addFromSource){
		   if(!assumeConstEvolution) EFlux += neutFluxSource.getDFDE(zMax,m,alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
		   else EFlux += neutFluxSource.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy)*neutrinoEnergy*neutrinoEnergy;
	       }

	       double logEFlux;
	       if(EFlux > 0.0){
		   logEFlux = Math.log(EFlux)/Math.log(10.0);
	       }else{
		   logEFlux = -14.0;
	       }
	       System.out.println("data " + logE + " 0.0 " + logEFlux + " 0.0");
	       logE += 0.1;
	   }

	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("endg");
	    
       }else{
	   System.out.println("titx Energy [GeV]");
	   System.out.println("tity dF/dLogE [cm^-2! sec^-1! sr^-1!])");
	   if(!addFromSource) System.out.println("scal 1.0e8 1.0e11 1.0e-20 1.0e-15");
	   else System.out.println("scal 1.0e4 1.0e11 1.0e-20 1.0e-15");
	   System.out.println("gwin 0.2 0.9 0.2 0.9");
	   System.out.println("logx");
	   System.out.println("logy");

	   double logE = 8.0;
	   if(addFromSource) logE = 4.0;
	   while(logE<11.0){
	       double neutrinoEnergy = Math.pow(10.0,logE);
	       double flux = 0.0;
	       if(!assumeConstEvolution)
		   flux = neutFlux.getDFDE(zMax,m,alpha,neutrinoEnergy)*neutrinoEnergy*Math.log(10.0);
	       else
		   flux = neutFlux.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy)*neutrinoEnergy*Math.log(10.0);

	       if(addFromSource){
		   if(!assumeConstEvolution) flux += neutFluxSource.getDFDE(zMax,m,alpha,neutrinoEnergy)*neutrinoEnergy*Math.log(10.0);
		   else flux += neutFluxSource.getDFDE(zMax,zConst,m,alpha,neutrinoEnergy)*neutrinoEnergy*Math.log(10.0);
	       }

	       if(flux>0.0){
		   System.out.println("data " + neutrinoEnergy + " 0.0 " + flux + " 0.0");
	       }
	       logE += 0.1;
	   }

	   System.out.println("join");
	   System.out.println("disp");
	   System.out.println("endg");
       }
   }
}
