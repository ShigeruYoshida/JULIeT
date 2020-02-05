package iceCube.uhe.neutrinoModel;

import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawCRFlux {

   public static void main(String[] args) throws IOException {

       double zMax = 2.0;
       double zConst = 2.0;
       double alpha = 2.5;
       double m = 2.0;
       double I_10PeV = 0.0;
       boolean plotDifferential = false;
       boolean assumeConstEvolution = false;
 
       if(args.length<3){
            System.out.println("Usage: DrawNeutrinoFluxFunction zMax powerIndex m (log differential? yes 1 no 0) (zconst)");
            System.exit(0);
       }else if (args.length>=3) {
	   zMax = Double.valueOf(args[0]).doubleValue();
	   alpha = Double.valueOf(args[1]).doubleValue();
	   m = Double.valueOf(args[2]).doubleValue();
       }if(args.length == 4){
	   zConst = Double.valueOf(args[3]).doubleValue();
	   System.err.println(" evolution costant above" + zConst);
	   assumeConstEvolution = true;
       }

       NeutrinoFluxFromSource neutFluxSource = null;
       neutFluxSource = new NeutrinoFluxFromSource();

       DataInputStream input = new DataInputStream(System.in); 
       BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
       String buffer; 

       System.err.print("maximum neutrino energy from source [GeV] ->"); 
       buffer   = d.readLine(); 
       double neutrinoEnergyMax = Double.valueOf(buffer).doubleValue();
       neutFluxSource.setMaximumNeutrinoEnergy(neutrinoEnergyMax);
       System.err.println("minimum target photon energy [GeV] " + neutFluxSource.getTargetPhotonEnergyMinimum()); 
       System.err.print("power law index of the target photon field spectrum ->"); 
       buffer   = d.readLine(); 
       double gamma = Double.valueOf(buffer).doubleValue();
       neutFluxSource.setPowerLawIndexOfRadiation(gamma);
       System.err.print("set the CR flux? (yes 1 no 0) ->"); 
       buffer   = d.readLine(); 
       if(Integer.valueOf(buffer).intValue()==1){
	   System.err.print("set the CR reference energy ->"); 
	   buffer   = d.readLine(); 
	   double energyRef = Double.valueOf(buffer).doubleValue();
	   neutFluxSource.setCREnergyReference(energyRef);
	   System.err.print("set the CR energy flux [GeV/cm2 sec sr] ->"); 
	   buffer   = d.readLine(); 
	   double eFlux = Double.valueOf(buffer).doubleValue();
	   if(!assumeConstEvolution){
	       neutFluxSource.setCRFluxAtReferenceEnergy(zMax, m, alpha, eFlux);
	   }else{
	       neutFluxSource.setCRFluxAtReferenceEnergy(zMax, zConst, m, alpha, eFlux);
	   }
       }
	   


       System.out.println("titx E[GeV]");
       System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
       System.out.println("scal 1.0e6 1.0e11 1.0e-10 3.0e-4");
       System.out.println("gwin 0.4 0.9 0.2 0.9");
       System.out.println("logx");
       System.out.println("logy");

       double logE = 6.0;
       while(logE<11.0){
	   double crEnergy = Math.pow(10.0,logE);
	   double EFlux = 0.0;
	   if(!assumeConstEvolution)
	       EFlux = neutFluxSource.getCRDFDE(zMax,m,alpha,crEnergy)*crEnergy*crEnergy;
	   else
	       EFlux = neutFluxSource.getCRDFDE(zMax,zConst,m,alpha,crEnergy)*crEnergy*crEnergy;

	       System.out.println("data " + crEnergy + " 0.0 " + EFlux + " 0.0");
	       logE += 0.1;
	   }

       System.out.println("join");
       System.out.println("disp");
       System.out.println("endg");
	    

   }
}
