package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class DrawIntegralNeutrinoFlux {

   public static void main(String[] args) throws IOException {

       double zMax = 2.0;
       double zConst = 1.0;
       double alpha = 2.5;
       double m = 2.0;
       double opticalDepth = 0.1;
       int modelID = 1;
       boolean drawAnalyticalFunction = false;
       boolean assumeConstEvolution = false;
       boolean addFromSource = false;
       NeutrinoFluxFunction neutFluxFunction = null;
       //NeutrinoFluxFunctionEnergetics neutFluxFunction = null;
       NeutrinoFluxFromSource neutFluxSource = null;
       NeutrinoFlux neutFlux = null;
 
       if(args.length!=3 && args.length!=4 && args.length!=1){
            System.out.println("Usage: DrawNeutrinoFluxFunction zMax powerIndex m");
            System.out.println("Or Usage: DrawNeutrinoFluxFunction zMax Zconst powerIndex m");
            System.out.println("or Usage: DrawNeutrinoFluxFunction modelID");
            System.exit(0);
       }else if (args.length == 3){
	   zMax = Double.valueOf(args[0]).doubleValue();
	   alpha = Double.valueOf(args[1]).doubleValue();
	   m = Double.valueOf(args[2]).doubleValue();
	   drawAnalyticalFunction = true;
       }else if (args.length == 4){
	   zMax = Double.valueOf(args[0]).doubleValue();
	   zConst = Double.valueOf(args[1]).doubleValue();
	   alpha = Double.valueOf(args[2]).doubleValue();
	   m = Double.valueOf(args[3]).doubleValue();
	   drawAnalyticalFunction = true;
	   assumeConstEvolution = true;
       }else if (args.length == 1){
	   modelID = Integer.valueOf(args[0]).intValue();
       }

       // set up the parameters
       double[] parameters = new double[5];
       if(drawAnalyticalFunction){
	   neutFluxFunction = new NeutrinoFluxFunction();
	   DataInputStream input = new DataInputStream(System.in); 
	   BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
	   String buffer; 
	   System.err.print("add from source ? (yes 0) ->"); 
	   buffer   = d.readLine(); 
	   if(Integer.valueOf(buffer).intValue()==0) addFromSource = true;
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
	       System.err.print("optical depth ->"); 
	       buffer   = d.readLine(); 
	       opticalDepth = Double.valueOf(buffer).doubleValue();
	       neutFluxSource.setOpticalDepth(opticalDepth);
	   }
	   parameters[0] = alpha;
	   parameters[1] = zMax;
	   if(!assumeConstEvolution) parameters[2] = m;
	   else{
	       parameters[2]= zConst;
	       parameters[3] = m;
	   }
       }else{
	   neutFlux = new NeutrinoFlux(modelID);
       }

       System.out.println("titx Energy [GeV]");
       System.out.println("tity Integral [GeV cm^-2 sec^-1 sr^-1])");
       if(!addFromSource) System.out.println("scal 1.0e7 1.0e10 1.0e-20 1.0e-15");
       else System.out.println("scal 1.0e5 1.0e10 1.0e-20 1.0e-13");
       System.out.println("gwin 0.4 0.9 0.2 0.9");

	double logE = 7.0;
	if(addFromSource) logE = 5.0;
	double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	double maxLogNeutrinoEnergy = 11.0;   // 1.0e11 GeV
	while(logE<10.0){
	    double neutrinoEnergy = Math.pow(10.0,logE);
	    double integralFlux = 0.0;
	    if(drawAnalyticalFunction){
		if(!assumeConstEvolution)
		    integralFlux = Integration.RombergIntegral(neutFluxFunction, 
				       0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
		else
		    integralFlux = Integration.RombergIntegral(neutFluxFunction, 
				       1, parameters, neutrinoEnergy, maxNeutrinoEnergy);
		if(addFromSource){
		    double maxNeutrinoEnergyFromSource = neutFluxSource.neutrinoEnergyMax;
		    if(neutrinoEnergy<maxNeutrinoEnergyFromSource){
			if(!assumeConstEvolution)
			    integralFlux += Integration.RombergIntegral(neutFluxSource, 
					       0, parameters, neutrinoEnergy, maxNeutrinoEnergyFromSource);
			else
			    integralFlux += Integration.RombergIntegral(neutFluxSource, 
					       1, parameters, neutrinoEnergy, maxNeutrinoEnergyFromSource);
		    }
		}
	    }else{
		for(int pID = 1; pID<=3; pID++){
		    parameters[0] = (double )pID;
		    integralFlux += Integration.RombergIntegral(neutFlux, 
							  1, parameters, logE, maxLogNeutrinoEnergy);
		}
	    }
	    System.out.println("data " + neutrinoEnergy + " 0.0 " + integralFlux + " 0.0");
	    logE += 0.1;
	}

	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
	    
   }
}
