package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.particles.*;
import java.io.*;

public class SearchEvolution {

   public static void main(String[] args) throws IOException {

       double zMax = 2.0;
       double zConst = 1.0;
       double alpha = 2.5;
       double m = 2.0;
       double fluxRatio=1.126; // 20% more to the original flux
       int functionIndex = 0;
       NeutrinoFluxFunction neutFluxFunction = null;
 
       if(args.length<2 ){
            System.out.println("Usage: DrawNeutrinoFluxFunction m zMax (zConst)");
            System.exit(0);
       }else{
	   m = Double.valueOf(args[0]).doubleValue();
	   zMax = Double.valueOf(args[1]).doubleValue();
	   if(args.length == 3){
	       zConst = Double.valueOf(args[2]).doubleValue();
	       functionIndex = 1;
	       System.err.format(" evolution constant above %f\n",zConst);
	   }
       }
       // set up the parameters
       double[] parameters = new double[5];
       neutFluxFunction = new NeutrinoFluxFunction();
       if(functionIndex == 0){
	   parameters[0] = alpha;
	   parameters[1] = zMax;
	   parameters[2] = m;
       }else{
	   parameters[0] = alpha;
	   parameters[1] = zMax;
	   parameters[2] = zConst;
	   parameters[3] = m;
       }

       double neutrinoEnergy = 1.0e7;     // 10^7 GeV
       double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	
       double integralFlux = Integration.RombergIntegral(neutFluxFunction, 
				       0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
       System.out.format("Integral flux (m Zmax) (%f %f)=%e\n",m,zMax,integralFlux);

       //
       // Now search for the evolution to give the flux exceeding the original by fluxRatio;
       //
       double integralFluxUpperLimit = integralFlux*fluxRatio;
       double delta_m = 1.0e-2;
       double m_search = m ;
       double m_max = 10.0; // set the boundary for fail safe
       do{
	   m_search += delta_m;
	   if(m_search >= m_max) break;

	   if(functionIndex == 0) parameters[2] = m_search;
	   else parameters[3] = m_search;

	   integralFlux = Integration.RombergIntegral(neutFluxFunction, 
			    0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
       }while(integralFlux<integralFluxUpperLimit);

       System.out.format("m = %e\n",m_search);

	    
   }
}
