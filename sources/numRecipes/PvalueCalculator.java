package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

/** Display the SpecialFunctions */
public class PvalueCalculator {

    public static void main(String[] args) throws IOException{

        double bg = 0.0;
	double x = 0.0;
	long m;

	if(args.length!=2){
	    System.out.println("Usage: PvalueCalculator background observed#");
	    System.exit(1);
	}else{
	    bg = Double.valueOf(args[0]).doubleValue();
	    x = Double.valueOf(args[1]).doubleValue();
	}

	double prob = 0.0;
	double event = 0.0;
	// prpbability conistent with the bg hypothesis
	while(event<x){
	    prob += SpecialFunctions.poisson(bg,(long )event);
	    event += 1.0;
	}
	double pValue = 1.0-prob;

	System.out.format(" p-value (BG=%f observed=%d)=%e\n",bg,(long)x,pValue);

	// Feldman-Cousins approach
	FeldmanCousins.inclusive = false;
	double pValueFC = 1.0-FeldmanCousins.probabilityOfNobserved(0.0, bg, (long)x);
	System.out.format(" p-value FC (BG=%f observed=%d)=%e\n",bg,(long)x,pValueFC);

        DataInputStream input = new DataInputStream(System.in);  
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input));  
        String buffer; 
	System.err.print("calculate C.L.  with uncertain BG? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    System.err.print("relative plus BG error ->");  
	    buffer   = d.readLine(); 
	    double maxBGerr = Double.valueOf(buffer).doubleValue();
	    System.err.print("relative minus BG error ->");  
	    buffer   = d.readLine(); 
	    double minBGerr = Double.valueOf(buffer).doubleValue();
	    FeldmanCousins.setRelativeBackgroundUncertainty(minBGerr,maxBGerr);

	    pValueFC = 1.0-FeldmanCousins.probabilityOfNobserved(0.0, bg, (long)x);
	    System.out.format(" p-value with sys FC (BG=%f observed=%d)=%e\n",
			      bg,(long)x,pValueFC);
	}

	// calculate sigma
	double sigma=1.0;
	double sigmaUpperBound = 100.0;
	double dSigma=0.02;
	double gaussProb = 0.0;
	do{
	    sigma += 0.01;
	    gaussProb = SpecialFunctions.integrateGauss(0.0,1.0,sigma,sigmaUpperBound);
	}while(gaussProb>=pValueFC);
  	    
	System.out.format(" %f sigma equivalent\n",sigma);
    }

}

