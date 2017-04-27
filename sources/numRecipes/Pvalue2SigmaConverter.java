package numRecipes;

import numRecipes.*;

import java.io.*;
import java.util.*;

public class Pvalue2SigmaConverter {

    public static void main(String[] args) throws IOException{

        double pValue = 0.0;

	if(args.length!=1){
	    System.out.println("Usage: Pvalue2SigmaConverter p-value");
	    System.exit(1);
	}else{
	    pValue = Double.valueOf(args[0]).doubleValue();
	}

	// calculate sigma
	double sigma=1.0;
	double sigmaUpperBound = 100.0;
	double gaussProb = 0.0;
	do{
	    sigma += 0.01;
	    gaussProb = SpecialFunctions.integrateGauss(0.0,1.0,sigma,sigmaUpperBound);
	}while(gaussProb>=pValue);
  	    
	System.out.format("p-value(%e) is  %f sigma equivalent\n",pValue,sigma);
    }

}

