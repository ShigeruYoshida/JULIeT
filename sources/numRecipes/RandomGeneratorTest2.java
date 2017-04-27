package numRecipes;

import numRecipes.*;
import java.util.*;


/** Testing of the Gaussian probability with Monte Carlo method based on
    the RandomGenerator class */
public class RandomGeneratorTest2  {

    private static final long times = 1000000;


    public static void main(String[] args){

	RandomGenerator random1 = new RandomGenerator( );
	double sigma = 0.0;
	long pass = 0;
	double r;

	if(args.length!=1){
	    System.out.println("Usage: RandomDoubleTest2 sigma");
	    System.exit(0);
	}else{
	    sigma = Double.valueOf(args[0]).doubleValue();
	}


	double prob = SpecialFunctions.integrateGauss(0.0,1.0,sigma,100.0);
	pass = 0;
	for(long i=0;i<times;i++){ // Gaussian
	    r =  random1.GetGaussianDouble(0.0,1.0);
	    if(r>=sigma) pass++;
	}
	double probMC = (double )pass/(double )times;
        System.out.println("Monte Carlo prob " + probMC + " +- " + 
			   Math.sqrt((double)pass)/(double )times + 
			   " Integrated prob " + prob);

    }
}

