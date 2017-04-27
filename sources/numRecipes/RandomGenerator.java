package numRecipes;

import java.io.*;
import java.util.*;

/** Utilities for Random generator in Numerical Calculation. 
    Some of them are based on the Numerical Recipes Library */

public class RandomGenerator implements Serializable{

    private RandomDouble generator = null;
    // private Random generatorUtil = null;

    private boolean haveNextNextGaussian = false;
    private double nextNextGaussian;


    /** default constructor : using the current time as a seed. */
    public RandomGenerator( ) {

	if(generator == null){ // Initialization
	    // generatorUtil = new Random(System.currentTimeMillis());
	    generator = new RandomDouble(System.currentTimeMillis());
	}

    }

    /** constructor with a given seed. */
    public RandomGenerator(long seed) {

	if(generator == null){ // Initialization
	    // generatorUtil = new Random(System.currentTimeMillis());
	    generator = new RandomDouble(seed);
	}

    }


    public double GetRandomDouble( ) { // Get the preudorandom distributed 
	                               // uniformly.
	double r = generator.nextDouble( );
	return r;
    }

    /*
    public void ChangeSeed(long seed) { 
	generatorUtil.setSeed(seed);
    }
    */

    public double GetGaussianDouble(double mu, double sigma ){ 
	                                // Get the preudorandom distributed 
	                                // following Gaussian.
	double r;
	if (haveNextNextGaussian) {
	    haveNextNextGaussian = false;
	    r =  nextNextGaussian;
	} else {
	    double v1, v2, s;
	    do { 
		v1 = 2 * generator.nextDouble() - 1;   // between -1.0 and 1.0
		v2 = 2 * generator.nextDouble() - 1;   // between -1.0 and 1.0
		s = v1 * v1 + v2 * v2;
	    } while (s >= 1 || s == 0);
	    double multiplier = Math.sqrt(-2 * Math.log(s)/s);
	    nextNextGaussian = v2 * multiplier;
	    haveNextNextGaussian = true;
	    r = v1 * multiplier;
	}
	double gaussValue = mu + r*sigma;

	return gaussValue;

    }


    public static double muOld = -1.0;
    public static double expMu;
    public static double alogMu;
    public static double sqMu;
    public static double GammaFactor;
    public long GetPoissonian(double mu){ // Get the preudorandom distributed
					  // following Poissonian.

	double estimatedNumber = -1.0;
	double y, t;

	if(mu != muOld){ // Setup the initial value.
	    muOld = mu;
	    alogMu = Math.log(mu);
	    sqMu = Math.sqrt(2.0*mu);
	    if(mu<12.0){
		expMu = Math.exp(-mu);
	    }else{
		expMu = mu*alogMu-SpecialFunctions.alogGamma(mu+1.0);
	    }
	}


	if(mu<12.0){ // Use direct method.
	    t = 1.0;
	    do{
		estimatedNumber += 1.0;
		t *= generator.nextDouble();
	    }while(t > expMu);

	}else{ // Use reject method.
	    do{
		do{
		    y = Math.tan(Math.PI*generator.nextDouble());
		    estimatedNumber = sqMu*y + mu;
		}while(estimatedNumber < 0.0);
		estimatedNumber = Math.floor(estimatedNumber);
		t = 0.9*(1.0+y*y)*Math.exp(estimatedNumber*alogMu-
 	        SpecialFunctions.alogGamma(estimatedNumber+1.0)-expMu);
	    }while(generator.nextDouble() > t);
	}

	return (long )estimatedNumber;
    }
}
