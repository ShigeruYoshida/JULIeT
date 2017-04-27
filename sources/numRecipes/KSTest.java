package numRecipes;

import numRecipes.*;
import java.util.*;

/** Utilities for Kolmogorov-Smirnov Test. 
    Based on the Numerical Recipes Library */

public class KSTest {

    private static final double EPS1 = 0.001; // Relative accuracy.
    private static final double EPS2 = 1.0e-8; // Relative accuracy.
    // Number near the smallest representable floating-point number.

    public KSTest( ){

    }

    /**
       calculate the KS significance for a given data array and a user-supplied func describing probability 
       by a model prediction in the range [lowerBound upperBound]. 
       Its normalization is taken care internally so that you do not have to worry it.
       When isCumulative = true, then func must be a cummulative probability distribution.
       isCumulative = false, then this method cauclates the cumullative distribution internally.
     */
    public static double getKSSignificance(double[] data, Function func, int functionIndex, double[] parameters,
					   double lowerBound, double upperBound, 
					   boolean isCumulative){
	// sorting the order
	Arrays.sort(data);
	// normaization constant
	double normalization = Integration.RombergIntegral(func, functionIndex, parameters,lowerBound, upperBound);

	// now making comparison
	double numberOfData = (double )(data.length);
	double ksD = 0.0;  // KS test static D
	for(int i=0; i< data.length; i++){
	    double sNx1st = (double )(i)/numberOfData;
	    double sNx2nd = (double )(i+1)/numberOfData;
	    double px = func.getFunction(functionIndex,parameters,data[i])/normalization;
	    if(!isCumulative){ // this is not cumulative - we calculate its cumulative distriution now
		px = Integration.RombergIntegral(func, functionIndex, parameters,lowerBound, data[i])/normalization;
	    }
	    double d = Math.abs(px-sNx1st);
	    double d2nd = Math.abs(px-sNx2nd);
	    if(d2nd>d) d = d2nd;

	    if(d>ksD) ksD = d;
	}

	double sqrtN = Math.sqrt(numberOfData);
	double lambda = sqrtN + 0.12 + 0.11/sqrtN;
	double prob = probabilityFunction(lambda*ksD);

	return prob;

    }

    /**
       Return the Kolmogorov-Smirnov probability Q_ks(lambda)
     */
    public static double probabilityFunction(double lambda){

	double expLambda = -2.0*lambda*lambda;
	double fac = 2.0;

	double sum = 0.0; double termbf = 0.0;
	for(int j=1;j<=100;j++){
	    double dJ = (double )j;
	    double term = fac*Math.exp(dJ*dJ*expLambda);
	    sum += term;
	    if(Math.abs(term) <= EPS1*termbf || Math.abs(term) <= EPS2*sum) return sum;
	    fac = -fac;  // flip signs in sum
	    termbf = Math.abs(term);
	}
	return (1.0);  // get here only by failing to converge
    }


}
