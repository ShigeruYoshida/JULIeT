package numRecipes;

import java.util.*;

/** Utilities for Numerical Integration. Some of them
    are based on the Numerical Recipes Library. 
    The interface class Function is required to call
    the member in this class */

public class Integration {

    private static final int ITMAX = 100; // Max.allowed number of iterations.
    private static final int maxStep =100; // Max.allowed step number for
                                          // the Romberg Integration.
    private static final int polynominalOfDegree = 5;
    private static double EPS = 5.0e-5;   // Relative accuracy.
    private static final int functionMax = 20;
    private static double[] grandSum = new double[functionMax];

    public Integration( ){ };


    /** Trapzd integraion method. This is the fundamenthal method
	for more advanced integration.
        Integration blacket [lowerBound upperBound] */
    public static double trapzd(Function func, int functionIndex, 
				double[] parameters,
				double lowerBound, double upperBound, 
				int iterationTimes) {

	double sum = 0.0;
	int it = 1;
	double x;
	double tnm;
	double del;

	if(iterationTimes == 1){
	    grandSum[functionIndex] = 0.5*(upperBound-lowerBound)*
		(func.getFunction(functionIndex, parameters, lowerBound)+
		 func.getFunction(functionIndex, parameters, upperBound));
	}else{
	    for(int j=1;j<iterationTimes-1;j++) it <<= 1;
	    tnm = (double )it;
	    del = (upperBound - lowerBound)/tnm;
	    x = lowerBound + 0.5*del;
	    double[] funcValuesToSum = new double[it];
	    for(int j=1;j<=it;j++){
			funcValuesToSum[j-1] = func.getFunction(functionIndex, parameters, x);
			x += del;
	    }
    	sum = KahanSum(funcValuesToSum);
	    grandSum[functionIndex] = 0.5*(grandSum[functionIndex]+del*sum);
	}


	return grandSum[functionIndex];
    }

    public static double KahanSum(double[] valuesToSum){
    	double sum = 0.0;

    	// Summation error
    	double c = 0.0;
    	for(int i=0; i<valuesToSum.length; i++){
    		double val = valuesToSum[i];
    		double y = val - c;
    		double t = sum + y;

    		c = (t - sum) - y;

    		sum = t;
    	}
    	return sum;
    }
    /** Iterated Trapzd integraion method. This is the fundamenthal method
	for more advanced integration.
        Integration blacket [lowerBound upperBound] */
    public static double iterateTrapzd(Function func, int functionIndex, 
				double[] parameters,
				double lowerBound, double upperBound, 
				int iterationTimes) {
	double sum = 0.0;

	if(iterationTimes<1||ITMAX<iterationTimes){
	    System.err.println("iteration times out of bound");
	    System.exit(0);
	}


	for(int iteration = 1; iteration<iterationTimes;iteration++){
	    sum =
	    Integration.trapzd(func, functionIndex, parameters,
			       lowerBound, upperBound, iteration);

	    //    System.out.println("func( " + iterationTimes + ")=" + sum);
	   
	}
	return(sum);
    }



    /** The Romberg method. This is the fundamenthal method
	for more advanced integration./
        Integration blacket [lowerBound upperBound] */
    public static double RombergIntegral(Function func, int functionIndex, 
			 double[] parameters,
			 double lowerBound, double upperBound){

	double[] h = new double[maxStep+2];
	double[] s = new double[maxStep+2];
	h[0] = 1.0;
	double[] hpart = new double[polynominalOfDegree];
	double[] spart = new double[polynominalOfDegree];

	for(int iterationTimes = 1; iterationTimes<=maxStep; iterationTimes++){
	    s[iterationTimes-1] = 
	    Integration.trapzd(func, functionIndex, parameters,
			       lowerBound, upperBound, iterationTimes);
	    //System.out.println("Iteration Times " + iterationTimes + 
	    //	    	       " s= " + s[iterationTimes-1]);
	    if(iterationTimes>=polynominalOfDegree){


		System.arraycopy(h,iterationTimes-polynominalOfDegree,
				 hpart,0,polynominalOfDegree);
		System.arraycopy(s,iterationTimes-polynominalOfDegree,
				 spart,0,polynominalOfDegree);

		double sum = Interpolation.polynominalInterpolate(hpart,spart,0.0);
		double delta = Interpolation.getErrorInPolynominalInterpolate( );
//  	        System.out.println("Iteration Times" + iterationTimes + 
//	   " sum= " +  sum + " delta = " + delta);

		if(Math.abs(delta)<=EPS*Math.abs(sum)){
		    //System.err.println("Iteration Times " + iterationTimes);
		    return sum;
		}
	    }
	    s[iterationTimes]=s[iterationTimes-1];
	    h[iterationTimes]=0.25*h[iterationTimes-1];
	}

	return 0.0;
    }

    /** Initialize the static variables such as grandSum.
	This may have to be called when different classes
	utilize this class at same time.*/
    public static void initializeStaticVariables( ){
	for(int i=0;i<functionMax;i++)	grandSum[i] = 0.0;
    }



    /** Change the relative accuracy in the inegration */
    public static void setRelativeAccuracy(double epsilon){
	EPS = epsilon;
    }

}
