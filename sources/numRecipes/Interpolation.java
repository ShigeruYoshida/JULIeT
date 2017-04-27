package numRecipes;

import java.util.*;

/** Utilities to interpolate data points (x_i, y_i) (i=0,...).
    A function in form of y=f(x) is obtained by a given data points. */

public class Interpolation {

    private static double delta = 0.0;



    /** Search the array index such that x is between xx[index]
	and xx[index+1]. xx must be monotonic, either increasing
	or decreasing. index =0 or index = n-1 is returned
	to indicate x is out of range.  n gives number of the elements
    i.e. xx[0....n-1].*/
    public static int searchIndex(double[] xx, double x, int n){

        int jm;
	int jl=0;
        int ju = n;
	int index;
        boolean ascnd;

        ascnd=(xx[n-1] >= xx[0]);
        while (ju-jl > 1) {
                jm=(ju+jl) >> 1;
                if ((x >= xx[jm]) == ascnd)
                        jl=jm;
                else
                        ju=jm;
        }
        if (x == xx[0]) index=0;
        else if(x == xx[n-1]) index=n-1;
        else index=jl;

	return index;
    }




    /** Given arrays xa[0....n-1] and ya[0....n-1], and given a value x,
	this method returns a value y by the polynominal interpolation.
	Number of data ponts n is qeual to the plynominal of degree. */
    public static double polynominalInterpolate(double[] xa, 
						double[] ya, double x){

	int ns = 0;
	int n = xa.length;
	double[] c = new double[n];
	double[] d = new double[n];
	System.arraycopy(ya,0,c,0,n);
	System.arraycopy(ya,0,d,0,n);

	double dif = Math.abs(x-xa[0]);
	double difNew;
	int closestIndex = 0;

	for(int i=0;i<n;i++){ // Here we find the index of the closest table.
	    difNew = Math.abs(x-xa[i]);
	    if(difNew< dif){
		closestIndex = i;
		dif = difNew;
	    }
	}

	double yOut = ya[closestIndex--];
	for(int m=1;m<n;m++){
	    for(int i=0;i<n-m;i++){
		double ho = xa[i]-x;
		double hp = xa[i+m]-x;
		double w = c[i+1]-d[i];
		double den = ho - hp;
		if(den == 0.0){
		    System.err.println("Error in polynominalInterpolate");
		    System.exit(0);
		}
		den = w/den;
		c[i]=ho*den;
		d[i]=hp*den;
	    }

	    delta = ((2.0*closestIndex < (n-m)) ? c[closestIndex+1] : 
		     d[closestIndex--]);
	    yOut += delta;
	}

	return yOut;
    }

    public static double getErrorInPolynominalInterpolate( ){
	return delta;
    }


    /** Interpolate with mTh pylinominal function. The data arrays
	xa[0.....n-1] and ya[0....n-1] must be larger than m, 
	polynominal of degree. */
    public static double mThPolynominalInterpolate(double[] xa, 
						   double[] ya, int n,
						   double x, int m){

	int index = Interpolation.searchIndex(xa, x, n);
	int kMax = ((index-(m-1)/2 > 1) ? index-(m-1)/2 : 1);
	int kMin = ((kMax < n+1-m) ? kMax : n+1-m);
	//System.err.println("Min = " + kMin + " kMax =" + kMax);

	double[] xpart = new double[m];
	double[] ypart = new double[m];

	System.arraycopy(xa,kMin-1,xpart,0,m);
	System.arraycopy(ya,kMin-1,ypart,0,m);


	double yOut = polynominalInterpolate(xpart,ypart,x);

	return yOut;
    }


}
