package numRecipes;

import java.util.*;

/** Utilities for Spectial Functions. Some of them
    are based on the Numerical Recipes Library */

public class SpecialFunctions implements Function {

    private static final int ITMAX = 100; // Max.allowed number of iterations.
    private static final double EPS = 3.0e-7; // Relative accuracy.
    private static final double FPMIN = 1.0e-30;
    // Number near the smallest representable floating-point number.

    private static final double[ ] GAMMA_Coefficient =
    {76.18009172947146,-86.50532032941677,
     24.01409824083091,-1.231739572450155,
     0.1208650973866179e-2,-0.5395239384953e-5};

    private static final double gaussTerm = Math.sqrt(2.0*Math.PI);

    public SpecialFunctions( ){

    }


    /** log(Gamma(a)) function  */
    public static double alogGamma(double xx) {
	double x;
	double y;
	double tmp;
        double ser=1.000000000190015;

	y = xx; x = xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*Math.log(tmp);
        for(int j=0;j<=5;j++) ser += GAMMA_Coefficient[j]/++y;
        return -tmp+Math.log(2.5066282746310005*ser/x);
    }

    /** gamma(a,x)/Gamma(a) from 0 to x */
    public static double incompleteGammaUptoX(double a, double x) {
	double gamLog = alogGamma(a);
	double sum = 1.0/a;
	double del = 1.0/a;
	double ap = a;
	int n;

	if(x<=0.0){    // x must be grater than or equal to 0.
	    return 0.0;
	}else{
	    for(n=1;n<=ITMAX;n++){
		ap += 1.0;
		del *= x/ap;
		sum += del;
		if(Math.abs(del) < Math.abs(sum)*EPS){
		    double gamma = sum*Math.exp(-x+a*Math.log(x)-gamLog);
		    return gamma;
		}
	    }

	    return Double.NaN;
	}
    }


    /** Gamma(a,x)/Gamma(a) from x to infinity*/
    public static double incompleteGammaDowntoX(double a, double x) {
	double gamLog = alogGamma(a);
	double checkoutIndex = 1.0/FPMIN;
	double b = x+1.0-a;
	double checkoutIndexD = 1.0/b;
	double h = checkoutIndexD;
	double aTerm;
	double del;
	int i;

	for(i=1;i<=ITMAX;i++){
	    aTerm = -(double )i*((double )i-a);
	    b += 2.0;
	    checkoutIndexD = aTerm*checkoutIndexD + b;
	    if(Math.abs(checkoutIndexD) < FPMIN) checkoutIndexD = FPMIN;
	    checkoutIndex = aTerm/checkoutIndex + b;
	    if(Math.abs(checkoutIndex) < FPMIN) checkoutIndex = FPMIN;
	    checkoutIndexD = 1.0/checkoutIndexD;
	    del = checkoutIndex*checkoutIndexD;
	    h *= del;
	    if(Math.abs(del-1.0) <EPS) break;
	}
	if(i>ITMAX) return Double.NaN;
	double gamma = Math.exp(-x+a*Math.log(x)-gamLog)*h;
	return gamma;
    }



    /** P(a x)= gamma(a,x)/Gamma(a) from 0 to x. Using imcompleteGamma methods
	described above saves iterations and CPU time. */
    public static double gammaP(double a, double x) {
	if(x<0.0 || a<=0.0) return Double.NaN;
	if(x<(a+1.0)) return incompleteGammaUptoX(a,x);
	else return 1.0-incompleteGammaDowntoX(a,x);
    }


    /** Q(a x)= Gamma(a,x)/Gamma(a) from 0 to x. Using imcompleteGamma methods
	described above saves iterations and CPU time. */
    public static double gammaQ(double a, double x) {
	if(x<0.0 || a<=0.0) return Double.NaN;
	if(x<(a+1.0)) return 1.0-incompleteGammaUptoX(a,x);
	else return incompleteGammaDowntoX(a,x);
    }


    /** Gaussian function */
    public static double gauss(double mu, double sigma, double x) {
	double chi = (x-mu)*(x-mu)/(2.0*sigma*sigma);
	double prob = 1.0/(gaussTerm*sigma)*Math.exp(-chi);
	return prob;
    }


    /** Poisson distribution function */
    public static double poisson(double mu, long m) {
	double p;
	double prob;
	double r_chi=0.0;
	long j;

	for(j=1;j<=m;j++) r_chi += Math.log((double )j);
	p= -mu+((double )m)*Math.log(mu)-r_chi;
	prob = Math.exp(p);
	return prob;
    }


    /** Binominal distribution  nCm mu^m(1-mu)^(n-m) */
    public static double binominal(double mu, long n, long m){
	if(m < 0 || n < m  || n < 0) return 0.0;
	double nChi = 0.0;
	for(long j=1;j<=n;j++) nChi += Math.log((double )j);
	double mChi = 0.0;
	for(long j=1;j<=m;j++) mChi += Math.log((double )j);
	double nmChi = 0.0;
	for(long j=1;j<=(n-m);j++) nmChi += Math.log((double )j);
	double p = Math.pow(mu, (double )m);
	double pBar = Math.pow(1.0-mu, (double )(n-m));
	double expProb = nChi-mChi-nmChi;
	double prob = Math.exp(expProb)*p*pBar;

	return prob;
    }

    /** 
	<pre>
	Method for interface <Function> .
	Interface the special functions given here
	to the utility routiner such as the Romberg
	Integration code that is desinged for a genereal
	function in form of Func(x).

	FunctionIndex    1     Gaussian gauss(mu, sigma, double x)
	                       parameters[0]=mu
	                       parameters[1]=sigma 
       </pre>
    */

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){

	double interfacedFunction;

	switch (functionIndex) {
	case 1 : 
	    double mu = parameters[0];
	    double sigma = parameters[1];
	    interfacedFunction = gauss(mu,sigma,x);
	    break;
	default:
	    interfacedFunction = 0.0;
	    System.err.println("Illegal parameters! Index" + functionIndex);
	    System.exit(0);
	}

	return interfacedFunction;

    }


    /** Integration of the Gaussian function.
	It uses the interface <Function> and
	the class <Integration> written in Integration.java .
        Integration blacket [lowerBound upperBound] */
    public static double integrateGauss(double mu, double sigma,
				    double lowerBound, double upperBound){
	int functionIndex =1; // Specifies gaussian( )
	double[] parameters;
	SpecialFunctions integratedFunctions = new SpecialFunctions( );
	double sum = 0.0;
	int iterationTimes =10;

	parameters = new double[2];

	parameters[0] = mu;
	parameters[1] = sigma;

	/*
	sum = Integration.iterateTrapzd(integratedFunctions, 
					functionIndex, parameters,
				lowerBound, upperBound, iterationTimes);
	*/

	sum = Integration.RombergIntegral(integratedFunctions, 
					  functionIndex, parameters,
					  lowerBound, upperBound);
	return(sum);
    }

}
