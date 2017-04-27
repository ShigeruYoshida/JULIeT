package numRecipes;

import numRecipes.*;

/** Display the SpecialFunctions */
public class SpecialFunctionsDemo {

    public static void main(String[] args){

        double a = 0.0;
	double x = 0.0;
	long m;

	if(args.length!=2){
	    System.out.println("Usage: SpecialFunctionsDemo a x");
	}else{
	    a = Double.valueOf(args[0]).doubleValue();
	    x = Double.valueOf(args[1]).doubleValue();
	}

	/** Gamma Functions **/
	System.out.println("Gamma( " + a + ")=" + 
	   Math.exp(SpecialFunctions.alogGamma(a)));

	/** IncompleteGamma P(a,x) */
	System.out.println("GammaP( " + a + "," + x + ")=" + 
	   SpecialFunctions.gammaP(a,x));

	/** IncompleteGamma Q(a,x) */
	System.out.println("GammaQ( " + a + "," + x + ")=" + 
	   SpecialFunctions.gammaQ(a,x));

	/** Gaussian G(a,sqrt(a),x) */
	System.out.println("Gauss( " + a + "," + x + ")=" + 
	   SpecialFunctions.gauss(a,Math.sqrt(a),x));

	/** Poisson P(a,x) */
	System.out.println("GammaQ( " + a + "," + x + ")=" + 
	   SpecialFunctions.poisson(a,(long )x));
    }

}

