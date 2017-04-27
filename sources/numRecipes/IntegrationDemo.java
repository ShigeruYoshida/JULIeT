package numRecipes;

import numRecipes.*;

/** Demo program for the Integration.class */
public class IntegrationDemo {

    public static void main(String[] args){

        double a = 0.0;
	double x = 0.0;
	long m;

	if(args.length!=2){
	    System.out.println("Usage: IntegrationDemo a x");
	    System.exit(0);
	}else{
	    a = Double.valueOf(args[0]).doubleValue();
	    x = Double.valueOf(args[1]).doubleValue();
	}

	/** Gaussian G(a,sqrt(a),x) */
	System.out.println("Gauss( " + a + "," + x + ")=" + 
	   SpecialFunctions.gauss(a,Math.sqrt(a),x));

	/** Integration of Gaussian  */
	System.out.println("Gauss( " + a + "," + x + ")=" + 
	   SpecialFunctions.integrateGauss(a,Math.sqrt(a),-50.0,50.0));
    }

}

