package numRecipes;

import numRecipes.*;


/** Display the pseudorandom numbers */
public class RandomGeneratorDemo {

    public static void main(String[] args){

	RandomGenerator random1 = new RandomGenerator( );
        double r;
	long m;
        int times = 0;

	if(args.length!=1){
	    System.out.println("Usage: RandomDoubleDemo repeat-times");
	}else{
	    times = Integer.valueOf(args[0]).intValue();
	}

	for(int i=0;i<times;i++){ // 1st trial
	    r =  random1.GetRandomDouble( );
	    System.out.println("Pseudorandom (" + i + "):" + r);
	}
        System.out.println("");

	for(int i=0;i<times;i++){ // 1st trial Gaussian
	    r =  random1.GetGaussianDouble(0.0,1.0);
	    System.out.println("Gaussian (" + i + "):" + r);
	}

	for(int i=0;i<times;i++){ // 1st trial Possionian
	    m =  random1.GetPoissonian(10.0);
	    System.out.println("Poisson (" + i + "):" + m);
	}

	/** Gamma Functions **/
	System.out.println("Gamma( " + times + ")=" + 
	   Math.exp(SpecialFunctions.alogGamma((double )times)));

    }

}

