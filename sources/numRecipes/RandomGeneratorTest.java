package numRecipes;

import numRecipes.*;
import java.util.*;


/** Testing of the uniformity in the pseudorandom numbers */
public class RandomGeneratorTest {

    private static final double scaleFactor = 1.0e4;
    private static final double gaussFactor = 5.0e2;
    private static final int dim = 10000;
    private static int index;
    private static long times = 0;


    public static void main(String[] args){

	RandomGenerator random1 = new RandomGenerator( );
        double r;
	long[ ] randomBin;
	long[ ] gaussBin;
	long[ ] poissonBin;

	if(args.length!=1){
	    System.out.println("Usage: RandomDoubleTest repeat-times");
	    System.exit(0);
	}else{
	    times = Long.valueOf(args[0]).longValue();
	}
        System.err.println("Repeat times (" + times + ")");

	//Generate the data array
	randomBin = new long[dim];
	gaussBin = new long[dim];
	poissonBin = new long[dim];
	for(int i=0;i<dim;i++){ // initialization
	    randomBin[i]=0;
	    gaussBin[i]=0;
	    poissonBin[i]=0;
	}

        System.err.println("Now generating Pseudorandom....");
	for(long i=0;i<times;i++){ // Pseudorandom
	    r =  scaleFactor*random1.GetRandomDouble( );
	    index = (int )Math.floor(r);
	    indexCheck(); //Check the range of index
	    randomBin[index]++;
	}
        System.out.println("done");

        System.err.println("Now generating Gaussian....");
	for(long i=0;i<times;i++){ // Gaussian
	    r =  gaussFactor*random1.GetGaussianDouble(0.0,1.0);
	    index =  (int )Math.floor(0.5*scaleFactor+r);
	    indexCheck(); //Check the range of index
	    gaussBin[index]++;
	}
        System.out.println("done");

        System.err.println("Now generating Poisson....");
	for(long i=0;i<times;i++){ // 1st trial Possionian
	    index =  (int )random1.GetPoissonian(10.0);
	    indexCheck(); //Check the range of index
	    poissonBin[index]++;
	}
        System.out.println("done");



	//Writing the data
	System.out.println("zone 2 2");
	System.out.println("tith Pseudorandom");
	for(int i=0;i<dim;i++){ // Pseudorandom
	    System.out.println("data " + i + " 0.0 " + randomBin[i] + " 0.0");
	}
	System.out.println("hist");
	System.out.println("disp");
	System.out.println("endg");

	System.out.println("tith Gaussian");
	for(int i=0;i<dim;i++){ // Gaussian
	    System.out.println("data " + i + " 0.0 " + gaussBin[i] + " 0.0");
	}
	System.out.println("hist");
	System.out.println("disp");
	System.out.println("endg");

	System.out.println("tith Poisson");
	for(int i=0;i<100;i++){ // Poisson
	    System.out.println("data " + i + " 0.0 " + poissonBin[i] + " 0.0");
	}
	System.out.println("hist");
	System.out.println("disp");
	System.out.println("endg");

    }

    public static void indexCheck(){

	if(index<0){
	    System.err.println("Iliegal number "+ index);
	    index = 0;
	}
	if(index>=dim){
	    System.err.println("Iliegal number "+ index);
	    index = dim-1;
	}

    }
}

