package numRecipes;

import numRecipes.*;


/** Display the pseudorandom numbers */
public class RandomGeneratorState{

    public static void main(String[] args){

		RandomGenerator random1 = new RandomGenerator(123);
        double r;
        int times = 0;

		if(args.length!=1){
		    System.out.println("Usage: RandomDoubleDemo repeat-times");
		}
		else{
		    times = Integer.valueOf(args[0]).intValue();
		}

		for(int i=0;i<times;i++){
		    r = random1.GetRandomDouble();
		    System.out.println("Pseudorandom (" + i + "):" + r);
		}

		long[] state = random1.GetState();

		for(int i=0;i<5;i++){
		    r = random1.GetRandomDouble();
		    System.out.println("Pseudorandom1 (" + i + "):" + r);
		}

		// for(int i=0; i<state.length; i++){
		//	System.out.println(state[i]);
		//}

		System.out.println("Before init");
		RandomGenerator random2 = new RandomGenerator(state);
		System.out.println("After init");

		for(int i=0;i<5;i++){
	        r = random2.GetRandomDouble();
	        System.out.println("Pseudorandom2 (" + i + "):" + r);
		}
    }

}

