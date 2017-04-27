package numRecipes;

/** Interface to get value of a given function in a given class. 
    This will be  used for numerical integration etc. */

public interface Function {

    double getFunction(int functionIndex, double[] parameters, double x);

}
