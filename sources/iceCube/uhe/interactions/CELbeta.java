package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import numRecipes.*;

/**
   This static class provides the inelasticity coefiicient 
   beta [cm^2/g] as a function of muon energy
*/

public class CELbeta {

    public final static double[] logEarray =
    {
	4.5,
	5.0, 5.5, 6.0, 6.5, 7.0, 7.5,
	8.0, 8.5, 9.0, 9.5, 10.0, 10.5,
	11.0, 11.5, 12.0
    };

    public final static double[] betaArray ={
	3.566625572457458E-6,
	3.566625572457458E-6,
	3.7785955974147316E-6,
	3.963023042283691E-6,
	4.218591206891876E-6,
	4.4619776009127244E-6,
	4.8267428368612324E-6,
	5.198886183075413E-6,
	5.636788707460893E-6,
	6.067432912765281E-6,
	6.52013071944048E-6,
	6.9571529244992115E-6,
	7.423794429292434E-6,
	7.886457778180616E-6,
	8.339887556707196E-6,
	8.821808646101188E-6
    };


    public static double getBeta(double logEnergy) {

	int degOfPol = 6;
	//if(logEnergy < 5.5) degOfPol = 4; // if E < 3x10^5 GeV

	double beta =  
	    Interpolation.mThPolynominalInterpolate(logEarray,betaArray,16,
						    logEnergy,degOfPol);
	return beta;
    }

    public static void main(String[] args) {

        System.out.println("titx Log E[GeV]");
        System.out.println("tity beta");
        System.out.println("scal 5.0 10.0 0.0 10.0");

	int dimension = Particle.getDimensionOfLogEnergyMatrix();
	for(int iLogE=0;iLogE<dimension;iLogE++){
	    double logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    double beta = CELbeta.getBeta(logE)*1.0e6;

	    System.out.println("data " + logE + " 0.0 " + beta + " 0.0");
	}

	//for(int i=0;i<15;i++){
	//    System.err.println(i + " logE=" + logEarray[i] + 
	//		       " beta=" + betaArray[i]);
	//}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
    }

}
