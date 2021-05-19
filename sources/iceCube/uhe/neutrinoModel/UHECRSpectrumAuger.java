package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.particles.*;

import java.io.*;

/** Calculate the UHECR spectrum measured by Auger collaboration.
    The analytical function fitted with the data by Auger
    is employed to calculate the flux. 
    Based on the ICRC 2017 publication (arXiv 1708.06592).


*/
public class UHECRSpectrumAuger implements Function{

    protected static double energy_ankle = 5.08e9;  // The "angle" energy [GeV]
    protected static double energy_suppression = 3.9e10; // The suppression energy [GeV]
    protected static double spectrum_norm = 8.97e-28;  // Normaliztion [/GeV cm2 sec str]
    protected static double alpha1 = 3.293;  // power law index for energies below ankle
    protected static double alpha2 = 2.53;  // power law index for energies above ankle
    protected static double delta_alpha = 2.5;

    /** constructor: do nothing */
    public UHECRSpectrumAuger(){}


    /** Calculate the differential spectrum
      [/GeV cm2 sec str]

      <pre>
       double cosmicRayEnergy :  UHECR energy [GeV]
      </pre>
    */
    public static double getCRDFDE(double cosmicRayEnergy){
	double xE = cosmicRayEnergy/energy_ankle;
	double flux = 0.0;
	if(xE<=1.0){
	    flux = spectrum_norm*Math.pow(xE,-alpha1);
	}else{
	    double xEsuppression = 1.0+Math.pow(cosmicRayEnergy/energy_suppression,delta_alpha);
	    double xEankle = 1.0+Math.pow(energy_ankle/energy_suppression,delta_alpha);
	    double suppressionFactor = xEankle/xEsuppression;
	    flux = spectrum_norm*Math.pow(xE,-alpha2)*suppressionFactor;
	}
	return flux;
    }

    public double getFunction(int functionIndex, double[] parameters, 
			      double x){
	return getCRDFDE(x);
    }


    //** Main function: Drawing the spectrum */
    public static void main(String[] args){
        System.out.println("titx E[GeV]");
       System.out.println("tity log (Flux E^2 [GeV cm^-2 sec^-1 sr^-1])");
       System.out.println("scal 1.0e5 1.0e11 1.0e-10 3.0e-4");
       System.out.println("gwin 0.2 1.0 0.2 0.7");
       System.out.println("logx");
       System.out.println("logy");
       System.out.println("lnth 2");
       System.out.println("lncl 12");

       double logE = 9.0;
       while(logE<11.0){
	   double crEnergy = Math.pow(10.0,logE);
	   double EFlux = UHECRSpectrumAuger.getCRDFDE(crEnergy)*crEnergy*crEnergy;
	   System.out.println("data " + crEnergy + " 0.0 " + EFlux + " 0.0");
	   logE += 0.1;
       }
       
       UHECRSpectrumAuger crFlux = new  UHECRSpectrumAuger();
       double[] parameters = new double[3];
       double cosmicRayEnergy = 1.0e10;
       double crEndEnergy = 1.0e11;
       double intensity = Integration.RombergIntegral(crFlux,0,parameters,cosmicRayEnergy,crEndEnergy);
       System.out.format("mssg Intensity at E > %e (%e)\n",cosmicRayEnergy,intensity);


       System.out.println("join");
       System.out.println("disp");
       System.out.println("endg");
   }

}
