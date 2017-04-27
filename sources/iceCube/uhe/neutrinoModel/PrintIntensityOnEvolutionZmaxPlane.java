package iceCube.uhe.neutrinoModel;

import numRecipes.*;
import iceCube.uhe.neutrinoModel.*;
import java.io.*;

public class PrintIntensityOnEvolutionZmaxPlane {

    public static void main(String[] args) throws IOException {
	double opticalDepth = 0.1;
	double neutrinoEnergy = 1.0e8; // [GeV]
	boolean drawSourceFunction = false;
	boolean basedOnIntegralFlux = true;
	NeutrinoFluxFunction neutFluxFunction = null;
	NeutrinoFluxFromSource neutFluxSource = null;

	// Interactive session to set the parameters
	DataInputStream input = new DataInputStream(System.in); 
	BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
	String buffer; 

	System.err.print("Threshold Neutrino Energy [GeV] ->"); 
	buffer   = d.readLine(); 
	neutrinoEnergy = Double.valueOf(buffer).doubleValue();

	System.err.print("Add neutrinos from sources? [yes(1)/no(0)] ->"); 
	buffer   = d.readLine(); 
	if(Integer.valueOf(buffer).intValue()==1) drawSourceFunction = true;

	if(drawSourceFunction){
	    System.err.print("Optical Depth ->"); 
	    buffer   = d.readLine(); 
	    opticalDepth = Double.valueOf(buffer).doubleValue();
	    neutFluxSource = new NeutrinoFluxFromSource();
	    neutFluxSource.setOpticalDepth(opticalDepth);
	}

	neutFluxFunction = new NeutrinoFluxFunction();
	double maxNeutrinoEnergy = 1.0e11;   // 1.0e11 GeV
	double alpha = 2.5;
	double mMax = 5.0;
	double zMax = 1.0;
	double zMaxBound = 5.0;
	double[] parameters = new double[5];
	while(zMax<=zMaxBound){
	    double m = 0.0;
	    while(m<=mMax){
		parameters[0] = alpha;
		parameters[1] = zMax;
		parameters[2] = m;
		double integralFlux = Integration.RombergIntegral(neutFluxFunction, 
							      0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
		if(drawSourceFunction){
		    integralFlux += 
			Integration.RombergIntegral(neutFluxSource, 
						    0, parameters, neutrinoEnergy, maxNeutrinoEnergy);
		}

		System.out.format("m %4.2f zMax %4.2f intensity %8.3e\n",m,zMax,integralFlux);

		m += 0.01;
	    }

	    zMax += 0.1;
	}

    }

}
	
    


 
