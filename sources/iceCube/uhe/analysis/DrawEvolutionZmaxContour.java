package iceCube.uhe.analysis;

import numRecipes.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;
import java.io.*;

public class DrawEvolutionZmaxContour {

    public static void main(String[] args) throws IOException {
	int mode = 1;
	double neutrinoEnergy = 1.0e8; // [GeV]
	double number = 1.0; 
	double confLevel = 0.95;
	double alpha = 2.5;
	NeutrinoFluxIceCube flux = null;

	if(args.length<2){
	    System.out.println("Usage:DrawEvolutionZmaxContour mode(1 analytical GZK 2 plus in-situ 3 numerical model) event_rate  (alpha default 2.5)");
	    System.exit(0);
	}else if(args.length==2){
	    mode = Integer.valueOf(args[0]).intValue();
	    number = Double.valueOf(args[1]).doubleValue();
	}else{
	    mode = Integer.valueOf(args[0]).intValue();
	    number = Double.valueOf(args[1]).doubleValue();
	    alpha =  Double.valueOf(args[2]).doubleValue();
	    System.err.format(" alpha set by %f\n",alpha);
	}

	//FeldmanCousins.setConfidenceLevel(confLevel);
	//number = FeldmanCousins.getUpperLimit(0.14,0);
	//System.err.format(" nEvents(%f conf level)=%f\n",confLevel,number);

	flux = new NeutrinoFluxIceCube(mode);

	// Draw Contour map on evolution-Zmax plane
        System.out.println("titx m");
        System.out.println("tity Zmax");
        System.out.println("scal 2.0 5.0 1.0 5.0");
        System.out.println("gwin 0.2 0.9 0.2 0.9");

	double logEmax = 11.98;   // 1.0e11.9 GeV
	//double logEmin = 8.0;   // 1.0e8 GeV
	double logEmin = 5.0;   // 1.0e5 GeV
	if(mode==2) logEmin=4.0;
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
		double numberOfEvents = Integration.RombergIntegral(flux, 0, parameters, logEmin,logEmax);

		//System.err.format("number of events(m= %4.2f zMax=%4.2f)=%f\n",m,zMax,numberOfEvents);
		if(numberOfEvents >= number){
		    System.out.format("data %4.2f 0.0 %4.2f\n",m,zMax);
		    break;
		}

		m += 0.01;
	    }

	    zMax += 0.1;
	}
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");


    }

}
	
    


 
