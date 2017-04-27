package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** Calculate the detectable Neutrino flux [GeV/cm^2 sec sr] at the surface
by running ProrpagationMatrixFlux. Also it can calculate the theretical
fluxes of cosmic neutrios with NeutrinoFlux class
in the neutrinoModel package.

Written by S.Yoshida 2007 April 11th
*/

public class RunPropagationMatrixFlux {

    static final double ln10 = Math.log(10.0);
    static final int[] model = {
	4,7,8,11
    }; // GZK, Z-Burst, TD-SUSY

    public static void main(String[] args) throws IOException{

	boolean calTheory = false;
	double numberOfEvents = 0.0;

	if(args.length!=2){
            System.out.println("Usage: RunI3ParticleFlux numEvents calTheory(1 true)");
            System.exit(0);
	}else{
	    numberOfEvents = Double.valueOf(args[0]).doubleValue();
	    if(Integer.valueOf(args[1]).intValue()==1) calTheory = true;
	}


	PropagationMatrixFlux iceFlux = new PropagationMatrixFlux();


	// Interactive session to set the inIceParticle flavor and doublets
        DataInputStream input = new DataInputStream(System.in); 
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input)); 
        String buffer; 

	boolean interactionEnd = false;
	do{
	    System.err.print("inIceParticle flavor ->"); 
	    buffer   = d.readLine(); 
	    int flavor = Integer.valueOf(buffer).intValue();

	    System.err.print("inIceParticle doublet ->"); 
	    buffer   = d.readLine(); 
	    int doublet = Integer.valueOf(buffer).intValue();

	    iceFlux.addInIceParticle(flavor,doublet);

	    System.err.print("No more inIce particles [yes(1)/no(0)] ->"); 
	    buffer   = d.readLine(); 
	    if(Integer.valueOf(buffer).intValue()==1) interactionEnd = true; 
	}while(!interactionEnd);

	//double time =  365.0*24.0*3600.0*5.0; // 5year
	double time =  2.091744e7; //  242.1 days in [sec]
	iceFlux.setObservationTime(time);
	iceFlux.calculateYield();

	    
	System.out.println("IceCube bound flux for " + numberOfEvents +
			   " detection");
        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2! [GeV cm^-2! sec^-1! sr^-1!])");
        System.out.println("scal 6.0 11.0 -10.0 -4.0");



	double epsilon = 1.0e-4; // round-off margin for binning
	double logNeutrinoEnergy = Particle.getLogEnergyMinimum() + epsilon;
	double logNeutrinoMaxEnergy = Particle.getLogEnergyMinimum()
	+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	while(logNeutrinoEnergy < logNeutrinoMaxEnergy){ // E < 10^12 GeV
	    double flux = iceFlux.getDFDLogE(logNeutrinoEnergy,numberOfEvents,false);
	    if(flux>0.0){
		double logEFlux = Math.log(flux)/ln10 + 
		    logNeutrinoEnergy -Math.log(ln10)/ln10;
		double logE = logNeutrinoEnergy-epsilon;
		System.out.println("data " + logE + " 0.0 " +
				   logEFlux  + " 0.0");
	    }
	    logNeutrinoEnergy += Particle.getDeltaLogEnergy();
	}


	System.out.println("join");
	System.out.println("disp");
	if(calTheory) System.out.println("cont");
	else System.out.println("endg");

	// Draw Theoretical model fluxes
	if(calTheory) {
	    for(int i=0;i<model.length;i++){
		System.out.println("Neutrino Model Flux " + model[i]);
		NeutrinoFlux neutFlux = new NeutrinoFlux(model[i]);

		logNeutrinoEnergy = 6.0;
		while(logNeutrinoEnergy<12.0){
		    double EFlux = 
			neutFlux.getEFlux(logNeutrinoEnergy,1)+ // nu_e
			neutFlux.getEFlux(logNeutrinoEnergy,2)+ // nu_mu
			neutFlux.getEFlux(logNeutrinoEnergy,3); //nu_tau
		    double logEFlux;
		    if(EFlux > 0.0){
			logEFlux = Math.log(EFlux)/ln10;
		    }else{
			logEFlux = -14.0;
		    }
		    System.out.println("data " + logNeutrinoEnergy + 
				       " 0.0 " + logEFlux + " 0.0");
		    logNeutrinoEnergy += 0.1;
		}
		System.out.println("join");
		System.out.println("disp");
		System.out.println("cont");
	    }
	    System.out.println("endg");
	}


    }


}
