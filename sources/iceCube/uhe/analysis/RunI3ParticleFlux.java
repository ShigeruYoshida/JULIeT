package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** Calculate the detectable Neutrino flux [GeV/cm^2 sec sr] at the surface
by running I3ParticleFlux. Also it can calculate the theretical
fluxes of cosmic neutrios with NeutrinoFlux class
in the neutrinoModel package.

Written by S.Yoshida 2007 April 9th
*/

public class RunI3ParticleFlux {

    static final double ln10 = Math.log(10.0);
    static final int[] model = {
	4,7,8,11
    }; // GZK, Z-Burst, TD-SUSY

    public static void main(String[] args) throws IOException{

	boolean calTheory = false;
	double numberOfEvents = 0.0;
	String muonDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1_mtx_I3Particles";
	    //"iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
	String tauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1TAU_mtx_I3Particles";
	    //"iceCube/uhe/analysis/EHEMCIC80DatSet499-641Slope1TAU_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1TAU_mtx_flux_I3Particles";
	String nueDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuE_mtx_I3Particles"
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuE_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuE_mtx_flux_I3Particles";
	String numuDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuMu_mtx_I3Particles";
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuMu_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuMu_mtx_flux_I3Particles";
	String nutauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuTau_mtx_I3Particles";
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuTau_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuTau_mtx_flux_I3Particles";

	if(args.length!=2){
            System.out.println("Usage: RunI3ParticleFlux numEvents calTheory(1 true)");
            System.exit(0);
	}else{
	    numberOfEvents = Double.valueOf(args[0]).doubleValue();
	    if(Integer.valueOf(args[1]).intValue()==1) calTheory = true;
	}

	//Criteria criteria = new Criteria();
	CriteriaIC22 criteria = new CriteriaIC22();
	//criteria.setThresholdOfLogNpe(4.6);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteria.setMinimumBound(0.1,4.6);
	//criteria.setThresholdOfNDOMs(80);

	//criteria.setThresholdOfNDOMs(80);
	//criteria.setThresholdOfLogNpe(5.2);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,6.0,6.5);
	//criteria.setMinimumBound(0.1,5.2);

	//criteria.setNPEScalingFactor(0.63);
	criteria.setCOBZcut();
	criteria.setEHESuperCut();

	InputStream in = ClassLoader.getSystemResourceAsStream(muonDataFileName);
	I3ParticleFlux iceFlux = new I3ParticleFlux(in,criteria,false); // does NOT use MCtruth
	in.close();

	in = ClassLoader.getSystemResourceAsStream(tauDataFileName);
	iceFlux.readI3Particles(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(nueDataFileName);
	iceFlux.readI3Particles(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(numuDataFileName);
	iceFlux.readI3Particles(in);
	in.close();

	in = ClassLoader.getSystemResourceAsStream(nutauDataFileName);
	iceFlux.readI3Particles(in);
	in.close();

	//iceFlux.observationTime = 365.0*24.0*3600.0*5.0; // 5year
	iceFlux.observationTime = 2.091744e7; //  242.1 days in [sec]


	//iceFlux.switchToMCTruth();

	System.out.println("IceCube bound flux for " + numberOfEvents +
			   " detection");
        System.out.println("titx Log E[GeV]");
        System.out.println("tity log (Flux E^2! [GeV cm^-2! sec^-1! sr^-1!])");
        System.out.println("scal 6.0 11.0 -10.0 -4.0");



	for(int iLogE = 0; iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
	    double logNeutrinoEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )iLogE;
	    double flux = iceFlux.getDFDLogE(logNeutrinoEnergy,numberOfEvents,false);
	    if(flux>0.0){
		double logEFlux = Math.log(flux)/ln10 + 
		    logNeutrinoEnergy -Math.log(ln10)/ln10;
		System.out.println("data " + logNeutrinoEnergy + " 0.0 " +
				   logEFlux  + " 0.0");
	    }
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

		double logNeutrinoEnergy = 6.0;
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
