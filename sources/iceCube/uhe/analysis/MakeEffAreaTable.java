package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** 
Calculate the in-ice effective area table
by running I3ParticleFlux.getInIceEffectiveArea(double logEnergy, 
double cosZenith, int flavor, int doublet). The calculated
table is outputs to std-out in the readable format 
for EffAreaTable class.

Written by S.Yoshida 2008 February 10th
*/

public class MakeEffAreaTable {

    static final double ln10 = Math.log(10.0);

    public static void main(String[] args) throws IOException{

	int flavor = 0;
	int doublet =0;

	String muonDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1Mu_mtx_allfluxes_I3Particles";
	String tauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCIC80DatSet499-641Slope1TAU_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1TAU_mtx_flux_I3Particles";
	String nueDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuE_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuE_mtx_flux_I3Particles";
	String numuDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuMu_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuMu_mtx_flux_I3Particles";
	String nutauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCIC80Slope1NuTau_mtx_flux_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC22Slope1NuTau_mtx_flux_I3Particles";

	if(args.length!=2){
            System.out.println("Usage:MakeEffAreaTable inice-flavor inice-doublet");
            System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    doublet = Integer.valueOf(args[1]).intValue();
	}

	System.err.println("input particle " + Particle.particleName(flavor,doublet));

	//Criteria criteria = new Criteria();
	//criteria.setEHESuperCut();
	CriteriaIC22 criteria = new CriteriaIC22();
        criteria.setCOBZcut(); 
	criteria.setEHESuperCut();

	// reading the I3Particle's
	InputStream in;
	I3ParticleFlux iceFlux = null;
	if(flavor==1 && doublet ==1){ // muon
	    in = ClassLoader.getSystemResourceAsStream(muonDataFileName);
	    iceFlux = new I3ParticleFlux(in,criteria,false);
	    in.close();
	}
	else if(flavor==2 && doublet ==1){ // tau
	    in = ClassLoader.getSystemResourceAsStream(tauDataFileName);
	    iceFlux = new I3ParticleFlux(in,criteria,false);
	    in.close();
	}
	else if(flavor==0 && doublet ==0){ // nu-e
	    in = ClassLoader.getSystemResourceAsStream(nueDataFileName);
	    iceFlux = new I3ParticleFlux(in,criteria,false);
	    in.close();
	}
	else if(flavor==1 && doublet ==0){ // nu-mu
	    in = ClassLoader.getSystemResourceAsStream(numuDataFileName);
	    iceFlux = new I3ParticleFlux(in,criteria,false);
	    in.close();
	}
	else if(flavor==2 && doublet ==0){ // nu-mu
	    in = ClassLoader.getSystemResourceAsStream(nutauDataFileName);
	    iceFlux = new I3ParticleFlux(in,criteria,false);
	    in.close();
	}
	else{ // wrong parameters!
	    System.err.println(" no such in-ice particles!");
	    System.exit(0);
	}

	//iceFlux.observationTime = 365.0*24.0*3600.0*5.0; // 5year
        iceFlux.observationTime = 2.091744e7; //  242.1 days in [sec] 
 
	iceFlux.switchToReco();

	int numCosZenithBin = 40;
	int numLogEbin = 30; // 30 data points
	for(int ie=1;ie<=numLogEbin;ie++){
	    double logInIceEnergy= 5.1+0.2*(double )(ie-1);
	    for(int j=0;j<numCosZenithBin;j++){
		double cosZenith = -0.975+0.05*(double )j;
		// area in [km^2]
		double area = iceFlux.getInIceEffectiveArea(logInIceEnergy,cosZenith, 
							    flavor, doublet)/1.0e10; 
		System.out.println(" " + ie + " " + area + " ");
	    }
	}

    }


}
