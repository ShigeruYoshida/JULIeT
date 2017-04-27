package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** 
Draw the in-ice effective area table
by running I3ParticleFlux.getInIceEffectiveArea(double logEnergy, 
double cosZenith, int flavor, int doublet). 

Written by S.Yoshida 2008 December 9th
*/

public class DrawInIceEffArea {

    static final double ln10 = Math.log(10.0);

    public static void main(String[] args) throws IOException{

	int flavor = 0;
	int doublet =0;
	double cosWidth = 1.0;double cosZenith = 0.0;

	String muonDataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";
	String tauDataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC80DatSet499-641Slope1TAU_mtx_flux_I3Particles";
	String nueDataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuE_mtx_flux_I3Particles";
	String numuDataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuMu_mtx_flux_I3Particles";
	String nutauDataFileName = 
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuTau_mtx_flux_I3Particles";

	if(args.length!=4){
            System.out.println("Usage:DrawInIceEffArea inice-flavor inice-doublet cosZenith cosZenithWidth");
            System.exit(0);
	}else{
	    flavor = Integer.valueOf(args[0]).intValue();
	    doublet = Integer.valueOf(args[1]).intValue();
	    cosZenith = Double.valueOf(args[2]).doubleValue();
	    cosWidth = Double.valueOf(args[3]).doubleValue();
	}

	System.err.println("input particle " + Particle.particleName(flavor,doublet));

	Criteria criteria = new Criteria();
	criteria.setEHESuperCut();

	// reading the I3Particle's
	InputStream in;
	I3ParticleFlux iceFlux = null;
	if(flavor==1 && doublet ==1){ // muon
	    in = ClassLoader.getSystemResourceAsStream(muonDataFileName);
	    iceFlux = new I3ParticleFlux(in);
	    in.close();
	}
	else if(flavor==2 && doublet ==1){ // tau
	    in = ClassLoader.getSystemResourceAsStream(tauDataFileName);
	    iceFlux = new I3ParticleFlux(in);
	    in.close();
	}
	else if(flavor==0 && doublet ==0){ // nu-e
	    in = ClassLoader.getSystemResourceAsStream(nueDataFileName);
	    iceFlux = new I3ParticleFlux(in);
	    in.close();
	}
	else if(flavor==1 && doublet ==0){ // nu-mu
	    in = ClassLoader.getSystemResourceAsStream(numuDataFileName);
	    iceFlux = new I3ParticleFlux(in);
	    in.close();
	}
	else if(flavor==2 && doublet ==0){ // nu-mu
	    in = ClassLoader.getSystemResourceAsStream(nutauDataFileName);
	    iceFlux = new I3ParticleFlux(in);
	    in.close();
	}
	else{ // wrong parameters!
	    System.err.println(" no such in-ice particles!");
	    System.exit(0);
	}

	iceFlux.setCriteria(criteria);
	iceFlux.observationTime = 365.0*24.0*3600.0*5.0; // 5year
	iceFlux.cosZenithBinWidth = cosWidth;

	iceFlux.switchToReco();

        System.out.println("titx Energy [GeV]");
        System.out.println("tity Area [km^2!])");
        System.out.println("scal 1.0e6 1.0e11 3.0e-2 3.0");

	int numCosZenithBin = (int)(2.0/cosWidth+0.01);
	int numLogEbin = 60; // 30 data points
	for(int ie=1;ie<=numLogEbin;ie++){
	    double logInIceEnergy= 5.1+0.1*(double )(ie-1);
	    // area in [km^2]
	    double area = iceFlux.getInIceEffectiveArea(logInIceEnergy,cosZenith, 
							flavor, doublet)/1.0e10; 
	    double energy = Math.pow(10.0,logInIceEnergy);
	    System.out.println("data " + energy + " 0.0 " +
			       area  + " 0.0");

	}
	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

    }


}
