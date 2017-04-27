package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;

/** Calculate the neutrino detection effective area [km^2 sr] 
(neutrino interaction probability convoluted)
by running I3ParticleFlux.

Written by S.Yoshida 2007 April 13th
*/

public class DrawNuEffAreaByI3Particle {

    static final double ln10 = Math.log(10.0);
    public static void main(String[] args) throws IOException{

	String muonDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1_mtx_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC80DatSet498-594Slope1Mu_mtx_allfluxes_I3Particles";
	String tauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1TAU_mtx_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC80DatSet499-641Slope1TAU_mtx_flux_I3Particles";
	String nueDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuE_mtx_I3Particles"
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuE_mtx_flux_I3Particles";
	String numuDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuMu_mtx_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuMu_mtx_flux_I3Particles";
	String nutauDataFileName = 
	    //"iceCube/uhe/analysis/EHEMCReducedIC9Slope1NuTau_mtx_I3Particles";
	    "iceCube/uhe/analysis/EHEMCIC80Slope1NuTau_mtx_flux_I3Particles";

	Criteria criteria = new Criteria();
	//criteria.setThresholdOfLogNpe(4.6);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,4.7,5.8);
	//criteria.setMinimumBound(0.1,4.6);
	//criteria.setThresholdOfNDOMs(80);

	//criteria.setThresholdOfNDOMs(80);
	//criteria.setThresholdOfLogNpe(5.2);
	//criteria.setRangeOfCosineOfZenith(0.1,1.0,6.0,6.5);
	//criteria.setMinimumBound(0.1,5.2);

	criteria.setEHESuperCut();

	double solidAngle = 4.0*Math.PI;


        // Interactive session to set the IceParticle spieces for
	// I3ParticleFlux to read out I3Particle data filled with propagation matrices.
	I3ParticleFlux iceFlux = null;
	InputStream in = null; 

        DataInputStream input = new DataInputStream(System.in);  
        BufferedReader  d     = new BufferedReader(new InputStreamReader(input));  
        String buffer; 

	// muon?
	System.err.print("count muon as inIce particles? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    in = ClassLoader.getSystemResourceAsStream(muonDataFileName);
	    if(iceFlux == null) iceFlux = new I3ParticleFlux(in);
	    else iceFlux.readI3Particles(in);
	    in.close();
	}
	// tau?
	System.err.print("count tau as inIce particles? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    in = ClassLoader.getSystemResourceAsStream(tauDataFileName);
	    if(iceFlux == null) iceFlux = new I3ParticleFlux(in);
	    else iceFlux.readI3Particles(in);
	    in.close();
	}
	// nu-e?
	System.err.print("count nu-e as inIce particles? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    in = ClassLoader.getSystemResourceAsStream(nueDataFileName);
	    if(iceFlux == null) iceFlux = new I3ParticleFlux(in);
	    else iceFlux.readI3Particles(in);
	    in.close();
	}
	// nu-mu?
	System.err.print("count nu-mu as inIce particles? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    in = ClassLoader.getSystemResourceAsStream(numuDataFileName);
	    if(iceFlux == null) iceFlux = new I3ParticleFlux(in);
	    else iceFlux.readI3Particles(in);
	    in.close();
	}
	// nu-tau?
	System.err.print("count nu-tau as inIce particles? [yes(1)/no(0)] ->");  
	buffer   = d.readLine();  
	if(Integer.valueOf(buffer).intValue()==1){
	    in = ClassLoader.getSystemResourceAsStream(nutauDataFileName);
	    if(iceFlux == null) iceFlux = new I3ParticleFlux(in);
	    else iceFlux.readI3Particles(in);
	    in.close();
	}



	iceFlux.setCriteria(criteria);
	//iceFlux.switchToMCTruth();
	iceFlux.switchToReco();

        System.out.println("titx Energy [GeV]");
        System.out.println("tity Area [km^2!])");
        System.out.println("scal 1.0e6 1.0e11 3.0e-7 0.1");



	for(int iLogE = 0; iLogE<Particle.getDimensionOfLogEnergyMatrix();iLogE++){
	    double logNeutrinoEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )iLogE;
	    double area = iceFlux.getYield(logNeutrinoEnergy)*3.0/
		(iceFlux.observationTime*1.0e10*solidAngle);
	    double nuEnergy = Math.pow(10.0,logNeutrinoEnergy);
	    System.out.println("data " + nuEnergy + " 0.0 " +
			       area  + " 0.0");
	}


	System.out.println("logx");
	System.out.println("logy");
	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");

    }


}
