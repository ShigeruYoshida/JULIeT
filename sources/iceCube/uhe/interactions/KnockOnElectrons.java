package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The  Knock-on electrons energy loss involving UHE leptons
    propagating undergound rock and ice are calculated in this subclass.

    The class variables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

*/

public class KnockOnElectrons extends Interactions implements Function{

    private double massRatio;
    private double[] para = new double[2];



    /** Constructor: Register the Particle and ParticlePoint classes */
    public KnockOnElectrons(Particle p, ParticlePoint s){
	super(p, s, 0);
	massRatio = mass/producedMass;
    }




    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	double momentum = Math.sqrt(energy*energy-mass*mass);
	double beta = momentum/energy;
	double factor = 2.0*Math.PI*Re*Re*producedMass/(beta*beta*energy)*
	    (1.0/(y*y)-1.0/y*(beta*beta/getYmax())+0.5);

	double delta = Alpha/(2.0*Math.PI)*Math.log(1.0+2.0*y*energy/producedMass)*
	    (Math.log(4.0*energy*energy*(1.0-y)/(mass*mass))-
	     Math.log(1.0+2.0*y*energy/producedMass));

	double chargeFactor = 0.0;
	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    chargeFactor += s.getCharge(i);
	}
	double dSigmaDy = factor*chargeFactor*(1.0+delta);
	//System.err.println("y= " + y + " dSigma/dy = " + dSigmaDy);
	return dSigmaDy;
    }



    /** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
    public boolean isValidInelasticity(double y){ 
	if((getYmin()+roundOffError)<= y && 
	   y <= (getYmax()-roundOffError)) return true;
	else return false;
    }



    /** Getting the range of allowed y for a given interaction */
    public double getYmin( ){
	if(energyCut == 0.0){
	    double yMin = s.getIonizationE()*1.0e-9/energy;
	    return yMin;
	}
	else {
	    return energyCut/energy;
	}
    }
    public double getYmax( ){
	double momentum = Math.sqrt(energy*energy-mass*mass);
	double beta = momentum/energy;
        double yMax = 2.0*producedMass*beta*beta*energy/
	    (producedMass*producedMass+mass*mass+2.0*producedMass*energy);
	return yMax;
    }


    /** Checking the particle kind involved with
	a given interaction. Only muon and tauon
	but not neutrinos and pions
	are allowed to be involved with the pair creation
	in the medium.
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()>=0 && 3>p.getFlavor())
	return true;
	else return false;
    }




    /** Show Name of the Interaction */
    public String interactionName(){
	String channel = "KnockOnElctrons ";
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("from ").concat(incidentParticle);
	return name;
    }
}


