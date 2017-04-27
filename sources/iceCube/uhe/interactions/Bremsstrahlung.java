package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The  Bremsstrahlung involving UHE leptons
    propagating undergound rock and ice are calculated in this subclass.

    The class variables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

*/
public class  Bremsstrahlung extends Interactions implements Function{

    private double massRatio;
    private double[] para = new double[2];



    /** Constructor: Register the Particle and ParticlePoint classes */
    public Bremsstrahlung(Particle p, ParticlePoint s){
	super(p, s, 0);
	massRatio = mass/producedMass;
    }




    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	double factor = 4.0*Alpha*Re*Re/(massRatio*massRatio)/y;
	double chargeFactor = 0.0;

	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    if(y<getYmaxCharge(i))
	    chargeFactor += s.getCharge(i)*s.getCharge(i)*
		(double )(s.getNumberOfAtoms(i))*getChargeFactor(y,i);
	}
	double dSigmaDy = factor*chargeFactor;
	//System.err.println("y= " + y + " dSigma/dy = " + dSigmaDy);
	return dSigmaDy;
    }




    /** Calculate the term which depends on the charge of the propation medium */
    public double getChargeFactor(double y, int ithSpecies){
	double yFactor1 = 2.0-2.0*y+y*y;
	double yFactor2 = 2.0/3.0*(1.0-y);
	double qMin = mass*mass*y/(2.0*energy*(1.0-y));
	double chargeTerm = Math.pow(s.getCharge(ithSpecies),1.0/3.0);
	double a1 = 111.7/(chargeTerm*producedMass);
	double a2 = 724.2/(chargeTerm*chargeTerm*producedMass);
	double x1 = a1*qMin;
	double x2 = a2*qMin;
	double muonMass = Particle.particleMasses[1][1];
	double qC = 1.9*muonMass/chargeTerm;
	double qsi = Math.sqrt(1.0+4.0*mass*mass/(qC*qC));
	double delta1;
	double delta2;
	if(s.getCharge(ithSpecies)==1.0){ // H
	    delta1 = 0.0;
	    delta2 = 0.0;
	}else { // Otherwise
	    delta1 = Math.log(mass/qC)+qsi/2.0*
		Math.log((qsi+1.0)/(qsi-1.0));
	    delta2 = Math.log(mass/qC)+qsi*(3.0-qsi*qsi)/4.0*
		Math.log((qsi+1.0)/(qsi-1.0)) + 2.0*mass*mass/(qC*qC);
	}
	double arcTan1=Math.atan(1.0/x1);
	double arcTan2=Math.atan(1.0/x2);

	double psi1 = 0.5*(1.0+Math.log(mass*mass*a1*a1/(1.0+x1*x1)))-x1*arcTan1+
	    (0.5*(1.0+Math.log(mass*mass*a2*a2/(1.0+x2*x2)))-x2*arcTan2)/
	    s.getCharge(ithSpecies) - delta1;
	double psi2 = 0.5*(2.0/3.0+Math.log(mass*mass*a1*a1/(1.0+x1*x1)))
	    +2.0*x1*x1*(1.0-x1*arcTan1+3.0/4.0*Math.log(x1*x1/(1.0+x1*x1)))+
	    (0.5*(2.0/3.0+Math.log(mass*mass*a2*a2/(1.0+x2*x2)))
	     +2.0*x2*x2*(1.0-x2*arcTan2+3.0/4.0*Math.log(x2*x2/(1.0+x2*x2))))/
	    s.getCharge(ithSpecies) - delta2;

	double term = yFactor1*psi1-yFactor2*psi2;
	return term;
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
	    double yMin = 0.0;
	    return yMin;
	}
	else {
	    return energyCut/energy;
	}
    }
    public double getYmax( ){
	double chargeTerm =1.0; // Pure Hydrogen as an upperbound
        double yMax = 1.0-3.0/4.0*Math.sqrt(E)*mass/energy*chargeTerm;
	return yMax;
    }
    public double getYmaxCharge(int ithSpecies){
	double chargeTerm = Math.pow(s.getCharge(ithSpecies),1.0/3.0);
	double yMax = 1.0-3.0/4.0*Math.sqrt(E)*mass/energy*chargeTerm;
	return yMax;
    }


    /** Checking the particle kind involved with
	a given interaction. Only muon and tauon
	but not neutrinos and pions
	are allowed to be involved with the Bremsstrahlung
	in the medium.
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()>=0 && 3>p.getFlavor())
	return true;
	else return false;
    }




    /** Show Name of the Interaction. Returns the Sting where 
	the interation name is written. */
    public String interactionName(){
	String channel = "Bremsstrahlung ";
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("from ").concat(incidentParticle);
	return name;
    }
}


