package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The  Dummy involving UHE leptons
    propagating undergound rock and ice are calculated in this subclass.

    This describes *unphysical* hypothetical interactions for
    debugging purposes.

*/
public class  DummyInteractions extends Interactions implements Function{

    private double massRatio;
    private double[] para = new double[2];



    /** Constructor: Register the Particle and ParticlePoint classes */
    public DummyInteractions(Particle p, ParticlePoint s){
	super(p, s, 0);
	massRatio = mass/producedMass;
    }




    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	double factor = 1.0e-26;     // [cm^2] somewhere between pair creation and photo-nucl
	double dSigmaDy = factor/(getYmax()-getYmin());
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
	    double yMin = 1.0e-8;
	    return yMin;

    }
    public double getYmax( ){
        double yMax = 1.0;
	return yMax;
    }

    /** Checking the particle kind involved with
	a given interaction. Only muon and tauon
	but not neutrinos and pions
	are allowed to be involved with the dummy interaction
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
	String channel = "Dummy ";
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("from ").concat(incidentParticle);
	return name;
    }
}


