package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** 
  The Glashow Resonance reaction with W into the hadronic decay
   <pre>
   \bar{nu_e} + e^{-1} -> hadrons
   </pre>
   The super class GlashowResonanceLeptonic class is used 
   in most of the calculation

   The inelasiticity parameter y is fixed to be always 0
   because all the final states are hadrons that generates cascades at once.
*/
public class  GlashowResonanceHadronic extends GlashowResonanceLeptonic implements Function{


    /** Ration Gamma_hadron/Gamma_muon = Gamma_hadron/Gamma_tau */
    protected final double ratioOfDecayWidth = 6.523;

    public GlashowResonanceHadronic(ParticlePoint s){
	super(s, 3); 
	// Progagating particle is hadron
    }


    /**
       return dSigma/dy [cm^2] where y = 1 - - E_{l^{-1}}/E_{\bar{\nu_e}}
    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	if((getYmax()-roundOffError)<= y && y <= getYmax()){
	    return getSigma()/roundOffError;
	}else{
	    return 0.0;
	}
    }

    /** return total cross section [cm^2] */
    public double getSigma(){
	return super.getSigma()*ratioOfDecayWidth;
    }

    public String interactionName(){
	String channel = "Glashow Resonance Hadronic ";
	return channel;
    }
}
