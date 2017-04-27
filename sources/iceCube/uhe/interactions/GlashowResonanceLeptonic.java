package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** 
  The Glashow Resonance reaction with W into the leptonic decay
   <pre>
   \bar{nu_e} + e^{-1} -> \bar{nu_l} + l ^{-1}
   </pre>
  The approximation that
   <pre>
      1. Masses of the produced leptons are negligible
      2. Energy of incoming anti electron neutrino is by far higher
         than the electron mass.
   </pre>
   has been introduced to calculate the differential cross section

   The inelasiticity parameter y is here difined as
   <pre>
     y = 1 - E_{l^{-1}}/E_{\bar{\nu_e}}
   </pre>

   Written by S. Yoshida November 20 2007

*/
public class  GlashowResonanceLeptonic extends Interactions implements Function{

    protected boolean isPerNucleon = true; // Cross secion is given as per target nucleon
    protected double chargePerNucleon = 0.0;  // number of electrons per nucleon
                                              // calculated in the constructor.

    /** Fermi Coupling Constant  [GeV ^{-2}] */
    public static final double G_F = 1.16639e-5;
    /** conversion constant in the natural unit [cm GeV] */
    public static final double hbar_c = 1.97327e-14;
    /** mass of W [GeV] */
    public static final double massW = 80.22;
    /** decay width of W [GeV] */
    public static final double gammaW = 2.12;
    /** mass of Electron [GeV] */
    public static final double massE = Particle.particleMasses[0][1]; // flavor 0 doublet 1

    /** Constructor: Register the ParticlePoint classes 
         and the prodocued flavor - 0 for e, 1 for mu, 2 for tau
    */
    public GlashowResonanceLeptonic(ParticlePoint s, int flavor){
	super(new Particle(0,0), s, flavor); 
	// Progagating particle is an electron neutrino

	// calculate the number of electrons per nucleon
	double massNumber = 0.0;
	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    chargePerNucleon += (double )(s.getNumberOfAtoms(i))*s.getCharge(i);
	    massNumber += (double )(s.getNumberOfAtoms(i))*s.getAtomicNumber(i);
	}
	chargePerNucleon = chargePerNucleon/massNumber;
    }

    /** Checking the particle kind involved with
	a given interaction. Only electron neutrino
	is allowed to be involved with the Glashow resonance
	in the medium. This is the overriden method from Interactions.java
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==0 && p.getFlavor()==0) return true;
	else return false;
    }



    /**
       Calculate the differential cross section as the one per nucleon.
       for instance in the case of ice, the number of electrons per nuclei
       is approximately (1 x 8 + 2 x1) - corresponding
       to H_2 0 !! -- is multiplied to the cross section.
    */
    public void calculateCrossSectionAsPerNucleon(){
	isPerNucleon = true;
    }


    /**
       Calculate the differential cross section as the one per electron.
       This is simple - because the Glashow resonance is subject to
       electron.
    */
    public void calculateCrossSectionAsPerElectron(){
	isPerNucleon = false;
    }


    /**
       return dSigma/dy [cm^2] where y = 1 - - E_{l^{-1}}/E_{\bar{\nu_e}}
    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;
	double invariant_s = 2.0*massE*energy;
	if(!isPerNucleon)
	    return areaFactorByWeakCoupling(invariant_s)*wResonance(invariant_s)*y*y;
	else
	    return chargePerNucleon*areaFactorByWeakCoupling(invariant_s)*
		wResonance(invariant_s)*y*y;
    }

    /** return total cross section [cm^2] */
    public double getSigma(){
	double invariant_s = 2.0*massE*energy;
	if(!isPerNucleon)
	    return areaFactorByWeakCoupling(invariant_s)*wResonance(invariant_s)/3.0;
	else
	    return chargePerNucleon*areaFactorByWeakCoupling(invariant_s)*
		wResonance(invariant_s)/3.0;
    }


    /** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
    public boolean isValidInelasticity(double y){ 
        if((getYmin()+roundOffError)<= y && 
           y <= getYmax()) return true;
        else return false;
    }

    public double getYmin(){
	return 0.0;
    }

    public double getYmax(){
	return 1.0;
    }

    /** the area given by the Weak coupling constant 
	<pre>
	A = G_F^2 s / pi    s: Lorentz invariant energy squared
        </pre>
    */
    protected double areaFactorByWeakCoupling(double invariant_s){
	return G_F*G_F*invariant_s*hbar_c*hbar_c/Math.PI;
    }

    /** W resonance function (dimension less)
        described by the Lorentzian */
    protected double wResonance(double invariant_s){
	double width = (invariant_s-massW*massW)*(invariant_s-massW*massW)
	    + massW*massW*gammaW*gammaW;
	return massW*massW*massW*massW/width;
    }

    public String interactionName(){
	String channel = "Glashow Resonance Leptonic ";
	return channel;
    }
}
