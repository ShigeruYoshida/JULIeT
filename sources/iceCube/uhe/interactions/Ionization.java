package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** 
    
    The Ionization loss formula for JULIeT application.
    
    This code introduces some tweaks to handle the coinuous processes
    like the ionization as the Interactions class has been originaly 
    designed for stochastic interaction process which dominates anyway
    in ultra-high energy region.


    The class variables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

*/
public class Ionization extends Interactions {

    private static final double ln10 = Math.log(10.0);

    /** Ion potential table */
    protected final static double[] IonizedPotential = 
    {7.50e-8, // ice [GeV]
     1.364e-7 // rock [GeV]
    };

    /** The coefficinets in the screening effects */
    protected final static double[] c_screen = {-3.5017, -3.7738}; // ice, rock
    protected final static double[] a_screen = {0.09116, 0.08301};
    protected final static double[] m_screen = {3.4773, 3.4120};
    protected final static double log10_screen = 2.0*Math.log(10.0);
    protected final static double[] x0_screen = {0.24, 0.0492};
    protected final static double[] x1_screen = {2.8004, 3.0549};

    /** The dimension constatnt in the dE/dX. */
    protected final static double Kion = 0.307075e-3; // [GeV/g/cm^2]

    /** Lorentz factor of the propagating particle */
    private double gamma = 1.0;

    /** Particle velocity beta */
    private double beta = 1.0;
 
    /** electron mass */
    private double electronMass = Particle.particleMasses[0][1]; // flavor 0 doublet 1


    /** Constructor. 
	Setting the coefficients depending on the propagation medium.
     */
    public Ionization(Particle p, ParticlePoint s){
	super(p,s,0); // "produced flavor" = 0 i.e.,photons(or electrons - does not matter)
	energyCut = 1.0;  // 1 GeV energy cut
	gamma = energy/mass; // The Lorentz factor
	beta = getBeta(gamma);
    }

    public void setIncidentParticleEnergy(double energy){
	super.setIncidentParticleEnergy(energy);
	gamma = energy/mass;
	beta = getBeta(gamma);
	//System.err.println("energy=" + energy + " gamma=" + gamma + " beta=" + beta);
	//System.err.println("  yMax=" + getYmax() + " yUp=" + getYupper());
    }
    public void setIncidentParticleEnergy(int iLogE){
	super.setIncidentParticleEnergy(iLogE);
	gamma = energy/mass;
	beta = getBeta(gamma);
    }

    public double getYcut(){
	return energyCut/energy;
    }

    public double getYmin(){
	return IonizedPotential[s.getMaterialNumber()]/energy;
    }

    public double getYmax(){
	double massRatio = electronMass/mass;
	double yMax = 2.0*electronMass*(gamma*gamma-1.0)/
	    (energy*(1.0+2.0*gamma*massRatio+massRatio*massRatio));
	double upperBound = 1.0-1.0/gamma;
	double y =  Math.min(yMax,upperBound);
	if(yMax<getYmin()) y = getYmin();
	return y;
    }

    public double getYupper(){
	return Math.min(getYmax(),getYcut());
    }

    /** calculate velocity from the Loarentz factor gamma */
    private double getBeta(double gamma){
	return Math.sqrt(1.0-1.0/(gamma*gamma));
    }


    /** dE/dX - A modified Bethe-Bloch formula [GeV/g/cm^2] */
    public double getDEDX(){
	double charge =0.0;
	double massNumber = 0.0;
	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    charge +=  s.getNumberOfAtoms(i)*s.getCharge(i);
	    massNumber +=  s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	}
	double normFactor = Kion*(charge/massNumber)/(2.0*beta*beta);
	double yMax = getYmax();
	double yUpper = getYupper();
	double ion = IonizedPotential[s.getMaterialNumber()];
	double logTerm = Math.log(2.0*electronMass*beta*beta*gamma*gamma*yUpper*energy/(ion*ion));
	double betaTerm = beta*beta*(1.0+yUpper/yMax);
	double gammaTerm = yUpper*yUpper/(2.0*2.0*(1.0+1.0/gamma)*(1.0+1.0/gamma));
	double lowerX = Math.log(getYmin());
	double upperX = Math.log(yUpper);
	//System.err.print(" now integration from " + lowerX + " to " + upperX);
	//double yFactor = Integration.RombergIntegral(this, 5, parameters,lowerX,upperX);
	double mainFactor =  normFactor*(logTerm-betaTerm+gammaTerm-getScreenFactor());
	//System.err.println(" ..done  main=" + mainFactor + " yFactor= " + energy*yFactor);
	//return mainFactor+yFactor*energy;
	return mainFactor;
    }

    /** the ionization y distribution **/
    public double getDNDyDX(double y){
        if(getYmin()> y || y > getYupper()) return 0.0;
	else{
	    double charge =0.0;
	    double massNumber = 0.0;
	    for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
		charge +=  s.getNumberOfAtoms(i)*s.getCharge(i);
		massNumber +=  s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	    }
	    double normFactor = Kion*(charge/massNumber)/(2.0*beta*beta*y*y*energy);
	    double yMax = getYmax();
	    double gammaTerm = y*y/(2.0*(1.0+1.0/gamma)*(1.0+1.0/gamma));
	    return normFactor*(1.0-beta*beta*y/yMax+gammaTerm);
	}
    }

    /** inelasticity correction */
    public double inelasticCorrection(double y){
        if(getYmin()> y || y > getYupper()) return 0.0;
	else{
	    double electronTerm = Math.log(1.0+2.0*y*energy/electronMass);
	    double yMax = getYmax();
	    double yMaxTerm = Math.log((1.0-y/yMax)/(1.0-y));
	    double gammaTerm = Math.log((2.0*gamma*(1.0-y)*electronMass)/(mass*y));
	    return (Alpha/(2.0*Math.PI))*(electronTerm*(2.0*yMaxTerm+gammaTerm)-yMaxTerm*yMaxTerm);
	}
    }

    /** 
	<pre>
	Method for interface <Function>. 
        Interface the differntial cross sections given here
        to the utility methods such as the Romberg
        Integration code that is desinged for a genereal
        function in form of Func(x).

	functionIndex     1     dsigma/dy
	functionIndex     2     dsigma/dz  z = y-1
	functionIndex     3     y x dsigma/dy 
	functionIndex     4     z x dsigma/dz 
	functionIndex     5    dN/dYdX(y=exp x) * coorection(y=exp x) y 
	</pre>
    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){ 
	if(functionIndex<=4) return super.getFunction(functionIndex, parameters,x);
	else{
	    double dSigma;

	    switch (functionIndex) {
	    case 5:
		double y = Math.exp(x);
		dSigma = getDNDyDX(y)*inelasticCorrection(y)*y;
		//dSigma = getDNDyDX(y)*y;
		//		System.err.println(" dN = " + dSigma + " " + x);
		break;
	    default:
		dSigma = 0.0;
		System.err.println("Illegal parameters! Index" + functionIndex);
		System.exit(0);
	    }

	    return dSigma;
	}

    }




    protected double getScreenFactor(){
	double x = Math.log(beta*gamma)/ln10;
	if(x<x0_screen[s.getMaterialNumber()]) return 0.0;
	else if (x<x1_screen[s.getMaterialNumber()]) 
	    return log10_screen*x+c_screen[s.getMaterialNumber()]+
		a_screen[s.getMaterialNumber()]*
		Math.pow(x1_screen[s.getMaterialNumber()]-x,m_screen[s.getMaterialNumber()]);
	else return log10_screen*x+c_screen[s.getMaterialNumber()];
    }

    /** Differential Cross section. This is not a true cross section but
        defined to give dE/dX, because the ionization is continuous process.
    */
    public double getDSigmaDy(double y){
	if(isValidInelasticity(y)){
	    double massNumber = 0.0;
	    for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
		massNumber +=  s.getNumberOfAtoms(i)*s.getAtomicNumber(i);
	    }
	    double dY = getYmax()-getYmin();
	    if(dY<0.0) dY = roundOffError;
	    return (massNumber/(s.NA*energy))*(2.0/(dY*dY))*getDEDX();
	}else{
	    return 0.0;
	}
    }



    /** Checking the range of the given inelasticity y 
	that is determined in an individual interaction channel.
    */
    public boolean isValidInelasticity(double y){ 
        if(getYmin()<= y && y <= getYmax()) return true;
        else return false;
    }

    /** Checking the particle kind involved with
	a given interaction. Only charged particles
	but not neutrinos
	are allowed to be involved with the ionization
	in the medium.
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==1) return true;
	else return false;
    }

    /** Show Name of the Interaction. Returns the Sting where 
	the interation name is written. */
    public String interactionName(){
	String channel = "Ionization ";
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("of ").concat(incidentParticle);
	return name;
    }


}
