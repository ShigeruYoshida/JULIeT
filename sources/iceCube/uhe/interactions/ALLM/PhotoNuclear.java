package iceCube.uhe.interactions;

import java.io.*;
import java.util.*;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import numRecipes.*;

/** The  PhotoNuclear interactions
    based on the deep-inelastic scattering formalism
    with the ALLM parameterization of the structure function.

    The class variables  Particle, and Point are necessary for the methods
    described here because the cross section depends on the particle propaty
    and the medium like Z and A(atomic number).

*/
public class  PhotoNuclear extends Interactions implements Function{

    private double[] para = new double[2];
    private final double ProtonMass = 938.28e-3; //Proton Mass [GeV]
    private final double GeVcm = 1.9732e-14; //GeV x cm
    private double nuclearMass = ProtonMass;



    /** Constructor: Register the Particle and ParticlePoint classes */
    public PhotoNuclear(Particle p, ParticlePoint s){
	super(p, s, 3); // produced particle is pion (flavor = 3)
    }




    /** 
	<pre>
	Differential cross section dsigma/dxdQ^2 [cm^2]
	l(kIncoming) N(p) --to l(kRecoiling) X 
	Qsquared = -(kIncoming-kRecoiling)^2 (negative four momentum squared) 
	x = Qsquared /(2 pq) 
	</pre>
    */
    public double getDSigmaDxDQsquared(double x, double Qsquared, int ithSpecies){
	double factor = GeVcm*GeVcm*
	    4.0*Math.PI*Alpha*Alpha/(Qsquared*Qsquared*x);
	double chargeFactor = (double )(s.getNumberOfAtoms(ithSpecies))*
		getDISterm(x,Qsquared,ithSpecies);
	double dSigmaDxDQsquared = 
	    factor*chargeFactor*getStructureFunction(x, Qsquared, ithSpecies);
	return dSigmaDxDQsquared;
    }


    /** The Deep-inelastic scattering term */
    public double getDISterm(double x, double Qsquared, int ithSpecies){
	//nuclearMass = s.getAtomicNumber(ithSpecies)*ProtonMass; 
	                                        // Nuclear Mass [GeV]
	double DISterm = 0.0;
	double y = Qsquared/(2.0*energy*nuclearMass*x);
	if(isValidInelasticity(y)&& isValidQsquared(Qsquared,y)){
	    DISterm = 1.0-y-Qsquared/(4.0*energy*energy)+
		(1.0-2.0*mass*mass/Qsquared)*y*y/2.0*
		(1.0+4.0*nuclearMass*nuclearMass*x*x/Qsquared);
	}

	return DISterm;
    }


    /** The ratio of the neutron structure function to that of proton */
    public double npRatio(double x){
	double px = 1.0 -1.85*x + 2.45*x*x - 2.35*x*x*x + x*x*x*x;
	return px;
    }

    /** The ratio of the nuclear structure function to that of nucleon */
    public double anRatio(double x, int ithSpecies){
	double an;
	if(x<=0.0014){
	    an = Math.pow(s.getAtomicNumber(ithSpecies), -0.1);
	}else if (x<0.04) {
	    double logX = Math.log(x)/Math.log(10.0);
	    an = Math.pow(s.getAtomicNumber(ithSpecies),0.069*logX+0.097);
	}else{
	    an = 1.0;
	}
	return an;
    }




    /** Structure Function with the ALLM parametarization */
    public double getStructureFunction(double x, double Qsquared, int ithSpecies){
	double Q0Squared = 0.46017; // [GeV^2]
	double LambdaSq = 0.06527;  // [GeV^2]
        double LnQoverLambda = 1.9530634738;

	double t = 
	    Math.log((Math.log((Qsquared+Q0Squared)/LambdaSq)/LnQoverLambda));

	double WSq = Qsquared/x+nuclearMass*nuclearMass - Qsquared; // CM enerrgy [GeV^2]

	double mPomeronSq = 49.457; // [GeV^2]
	double xPomeron = (Qsquared+mPomeronSq)/
	    (Qsquared+mPomeronSq+WSq-nuclearMass*nuclearMass);
	double fPomeron = getcPomeron(t)*Math.pow(xPomeron,getaPomeron(t))*
	    Math.pow((1.0-x), getbPomeron(t));

	double mReggeonSq = 0.15052; // [GeV^2]
	double xReggeon = (Qsquared+mReggeonSq)/
	    (Qsquared+mReggeonSq+WSq-nuclearMass*nuclearMass);
	double fReggeon = getcReggeon(t)*Math.pow(xReggeon,getaReggeon(t))*
	    Math.pow((1.0-x), getbReggeon(t));


	double m0Squared = 0.31985; // [GeV^2]
	double allmTerm =Qsquared/(Qsquared+m0Squared)*(fPomeron+fReggeon);

	double nuclearTerm = anRatio(x, ithSpecies)*allmTerm*
	    (s.getCharge(ithSpecies)+
            (s.getAtomicNumber(ithSpecies)-s.getCharge(ithSpecies))*npRatio(x));

	return nuclearTerm;

    }

    public double getcPomeron(double t){
	double cp1 = 0.28067;
	double cp2 = 0.22291;
	double cp3 = 2.1979;

	double tTerm = 1.0+Math.pow(t,cp3);

	double cp = cp1+(cp1-cp2)*((1.0/tTerm)-1.0);
	return cp;
    }
    public double getaPomeron(double t){
	double ap1 = -0.0808;
	double ap2 = -0.44812;
	double ap3 = 1.1709;

	double tTerm = 1.0+Math.pow(t,ap3);

	double ap = ap1+(ap1-ap2)*((1.0/tTerm)-1.0);
	return ap;
    }

    public double getcReggeon(double t){
	double cr1 = 0.80107;
	double cr2 = 0.97307;
	double cr3 = 3.4942;

	double tTerm = Math.pow(t,cr3);

	double cr = cr1+cr2*tTerm;
	return cr;
    }
    public double getaReggeon(double t){
	double ar1 = 0.58400;
	double ar2 = 0.37888;
	double ar3 = 2.6063;

	double tTerm = Math.pow(t,ar3);

	double ar = ar1+ar2*tTerm;
	return ar;
    }

    public double getbPomeron(double t){
	double bp1 = 0.60243;
	double bp2 = 1.3754;
	double bp3 = 1.8439;

	double tTerm = Math.pow(t,bp3);

	double bp = bp1+bp2*tTerm;
	return bp;
    }
    public double getbReggeon(double t){
	double br1 = 0.10711;
	double br2 = 1.9386;
	double br3 = 0.49338;

	double tTerm = Math.pow(t,br3);

	double br = br1+br2*tTerm;
	return br;
    }

    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter    */
    public double getDSigmaDy(double y){
	if(!isValidInelasticity(y)) return 0.0;

	double dSigmaDy = 0.0;
	for(int i=0;i<s.NumberOfSpecies[s.getMaterialNumber( )];i++){
	    //nuclearMass = s.getAtomicNumber(i)*ProtonMass; 
	                                        // Nuclear Mass [GeV]
	    double QsquaredMin = mass*mass*y*y/(1.0-y);
	    double QsquaredMax = 2.0*nuclearMass*energy*y-
		((nuclearMass+producedMass)*(nuclearMass+producedMass)-
		 nuclearMass*nuclearMass);

	    para[0] = y; para[1] = (double )i;	    
	    PhotoNuclear photonuclear = this;
	    dSigmaDy += Integration.RombergIntegral(photonuclear,
				5, para, QsquaredMin, QsquaredMax);
	}

	return dSigmaDy;
    }




    /** Allowed range of Qsquared which is defined by
	y, the nuclear mass, the lepton mass and energy.*/
    public boolean isValidQsquared(double Qsquared, double y){
	double QsquaredMin = mass*mass*y*y/(1.0-y);
	double QsquaredMax = 2.0*nuclearMass*energy*y-
	    ((nuclearMass+producedMass)*(nuclearMass+producedMass)-
	     nuclearMass*nuclearMass);
	if(QsquaredMin<=Qsquared && Qsquared<=QsquaredMax) return true;
	else return false;
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
	    double yMin = (2.0*nuclearMass*producedMass+producedMass*producedMass)/
		(2.0*nuclearMass*energy);
	    return yMin;
	}
	else {
	    return energyCut/energy;
	}
    }
    public double getYmax( ){
        double yMax = 1.0-mass/energy;
	return yMax;
    }



    /** 
	<pre>
	Override Method for interface <Function>. 
        Interface the differntial cross sections given here
        to the utility routiner such as the Romberg
        Integration code that is desinged for a genereal
        function in form of Func(x).

	functionIndex     1     dsigma/dy
	functionIndex     2     dsigma/dz  z = y-1
	functionIndex     3     y x dsigma/dy 
	functionIndex     4     z x dsigma/dz 
	functionIndex     5     dsigma/dxdQsquared * Qsquared/2MEy^2 
	</pre>
    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){ 

	double dSigma;
	double y;

        switch (functionIndex) {
        case 1 : 
	    dSigma = getDSigmaDy(x);
	    break;
	case 2 :
	    dSigma = getDSigmaDz(x);
	    break;
	case 3:
	    dSigma = getDSigmaDy(x)*x;
	    break;
	case 4:
	    dSigma = getDSigmaDz(x)*x;
	    break;
	case 5:
	    y = para[0];
	    int ithSpecies = (int )para[1];
	    //nuclearMass = s.getAtomicNumber(ithSpecies)*ProtonMass; 
	    double Qsquared = x;
	    double p = Qsquared/(2.0*nuclearMass*energy);
	    double xB = p/y;
	    dSigma = getDSigmaDxDQsquared(xB,Qsquared,ithSpecies)*p/(y*y);
	    break;
        default:
            dSigma = 0.0;
            System.err.println("Illegal parameters! Index" + functionIndex);
            System.exit(0);
        }

        return dSigma;

    }

    /** Checking the particle kind involved with
	a given interaction. Only muon and tauon
	but not neutrinos and pions
	are allowed to be involved with the photo-nuclear interaction
	in the medium.
    */
    public boolean isValidParticle(Particle p){
	if(p.getDoublet()==1 && p.getFlavor()>=0 && 3>p.getFlavor())
	return true;
	else return false;
    }




    /** Show Name of the Interaction */
    public String interactionName(){
	String channel = "Photonuclear interaction ";
	String incidentParticle = 
	    p.particleName(p.getFlavor(), p.getDoublet());
	String name = channel.concat("from ").concat(incidentParticle);
	return name;
    }
}


