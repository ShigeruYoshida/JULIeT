package iceCube.uhe.interactions;

import iceCube.uhe.particles.Particle;
import iceCube.uhe.points.ParticlePoint;
import java.io.*;
import numRecipes.Function;
import numRecipes.Interpolation;

public class PhotoNuclear extends Interactions implements Function {

    static final double ln10 = Math.log(10.0);
    private String photoNuFile    = null;
    private String photoNuFileMu  = "iceCube/uhe/interactions/BB/muPhoto.dat";
    private String photoNuFileTau = "iceCube/uhe/interactions/BB/tauPhoto.dat";
    private double[]   logEArray    =new double[7];
    private double[][] parameter    = new double[7][8];
    private double[]   hardTermTable= new double[7];
    private final double ProtonMass = 938.28e-3;  //Proton Mass
    private double nuclearMass = ProtonMass;

    /** Constructor */
    public PhotoNuclear(Particle p, ParticlePoint s) throws IOException{

        super(p, s, 3); // produced particle is pion

        if(p.getFlavor() == 1 && p.getDoublet() == 1) // Mu
            photoNuFile = photoNuFileMu;
        else if(p.getFlavor() == 2 && p.getDoublet() == 1) // Tau
            photoNuFile = photoNuFileTau;

        File file = new File(photoNuFile);
        BufferedReader in = new BufferedReader(new FileReader(file));
	String buffer;
	char   separator = ' '; int sepstart = 0; int sep = 0;

        for(int e=0; e<7; e++){
	    buffer = in.readLine();

	    // Read log energy [GeV]
	    sepstart = 0;
	    sep = buffer.indexOf(separator,sepstart+1);
            logEArray[e] = Double.valueOf(buffer.substring(sepstart, sep)).doubleValue();

	    // Read parameters (a0->a7)
            for(int k=0; k<8; k++){
		sepstart = sep;
                sep = buffer.indexOf(separator, sepstart+1);
                parameter[e][k] = Double.valueOf(buffer.substring(sepstart,sep)).doubleValue();
            }

        }

        in.close();
    }

    /** Differential cross section dsigma/dy [cm^2]
	y = 1 - Erecoiling/Eincoming -- inelasticity parameter */

    public double getDSigmaDy(double y){
        if(!isValidInelasticity(y)) return 0.0;
        double dSigmaDy = 0.0;
        double softTerm = 0.0;
        double atomicWeight = 0.0;
        double factor = Alpha/(8.0*Math.PI);

	for(int i=0; i<s.NumberOfSpecies[s.getMaterialNumber()];i++){
	    softTerm += s.getAtomicNumber(i)*(double)s.getNumberOfAtoms(i)*getSoftTerm(y, i);
	    atomicWeight += s.getAtomicNumber(i)*(double)s.getNumberOfAtoms(i);
        } 

        dSigmaDy = factor*softTerm*getAbsorptionTerm(y)*1.0e-30*y + atomicWeight*getHardTerm(y);
        return dSigmaDy;
    }

    private double getSoftTerm(double y, int i){
        double soft = 0.0;
        double h    = 1.0 - 2.0/y + 2.0/y/y;
        double z    = 0.00282*Math.pow(s.getAtomicNumber(i), 1.0/3.0)*getAbsorptionTerm(y);
        double t    = (mass*mass*y*y)/(1.0-y);
        double mSq  = mass*mass;
        double m1Sq = 0.54*0.54; //[GeV^2]
        double m2Sq = 1.80*1.80; //[GeV^2]
        soft = (h+2.0*mSq/m2Sq)*Math.log(1.0+m2Sq/t)
	    - 2.0*mSq/t*(1.0-0.25*m2Sq/t*Math.log(1.0+t/m2Sq))
	    + getG(z, s.getCharge(i))*(h*(Math.log(1.0+m1Sq/t)-m1Sq/(m1Sq + t))
				       + 4.0*mSq/m1Sq*Math.log(1.0+m1Sq/t)
				       - 2.0*mSq/t*(1.0-(0.25*m1Sq-t)/(m1Sq + t)));
        return soft;
    }

    private double getAbsorptionTerm(double y){
        double absorptionTerm = 114.3+1.647*Math.pow(Math.log(0.0213*y*energy),2.0);
        return absorptionTerm;
    }

    private double getG(double z, double charge){
        double gz = 0.0;
        if(charge==1.0)
            gz = 3.0;
        else
            gz = 9.0/z*(0.5+((1.0+z)* Math.exp(-z)-1.0)/z/z);
        return gz;
    }

    private double getHardTerm(double y){
        double hardPart  = 0.0;
        double logEnergy = Math.log(energy)/ln10;

        if(isValidInelasticity(y)){
            int index = Interpolation.searchIndex(logEArray, logEnergy, 7);
            if(index >= 5) index = 5;
            double logELow = logEArray[index];
            double logEUp = logEArray[index+1];
            double hardPartLow = 0.0;
            double hardPartUp = 0.0;
            double logTerm = Math.log(y) / ln10;

            for(int k=0; k<8; k++){
                hardPartLow += 1.0/y*parameter[index][k]*Math.pow(logTerm,(double)k)*1.0e-30;
                hardPartUp  += 1.0/y*parameter[index+1][k]*Math.pow(logTerm,(double)k)*1.0e-30;
            }

            hardPart = hardPartLow + (hardPartUp-hardPartLow)*(logEnergy-logELow)/(logEUp - logELow);
        }
        return hardPart;
    }

    /** Checking the range of the given inelasticity y
	that is determined in an individual interaction channel. */
    public boolean isValidInelasticity(double y){
	if((getYmin()+roundOffError) <= y &&
	   y <= (getYmax() - roundOffError)) return true;
	else return false;
    }

    /** Getting the range of allowed y for a given interaction */
    public double getYmin(){
        if(energyCut == 0.0){
            double yMin = 1.0e-6;
            return yMin;
        }else{
            return energyCut/energy;
        }
    }

    public double getYmax(){
        double yMax = 1.0-mass/energy;
        return yMax;
    }

    
    /** Checking the particle kind involved with a given interaction.
	Only muon and tauon but no neutrinos and pions are allowed to be
	involved with the photo-nuclear interaction in the medium. */
    public boolean isValidParticle(Particle p){
        if(p.getDoublet()==1 && p.getFlavor()>=0 && 3>p.getFlavor())
	    return true;
	else return false;
    }

    /** Show the name of interaction */
    public String interactionName(){
        String channel = "Photonuclear interaction ";
        String incidentParticle = p.particleName(p.getFlavor(), p.getDoublet());
        String name = channel.concat("from ").concat(incidentParticle);
        return name;
    }


}

