package iceCube.uhe.interactions;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import java.io.*;

/** Make the differential cross section table of Photo-Nuclear Interactions.
    The generated data will be used for PhotoNuclearFit object. */
public class MakePhotoNuclearTable {

    public static void main(String[] args) throws IOException {

        String fileName = null;
	int indexE;
	int dim = 32;
	int flavor = 1;
        int doublet =1;        // Charged Lepton
	int material = 0;      // Ice->0, Rock->1
	double energy = 1.0e9; // 1EeV= 10^9 GeV
	double epsilon = 0.01;
	double ln10 = Math.log(10.0);
	double[ ] logEArray = new double[dim];
	double[ ] logDyArray = new double[14];
	double[ ][ ] yDsigmaArray = new double[dim][14];



        if(args.length!=2){
            System.out.println("Usage: MakePhotoNuclearTable file-name flavor");
	    System.exit(0);
        }else{
            fileName = args[0];
            flavor = Integer.valueOf(args[1]).intValue();
        }


        // Generate the particle class.
        Particle lepton = 
	    new Particle(flavor, doublet, energy); // Charged Leptons


        System.err.println("The Particle Name is " + 
        lepton.particleName(lepton.getFlavor(), lepton.getDoublet()));

	// Generate the ParticlePoint class.
	ParticlePoint s = new ParticlePoint(0.0, 5.0*Math.PI/180.0,material);

	for(int i=0;i<s.NumberOfSpecies[material];i++){
	    System.out.println("Charge " + s.getCharge(i));
	    System.out.println("Atomic Number " + s.getAtomicNumber(i));
	}

	//Generate object of the Photonuclear interaction.
	PhotoNuclear photoNucl = new PhotoNuclear(lepton, s);
	System.err.println(photoNucl.interactionName( ));

	numRecipes.Integration.setRelativeAccuracy(epsilon);
	// integration accuracy for save CPUtime

	photoNucl.setIncidentParticleEnergy(0);
	double Ymax = photoNucl.getYmax();

	int iLogE;
	for(iLogE=0;iLogE<=lepton.getDimensionOfLogEnergyMatrix();iLogE+=20){
	    if(iLogE>600){
		iLogE = lepton.getDimensionOfLogEnergyMatrix();
		indexE = dim-1;
	    }
	    else{
		indexE = iLogE/20;
	    }
	    photoNucl.setIncidentParticleEnergy(iLogE);
	    energy = photoNucl.getIncidentParticleEnergy();
	    System.err.println("The Incident energy " + 
			   energy + " GeV");
	    logEArray[indexE] = 
		Math.log(energy)/ln10;

	    // differential cross sections 
	    int jLogY;
	    if(indexE%5!=0){
		System.err.println("iLogE " + indexE + " jLogE " + (indexE/5)*5);
		for(jLogY=0;jLogY<14;jLogY++){
		    yDsigmaArray[indexE][jLogY] = yDsigmaArray[(indexE/5)*5][jLogY];
		}
	    }else{		
		for(jLogY=0;jLogY<14;jLogY++){
		    double Y;
		    if(jLogY<3){
			Y = Math.pow(10.0,-0.05*(double )jLogY)*Ymax;
		    }else if(jLogY<6){
			Y = Math.pow(10.0,-0.33*(double )(jLogY-2))*Ymax;
		    }else{
			Y = Math.pow(10.0,-(0.5*(double )(jLogY-5)+0.99))*Ymax;
		    }
		    logDyArray[jLogY] = Math.log(Y)/ln10;
		    double dSigma = photoNucl.getDSigmaDy(Y);
		    if(dSigma>0.0){
			yDsigmaArray[indexE][jLogY]= 
			    Math.log(dSigma*Y/(energy*1.0e-32))/ln10;
			System.err.println("  Transfer Matrix " + jLogY + 
					   " Y " + Y + " " +
					   yDsigmaArray[indexE][jLogY]);
		    }else{
			yDsigmaArray[indexE][jLogY]= -50.0;
		    }
		}

	    }
	}

        DataOutputStream out =  
            new DataOutputStream(new FileOutputStream(fileName));

	int jLogY;
	for(jLogY=0;jLogY<14;jLogY++) out.writeDouble(logDyArray[jLogY]);
	for(iLogE=0;iLogE<dim;iLogE++){
	    out.writeDouble(logEArray[iLogE]);
	    for(jLogY=0;jLogY<14;jLogY++){
		out.writeDouble(yDsigmaArray[iLogE][jLogY]);
	    }
	}

        out.close( );

    }
}
