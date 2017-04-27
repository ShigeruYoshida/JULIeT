package iceCube.uhe.particles;


/** 
    <pre>
    This class creates the objects of array "logEnergyMatrix"
    for all leptons (+ pions)  defined by Particle.class. 
    The created objects of particles are accomodated in the array
             flavor  0      1      2        3
doublet 
    0              e-nu   mu-nu   tau-nu    hadron(pi0)
    1              e      muon    tauon     hadron(pi+)

    The Particle object array particles[flavor][doublet][ilogE]
    is generated. iLogE is related to logEnergy by
    logEnergy = Particle.getLogEnergyMinimum()+
    Particle.getDeltaLogEnergy()*(double )iLogE;

    Each element of particles[flavor][doublet][ilogE]
    has itws own array to accomodate dn/dLogE type of distribution.
    Its binwidth and origin can be again defined by
    Particle.getLogEnergyMinimum() and Particle.getDeltaLogEnergy().
    Its dimension is defined by Particle.getDimensionOfLogEnergyMatrix().
    </pre>
*/


public class ParticleArray {

    private Particle[][][] particles;

    public ParticleArray( ){
	particles = new 
	    Particle[Particle.NumberOfFlavor][Particle.NumberOfDoublet]
	    [Particle.getDimensionOfLogEnergyMatrix()];

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    for(int doublet=0;doublet<Particle.NumberOfDoublet;doublet++){
		for(int iLogE=0;iLogE<Particle.getDimensionOfLogEnergyMatrix();
		    iLogE++){
		    double logEnergy = Particle.getLogEnergyMinimum()+
			Particle.getDeltaLogEnergy()*(double )iLogE;
		    double energy = Math.pow(10.0,logEnergy);
		
		    particles[flavor][doublet][iLogE] = 
			new Particle(flavor,doublet,energy);

	    particles[flavor][doublet][iLogE].generateLogEnergyMatrix( );
		}
	    }
	}
    }

    public Particle getParticle(int flavor, int doublet, int ilogE){
	return particles[flavor][doublet][ilogE];
    }
}


