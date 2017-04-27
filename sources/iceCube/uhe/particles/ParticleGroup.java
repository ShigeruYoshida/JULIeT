package iceCube.uhe.particles;


/** 
    <pre>
    This class creates the objects for all leptons (+ pions)
    defined by Particle.class. 
    The created objects of particles are accomodated in the array
             flavor  0      1      2        3
doublet 
    0              e-nu   mu-nu   tau-nu    hadron(pi0)
    1              e      muon    tauon     hadron(pi+)
    </pre>
*/

public class ParticleGroup {

    private Particle[][] particles;

    public ParticleGroup( ){
	particles = 
	new Particle[Particle.NumberOfFlavor][Particle.NumberOfDoublet];

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    for(int doublet=0;doublet<Particle.NumberOfDoublet;doublet++){
		particles[flavor][doublet] = new Particle(flavor, doublet);
		particles[flavor][doublet].generateLogEnergyMatrix( );
	    }
	}
    }

    public ParticleGroup(double energy){
	particles = 
	new Particle[Particle.NumberOfFlavor][Particle.NumberOfDoublet];

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    for(int doublet=0;doublet<Particle.NumberOfDoublet;doublet++){
		particles[flavor][doublet] = 
		new Particle(flavor,doublet,energy);
		particles[flavor][doublet].generateLogEnergyMatrix( );
	    }
	}
    }

    public Particle getParticle(int flavor, int doublet){
	return particles[flavor][doublet];
    }
}


