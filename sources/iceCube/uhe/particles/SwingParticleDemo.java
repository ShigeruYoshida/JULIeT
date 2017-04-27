package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.net.URL;


/** Particle.class Demo program based on Swing User Interface */
public class SwingParticleDemo implements ActionListener {

    final static int START_INDEX = 0;
    Particle uheParticle;
    private int flavor = 0;
    private int doublet = 0;
    private double energy = 0.0;

    // The Generator Panel
    JPanel particleGenerator;
    JLabel javaIconLabel;
    JLabel energyLabel;
    JTextField inputEnergy;
    JComboBox particleFlavor = null;
    JComboBox particleDoublet = null;
    JButton generateParticle;

    // The Display Panel
    JPanel particleDisplay;
    JLabel particleMass, particleLifeTime, particleName, particleEnergy;

    // Main Panel
    JPanel mainPanel;

    //constructor
    public SwingParticleDemo( ) {

	// Create the particle generator and display panels.
	particleGenerator = new JPanel();
	particleDisplay = new JPanel();

        // Add various widgets to the sub panels.
        addWidgets();
        // Create the main panel to contain the two sub panels.
        mainPanel = new JPanel();
        mainPanel.setLayout(new GridLayout(2,1,5,5));
        mainPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));

        // Add the select and display panels to the main panel.
        mainPanel.add(particleGenerator);
        mainPanel.add(particleDisplay);
    }


    // Create and add the widgets to select particle spieces and
    // its propaty.
    private void addWidgets(){
	// Get the images and put them into ImageIcon.
	String imageJava = "images/javalogo52x88.gif";
	URL javaIconURL =ClassLoader.getSystemResource(imageJava);
	ImageIcon javaIcon = new ImageIcon(javaIconURL);

	// The widgets to select and generate particle.
	javaIconLabel = new JLabel();
        javaIconLabel.setHorizontalAlignment(JLabel.CENTER);
        javaIconLabel.setVerticalAlignment(JLabel.CENTER);
        javaIconLabel.setVerticalTextPosition(JLabel.CENTER);
        javaIconLabel.setHorizontalTextPosition(JLabel.CENTER);
        javaIconLabel.setBorder(BorderFactory.createCompoundBorder(
                            BorderFactory.createLoweredBevelBorder(),
                            BorderFactory.createEmptyBorder(5,5,5,5)));

        javaIconLabel.setBorder(BorderFactory.createCompoundBorder(
                            BorderFactory.createEmptyBorder(0,0,10,0),
                            javaIconLabel.getBorder()));
	javaIconLabel.setIcon(javaIcon);
	javaIconLabel.setText("");

	String[] flavorString = 
	{"e   Flavor","mu  Flavor","tau Flavor","Hadron"};
	particleFlavor = new JComboBox(flavorString);
	particleFlavor.setSelectedIndex(START_INDEX);

	String[] doubletString = {"Doblet 1","Boublet 2"};
	particleDoublet = new JComboBox(doubletString);
	particleDoublet.setSelectedIndex(START_INDEX);

	energyLabel = new JLabel("log (Eenergy [GeV])");
	inputEnergy = new JTextField(4);

	generateParticle = new JButton("Generate...");

	//Add the labels to the panal.
        particleGenerator.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Select Particle"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
	particleGenerator.setLayout(gridbag);
	c.fill = GridBagConstraints.HORIZONTAL; 
	c.weightx = 0.5;

        c.gridx = 0;
        c.gridy = 0;
        c.gridheight = 4;
        gridbag.setConstraints(javaIconLabel, c);
	particleGenerator.add(javaIconLabel);

        c.gridx = 1;
        c.gridy = 0;
        c.gridheight = 1;
        gridbag.setConstraints(particleFlavor, c);
	particleGenerator.add(particleFlavor);

        c.gridx = 2;
        c.gridy = 0;
        gridbag.setConstraints(particleDoublet, c);
	particleGenerator.add(particleDoublet);

        c.gridx = 1;
        c.gridy = 1;
        gridbag.setConstraints(energyLabel, c);
	particleGenerator.add(energyLabel);

        c.gridx = 2;
        c.gridy = 1;
        gridbag.setConstraints(inputEnergy, c);
	particleGenerator.add(inputEnergy);

        c.gridx = 2;
        c.gridy = 3;
        gridbag.setConstraints(generateParticle, c);
	particleGenerator.add(generateParticle);

	// The widgets to display particle propaties.
	particleName = new JLabel("Name",SwingConstants.LEFT);
	particleMass = new JLabel("Mass",SwingConstants.LEFT);
	particleLifeTime = new JLabel("LifeTime",SwingConstants.LEFT);
	particleEnergy = new JLabel("Energy",SwingConstants.LEFT);

	//Add the labels to the panal.
        particleDisplay.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Particle Propaty"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	particleDisplay.setLayout(new GridLayout(2,2));
	particleDisplay.add(particleName);
	particleDisplay.add(particleMass);
	particleDisplay.add(particleLifeTime);
	particleDisplay.add(particleEnergy);

	// Listen to events from Generatorbuttun.
	generateParticle.addActionListener(this);

    }

    // Implementation of ActionListener interface.
    public void actionPerformed(ActionEvent event) {
	flavor = particleFlavor.getSelectedIndex();
	doublet = particleDoublet.getSelectedIndex();
	if(Particle.isValidFlavor(flavor) && Particle.isValidDoublet(doublet)){
	    double logE = Double.parseDouble(inputEnergy.getText());
	    energy = Math.pow(10.0,logE);
	    uheParticle = new Particle(flavor,doublet,energy);
	    particleMass.setText("Mass "+ uheParticle.getMass() + "[GeV]");
	    particleLifeTime.setText(uheParticle.getLifeTime() + "[sec]");
	    particleEnergy.setText("Energy "+ 
				   uheParticle.getEnergy() + "[GeV]");
	    particleName.setText(
		uheParticle.particleName(uheParticle.getFlavor(),
					 uheParticle.getDoublet()));
	}else{
	    particleName.setText("Invalid flavor/doublet");
	}
    }



    // main method
    public static void main(String[] args) {
	SwingParticleDemo swingParticle = new SwingParticleDemo();
	JFrame particleFrame = new JFrame("Particle Generator");

        // Set the look and feel.
        try {
            UIManager.setLookAndFeel(
                UIManager.getCrossPlatformLookAndFeelClassName());
        } catch(Exception e) {}

	particleFrame.setContentPane(swingParticle.mainPanel);
	particleFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        // Show the generator.
        particleFrame.pack();
        particleFrame.setVisible(true);
    }
}
