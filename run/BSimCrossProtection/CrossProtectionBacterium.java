package BSimCrossProtection;

import bsim.BSim; 

import bsim.BSimChemicalField;
import bsim.winter2021.Bacterium;
import javax.vecmath.Vector3d;

import java.util.*;
/**
 * This class represents a bacteria for a cross-protection simulation. 
 * The growth rate of the bacteria depends on it's initial growth rate and the surrounding
 * concentration of antibiotics. 
 * If a cell is surrounded by an antibiotic it isn't resistant to, the growth rate will decrease. 
 * If the concentration of antibiotic around a bacteria is greater than a certain threshold, the cell 
 * will start to shrink and die (turns red in the simulation). 
 */
public class CrossProtectionBacterium extends Bacterium {
	public int lifetime;
	
    /** Scales randomVec during division **/
    static double twist;
    /** Scales longVec during division **/
    static double push;

	/** Chemical field representing the antibiotic the bacteria is resistant to. */
	protected BSimChemicalField resistant_field;
	/** Chemical field representing the antibiotic the bacteria is not resistant to. */
	protected BSimChemicalField antibiotic_field;
	
	/** Initial growth rate of the bacteria. */
	public static double initial_el_mean = 0.2;
	public static double initial_el_stdv = 0.05;
	final private double initial_el_rate;
	
	/** Rate at which the bacteria shrinks when dying. **/
	final private double shrink_stdv = 0.0434;			
	final private double shrink_mean = -0.117;
	final private double shrink_rate;
	
	/** Scales the growth rate of the bacteria. */
	final private double scale_1 = 1;//2;
	/** Scales the concentration for the growth rate of the bacteria. */
	final private double scale_2 = 0.1;
	
	/** Threshold of antibiotics for bacterium growth to slow [Molecules/(micron)^3]. */
	final private double growth_threshold;  
	/** Threshold of antibiotics for bacterium to start dying [Molecules/(micron)^3]. */
	final private double conc_threshold;  
	/** Flag for enzyme production. */
	private boolean production;			
	/** Flag to indicate toxin levels are above threshold. */
	private boolean aboveConcThreshold;
	/** Amount of enzymes produced to decay the antibiotic field. */
	private double enzymeNum;

	
    // Function you call when you want to make a new bacterium object
    public CrossProtectionBacterium(BSim sim, BSimChemicalField antibiotic, BSimChemicalField resistant, Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        
        lifetime = 0;
        
        this.antibiotic_field = antibiotic;	
        this.resistant_field = resistant;
        
        growth_threshold = 1e2;
        conc_threshold = 3e2;	
        production = true;
        aboveConcThreshold = false;
        enzymeNum = 3e3;
        
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        
        // Calculate the initial elongation rate 
        initial_el_rate = initial_el_stdv * bacRng.nextGaussian() + initial_el_mean;
        setK_growth(initial_el_rate);
        // Calculate the rate at which the bacteria will shrink
        shrink_rate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();
        
        // If toxin levels are above the threshold, the bacteria starts to die
        if ( antibiotic_field.getConc(position) > getConcThreshold() ) {
        	aboveConcThreshold = true;
        	
        	// Bacteria starts shrinking 
        	Random bacRng = new Random();	
	        setK_growth( shrink_rate );
        	//setK_growth(0.0);		// stops growing
	        //setProduction(false); 	// stops producing enzymes when dying
        }
        // If toxin levels are below the threshold, 
        // cell growth rate is lowered depending on surrounding antibiotic concentration 
        // (values and scaling are arbitrary for now)
        else {
        	aboveConcThreshold = false;
            if ( antibiotic_field.getConc(position) > growth_threshold ) {
            	setK_growth( initial_el_rate * scale_1/(scale_1 + antibiotic_field.getConc(position) * scale_2 * sim.getDt()));
            } 
        }
        
		// Production of enzymes to decay the antibiotic the bacteris is resistant to
		if ( isProducing() ) {						// Steady production of enzymes (enzyme/hr)
			resistant_field.addQuantity(position, -enzymeNum * sim.getDt() );	
		}
		
    }
	
    /** Returns the antibiotic threshold for a bacterium. */
    public double getConcThreshold() { return this.conc_threshold; }
    /** Returns the flag for the internal enzyme production of a bacterium. */
    public boolean isProducing() { return this.production; }
    /** Sets the flag for the internal enzyme production of a bacterium. */
    public void setProduction( boolean production) { this.production = production; }
    /** Returns the flag that indicates toxin levels above threshold. */
    public boolean isAboveThreshold() { return aboveConcThreshold; }
    /** Sets the elongation mean. **/
    public void set_elMean(double el_mean) { this.initial_el_mean = el_mean; }
    /** Sets the elongation stdv. **/
    public void set_elStdv(double el_stdv) { this.initial_el_stdv = el_stdv; }
    
    /** Sets the value of the twist during division. **/
    public void setTwist(double t) {twist = t;}
    /** Sets the value of the push during division. **/
    public void setPush(double p) {push = p;}
    
    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public CrossProtectionBacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble()/twist,rng.nextDouble()/twist,rng.nextDouble()/twist);
        System.out.println("Bacterium " + this.id + " is dividing...");

        Vector3d u = new Vector3d(); 
        u.sub(this.x2, this.x1);

        // Decides where to split bacterium
        double divPert = 0.1*L*(rng.nextDouble() - 0.5); // 0.1 is arbitrary scaling here.

        // Get length of bacterium
        double L_actual = u.length();

        // Computes lengths of child bacteria
        double L1 = L_actual*0.5*(1 + divPert) - radius;
        double L2 = L_actual*0.5*(1 - divPert) - radius;

        // Finds new endpoints of mother and daughter
        Vector3d x2_new = new Vector3d();
        x2_new.scaleAdd(L1/L_actual, u, this.x1);
        Vector3d longVec = new Vector3d();
        longVec.scaleAdd(-1,this.x2,this.x1); 			// Push along bacterium length
        longVec.scale(push*rng.nextDouble()); 			// Push is applied to bacterium
        
        // Impulse, not a force.
        longVec.add(randomVec);
        x2_new.add(longVec);

        Vector3d x1_child = new Vector3d();
        x1_child.scaleAdd(-(L2/L_actual), u, this.x2);
        Vector3d longVec2=new Vector3d();
        longVec2.scaleAdd(-1,longVec);
        x1_child.add(longVec2); 						// Push applied to child cell

        // Set the child cell.
        // Creates new bacterium called child and adds it to the lists, gives posns, infected status and chemical field status
        CrossProtectionBacterium child = new CrossProtectionBacterium(sim, antibiotic_field, resistant_field, x1_child, new Vector3d(this.x2));
        lifetime = 0;
		// Asymmetrical growth occurs at division node
        this.initialise(L1, x2_new, this.x1);			// Swap x1 and x2 for the mother after division for asymmetrical elongation
		child.L = L2;
        
        // add child to list of children - Sohaib Nadeem
        addChild(child);
                
        // Calculate angle between daughter cells at division - Sheng Fang
        angle_initial = coordinate(child);

        // Prints a line whenever a new bacterium is made
        System.out.println("Child ID id " + child.id);
        return child;
    }


}
