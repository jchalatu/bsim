package Winter2021;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import javax.vecmath.Vector3d;
import java.awt.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import java.lang.Math;
/**
 */
public class Bacterium extends BSimCapsuleBacterium {

	/** Chemical field of the simulation. */
	private BSimChemicalField field;
	
	/** Concentration of phage produced by a bacterium. */
	final private double productionConc;				
	/** Phage production delay in a bacterium (In milliseconds). */
	final private long productionDelay;				 
	/** Bacterium life span after infection (In milliseconds). */
	final private long lifeSpan;					
	/** Threshold of phage density for bacterium to become infected [Molecules/(micron)^3]. */
	final private double threshold;  
	/** Rate the cell shrinks when undergoing cell death. */
	final private double shrinkRate;		
	
	/** Status of phage occupation within bacterium. */
	private boolean infected;
	/** Flag for phage production. */
	private boolean production;			
	/** Flag for cell death. */
	private boolean isDying;	
	/** Timer to schedule internal phage production. */
    protected Timer timer;
    
    /** Sets the status of phage infection for a bacterium. */
    public void setInfected( boolean infected) {
        this.infected = infected;
    }
    
    /** Returns the status of phage infection for a bacterium. */
    public boolean isInfected() {
        return this.infected;
    }
	
    /** Returns the phage threshold for a bacterium. */
    public double getThreshold() {
        return this.threshold;
    }
    
    /** Returns the flag for the internal phage production of a bacterium. */
    public boolean isProducing() {
        return this.production;
    }
    
    /** Sets the flag for the internal phage production of a bacterium. */
    public void setProduction( boolean production) {
        this.production = production;
    }
    
    /** Returns the flag for cell death. */
    public boolean isDying() {
        return this.isDying;
    }
    
    /** Sets the flag for cell death. */
    public void setDying( boolean d) {
        this.isDying = d;
    }
	
    // Function you call when you want to make a new bacterium object
    public Bacterium(BSim sim, BSimChemicalField field, Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        
        this.field = field;						
        productionConc = 1e6;
        productionDelay = 2000;
        lifeSpan = 4000;
        threshold = 1e5;
        shrinkRate = -1;
        
        infected = false;
        production = false;
        isDying = false;  

        timer = new Timer();
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();
        
        /** Cell death. */
        
	    TimerTask cellDeath = new TimerTask() {
	    	public void run() {
	    		setDying(true);					// Cell death begins
	    		setK_growth(shrinkRate);		// Cell starts shrinking
	    		setProduction(false); 			// Cell stops producing phage
			}
	    };
	    timer.schedule(cellDeath, lifeSpan);
        
        /** Updated to allow phage infection and production. */
        
		// Infection occurs if current phage concentration exceeds a certain threshold
		if( field.getConc(position) > getThreshold() ) {
			setInfected( true ); 				// Update infection status
		}
		
		// Accounts for the delay regarding internal infection dynamics 
		// when first being infected
		if ( isInfected() ) {
			
			if ( !isDying ) {
				double infectedGrowthRate = 1;		// 0.05
	    		setK_growth(infectedGrowthRate);	// Update growth rate
			}
			else {
				setK_growth(shrinkRate);			// Update growth rate
			}
		    
			if ( isProducing() ) {
				field.addQuantity(position, productionConc);	
			}
			else {
				// Internal production delay when first infected
			    TimerTask task = new TimerTask() {
			    	public void run() {
			    		field.addQuantity(position, productionConc);
			    		setProduction( true );
					}
			    };
			    timer.schedule(task, productionDelay);
			}
			
		}
			
    }

    // Allows us to change growth rate mechanics for individual cells
    @Override
    public void setK_growth(double k_growth) {
        super.setK_growth(k_growth);
    }

    // Allows us to change the division threshold of individual cells
    public void setElongationThreshold(double len) { this.L_th=len; }

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public Bacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble()/100,rng.nextDouble()/100,rng.nextDouble()/100);
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
        longVec.scale(0.05*rng.nextDouble()); 			// Push is applied to bacterium
        
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
        Bacterium child = new Bacterium(sim, field, x1_child, new Vector3d(this.x2));
        this.initialise(L1, this.x1, x2_new);
        child.L = L2;

        // Prints a line whenever a new bacterium is made
        System.out.println("Child ID id " + child.id);
        return child;
    }


}
