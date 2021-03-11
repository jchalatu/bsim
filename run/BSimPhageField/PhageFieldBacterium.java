package BSimPhageField;

import bsim.BSim; 

import bsim.BSimChemicalField;
import bsim.winter2021.Bacterium;

import javax.vecmath.Vector3d;
import java.util.*;
/**
 */
public class PhageFieldBacterium extends Bacterium {
	
	/** Chemical field of the simulation. */
	protected BSimChemicalField field;
	
	/** Concentration of phage produced by a bacterium. */
	final private double phageNum;				
	/** Phage production delay in a bacterium (In hours). */
	final private long productionDelay;		
	/** Bacterium life span after infection (In hours). */
	private int lifeSpan;					
	/** Threshold of phage concentration for bacterium to become infected [Molecules/(micron)^3]. */
	final private double threshold;  
	
	/** Phage production delay counter in a bacterium. */
	private int pDelayCount;	
	/** Life span counter in a bacterium. */
	private int lifeCount;	
	
	/** Status of phage occupation within bacterium. */
	private boolean infected;
	/** Flag for phage production. */
	private boolean production;			
	/** Flag for cell death. */
	private boolean isDying;
	
	/** Time of infection. */
	private double infection_time = -1;
    
    /** Elongation rate of infected bacteria. */
    private double infectedGrowthRate;
    private double inf_growth_mean = 0.217;			// From stork paper (0.217 +-0.0434/hr)
    private double inf_growth_stdv = 0.0434;		
    
	/** Elongation rate of dying cells. */
    private double shrinkRate;
    private double shrink_mean = -0.217;			// From stork paper (-0.217 +-0.0434/hr)
    private double shrink_stdv = 0.0434;	
    
	// Function you call when you want to make a new bacterium object
    public PhageFieldBacterium(BSim sim, BSimChemicalField field, Vector3d px1, Vector3d px2) {
		super(sim, px1, px2);
		
        this.field = field;						
        phageNum = 1000;				// 1000 phages/hr
        productionDelay = 2; 			// hrs 
        lifeSpan = 20;					// hrs
        threshold = 1e2;	
        
        pDelayCount = 0;
        lifeCount = 0;
        
        infected = false;
        production = false;
        isDying = false;  
        
        Random bacRng = new Random(); 		// Random number generator
        
        // Calculated once
        infectedGrowthRate = inf_growth_stdv * bacRng.nextGaussian() + inf_growth_mean;
        shrinkRate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
	}	

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();
        
		// Infection occurs if phage concentration around cell exceeds a certain threshold
		if( field.getConc(position) > getThreshold() ) {
			setInfected( true ); 				// Update infection status
		}
		
        // Save time of infection
        if ( getInfectionTime() == -1 && isInfected() ) {
        	setInfectionTime( sim.getTimestep() * sim.getDt() );
        }
        
		// Cell behavior is changed when infected
		if ( isInfected() ) {
			
	        // Assigns a growth rate to infected bacterium according to a normal distribution
			// Infected cells have lower growth rate
	        setK_growth(infectedGrowthRate);
			
			// Cell death
			if ( lifeCount == lifeSpan / sim.getDt() ) {
	    		setDying(true);						// Cell death begins
	    		//elongation_rate = shrinkRate;
	    		setK_growth(shrinkRate);		// Cell starts shrinking
	    		setProduction(false); 				// Cell stops producing phage
			}
			else {
				lifeCount ++;						// Update counter for life span
			}
		    
			if ( isProducing() ) {					// Steady production of phage (phage/hr)
				field.addQuantity(position, phageNum * sim.getDt() );	
			}
			
			// Accounts for the delay regarding internal infection dynamics 
			// when first being infected
			else {									
				if ( pDelayCount == productionDelay / sim.getDt() ) {	// Delay in hours
					setProduction( true );
				}
				pDelayCount ++;
			}
		}
    }
    
    /** Sets the life span of a bacterium. */
    public void setLifeSpan( int lifeSpan ) {this.lifeSpan = lifeSpan;}
    
    /** Sets the status of phage infection for a bacterium. */
    public void setInfected( boolean infected) {this.infected = infected;}
    /** Returns the status of phage infection for a bacterium. */
    public boolean isInfected() {return this.infected;}
	
    /** Returns the phage threshold for a bacterium. */
    public double getThreshold() {return this.threshold;}
    
    /** Returns the flag for the internal phage production of a bacterium. */
    public boolean isProducing() {return this.production;}
    /** Sets the flag for the internal phage production of a bacterium. */
    public void setProduction( boolean production) {this.production = production;}
    
    /** Returns the flag for cell death. */
    public boolean isDying() {return this.isDying;}
    /** Sets the flag for cell death. */
    public void setDying( boolean d ) {this.isDying = d;}
    
    /** Returns the time of infection. */
    public double getInfectionTime() {return infection_time;}
    /** Sets the time of infection. */
    public void setInfectionTime( double t ) {infection_time = t;}
    
    // Allows us to change growth rate mechanics for individual cells
    @Override
    public void setK_growth(double k_growth) {
        super.setK_growth(k_growth);
    }

    // Allows us to change the division threshold of individual cells
    public void setElongationThreshold(double len) { this.L_th=len; }

	@Override
    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public PhageFieldBacterium divide() {
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
        PhageFieldBacterium child = new PhageFieldBacterium(sim, field, x1_child, new Vector3d(this.x2));
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
