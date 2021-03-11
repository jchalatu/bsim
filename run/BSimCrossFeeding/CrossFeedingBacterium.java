package BSimCrossFeeding;

import bsim.BSim;

import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.winter2021.Bacterium;
import javax.vecmath.Vector3d;

import com.beust.jcommander.Parameter;

import java.awt.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import java.lang.Math;
/**
 */
public class CrossFeedingBacterium extends Bacterium {

	/** Chemical field representing the amino acid the bacteria produces. */
	protected BSimChemicalField production_field;
	/** Chemical field representing the antibiotic the bacteria consumes. */
	protected BSimChemicalField consumption_field;
	
	/** Initial growth rate of the bacteria. */
	final private double initial_growth_mean = 0.01;	
	final private double initial_growth_stdv = 0.001;
	final private double initial_growth_rate;
	
	/** Flag for amino acid production (/hr). */
	private boolean production;			
	/** Amount of amino acid produced (/hr). */
	private double productionNum;
	
	/** Amount of amino acid consumed (/hr). */
	private double consumptionRate;
	/** Max amount of amino acid consumed by a bacterium (/hr). */
	private double consumptionMax;
	
	// Constants for Monod Equation
	
    /** Maximum growth rate of the cell (um/hr). */
    final private double mu_max = 1.3;
    /** The "half-velocity constant"; the value of [S] when mu/mu_max = 0.5 (g/dm^3). */
    final private double K_s = 2.2e-1;//2.2e-5;
    
    // Constants for yield conversion
    final double yield_coefficient = 1.0;			// Arbitrary value between 0 and 1
    double total_bio_mass = 1.0;			
    
    /** Function you call when you want to make a new bacterium object. */
    public CrossFeedingBacterium(BSim sim, BSimChemicalField production_field, BSimChemicalField consumption_field, Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        
        this.production_field = production_field;	
        this.consumption_field = consumption_field;
        
        production = true;
        productionNum = 1e3;
        consumptionRate = 1.0;//0.8;
        consumptionMax = 1e3;
        
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        initial_growth_rate = initial_growth_stdv * bacRng.nextGaussian() + initial_growth_mean;
        setK_growth(initial_growth_rate);
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();
        
        // Consumption field decays as it is consumed
        if ( consumption_field.getConc(position) > 0 ) {
        	double consumptionNum = consumption_field.getConc(position) * consumptionRate * sim.getDt();
        	
        	// The maximum amount of amino acids able to be consumed by a bacterium
        	if ( consumptionNum > consumptionMax ) {
        		consumption_field.addQuantity( position, -consumptionMax );
        	}
        	else {
        		consumption_field.addQuantity( position, -consumptionNum );
        	}
        	
    		// Cell growth rate is dependent on consumption 
            double growth = (mu_max * sim.getDt()) * ( consumptionNum / ( K_s + consumptionNum ) );
        	//double growth = (mu_max) * ( consumptionNum / ( K_s + consumptionNum ) );
            
            total_bio_mass += growth;
            
            // Yield conversion
            double growth_rate = (growth * total_bio_mass) / yield_coefficient;
            setK_growth(growth_rate);
            
            //System.out.println(total_bio_mass);
            System.out.println(consumptionNum + " " + growth_rate);
        }
		
		// Production of amino acid is dependent on nutrient field ( in access )
        // Production in this simulation is constant
		if ( isProducing() ) {						
			production_field.addQuantity(position, productionNum * sim.getDt() );	
		}
		
    }
    
    /** Returns the flag for the internal amino acid production of a bacterium. */
    public boolean isProducing() {
        return this.production;
    }
    
    /** Sets the flag for the internal amino acid production of a bacterium. */
    public void setProduction( boolean production) {
        this.production = production;
    }

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public CrossFeedingBacterium divide() {
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
        CrossFeedingBacterium child = new CrossFeedingBacterium(sim, production_field, consumption_field, x1_child, new Vector3d(this.x2));
														// Asymmetrical growth occurs at division node
        this.initialise(L1, x2_new, this.x1);			// Swap x1 and x2 for the mother after division for asymmetrical elongation
		child.L = L2;
        
        // add child to list of children - Sohaib Nadeem
        addChild(child);
                
        // Calculate angle between daughter cells at division 
        angle_initial = coordinate(child);

        // Prints a line whenever a new bacterium is made
        System.out.println("Child ID id " + child.id);
        return child;
    }


}
