package bsim.winter2021;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;

import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.List;

/**
 */
public class Bacterium extends BSimCapsuleBacterium {
    public final long origin_id;
    public long parent_id;
    public int lifetime;

    /** These fields are used to build a tree structure to keep track of lineage - Sohaib Nadeem */
    protected List<Bacterium> children;
    protected Bacterium parent;

    /** Scales randomVec during division **/
    public static double twist = 0.1; // 0.01
    /** Scales longVec during division **/
    public static double push = 0.05;

    /** Enables asymmetric growth. */
    static boolean asymmetric_growth = false;
    /** Length threshold for asymmetric growth (um). */
	static double L_asym = 3.75;
	/** Allows the cell to grow asymmetrically. */
	static double asymmetry = 0.1;
	/** The amount of force added back to achieve symmetric growth. */
	static double symmetry_adjustment = 0.05;

    //function you call when you want to make a new bacterium object
    public Bacterium(BSim sim, Vector3d px1, Vector3d px2){
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        this.origin_id = id;
        this.parent_id = -1;
        lifetime = 0;

        // initialize fields for lineage tree - Sohaib Nadeem
        children = new ArrayList<>();
        parent = null;
    }

    //function you call when you want to make a bacterium object that is a child of another
    public Bacterium(BSim sim, Vector3d px1, Vector3d px2, long origin_id, long parent_id){
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        this.origin_id = origin_id;
        this.parent_id = parent_id;
        lifetime = 0;

        // initialize fields for lineage tree - Sohaib Nadeem
        children = new ArrayList<>();
        parent = null;
    }
    
    /** Sets the value for asymmetric growth threshold. **/
    public void setLAsym(double length) {L_asym = length;}
    /** Sets the value of asymmetry. **/
    public void setAsym(double a) { asymmetry = a; }
    /** Sets the value of symmetric growth. **/
    public void setSym(double s) { symmetry_adjustment = s; }
    
    /** Sets the value of the twist during division. **/
    public void setTwist(double t) {this.twist = t;}
    /** Sets the value of the push during division. **/
    public void setPush(double p) {this.push = p;}
    
	@Override
    // Function which computes the internal spring force acting on the endpoints of the cell
    // We need the internal force to prevent the cell from lengthening as a result of forces acting individually
    // on the endpoints of the cell.
    // Responsible for the growth of the bacteria
	/** Updated to implement asymmetrical elongation. */
    public void computeSelfForce() {
        if (asymmetric_growth) {
            // create new number representing strength of the internal force -> initially set to zero
            double internalPotential = 0;

            // vector from x2 to x1
            Vector3d seg = new Vector3d();
            seg.sub(x2, x1);

            // checks whether there is a discrepancy between endpoint distance and "actual length"
            // if there is, this means that the internal force is needed to pull the endpoints back into position
            double lengthDiff = seg.length() - L;
            seg.normalize();

            if (lengthDiff < 0) {
                internalPotential = 0.5 * k_int * Math.pow(lengthDiff, 2);
            } else {
                internalPotential = -0.5 * k_int * Math.pow(lengthDiff, 2);
            }

            // Cell gradually starts growing symmetrically after the length threshold is met
            // and  as asym approaches 1.0
            if (L >= L_asym && asymmetry <= 1.0) {
                asymmetry += symmetry_adjustment * sim.getDt();
            }

            // Elongate asymmetrically until the length threshold is met
            // Internal potential is doubled to account for the x1force
            this.x1force.scaleAdd(-internalPotential * (2 - asymmetry), seg, this.x1force);
            this.x2force.scaleAdd(internalPotential * asymmetry, seg, this.x2force);
        } else {
            super.computeSelfForce();
        }
    }

    // in case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { //runs at every time step
    	//System.out.println("twist: " + twist);
    	//System.out.println("push: " + push);
        super.action();
    }

    // allows us to change growthrate mechanics for individual cells
    @Override
    public void setK_growth(double k_growth) {
        super.setK_growth(k_growth);
    }

    // allows us to change the division threshold of individual cells
    public void setElongationThreshold(double len) { this.L_th=len; }

    // for lineage tree - Sohaib Nadeem
    public void addChild(Bacterium child) {
        children.add(child);
        child.parent = this;
    }

	@Override
    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public Bacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble(), rng.nextDouble(), rng.nextDouble());
        randomVec.scale(twist);
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
        Bacterium child = new Bacterium(sim, x1_child, new Vector3d(this.x2), this.origin_id, this.id);
        this.parent_id = this.id;
        this.lifetime = 0;
        // this.initialise(L1, this.x1, x2_new); // for symmetric growth
        // Asymmetrical growth occurs at division node
        // so we need to swap x1 and x2 for the mother after division for asymmetrical elongation
        // This does not affect symmetric growth
        this.initialise(L1, x2_new, this.x1);
        child.L = L2;

        // add child to list of children - Sohaib Nadeem
        addChild(child);

        // Calculate angle between daughter cells at division
        angle_initial = coordinate(child);

        // Prints a line whenever a new bacterium is made
        System.out.println("Child ID is " + child.id);
        return child;
    }

}
