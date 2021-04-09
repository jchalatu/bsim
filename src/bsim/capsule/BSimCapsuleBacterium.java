package bsim.capsule;

import bsim.BSim;

import bsim.ode.BSimOdeSystem;
import bsim.winter2021.RectangleIntersection;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.Random;
import java.util.Timer;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicLong;
import java.lang.Math;

public class BSimCapsuleBacterium {

    // An "atomic long" is an object which produces a list of long integers
    // When you want the next integer in the list, you call .getAndIncrement()
    // This allows you to label each bacterium with an ID without having to know how many
    // bacteria there will be in total.
    // This variable is static. This means that it does not belong to an instance of
    // BSimCapsuleBacterium but it instead belongs to the class as a whole. Therefore
    // you do not need to call this.NEXT_ID.getAndIncrement(). You can just call
    // BSimCapsuleBacterium.NEXT_ID.getAndIncrement() without referring to a particular
    // bacterium.
    static final AtomicLong NEXT_ID = new AtomicLong(0);

    // If this is the first bacterium that is created, we will get id=0
    // If it's the second, you get id=1, and so on...
    public final long id = NEXT_ID.getAndIncrement();

    // used for generating random numbers using java's random number generator
    protected Random rng = new Random();

    // coordinate vectors for the endpoints of the bacterium, initilaized at 0,0,0
    public Vector3d x1 = new Vector3d(0,0,0);
    public Vector3d x2 = new Vector3d(0,0,0);

    // center of the bacterium. note that position = (x2-x1)/2, initialize
    public Vector3d position = new Vector3d(0,0,0);

    // force exerted on each endpoint of the bacterium, initialize
    public Vector3d x1force = new Vector3d(0,0,0);
    public Vector3d x2force = new Vector3d(0,0,0);

    // width of the bacterium (called radius because it also represents
    // how large the orbs are at each endpoint), all equal right now.
    public double radius = 1.0/2.0;

    // if no length is specified, this will set the initial length to 1 micrometer
    public double L_initial = 1;

    // actual length of bacterium
    public double L = L_initial;

    // maximum possible length of bacterium - it will never grow larger than this
    public double L_max = 15*L_initial;

    // how long the bacterium must be before there it divides
    public double L_th = 7; //make sure L_max is greater than L_th, default is 7
    
	/** Length threshold for asymmetric growth (um). */
	final public double L_asym = 5;

    // growth rate of bacterium -> increases length each second
    public double k_growth = 0.02;

    // true/false value representing whether the bacterium is infected or not
    // this can be renamed to "status" or anything else to represent a true/false property of a bacterium
    //public boolean infected;

    // function which allows you to set the infected status of a bacterium from another class
    /*public void setInfected(Boolean infected){
        this.infected = infected;
    }*/
	
    /** Angle between daughter cells at division. */
    public double angle_initial = 0.0; 
    
	/** Status of phage occupation within bacterium. */
	private boolean infected = false;

    // constructor for BSimCapsuleBacterium
    public BSimCapsuleBacterium(BSim _sim, Vector3d _x1, Vector3d _x2) {
        x1 = _x1;
        x2 = _x2;
        sim = _sim;
        setBrownianForceMagnitude(); //currently unused

        Vector3d u = new Vector3d();
        u.sub(this.x2, this.x1);
        this.position.scaleAdd(0.5, u, this.x1);
    }

    // function to allow you to set individual growth rates
    public void setK_growth(double k_growth) {
        this.k_growth = k_growth;
    }

    // returns growth rate of this cell (redundant)
    public double getK_growth() {
        return k_growth;
    }

    // checks whether cell is infected and returns 1 (infected) or 0 (not infected)
    // querying infection status
    public int getInfected() {
        if( infected == true ){
            return 1;
        }
        else
            return 0;
    }

    // Spring constants
    public double k_int = 50.0;     		// internal force (5e-5 from Storck)
    public double k_wall = 50.0 * 10;   			// overlap with wall
    public double k_cell = 50.0 * 10;   			// overlapping cells (1e-4 from Storck)
    public double k_filial = 5e-7;			// end to end attraction
    public double k_longfilial = 5e-7;		// opposite end repulsion when filial is active
    public double k_sticking = 0.01/*10.0*/;   		// side to side attraction
    public double k_wallstick = 0.00;  		// cell to wall attraction (0 right now since walls are leaky)

    // endpoint-to-endpoint ranges within which which these forces act, or (endpoint-to-wall distance)
    public double range_filial = 0.01; 		// maximum range at which filial forces activate
    public double range_sticking = 5.0/*0.6, 2.0*/; 	// max range ... sticking forces
    public double range_wallstick = 0.25; 	// max range ... attractive wall forces

    // contact range extension (pilus length) and contact threshold for sticking force (side to side attraction)
    public double contact_range_extension = 0.25;
    public double contact_threshold = 4.0;
    
    private double EPS = 1e-12; // machine precision value for use in checking whether quantity is zero
    // You can use it like this: value x is effectively zero if Math.abs(x)<EPS
    // so instead of writing if(x==0) you should write if(Math.abs(x)<=EPS) to ensure precision

    // stores data about the simulation so that each cell can use that information
    protected BSim sim;

    // function which is called at each time step. currently just updates the position of the cell
    // just defines the midpoint
    public void action(){
        Vector3d u = new Vector3d();
        u.sub(this.x2, this.x1);
        this.position.scaleAdd(0.5, u, this.x1);
    }

    // if this is set to a positive number, the bacteria will wiggle around due to brownian forces
    // see Brownian Force on youtube
    public double brownianForceMagnitude;

    // Function which is called at each time step. This updates the length of the bacterium according
    // to a logistic growth model(how much it grows each time step)
    // dL/dt = kL(1-L/Lmax)
    // explicit Euler step for logistic: could be linear
    public void grow() {
        L = L + sim.getDt() * this.k_growth * L * (1 - (L / L_max));
    }

    // function called at the very beginning of the simulation.
    // This sets the initial position and orientation of the cell.
    public void initialise(double _L, Vector3d _x1, Vector3d _x2){
        this.L = _L;
        this.x1 = new Vector3d(_x1);
        this.x2 = new Vector3d(_x2);
        Vector3d u = new Vector3d();
        u.sub(this.x2, this.x1);
        Vector3d _pos = new Vector3d();
        _pos.sub(_x2,_x1);
        _pos.scaleAdd(0.5,u,_x1);
        this.position=_pos;

        //this.position.scaleAdd(0.5, u, this.x1);
    }

    // we dont use this divide function, we use the divide function in bacterium
    public BSimCapsuleBacterium divide() {

        Vector3d u = new Vector3d(); u.sub(this.x2, this.x1); // vector pointing from x1 to x2

        // picks a random point on the cell (near the middle) where the division will start
        double divPert = 0.1*L*(rng.nextDouble() - 0.5);

        // checks length of cell
        double L_actual = u.length();

        // computes the length of the new daughter cells.
        double L1 = L_actual*0.5*(1 + divPert) - radius;
        double L2 = L_actual*0.5*(1 - divPert) - radius;

        // computes the positions of the daughter cells
        Vector3d x2_new = new Vector3d();
        x2_new.scaleAdd(L1/L_actual, u, this.x1);
        x2_new.add(new Vector3d(0.05*L_initial*(rng.nextDouble() - 0.5),
                                0.05*L_initial*(rng.nextDouble() - 0.5),
                                0.05*L_initial*(rng.nextDouble() - 0.5)));
        // the 0.05*L_initial part here just causes the new cells to be pushed a
        // little bit away from each other to avoid collisions

        Vector3d x1_child = new Vector3d();
        x1_child.scaleAdd(-(L2/L_actual), u, this.x2);
        x1_child.add(new Vector3d(0.05*L_initial*(rng.nextDouble() - 0.5),
                                  0.05*L_initial*(rng.nextDouble() - 0.5),
                                  0.05*L_initial*(rng.nextDouble() - 0.5)));

        // creates one new bacterium object (the "child" cell)
        BSimCapsuleBacterium child = new BSimCapsuleBacterium(sim, x1_child, new Vector3d(this.x2));

        // resets current cell (warning from original author: this MUST be done AFTER creating the child)
        this.initialise(L1, this.x1, x2_new);

        // sets length of child cell and returns it so that it may be added to the list of cells
        child.L = L2;
        return child;
    }

    // function which computes the internal spring force acting on the endpoints of the cell
    // We need the internal force to prevent the cell from lengthening as a result of forces acting individually
    // on the endpoints of the cell.
    // This section has not been modified from the original bsim package
    public void computeSelfForce() {
        // create new number representing strength of the internal force -> initially set to zero
        double internalPotential = 0;

        // vector from x2 to x1
        Vector3d seg = new Vector3d();
        seg.sub(x2, x1);

        // checks whether there is a discrepancy between endpoint distance and "actual length"
        // if there is, this means that the internal force is needed to pull the endpoints back into position
        double lengthDiff = seg.length() - L;
        seg.normalize();

        if(lengthDiff < 0) {
            internalPotential = 0.5 * k_int * Math.pow(lengthDiff, 2);
        }
        else {
            internalPotential = -0.5 * k_int * Math.pow(lengthDiff, 2);
        }

        this.x1force.scaleAdd(-internalPotential, seg, this.x1force);
        this.x2force.scaleAdd(internalPotential, seg, this.x2force);
    }

    // computes the force acting on the cell from walls below the cell (wall overlap)
    // not sorted out yet
    public void wallBelow(double pointCoord, Vector3d theForce, Vector3d forceDir){
        if(radius - pointCoord > 0){
            theForce.scaleAdd(0.4*k_wall*Math.pow(radius - pointCoord, 2.5), forceDir, theForce);
        }
    }
    // computes force acting on cell from walls above
    public void wallAbove(double pointCoord, Vector3d theForce, Vector3d forceDir, double simBound){
        if(radius + pointCoord - simBound > 0){
            theForce.scaleAdd(0.4*k_wall*Math.pow(radius + pointCoord - simBound, 2.5), forceDir, theForce);
        }
    }

    // checks whether the wall below the cell is the boundary of the simulation, or if it is the wall of a channel
    // within a larger region
    public boolean flowBelow(double pointCoord, Vector3d theForce, Vector3d forceDir, double simBound){
        if(radius - pointCoord > 0){
            theForce.scaleAdd(Math.abs(radius - pointCoord), forceDir, theForce);
            return true;
        }
        else return false;
    }

    // similar to previous
    public boolean flowAbove(double pointCoord, Vector3d theForce, Vector3d forceDir, double simBound){
        // Ignore the radius for now... (COMMENT ADDED BY ORIGINAL AUTHOR - what does it mean?)
        if(radius + pointCoord - simBound > 0){
            theForce.scaleAdd(Math.abs(radius + pointCoord - simBound), forceDir, theForce);
            return true;
        }
        else return false;
    }

    // computes force on bacterium from background fluid motion -> seems like this code is unfinished
    // notice that later in the file, this force is not used. if background flows are considered later, this should be modified
    public void computeFlowForce() {
        // TEST - apply velocity on BOTTOM
        flowAbove(x1.y, x1force, new Vector3d(0.5, 0, 0), sim.getBound().y);
        flowAbove(x2.y, x2force, new Vector3d(0.5, 0, 0), sim.getBound().y);

        // TEST - apply velocity on TOP
        flowBelow(x1.y, x1force, new Vector3d(0.5, 0, 0), sim.getBound().y);
        flowBelow(x2.y, x2force, new Vector3d(0.5, 0, 0), sim.getBound().y);
    }
	
    /*
    COMMENT ADDED BY ORIGINAL AUTHOR
    TODO - need proper wall contact location computation to implement flow past.
    TODO - for each boundary, apply a BoundaryCondition interface with BoundaryCondition.resolve() or something.
    - need a point of application so that force acts at both x1 and x2, or for a torque to act on the cell
     */
    // It seems like the original authors did not correctly account for cases where the cell may come in contact
    // with the corner of a wall. In this case, the corner should apply a torque to the cell based on where it collides.
    // This seems complicated.
    // For now, this function should only be used to compute forces between flat walls and a bacterium
	/** Updated to allow for different wall boundary conditions. */
    public void computeWallForce(){
        // TODO::: Ideally, there should also be a bounds check on the side NEXT to the one from which bacs can exit
        /**
         * i.e.,
         *
         * open, flow - - - - - - - ->
         *            |            |  should have a bounds check here @ top so that bacs being pushed by the 'flow'
         *  closed    |            |  are allowed to continue moving right, above the RHS wall, rather than being
         *            .            .  *stopped* by the RHS bound check!
         *
         */
    	
    	if ( sim.wall_boundaries[0] ) {
    		wallBelow(x1.x, x1force, new Vector3d(1,0,0));
    		wallBelow(x2.x, x2force, new Vector3d(1,0,0));
    	}
    	if ( sim.wall_boundaries[1] ) {
    		wallBelow(x1.y, x1force, new Vector3d(0,1,0)); // TOP //
    		wallBelow(x2.y, x2force, new Vector3d(0,1,0)); // TOP //
    	}
    	if ( sim.wall_boundaries[2] ) {
    		wallBelow(x1.z, x1force, new Vector3d(0,0,1));
    		wallBelow(x2.z, x2force, new Vector3d(0,0,1));
    	}
        
        if ( sim.wall_boundaries[3] ) {
        	wallAbove(x1.x, x1force, new Vector3d(-1,0,0), sim.getBound().x);
        	wallAbove(x2.x, x2force, new Vector3d(-1,0,0), sim.getBound().x);
        }
        if ( sim.wall_boundaries[4] ) {
        	wallAbove(x1.y, x1force, new Vector3d(0, -1, 0), sim.getBound().y); // BOTTOM //
        	wallAbove(x2.y, x2force, new Vector3d(0, -1, 0), sim.getBound().y); // BOTTOM //
        }
        if ( sim.wall_boundaries[5] ) {
        	wallAbove(x1.z, x1force, new Vector3d(0, 0, -1), sim.getBound().z);
        	wallAbove(x2.z, x2force, new Vector3d(0, 0, -1), sim.getBound().z);
        }
        
        ///////////////////////////Catie Terrey
        // added forces from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4017289/

        // this accounts for attractive forces between bacteria and walls due to hydrostatic interactions
        if(this.x2.y < radius + range_wallstick){
            double wallForce=-k_wallstick*(this.x2.y-radius-0.01); //spring attaching cell to wall
            Vector3d wallForceVec = new Vector3d(0,1,0); //vector along which this force will act
            this.x2force.scaleAdd(wallForce,wallForceVec,this.x2force);
        }
        else if (sim.getBound().y-this.x2.y < radius + range_wallstick){
            double wallForce=k_wallstick*(sim.getBound().y-this.x2.y -radius+0.01);
            Vector3d wallForceVec = new Vector3d(0,1,0);
            this.x2force.scaleAdd(wallForce,wallForceVec,this.x2force);
        }
        if(this.x1.y < radius + range_wallstick){
            double wallForce=-k_wallstick*(this.x1.y-radius-0.01);
            Vector3d wallForceVec = new Vector3d(0,1,0);
            this.x1force.scaleAdd(wallForce,wallForceVec,this.x1force);
        }
        else if (sim.getBound().y-this.x1.y < radius + range_wallstick){
            double wallForce=k_wallstick*(sim.getBound().y-this.x1.y -radius+0.01);
            Vector3d wallForceVec = new Vector3d(0,1,0);
            this.x1force.scaleAdd(wallForce,wallForceVec,this.x1force);
        }
        if(this.x2.x < radius + range_wallstick){
            double wallForce=-k_wallstick*(this.x2.x-radius-0.01);
            Vector3d wallForceVec = new Vector3d(1,0,0);
            this.x2force.scaleAdd(wallForce,wallForceVec,this.x2force);
        }
        else if (sim.getBound().x-this.x2.x < radius + range_wallstick){
            double wallForce=k_wallstick*(sim.getBound().x-this.x2.x-radius+0.01);
            Vector3d wallForceVec = new Vector3d(1,0,0);
            this.x2force.scaleAdd(wallForce,wallForceVec,this.x2force);
        }
        if(this.x1.x < radius + range_wallstick){
            double wallForce=-k_wallstick*(this.x1.x-radius-0.01);
            Vector3d wallForceVec = new Vector3d(1,0,0);
            this.x1force.scaleAdd(wallForce,wallForceVec,this.x1force);
        }
        else if (sim.getBound().x-this.x1.x < radius + range_wallstick){
            double wallForce=k_wallstick*(sim.getBound().x-this.x1.x -radius+0.01);
            Vector3d wallForceVec = new Vector3d(1,0,0);
            this.x1force.scaleAdd(wallForce,wallForceVec,this.x1force);
        }
        // no forces on the z-axis?
    }

    private double getNeighbourDist(BSimCapsuleBacterium neighbour_bac) {
        double neighbourDist;

        Vector3d u = new Vector3d(); 					// define u as a vector
        u.sub(this.x2, this.x1); 						// u = x2 - x1 = vector pointing from one endpoint to the other
        Vector3d v = new Vector3d(); 					// define v as a vector
        v.sub(neighbour_bac.x2, neighbour_bac.x1); 		// v = x2 - x1 (neighbour cell)
        Vector3d w = new Vector3d(); 					// define w as a vector
        w.sub(this.x1, neighbour_bac.x1); 				// difference between the two x1 points

        double a = u.dot(u);         // always >= 0
        double b = u.dot(v);
        double c = v.dot(v);         // always >= 0
        double d = u.dot(w);
        double e = v.dot(w);
        double D = a * c - b * b; 	// figure out what D is by visualizing what it is     // always >= 0
        double sc = D;
        double sN = D;
        double sD = D;       		// sc = sN / sD, default sD = D >= 0
        double tc = D;
        double tN = D;
        double tD = D;       		// tc = tN / tD, default tD = D >= 0

        // compute the line parameters of the two closest points
        if (D < EPS) { 				// the lines are almost parallel
            sN = 0.0;         		// force using point P0 on segment S1
            sD = 1.0;         		// to prevent possible division by 0.0 later
            tN = e;
            tD = c;
        } else {                 	// get the closest points on the infinite lines
            sN = (b * e - c * d);
            tN = (a * e - b * d);
            if (sN < 0.0) {        	// sc < 0 => the s=0 edge is visible
                sN = 0.0;
                tN = e;
                tD = c;
            } else if (sN > sD) {  	// sc > 1  => the s=1 edge is visible
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN < 0.0) {            	// tc < 0 => the t=0 edge is visible
            tN = 0.0;
            // recompute sc for this edge
            if (-d < 0.0)
                sN = 0.0;
            else if (-d > a)
                sN = sD;
            else {
                sN = -d;
                sD = a;
            }
        } else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
            tN = tD;
            // recompute sc for this edge
            if ((-d + b) < 0.0)
                sN = 0;
            else if ((-d + b) > a)
                sN = sD;
            else {
                sN = (-d + b);
                sD = a;
            }
        }

        // finally do the division to get sc and tc
        sc = (Math.abs(sN) < EPS ? 0.0 : sN / sD);
        tc = (Math.abs(tN) < EPS ? 0.0 : tN / tD);

        Vector3d dP = new Vector3d(w);
        dP.scaleAdd(sc, u, dP);
        dP.scaleAdd(-tc, v, dP);

        neighbourDist = dP.length();

        return neighbourDist;
    }

    private double getBoundingBoxIntersectionArea(BSimCapsuleBacterium neighbour_bac, double r) {
        double bac_radius = radius + r;
        double neighbour_bac_radius = neighbour_bac.radius + r;
        /*
        double bac_x = Math.min(x1.x, x2.x) - bac_radius;
        double bac_y = Math.max(x1.y, x2.y) + bac_radius;
        double bac_w = Math.abs(x1.x - x2.x) + 2 * bac_radius;
        double bac_h = Math.abs(x1.y - x2.y) + 2 * bac_radius;
         */
        Rectangle2D.Double bac_bbox = new Rectangle2D.Double(Math.min(x1.x, x2.x) - bac_radius, Math.max(x1.y, x2.y) + bac_radius,
                Math.abs(x1.x - x2.x) + 2 * bac_radius, Math.abs(x1.y - x2.y) + 2 * bac_radius);
        Rectangle2D.Double neighbour_bac_bbox = new Rectangle2D.Double(Math.min(neighbour_bac.x1.x, neighbour_bac.x2.x) - neighbour_bac_radius,
                Math.max(neighbour_bac.x1.y, neighbour_bac.x2.y) + neighbour_bac_radius,
                Math.abs(neighbour_bac.x1.x - neighbour_bac.x2.x) + 2 * neighbour_bac_radius,
                Math.abs(neighbour_bac.x1.y - neighbour_bac.x2.y) + 2 * neighbour_bac_radius);
        Rectangle2D intersection = bac_bbox.createIntersection(neighbour_bac_bbox);
        double area = intersection.getHeight() * intersection.getWidth();
        return (area > 0 ? area : 0);
    }

    private boolean getBoundingBoxIntersection(BSimCapsuleBacterium neighbour_bac, double r) {
        double bac_radius = radius + r;
        double neighbour_bac_radius = neighbour_bac.radius + r;

        Rectangle2D.Double bac_bbox = new Rectangle2D.Double(Math.min(x1.x, x2.x) - bac_radius, Math.max(x1.y, x2.y) + bac_radius,
                Math.abs(x1.x - x2.x) + 2 * bac_radius, Math.abs(x1.y - x2.y) + 2 * bac_radius);
        Rectangle2D.Double neighbour_bac_bbox = new Rectangle2D.Double(Math.min(neighbour_bac.x1.x, neighbour_bac.x2.x) - neighbour_bac_radius,
                Math.max(neighbour_bac.x1.y, neighbour_bac.x2.y) + neighbour_bac_radius,
                Math.abs(neighbour_bac.x1.x - neighbour_bac.x2.x) + 2 * neighbour_bac_radius,
                Math.abs(neighbour_bac.x1.y - neighbour_bac.x2.y) + 2 * neighbour_bac_radius);
        return bac_bbox.intersects(neighbour_bac_bbox);
    }

    private double getOrientedBoundingBoxIntersectionPerimeter(BSimCapsuleBacterium neighbour_bac, double r) {
        double bac_length = Math.sqrt(Math.pow(x2.x - x1.x, 2) + Math.pow(x2.y - x1.y, 2));
        double neighbour_bac_length = Math.sqrt(Math.pow(neighbour_bac.x2.x - neighbour_bac.x1.x, 2)
                + Math.pow(neighbour_bac.x2.y - neighbour_bac.x1.y, 2));
        double bac_angle = Math.atan((x2.y - x1.y) / (x2.x - x1.x));
        double neighbour_bac_angle = Math.atan((neighbour_bac.x2.y - neighbour_bac.x1.y) / (neighbour_bac.x2.x - neighbour_bac.x1.x));

        Vector2d[] rectangle_bac = RectangleIntersection.rectangle_vertices(position.x, position.y,
                bac_length + 2 * (radius + r), 2 * (radius + r), bac_angle);
        Vector2d[] rectangle_neighbour = RectangleIntersection.rectangle_vertices(neighbour_bac.position.x,
                neighbour_bac.position.y, neighbour_bac_length + 2 * (neighbour_bac.radius + r),
                2 * (neighbour_bac.radius + r), neighbour_bac_angle);

        return RectangleIntersection.rectangle_intersection_perimeter(rectangle_bac, rectangle_neighbour);
    }

    // computes force between neighboring cells
    // originally htis was just overlap
    // Catie add filial and sticking forces
    public void computeNeighbourForce(BSimCapsuleBacterium neighbour_bac){

        double neighbourDist = 0.0;
        boolean neighbourDist_calculated = false;

        ////////////////overlap first////////////

        // this first section computes the force due to collisions between cells
        /*
        COMMENT FROM ORIGINAL AUTHOR
        Vector returned (dP) is from the second bac, heading to the first
        Geometric Tools for Computer Graphics book,
        Also,
        http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
        Thank you!
        */


    	// Calculates the distance between two lines (see link above)
        Vector3d dist = new Vector3d(); 
        dist.sub(this.position, neighbour_bac.position); // position = midpoint of cell

        double rDist = (this.L + neighbour_bac.L)*0.5 + (this.radius + neighbour_bac.radius);
        Vector3d u = new Vector3d(); 					// define u as a vector
        u.sub(this.x2, this.x1); 						// u = x2 - x1 = vector pointing from one endpoint to the other
        Vector3d v = new Vector3d(); 					// define v as a vector
        v.sub(neighbour_bac.x2, neighbour_bac.x1); 		// v = x2 - x1 (neighbour cell)
        Vector3d w = new Vector3d(); 					// define w as a vector
        w.sub(this.x1, neighbour_bac.x1); 				// difference between the two x1 points
        if(dist.dot(dist) < rDist*rDist) {
            double a = u.dot(u);         // always >= 0
            double b = u.dot(v);
            double c = v.dot(v);         // always >= 0
            double d = u.dot(w);
            double e = v.dot(w);
            double D = a * c - b * b; 	// figure out what D is by visualizing what it is     // always >= 0
            double sc = D;
            double sN = D;
            double sD = D;       		// sc = sN / sD, default sD = D >= 0
            double tc = D;
            double tN = D;
            double tD = D;       		// tc = tN / tD, default tD = D >= 0

            // compute the line parameters of the two closest points
            if (D < EPS) { 				// the lines are almost parallel
                sN = 0.0;         		// force using point P0 on segment S1
                sD = 1.0;         		// to prevent possible division by 0.0 later
                tN = e;
                tD = c;
            } else {                 	// get the closest points on the infinite lines
                sN = (b * e - c * d);
                tN = (a * e - b * d);
                if (sN < 0.0) {        	// sc < 0 => the s=0 edge is visible
                    sN = 0.0;
                    tN = e;
                    tD = c;
                } else if (sN > sD) {  	// sc > 1  => the s=1 edge is visible
                    sN = sD;
                    tN = e + b;
                    tD = c;
                }
            }

            if (tN < 0.0) {            	// tc < 0 => the t=0 edge is visible
                tN = 0.0;
                // recompute sc for this edge
                if (-d < 0.0)
                    sN = 0.0;
                else if (-d > a)
                    sN = sD;
                else {
                    sN = -d;
                    sD = a;
                }
            } else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
                tN = tD;
                // recompute sc for this edge
                if ((-d + b) < 0.0)
                    sN = 0;
                else if ((-d + b) > a)
                    sN = sD;
                else {
                    sN = (-d + b);
                    sD = a;
                }
            }

            // finally do the division to get sc and tc
            sc = (Math.abs(sN) < EPS ? 0.0 : sN / sD);
            tc = (Math.abs(tN) < EPS ? 0.0 : tN / tD);

            Vector3d dP = new Vector3d(w);			
            dP.scaleAdd(sc, u, dP);					
            dP.scaleAdd(-tc, v, dP);

            neighbourDist = dP.length();
            neighbourDist_calculated = true;

            // If a collision occurs
            if (neighbourDist < 2 * radius) {
                // Try this; if necessary we can simplify to a linear force
                double repulsionStrength = 0.4 * k_cell * Math.pow(2 * radius - neighbourDist, 2.5);
                // why is this to the 2.5 power?

/////////////////////////////////////////////////////////////////////////////////////
                    /*
                    OK, this section is not right.
                    We want the projection of dP (or 0.5*dP???)
                    onto the vector from intersection point to x1 (or x2)
                    ... Actually, not so sure. I think we could use the dP weighted by inverse distance
                    between the intersection point and x1/x2 as in Storck et al.
                    ... Or the single sphere approximation from Volfson.
                    Test + compare the alternatives if this doesn't work.
                    */
//                            Vector3d x1comp = new Vector3d(0,0,0);
//                            Vector3d x2comp = new Vector3d(0,0,0);

//                            x1comp.scale(0.5, dP);
//                            x2comp.scale(0.5, dP);
//
//                            x1comp.scaleAdd(-sc, u, x1comp);
//                            x2comp.scaleAdd(1-sc, u, x2comp);
/////////////////////////////////////////////////////////////////////////////////////

                // ***********
                // *********** NEW
                // ***********
//            Vector3d halfDp = new Vector3d(dP);
//            halfDp.scale(0.5);
//
//            Vector3d contactTox1 = new Vector3d();
//            contactTox1.scaleAdd(-sc, u, halfDp);
//
//            Vector3d contactTox2 = new Vector3d();
//            contactTox2.scaleAdd(1-sc, u, halfDp);
//
//            Vector3d contactToNeighbourx1 = new Vector3d(halfDp);
//            contactToNeighbourx1.negate();
//            contactToNeighbourx1.scaleAdd(-tc, v, contactToNeighbourx1);
//
//            Vector3d contactToNeighbourx2 = new Vector3d(halfDp);
//            contactToNeighbourx2.negate();
//            contactToNeighbourx2.scaleAdd(1-tc, v, contactToNeighbourx2);
//
////            double rsx1 = repulsionStrength/Math.pow(halfDp.length() + contactTox1.length(), 2);
////            double rsx2 = repulsionStrength/Math.pow(halfDp.length() + contactTox2.length(), 2);
////            double rsnx1 = repulsionStrength/Math.pow(halfDp.length() + contactToNeighbourx1.length(), 2);
////            double rsnx2 = repulsionStrength/Math.pow(halfDp.length() + contactToNeighbourx2.length(), 2);
//
//            double rsx1 = repulsionStrength*(u.length() + 2*radius - contactTox1.length())/(u.length() + 2*radius);
//            double rsx2 = repulsionStrength*(u.length() + 2*radius - contactTox2.length())/(u.length() + 2*radius);
//            double rsnx1 = repulsionStrength*(v.length() + 2*radius - contactToNeighbourx1.length())/(v.length() + 2*radius);
//            double rsnx2 = repulsionStrength*(v.length() + 2*radius - contactToNeighbourx2.length())/(v.length() + 2*radius);
//
////            System.out.format("%f, %f, %f, %f %n", rsx1, rsx2, rsnx1, rsnx2);
//            contactTox1.normalize();
//            contactTox2.normalize();
//            contactToNeighbourx1.normalize();
//            contactToNeighbourx2.normalize();
//
//            dP.normalize();
//
////            System.out.println(repulsionStrength/contactTox1.length() + "; " +
////                    repulsionStrength/contactTox2.length() + "; " +
////                    repulsionStrength/contactToNeighbourx1.length() + "; " +
////                    repulsionStrength/contactToNeighbourx2.length());
//
//            this.x1force.scaleAdd(rsx1, dP, this.x1force);
//            this.x2force.scaleAdd(rsx2, dP, this.x2force);
//
//            neighbour_bac.x1force.scaleAdd(-rsnx1, dP, neighbour_bac.x1force);
//            neighbour_bac.x2force.scaleAdd(-rsnx2, dP, neighbour_bac.x2force);
//
////            this.x1force.scaleAdd(rsx1, contactTox1, this.x1force);
////            this.x2force.scaleAdd(rsx2, contactTox2, this.x2force);
////
////            neighbour_bac.x1force.scaleAdd(rsnx1, contactToNeighbourx1, neighbour_bac.x1force);
////            neighbour_bac.x2force.scaleAdd(rsnx2, contactToNeighbourx2, neighbour_bac.x2force);
                
                // ***********
                // *********** /NEW
                // ***********

                // ***********
                // *********** OLD
                // ***********
                dP.normalize();

                this.x1force.scaleAdd((1.0 - sc) * repulsionStrength, dP, this.x1force);
                this.x2force.scaleAdd(sc * repulsionStrength, dP, this.x2force);

                //final calculation of overlap forces
                neighbour_bac.x1force.scaleAdd(-(1.0 - tc) * repulsionStrength, dP, neighbour_bac.x1force);
                neighbour_bac.x2force.scaleAdd(-tc * repulsionStrength, dP, neighbour_bac.x2force);

//                this.x1force.scaleAdd((1.1 - sc) * repulsionStrength, dP, this.x1force);
//                this.x2force.scaleAdd((sc + 0.1) * repulsionStrength, dP, this.x2force);
//
//                neighbour_bac.x1force.scaleAdd(-(1.1 - tc) * repulsionStrength, dP, neighbour_bac.x1force);
//                neighbour_bac.x2force.scaleAdd(-(tc + 0.1) * repulsionStrength, dP, neighbour_bac.x2force);

                // ***********
                // *********** /OLD
                // ***********
            }

        }
        // END OF SECTION

        /////////////Catie's new forces/////////
        // this next section computes the force between neighboring cells as a result of cell-cell attraction
        // this can happen due to membrane proteins as well as hydrostatic interactions or friction
        //distances between enpoints
        Vector3d dis11=new Vector3d();
        dis11.sub(this.x1,neighbour_bac.x1);
        Vector3d dis22=new Vector3d();
        dis22.sub(this.x2,neighbour_bac.x2);
        Vector3d dis12=new Vector3d();
        dis12.sub(this.x1,neighbour_bac.x2);
        Vector3d dis21=new Vector3d();
        dis21.sub(this.x2,neighbour_bac.x1); //displacement between endpoints of cell and neighbour

        double d11=dis11.length();
        double d22=dis22.length();
        double d21=dis21.length();
        double d12=dis12.length(); //distance between endpoints of cell and neighbour

        //double angleuv=Math.acos(u.dot(v)/(u.length()*v.length())); //angle between cell and its neighbour
        double angleuv = u.angle(v);
        //used as a condition on filial force extistence

        double thisLength=this.L; //length of current cell
        double neighborLength=neighbour_bac.L; //length of neighbour

        Vector3d filialAxis=new Vector3d(); // vector from end of cell to end of cell
        Vector3d longFilialAxis = new Vector3d();

        double minimumAngle = 7 * Math.PI/8; //cells must be oriented to within pi/8 rads (22.5 degrees) of each other to stick end to end
        //filial links (secondary structure -> causes chains to form)
        //might want to include condition on being sisters
        // how to tell if two cells are sisters??
        //4 conditions (pairs of endpoints) to see if any are sufficnetly close and right angle
        // comparing x1 to x1. checking if this is the smallest of all 4, and is lesx than filial range
        // and angle is less than 2.5 degrees
        
        if (d11 < d22 && d11 < d12 && d11 < d21 && d11 < 2 * radius + range_filial && angleuv >= minimumAngle){
            //if the closest two points are x11 and x12
            filialAxis=dis11; // set axis to x11-x12
            double stickyForceStrength=-k_filial*(d11-2 * radius-0.05*radius); //spring attaching endpoints
            filialAxis.normalize(); // normalize axis vector to get correct scaling
            //this.x1force.scaleAdd(stickyForceStrength,filialAxis,this.x1force); // adds the force to the endpoint
            //neighbour_bac.x1force.scaleAdd(-stickyForceStrength,filialAxis,neighbour_bac.x1force);

            //also a spring attaching furthest endpoints (but repulsive!) -> maybe this should be weaker?
            longFilialAxis=dis22;
            double longFilialForceStrength=k_longfilial*(d22-2*radius-thisLength-neighborLength-0.05*radius);
            longFilialAxis.normalize();
            
            //this.x2force.scaleAdd(longFilialForceStrength,longFilialAxis,this.x2force);
            //neighbour_bac.x2force.scaleAdd(-longFilialForceStrength,longFilialAxis,neighbour_bac.x2force);
        }
        else if(d22 < d11 && d22 < d12 && d22 < d21 && d22 < 2 * radius + range_filial && angleuv >= minimumAngle){
            //if the closest two points are x21 and x22
            filialAxis=dis22; // set axis to x21-x22
            double stickyForceStrength=-k_filial*(d22-2 * radius-0.05*radius);
            filialAxis.normalize();
            //this.x2force.scaleAdd(stickyForceStrength,filialAxis,this.x2force);
            //neighbour_bac.x2force.scaleAdd(-stickyForceStrength,filialAxis,neighbour_bac.x2force);

            longFilialAxis=dis11;
            double longFilialForceStrength=k_longfilial*(d11-2*radius-thisLength-neighborLength-0.05*radius);
            longFilialAxis.normalize();
            //this.x2force.scaleAdd(longFilialForceStrength,longFilialAxis,this.x2force);
            //neighbour_bac.x2force.scaleAdd(-longFilialForceStrength,longFilialAxis,neighbour_bac.x2force);
        }
        else if(d12 < d11 && d12 < d21 && d12 < d22 && d12 < 2 * radius + range_filial && (Math.PI - angleuv) >= minimumAngle){
            //if the closest two points are x11 and x22
            filialAxis=dis12; // set axis to x11-x22
            double stickyForceStrength=-k_filial*(d12-2 * radius -0.05*radius);
            filialAxis.normalize();
            //this.x1force.scaleAdd(stickyForceStrength,filialAxis,this.x1force);
            //neighbour_bac.x2force.scaleAdd(-stickyForceStrength,filialAxis,neighbour_bac.x2force);

            longFilialAxis=dis21;
            double longFilialForceStrength=k_longfilial*(d21-2*radius-thisLength-neighborLength-0.05*radius);
            longFilialAxis.normalize();
            //this.x2force.scaleAdd(longFilialForceStrength,longFilialAxis,this.x2force);
            //neighbour_bac.x2force.scaleAdd(-longFilialForceStrength,longFilialAxis,neighbour_bac.x2force);
        }
        else if(d21 < d11 && d21 < d12 && d21 < d22 && d21 < 2 * radius + range_filial && (Math.PI - angleuv) >= minimumAngle){
            //if the closest two points are x21 and x12
            filialAxis=dis21; // set axis to x21-x12
            double stickyForceStrength=-k_filial*(d21-2 * radius -0.05*radius);
            filialAxis.normalize();
            //this.x2force.scaleAdd(stickyForceStrength,filialAxis,this.x2force);
            //neighbour_bac.x1force.scaleAdd(-stickyForceStrength,filialAxis,neighbour_bac.x1force);

            longFilialAxis=dis12;
            double longFilialForceStrength=k_longfilial*(d12-2*radius-thisLength-neighborLength-0.05*radius);
            longFilialAxis.normalize();
            //this.x2force.scaleAdd(longFilialForceStrength,longFilialAxis,this.x2force);
            //neighbour_bac.x2force.scaleAdd(-longFilialForceStrength,longFilialAxis,neighbour_bac.x2force);
        }

        // sticking forces: end to end and across
        //sticking links (tertiary structure -> changes cell density)
        // these bonds happen between cells when they are side-by-side
        if (dist.length()<range_sticking) { //dist is distance from midpoint to midpoint
            //System.out.println(dist.length());
            double maxCellLength = Math.max(thisLength, neighborLength); // checks which cell is longer

            // computes rest length of diagonally oriented springs -> sqrt((2r)^2 + L^2)
            double stickingRestLong = Math.pow(4 * radius * radius + maxCellLength * maxCellLength, 0.5);

            // computes rest length of springs between closest endpoints
            double stickingRestShort = 2.1 * radius; //dont let this be too close to 2.0 or the simulation will chug due to collisions

            Vector3d axis11 = dis11;
            axis11.normalize();
            Vector3d axis22 = dis22;
            axis22.normalize();
            Vector3d axis12 = dis12;
            axis12.normalize();
            Vector3d axis21 = dis21;
            axis21.normalize(); // gets directions for forces to point -> normalized displacement L/|L| as in paper

            double dot_product = u.dot(v);

            //double neighbourDist = getNeighbourDist(neighbour_bac);
            //double contact_range = 2.0;
            //double orientation_factor = 1/*Math.abs(dot_product / (u.length() * v.length()*/;
            //double dampening = 1 - ( orientation_factor * Math.pow(2 * radius + contact_range - neighbourDist, 2) / Math.pow(contact_range, 2) );
            //double dampening = 1.0;

            double damping = 1.0;
            if (!neighbourDist_calculated) {
                neighbourDist = getNeighbourDist(neighbour_bac);
            }
            if (neighbourDist < 2 * (radius + contact_range_extension) /*getBoundingBoxIntersection(neighbour_bac, contact_range_extension)*/) {
                double contact_area = getOrientedBoundingBoxIntersectionPerimeter(neighbour_bac, contact_range_extension);
                //System.out.println(contact_area + " " + this.L + " " + neighbour_bac.L);
                damping = Math.min(1.0, 1 - contact_area / contact_threshold);
                //damping = 1.0;
            }

            if (dot_product > 0 /*d11<d12 && d11 < d21*/){
                // if cells are oriented same direction
                double strength11 = -k_sticking*(d11-stickingRestShort) * damping;
                this.x1force.scaleAdd(strength11,axis11,this.x1force);
                neighbour_bac.x1force.scaleAdd(-strength11,axis11,neighbour_bac.x1force);

                double strength22 = -k_sticking*(d22-stickingRestShort) * damping;
                this.x2force.scaleAdd(strength22,axis22,this.x2force);
                neighbour_bac.x2force.scaleAdd(-strength22,axis22,neighbour_bac.x2force);

                double strength12 = -k_sticking*(d12-stickingRestLong) * damping;
                this.x1force.scaleAdd(strength12,axis12,this.x1force);
                neighbour_bac.x2force.scaleAdd(-strength12,axis12,neighbour_bac.x2force);

                double strength21 = -k_sticking*(d21-stickingRestLong) * damping;
                this.x2force.scaleAdd(strength21,axis21,this.x2force);
                neighbour_bac.x1force.scaleAdd(-strength21,axis21,neighbour_bac.x1force);
            }
            else {
                // if cells are oriented opposite direction
                double strength11 = -k_sticking*(d11-stickingRestLong) * damping;
                this.x1force.scaleAdd(strength11,axis11,this.x1force);
                neighbour_bac.x1force.scaleAdd(-strength11,axis11,neighbour_bac.x1force);

                double strength22 = -k_sticking*(d22-stickingRestLong) * damping;
                this.x2force.scaleAdd(strength22,axis22,this.x2force);
                neighbour_bac.x2force.scaleAdd(-strength22,axis22,neighbour_bac.x2force);

                double strength12 = -k_sticking*(d12-stickingRestShort) * damping;
                this.x1force.scaleAdd(strength12,axis12,this.x1force);
                neighbour_bac.x2force.scaleAdd(-strength12,axis12,neighbour_bac.x2force);

                double strength21 = -k_sticking*(d21-stickingRestShort) * damping;
                this.x2force.scaleAdd(strength21,axis21,this.x2force);
                neighbour_bac.x1force.scaleAdd(-strength21,axis21,neighbour_bac.x1force);
            }
        }
    }

    // function which applies all relevant forces to the bacterium
    public void applyForces(BSimCapsuleBacterium neighbour_bac){
        computeSelfForce();
        computeWallForce();
        computeNeighbourForce(neighbour_bac);
//        brownianForce();
    }

    // object which computes how far the bacterium endpoints should be pushed due to the forces
    class ForceDisplacement implements BSimOdeSystem {
        private int numEq = 3;                // System of 7 equations

        private double[] force = new double[numEq];

        public void setForce(Vector3d appliedForce){
            force[0] = appliedForce.x;
            force[1] = appliedForce.y;
            force[2] = appliedForce.z;
        }

        // The equations
        public double[] derivativeSystem(double x, double[] y) {
            double[] dy = new double[numEq];

            dy[0] = force[0];
            dy[1] = force[1];
            dy[2] = force[2];

            return dy;
        }

        public int getNumEq() {
            return numEq;
        }

        // Initial conditions for the ODE system
        public double[] getICs() {
            double[] ics = {0,0,0};
            return ics;
        }
    }

//    private ForceDisplacement dfx1 = new ForceDisplacement();
//    private ForceDisplacement dfx2 = new ForceDisplacement();

    // sets all forces on the cell to zero after the forces have been applied and the displacement has been computed
    // this resets the simulation so that in the next timestep it is fresh
    public void setAllForcesZero(){
        this.x1force.set(0,0,0);
        this.x2force.set(0,0,0);
    }

    // updates the position of the cell according to forces
    public void updatePosition(){
//        System.out.println("Updating Position");
//        // NEW: update using RK45
//        dfx1.setForce(x1force);
//        dfx2.setForce(x2force);
//
//        // Solve the ode system
//        double[] x1New = BSimOdeSolver.rungeKutta45(dfx1, sim.getTime(), new double[]{x1.x, x1.y, x1.z}, sim.getDt());
//        double[] x2New = BSimOdeSolver.rungeKutta45(dfx2, sim.getTime(), new double[]{x2.x, x2.y, x2.z}, sim.getDt());
//
//        // Apply the force to x1 and x2
//        x1 = new Vector3d(x1New[0], x1New[1], x1New[2]);
//        x2 = new Vector3d(x2New[0], x2New[1], x2New[2]);

        // OLD: Euler
        this.x1.scaleAdd(sim.getDt(), this.x1force, this.x1);
        this.x2.scaleAdd(sim.getDt(), this.x2force, this.x2);

        this.setAllForcesZero();
    }

    public void setBrownianForceMagnitude() {
        brownianForceMagnitude = Math.sqrt(2*stokesCoefficient()*BSim.BOLTZMANN*sim.getTemperature()/sim.getDt())*Math.pow(10,9);
    }

    public double stokesCoefficient() { return 6.0*Math.PI*radius*sim.getVisc(); } // micrometers*Pa sec

    /**
     * Applies a Brownian force to the particle. The applied force is a function of
     * radius, viscosity and temperature; if viscosity or temperature is changed externally,
     * you should call setBrownianForceMagnitude() again
     */
    public void brownianForce() {
        Vector3d f = new Vector3d(rng.nextGaussian(), rng.nextGaussian(), rng.nextGaussian());
        f.scale(brownianForceMagnitude);
        x1force.add(f);
        x2force.add(f);
    }

    // return angle between cell and other cell
    public double coordinate(BSimCapsuleBacterium neighbour_bac){
        Vector3d d1= new Vector3d();
        Vector3d d2= new Vector3d();
        d1.sub(this.x1,this.x2);
        d2.sub(neighbour_bac.x1,neighbour_bac.x2);

        double angle = Math.acos((d1.x * d2.x+ d1.y*d2.y+
                    d1.z*d2.z)/(d1.length()*d2.length()));
        return angle;
    }

    // returns orientation of cell
    public double direction(){
        Vector3d d = new Vector3d();
        d.sub(this.x1 ,this.x2);
        double direction = d.angle(new Vector3d(1,0,0));
        return direction;
    }

    // angle from y-axis between -90 and 90 used by cell profiler - Sohaib Nadeem
    public double cell_profiler_angle() {
        Vector3d d = new Vector3d();
        d.sub(this.x1 ,this.x2);
        /*
        if (d.y < 0) {
            d.scale(-1);
        }
        return d.angle(new Vector3d(1,0,0)) - Math.PI / 2;
         */
        double angle = d.angle(new Vector3d(0,1,0));
        if (d.x < 0) {
            angle *= -1;
        }
        if (angle > Math.PI / 2) {
            angle -= Math.PI;
        }
        if (angle < -Math.PI / 2) {
            angle += Math.PI;
        }
        return angle;
    }

}
