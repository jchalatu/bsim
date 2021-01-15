package PhysModBsim;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;

import javax.vecmath.Vector3d;
import java.util.ArrayList;
import java.util.List;

/**
 */
public class Bacterium extends BSimCapsuleBacterium {

    // these fields are used to build a tree structure to keep track
    // of lineage - Sohaib Nadeem
    private List<Bacterium> children;
    private Bacterium parent;

    //function you call when you want to make a new bacterium object
    public Bacterium(BSim sim, Vector3d px1, Vector3d px2){
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);

        // initialize fields for lineage tree - Sohaib Nadeem
        children = new ArrayList<>();
        parent = null;
    }

    // in case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { //runs at every time step
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

    // this function is called when the bacterium has passed its' division threshold and is ready to divide.
    public Bacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble()/100,rng.nextDouble()/100,rng.nextDouble()/100);
        System.out.println("Bacterium " + this.id + " is dividing...");

        Vector3d u = new Vector3d(); u.sub(this.x2, this.x1);

        //decides where to split bacterium
        double divPert = 0.1*L*(rng.nextDouble() - 0.5); //0.1 is arbitrary scaling here.

        // get length of bacterium
        double L_actual = u.length();

        // computes lengths of child bacteria
        double L1 = L_actual*0.5*(1 + divPert) - radius;
        double L2 = L_actual*0.5*(1 - divPert) - radius;

        //finds new endpoints of mother and daughter
        Vector3d x2_new = new Vector3d();
        x2_new.scaleAdd(L1/L_actual, u, this.x1);
        Vector3d longVec = new Vector3d();
        longVec.scaleAdd(-1,this.x2,this.x1); //push along bacterium length
        longVec.scale(0.05*rng.nextDouble()); //push is applied to bacterium
        // impulse, not a force.
        longVec.add(randomVec);
        x2_new.add(longVec);

        Vector3d x1_child = new Vector3d();
        x1_child.scaleAdd(-(L2/L_actual), u, this.x2);
        Vector3d longVec2=new Vector3d();
        longVec2.scaleAdd(-1,longVec);
        //x1_child.add(longVec2); //push applied to child cell


        // Set the child cell.
        //creates new bacterium called child and adds it to the lists, gives posns, infected status and chemical field status
        Bacterium child = new Bacterium(sim, x1_child, new Vector3d(this.x2));
        this.initialise(L1, this.x1, x2_new);
        child.L = L2;

        // add child to list of children - Sohaib Nadeem
        addChild(child);

        // prints a line whenever a new bacterium is made
        System.out.println("Child ID is " + child.id);
        return child;
    }

    // seemingly unused code
//    int count=0;
//    public void set_count(int n){
//        this.count = n;
//    }
//    public int get_count(){
//        return count;
//    }

}




/*
    // computes force acting on cell from walls above
    public double computeVolumeExclusionForce(double k, double pen_depth) {
        if (pen_depth > 0) {
            return 0.4 * k * Math.pow(pen_depth, 2.5);
        }
        else {
            return 0;
        }
    }

    // computes force acting on cell from walls above
    public double computeWallAttractionForce(double dist) {
        if (dist < radius + range_wallstick) {
            return -1 * k_wallstick * (dist - radius);
        }
        else {
            return 0;
        }
    }

 */

// this accounts for attractive forces between bacteria and walls due to hydrostatic interactions
        /*
        x1force.scaleAdd(computeWallAttractionForce(x1.x), new Vector3d(1, 0, 0));
        x1force.scaleAdd(computeWallAttractionForce(sim.getBound().x - x1.x), new Vector3d(-1, 0, 0));
        x1force.scaleAdd(computeWallAttractionForce(x1.y), new Vector3d(0, 1, 0));
        x1force.scaleAdd(computeWallAttractionForce(sim.getBound().y - x1.y), new Vector3d(0, -1, 0));
        x1force.scaleAdd(computeWallAttractionForce(x1.z), new Vector3d(0, 0, 1));
        x1force.scaleAdd(computeWallAttractionForce(sim.getBound().z - x1.z), new Vector3d(0, 0, -1));

        x2force.scaleAdd(computeWallAttractionForce(x2.x), new Vector3d(1, 0, 0));
        x2force.scaleAdd(computeWallAttractionForce(sim.getBound().x - x2.x), new Vector3d(-1, 0, 0));
        x2force.scaleAdd(computeWallAttractionForce(x2.y), new Vector3d(0, 1, 0));
        x2force.scaleAdd(computeWallAttractionForce(sim.getBound().y - x2.y), new Vector3d(0, -1, 0));
        x2force.scaleAdd(computeWallAttractionForce(x2.z), new Vector3d(0, 0, 1));
        x2force.scaleAdd(computeWallAttractionForce(sim.getBound().z - x2.z), new Vector3d(0, 0, -1));
        */

/*
        // volume exclusion forces for all walls on x1
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x1.z), new Vector3d(0,0,1)); // bottom
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x1.z - sim.getBound().z), new Vector3d(0,0,-1)); // top
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x1.x), new Vector3d(1,0,0)); // left
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x1.x - sim.getBound().x), new Vector3d(-1,0,0)); //right
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x1.y), new Vector3d(0,1,0)); // behind
        x1force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x1.y - sim.getBound().y), new Vector3d(0,-1,0)); // in front
        // volume exclusion forces for all walls on x2
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x2.z), new Vector3d(0,0,1)); // bottom
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x2.z - sim.getBound().z), new Vector3d(0,0,-1)); // top
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x2.x), new Vector3d(1,0,0)); // left
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x2.x - sim.getBound().x), new Vector3d(-1,0,0)); //right
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius - x2.y), new Vector3d(0,1,0)); // behind
        x2force.scaleAdd(computeVolumeExclusionForce(k_wall, radius + x2.y - sim.getBound().y), new Vector3d(0,-1,0)); // in front
        */

/*
double minimumAngle = 7 /8 * Math.PI; //cells must be oriented to within pi/8 rads (22.5 degrees) of each other to stick end to end
        //filial links (secondary structure -> causes chains to form)
        //might want to include condition on being sisters
        // how to tell if two cells are sisters??
        //4 conditions (pairs of endpoints) to see if any are sufficnetly close and right angle
        // comparing x1 to x1. checking if this is the smallest of all 4, and is lesx than filial range
        // and angle is less than 2.5 degrees

        double short_distance, long_distance;
        Vector3d short_filial_force, long_filial_force,
                 short_filial_force_neighbor, long_filial_force_neighbor;

        if (d11 < d22 && d11 < d12 && d11 < d21 && (d11 < 2 * radius + range_filial) && angleuv >= minimumAngle) {
            filialAxis=dis11;
            short_distance = d11;
            short_filial_force = x1force;
            short_filial_force_neighbor = neighbour_bac.x1force;

            longFilialAxis=dis22;
            long_distance = d22;
            long_filial_force = x2force;
            long_filial_force_neighbor = neighbour_bac.x2force;
        }
        else if(d22 < d11 && d22 < d12 && d22 < d21) {
            filialAxis = dis22;
            short_distance = d22;
            short_filial_force = x2force;
            short_filial_force_neighbor = neighbour_bac.x2force;

            longFilialAxis=dis11;
            long_distance = d11;
            long_filial_force = x1force;
            long_filial_force_neighbor = neighbour_bac.x1force;
        }
        else if(d12 < d11 && d12 < d21 && d12 < d22) {
            filialAxis=dis12;
            short_distance = d12;
            short_filial_force = x1force;
            short_filial_force_neighbor = neighbour_bac.x2force;

            longFilialAxis=dis21;
            long_distance = d21;
            long_filial_force = x2force;
            long_filial_force_neighbor = neighbour_bac.x1force;

            angleuv = Math.PI - angleuv;
        }
        else {
            filialAxis=dis21;
            short_distance = d21;
            short_filial_force = x2force;
            short_filial_force_neighbor = neighbour_bac.x1force;

            longFilialAxis=dis12;
            long_distance = d12;
            long_filial_force = x1force;
            long_filial_force_neighbor = neighbour_bac.x2force;

            angleuv = Math.PI - angleuv;
        }

        if ((short_distance < 2 * radius + range_filial) && angleuv >= minimumAngle) {
            //if the closest two points are x11 and x12
            double stickyForceStrength= -k_filial * (short_distance - 2 * radius-0.05*radius); //spring attaching endpoints
            filialAxis.normalize(); // normalize axis vector to get correct scaling
            short_filial_force.scaleAdd(stickyForceStrength,filialAxis); // adds the force to the endpoint
            short_filial_force_neighbor.scaleAdd(-stickyForceStrength,filialAxis);
            //also a spring attaching furthest endpoints (but repulsive!) -> maybe this should be weaker?
            double longFilialForceStrength=k_longfilial*(long_distance-2*radius-thisLength-neighborLength-0.05*radius);
            longFilialAxis.normalize();
            long_filial_force.scaleAdd(longFilialForceStrength,longFilialAxis);
            long_filial_force_neighbor.scaleAdd(-longFilialForceStrength,longFilialAxis);
        }
 */