package BSimAsymGrowth;

import bsim.BSim; 
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import bsim.winter2021.Bacterium;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import javax.vecmath.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import java.io.File;
 
/**
 * This class simulates the growth of cells asymmetrically. The nodes of the bacteria 
 * grow at different rates until a certain length threshold where they grow symmetrically again. 
 */
public class BSimAsymGrowth {
    static final double pixel_to_um_ratio = 13.89;
    final int width_pixels = 800;
    final int height_pixels = 500;
    final double width_um = width_pixels / pixel_to_um_ratio; // should be kept as constants for cell prof.
    final double height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used
	
    /** Whether to enable growth in the ticker etc. or not... */
    private static final boolean WITH_GROWTH = true;
    
    /** Wall boundaries of the bacteria (+x, +y, +z, -x, -y, -z). */
    private final boolean [] wall_boundaries = {false, false, true, false, false, true};
    
    // Boundaries
    // Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = false;
    
    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;
    
    // Simulation setup parameters. Set dimensions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{width_um, height_um, 1.}));

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    // optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 1;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;
    
    //growth rate standard deviation
    @Parameter(names="-el_stdv",arity=1,description = "elongation rate standard deviation")
    public static double el_stdv = 0.277;	
    @Parameter(names="-el_mean",arity=1,description = "elongation rate mean")
    public static double el_mean = 1.23;

    //elongation threshold standard deviation
    @Parameter(names="-div_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double div_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names="-div_mean",arity=1,description = "elongation threshold mean")
    public static double div_mean = 7.0;
    
    // Simulation Time
    @Parameter(names="-simt",arity=1,description = "simulation time")
    public static double sim_time = 6.5;
    @Parameter(names="-simdt",arity=1,description = "simulation time step")
    public static double sim_dt = 0.05;
    @Parameter(names="-export_time",arity=1,description = "export time")
    public static double export_time = 0.5;// Previously was 10, and simulation time was 100
    
    // internal force
    @Parameter(names="-k_int",arity=1,description = "internal force")
    public static double k_int = 50.0;
    // cell-cell collision force
    @Parameter(names="-k_cell",arity=1,description = "cell-cell collision force")
    public static double k_cell = 50.0;
    // sticking force
    @Parameter(names="-k_stick",arity=1,description = "side-to-side attraction")
    public static double k_sticking = 10.0;
    
    // sticking range
    @Parameter(names="-rng_stick",arity=1,description = "max range side-to-side attraction")
    public static double range_sticking = 0.6;
    
    // twist
    @Parameter(names="-twist",arity=1,description = "twist")
    public static double twist = 0.1;
    // push
    @Parameter(names="-push",arity=1,description = "push")
    public static double push = 0.05;
    
    // asymmetric growth threshold
    @Parameter(names="-l_asym",arity=1,description = "asymmetric growth threshold")
    public static double L_asym = 3.75;
    // value of asymmetry
    @Parameter(names="-asym",arity=1,description = "asymmetry")
    public static double asymmetry = 0.1;
    // symmetric growth
    @Parameter(names="-sym",arity=1,description = "symmetric growth")
    public static double sym_growth = 0.05;

    /** Main Function. 
     * This is the very first function that runs in the simulation.
     * This runs the simulation. */ 
    public static void main(String[] args) {

        // Creates new simulation data object
        // Initializing storage unit for simulation
    	BSimAsymGrowth bsim_ex = new BSimAsymGrowth();

        // Starts up JCommander which allows you to read options from the command line more easily
        // Command line stuff
        new JCommander(bsim_ex, args);

        // Begins the simulation
        bsim_ex.run();
    }
    
    /** Creates a new Bacterium object. */
    public static Bacterium createBacterium( BSim sim ) {
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        
        // Random initial positions 
        Vector3d pos1 = new Vector3d(Math.random()*sim.getBound().x, 
				Math.random()*sim.getBound().y, 
				Math.random()*sim.getBound().z);
	
        double r = BSimCapsuleBacterium.L_th * Math.sqrt(Math.random());
        double theta = Math.random() * 2 * Math.PI;
        Vector3d pos2 = new Vector3d(pos1.x + r * Math.cos(theta), 
				pos1.y + r * Math.sin(theta), 
				pos1.z);
        
        // Check if the random coordinates are within bounds
        while(pos2.x >= sim.getBound().x || pos2.x <= 0 || pos2.y >= sim.getBound().y || pos2.y <= 0) {
        	pos1 = new Vector3d(Math.random()*sim.getBound().x, 
    				Math.random()*sim.getBound().y, 
    				Math.random()*sim.getBound().z);
        	pos2 = new Vector3d(pos1.x + r * Math.cos(theta), 
    				pos1.y + r * Math.sin(theta), 
    				pos1.z);
        }
        
        // Creates a new bacterium object whose endpoints correspond to the above data
        Bacterium bacterium = new Bacterium(sim, pos1, pos2);

        // Determine the vector between the endpoints
        // If the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initialized to length 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(pos2, pos1); 						// Sub is subtract
        double length = dispx1x2.length(); 				// Determined.
        
        if (length < bacterium.L_max) {
            bacterium.initialise(length, pos1, pos2); 	// Redundant to record length, but ok.
        }
        
        // Assigns a growth rate and a division length to bacterium according to a normal distribution
        double growthRate = el_stdv * bacRng.nextGaussian() + el_mean;
        bacterium.setK_growth(growthRate);
        
        double lengthThreshold = div_stdv * bacRng.nextGaussian() + div_mean;
        bacterium.setElongationThreshold(lengthThreshold);
        
        // Assigns the specified forces, range, and impulses
        bacterium.setIntForce(k_int);
        bacterium.setCellForce(k_cell);
        bacterium.setStickForce(k_sticking);
        bacterium.setStickingRange(range_sticking);
        bacterium.setTwist(twist);
        bacterium.setPush(push);
        
        bacterium.setLAsym(L_asym);
        bacterium.setAsym(asymmetry);
        bacterium.setSym(sym_growth);

        return bacterium;
    }
    
    // This function does a lot. let's break it down
    // 1. Initializes simulation (how long it runs, etc)
    // 2. Reads initial position data for bacteria and creates new bacterium objects accordingly
    // 3. Runs simulation loop until the simulation time is up
    // |----> Updates bacteria positions according to the forces acting on them
    // |----> Logs all data from the simulation into an excel sheet
    // |----> Saves images of simulation
    public void run() {
    	
		/*********************************************************
		 * Define simulation domain size and start time
		 */
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        // Saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
        final BSim sim = new BSim();
        sim.setDt(0.05);				    // Set simulation timestep in time units (0.01)
        									// Let the time units be in hours
        sim.setSimulationTime(100);       	// Specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");		    // Time Format for display on images
        sim.setBound(simX, simY, simZ);		// Simulation domain Boundaries

		// NOTE - solid = false sets a periodic boundary condition. 
		// This overrides leakiness!
		sim.setSolid(true, true, true);		// Solid (true) or wrapping (false) boundaries
        // Leaky -> bacteria can escape from sides of six faces of the box. 
		// Only usable if fixedbounds allows it
		
		// Set the wall boundaries
		sim.setWallBoundaries(wall_boundaries);
		
		// Set a closed or open boundary for phage diffusion
		if ( fixedBounds ) {
			sim.setLeaky(false, false, false, false, false, false);
			sim.setLeakyRate(0, 0, 0, 0, 0, 0);
		}
		else {
			sim.setLeaky(true, true, true, true, false, false);
			sim.setLeakyRate(5, 5, 5, 5, 0, 0);
		}
		
        /*********************************************************
         * Create the bacteria
         */
		
        // Separate lists of bacteria in case we want to manipulate the species individually
        // If multiple subpopulations, they'd be initialized separately, they'd be kept in different
        // need an array for each subpopulation, members can be repeated.
        final ArrayList<Bacterium> bac = new ArrayList();
        
        /** Track all of the bacteria in the simulation, for use of common methods etc.
        A general class, no sub-population specifics. */
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Gets the location of the file that is currently running
        // Specify output file path
        String systemPath = new File("").getAbsolutePath()+"\\SingleCellSims";
        
        /** Creating random bacterium objects to demonstrate phage field */

		// Create new phage sensing bacteria objects randomly in space
		while( bacteriaAll.size() < initialPopulation) {
            
            // Creates a new bacterium object whose endpoints correspond to the above data
			Bacterium bacterium = createBacterium(sim);
			
			// Adds the newly created bacterium to our lists for tracking purposes
			bac.add(bacterium); 			// For separate sub-populations
			bacteriaAll.add(bacterium);		// For all cells	
		}

        // Set up stuff for growth. Placeholders for the recently born and dead
        final ArrayList<Bacterium> bac_born = new ArrayList();
        final ArrayList<Bacterium> bac_dead = new ArrayList();

        // Internal machinery - dont worry about this line
        // Some kind of initialize of mover
        final Mover mover;
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
        
        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; // logs data every 100 timesteps
        BasicTicker ticker = new BasicTicker(sim, bac, bacteriaAll, LOG_INTERVAL, bacRng, el_stdv, el_mean,
                div_stdv, div_mean);
        ticker.setGrowth(WITH_GROWTH);			// enables bacteria growth
        sim.setTicker(ticker);
        
        /*********************************************************
         * Set up the drawer
         */
        BasicDrawer drawer = new BasicDrawer(sim, width_pixels, height_pixels, pixel_to_um_ratio, bac);
        sim.setDrawer(drawer);

//        /*********************************************************
//         * Set up the drawer
//         */
//        BSimDrawer drawer = new BSimP3DDrawer(sim, 800, 600) {	
//            /**
//             * Draw the default cuboid boundary of the simulation as a partially transparent box
//             * with a wireframe outline surrounding it.
//             */
//            @Override
//            public void boundaries() {
//                p3d.noFill();
//                p3d.stroke(128, 128, 255);
//                p3d.pushMatrix();
//                p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
//                p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
//                p3d.popMatrix();
//                p3d.noStroke();
//            }
//
//            @Override
//            public void draw(Graphics2D g) {
//                p3d.beginDraw();
//
//                if(!cameraIsInitialised){
//                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
//                    p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f,
//                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
//                            simX > simY ? (float)simX : (float)simY,
//                            (float)bound.x*0.5f, (float)bound.y*0.5f, 0,
//                            0,1,0);
//                    cameraIsInitialised = true;
//                }
//
//                p3d.textFont(font);
//                p3d.textMode(PConstants.SCREEN);
//
//                p3d.sphereDetail(10);
//                p3d.noStroke();
//                p3d.background(255, 255,255);
//
//                scene(p3d);
//                boundaries();
//                time();
//
//                p3d.endDraw();
//                g.drawImage(p3d.image, 0,0, null);
//            }
//
//            /**
//             * Draw the formatted simulation time to screen.
//             */
//            @Override
//            public void time() {
//                p3d.fill(0);
//                //p3d.text(sim.getFormattedTimeHours(), 50, 50);
//                p3d.text(sim.getFormattedTime(), 50, 50);
//            }
//
//            @Override
//            public void scene(PGraphics3D p3d) {
//                p3d.ambientLight(128, 128, 128);
//                p3d.directionalLight(128, 128, 128, 1, 1, -1);
//				
//				// Draw the infected bacteria in red and the non-infected bacteria in green
//				for(Bacterium b : bac) {
//					draw(b, Color.GREEN);
//				}		
//
//            }
//        };
//        sim.setDrawer(drawer);
        
        export = false;
        if(export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio;

            if(fixedBounds){
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            String filePath = BSimUtils.generateDirectoryPath(systemPath +"/" + simParameters + "/");
//            String filePath = BSimUtils.generateDirectoryPath("/home/am6465/tmp-results/" + simParameters + "/");

            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Catie Terrey Fall 2020."); //change name when new person :)
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: "+ initialPopulation);
                    write("Ratio " + populationRatio);

                    if(fixedBounds){
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }
                }

                @Override
                public void during() {

                }
            };
            metaLogger.setDt(10);//3600);			// Set export time step
            sim.addExporter(metaLogger);

            BSimLogger posLogger = new BSimLogger(sim, filePath + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z, growth_rate,directions");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "," + formatter.format(b.getK_growth())
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger.setDt(10);			// set export time step for csv file
            sim.addExporter(posLogger);


            BSimLogger sumLogger = new BSimLogger(sim, filePath + "summary.csv") {


                @Override
                public void before() {
                    super.before();
                    write("time,id, status, p1x, p1y, p1z, p2x, p2y, p2z, px, py, pz, growth_rate, directions");
                }

                @Override
                public void during() {
                    String buffer = new String();
                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += sim.getFormattedTime()+","+b.id
                                + "," + b.getInfected()
                                + "," + b.x1.x
                                + "," + b.x1.y
                                + "," + b.x1.z
                                + "," + b.x2.x
                                + "," + b.x2.y
                                + "," + b.x2.z
                                + "," + b.position.x
                                + "," + b.position.y
                                + "," + b.position.z
                                + "," + b.getK_growth()
                                + "," + b.direction()
                                + "\n";
                    }

                    write(buffer);
                }
            };
            sumLogger.setDt(10);			// Set export time step
            sim.addExporter(sumLogger);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath );
            imageExporter.setDt(10); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);

            sim.export();

            /**
             * Drawing a java plot once we're done?
             * See TwoCellsSplitGRNTest
             */

        } 
        
        // Run the simulation
        else {
            sim.preview();
        }
        
        long simulationEndTime = System.nanoTime();

        //System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");
        
    }
    
    
}


