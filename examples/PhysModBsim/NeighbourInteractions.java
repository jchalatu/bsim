package PhysModBsim;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.*;
import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.*;

// Fall 2020 Project
// Growth experiments form Fall 2020
public class NeighbourInteractions {

    // Sohaib Nadeem
    static final double pixel_to_um_ratio = 13.89;
    final double width_pixels = 700;
    final double height_pixels = 250;
    final double width_um = width_pixels / pixel_to_um_ratio; // should be kept as constants for cell prof.
    final double height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used

    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;

    //set dimnesions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[] {width_um, height_um, 1.}));

    // Boundaries
    //Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = true;

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    //optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 10;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;

    //growth rate standard deviation
    @Parameter(names="-gr_stdv",arity=1,description = "growth rate standard deviation")
    public static double growth_stdv=0.05;
    //growth rate mean
    @Parameter(names="-gr_mean",arity=1,description = "growth rate mean")
    public static double growth_mean=0.2;

    //elongation threshold standard deviation
    @Parameter(names="-len_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double length_stdv=0.1;
    //elongation threshold mean
    @Parameter(names="-len_mean",arity=1,description = "elongation threshold mean")
    public static double length_mean=7.0;

    /**
     * Whether to enable growth in the ticker etc. or not...
     */
    // if false nothing can grow
    private static final boolean WITH_GROWTH = true;

    static Random bacRng;


    // this is the very first function that runs in the simulation
    // this runs the simulation
    public static void main(String[] args) {

        // creates new simulation data object
        // intializing storage unit for simualation
        NeighbourInteractions bsim_ex = new NeighbourInteractions();

        // starts up JCommander which allows you to read options from the command line more easily
        // command line stuff
        new JCommander(bsim_ex, args);

        // begins simulation
        bsim_ex.run();
    }

    public static Bacterium createBacterium(BSim sim, Vector3d x1, Vector3d x2) {
        // creates a new bacterium object whose endpoints correspond to the above data
        Bacterium bacterium = new Bacterium(sim, x1, x2);

        // determine the vector between the endpoints
        // if the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initiatlied to lengthe 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(x2,x1); // sub is subtract
        double length = dispx1x2.length(); // determined.
        if (length < bacterium.L_max) {
            bacterium.initialise(length, x1, x2); // redundant to record length, but ok.
        }

        // assigns a growth rate and a division length to bacterium according to a normal distribution
        double growthRate = growth_stdv * bacRng.nextGaussian() + growth_mean;
        bacterium.setK_growth(growthRate);
        double lengthThreshold = length_stdv * bacRng.nextGaussian() + length_mean;
        bacterium.setElongationThreshold(lengthThreshold);

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
        // initializes simulation domain size
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        // saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

        // create the simulation object
        final BSim sim = new BSim();
        sim.setDt(0.05);				    // set Simulation Timestep in time units
        sim.setSimulationTime(100);       // specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");		    // Time Format for display on images
        sim.setBound(simX, simY, simZ);		// Simulation domain Boundaries


        // Boundaries periodicity: true means walls, false means periodic
        sim.setSolid(true, true, true);

        // Leaky -> bacteria can escape from sides of six faces of the box. only useable if fixedbounds allows it
        // leaky means completely open
        // leakyrate is for small molecules, not quite sure what this does
        if(!fixedBounds) {
            sim.setLeaky(true, true, true, true,false, false);
            sim.setLeakyRate(0.1/60.0, 0.1/60.0, 0.1/60.0, 0.1/60.0, 0, 0);
        }

        /*********************************************************
         * Create the bacteria
         */
        // Separate lists of bacteria in case we want to manipulate the species individually
        // if multiple subpopulations, they'd be initiallized separately, they'd be kept in different
        // need an array for each subpopulation, members can be repeated.
        final ArrayList<Bacterium> bac = new ArrayList();
        // Track all of the bacteria in the simulation, for use of common methods etc
        // A genneral class, no sub-population specifics
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

        bacRng = new Random(); //random number generator
        bacRng.setSeed(50); // initializes random number generator

        // empty list which will later contain the endpoint of rectangle positions 4=x1,y1,x2,y2
        double[][] initEndpoints = new double[4][];

        // gets the location of the file that is currently running
        //specify output file path
        String systemPath = new File("").getAbsolutePath()+"\\SingleCellSims";

        BufferedReader csvReader = null;
        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader("C:\\Users\\sohai\\IdeaProjects\\bsim\\examples\\PhysModBsim\\MyExpt_EditedObjects8.csv"));
        } catch (FileNotFoundException e) {
            // if that doesn't work, print out an error
            e.printStackTrace();
        }
        try {
            String row = csvReader.readLine();
            // check if fields match those of the required csv file from CellProfiler
            while( (row = csvReader.readLine()) != null ) {
                double cell_info[] = Arrays.stream(row.split(",")).mapToDouble(Double::parseDouble).toArray();

                if ( (int) cell_info[0] == 1 ) {
                    double cell_length = ( cell_info[16] - cell_info[22] ) / pixel_to_um_ratio;
                    double cell_orientation = cell_info[23] + (Math.PI / 2);
                    double cell_center_x = cell_info[8] / pixel_to_um_ratio;
                    double cell_center_y = cell_info[9] / pixel_to_um_ratio;

                    double axis_x = cell_length * Math.cos(cell_orientation); // use trig. identity to simplify?
                    double axis_y = cell_length * Math.sin(cell_orientation); // use trig. identity to simplify?

                    Vector3d x1 = new Vector3d(cell_center_x + (axis_x / 2), cell_center_y + (axis_y / 2), 0.5/*bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0*/);
                    Vector3d x2 = new Vector3d(cell_center_x - (axis_x / 2), cell_center_y - (axis_y / 2), 0.5/*bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0*/);

                    Bacterium bac0 = createBacterium(sim, x1, x2);
                    // adds the newly created bacterium to our lists for tracking purposes
                    bac.add(bac0); //for separate subpopulations
                    bacteriaAll.add(bac0);  // for all cells
                }
                else {
                    break;
                }
            }
        } catch(IOException e) {
            e.printStackTrace(); // if there is an error, this will just print out the message
        }

        /*
        // creates a new csvreader object which can extract data from .csv files
        // Catue added this, to read a CSV file, which specified endpoints of each cell, one per row
        BufferedReader csvReader = null;
        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader("C:\\Users\\sohai\\OneDrive\\Desktop\\Work\\bsim-master - Winter2021\\examples\\PhysModBsim\\onecell.csv"));
        } catch (FileNotFoundException e) {
            // if that doesn't work, print out an error
            e.printStackTrace();
        }

        //reading the data
        // if loading the file works, try reading the content
        try {
            // goes through each row of the excel sheet and pulls out the initial positions
            String row = csvReader.readLine();
            int i=0;
            while (row != null) {
                // row.split takes a single line of the excel sheet and chops it up into the columns
                // maptodouble then takes the values in those columns and converts them to Java double data format
                // toarray then finally converts the data into an array
                initEndpoints[i] = Arrays.stream(row.split(",")).mapToDouble(Double::parseDouble).toArray();

                row = csvReader.readLine();
                i++;
            }
            csvReader.close(); //finally, close the file once all data is extracted
        } catch(IOException e) {
            e.printStackTrace(); // if there is an error, this will just print out the message
        }

        // now that the data is extracted, we can create the bacterium objects
        for(int j = 0; j < initEndpoints[0].length; j++){
            // initializes the endpoints of each bacterium from the array of endpoints
            // z-dimension is a small number, randomly generated, not sure why.
            Vector3d x1 = new Vector3d(initEndpoints[0][j]/13.89,initEndpoints[1][j]/13.89,bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0);
            Vector3d x2 = new Vector3d(initEndpoints[2][j]/13.89,initEndpoints[3][j]/13.89,bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0);
            // note that the endpoint positions are scaled by 13.89, since the images are a bit more than 2000 pixels wide
            // while the simulation is rougly 200 micrometers. the conversion factor ends up being 13.89
            //pixel to um scaling

            Bacterium bac0 = createBacterium(sim, x1, x2);

            // adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); //for separate subpopulations
            bacteriaAll.add(bac0);  // for all cells
        }
         */

        // Set up stuff for growth. Placeholders for the recently born and dead
        final ArrayList<Bacterium> bac_born = new ArrayList();
        final ArrayList<Bacterium> bac_dead = new ArrayList();

        // internal machinery - dont worry about this line
        // some kind of initailaize of mover
        final Mover mover;
        mover = new RelaxationMoverGrid(bacteriaAll, sim);

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; // logs data every 100 timesteps

        // This one is a bit long too. Let's break it up
        // 1. Begins an "action" -> this represents one timestep
        // 2. Tells each bacterium to perform their action() function
        // 3. Updates each chemical field in the simulation
        // 4. Bacteria are then told to grow()
        // 5. bacteria which are longer than their threshold are told to divide()
        // 6. forces are applied and bacteria move around
        // 7. bacteria which are out of bounds are removed from the simulation
        BSimTicker ticker = new BSimTicker() {
            @Override
            public void tick() {
                // ********************************************** Action
                long startTimeAction = System.nanoTime(); //wall-clock time, for diagnosing timing

                for(BSimCapsuleBacterium b : bacteriaAll) {
                    b.action(); //bacteria do action at each time step
                    // calculates and stores the midpoint of the cell.
                }

                long endTimeAction = System.nanoTime();
                if((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }  //outputing how long each step took, once every log interval.

                // ********************************************** Chemical fields
                startTimeAction = System.nanoTime();

                //field activity would go here.

                endTimeAction = System.nanoTime();
                if((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }

                // ********************************************** Growth related activities if enabled.
                if(WITH_GROWTH) {

                    // ********************************************** Growth and division
                    startTimeAction = System.nanoTime(); //start action timer

                    for (Bacterium b : bac) { //loop over bac array
                        b.grow();

                        // Divide if grown past threshold
                        if (b.L >= b.L_th) {
                            bac_born.add(b.divide());  // add daughter to newborn class,  'mother' keeps her status
                        }
                    }
                    bac.addAll(bac_born); //adds all the newborn daughters
                    bacteriaAll.addAll(bac_born); //adds all the newborn daughters
                    // just added - Sohaib Nadeem
                    for(Bacterium b : bac_born) {
                        // assigns a growth rate and a division length to each bacterium according to a normal distribution
                        double growthRate=growth_stdv*bacRng.nextGaussian() + growth_mean;
                        b.setK_growth(growthRate);

                        double lengthThreshold=length_stdv*bacRng.nextGaussian()+length_mean;
                        b.setElongationThreshold(lengthThreshold);
                    }
                    bac_born.clear(); // cleared for next time-step

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }
                    //above: prints out information abt bacteria when u want it to
                    // ********************************************** Neighbour interactions
                    startTimeAction = System.nanoTime();

                    mover.move();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }

                    // ********************************************** Boundaries/removal
                    startTimeAction = System.nanoTime();
                    // Removal
                    for (Bacterium b : bac) {
//                         Kick out if past the top or bottom boundaries
//                        if ((b.x1.y < 0) && (b.x2.y < 0)) {
//                            act_dead.add(b);
//                        }
//                        if ((b.x1.y > sim.getBound().y) && (b.x2.y > sim.getBound().y)) {
//                            act_dead.add(b);
//                        }
                        // kick out if past any boundary
                        if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
                            bac_dead.add(b);
                        } //bacteria out of bounds = dead
                    }
                    bac.removeAll(bac_dead);
                    bacteriaAll.removeAll(bac_dead);
                    bac_dead.clear();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }
                }

                startTimeAction = System.nanoTime();

                endTimeAction = System.nanoTime();
                if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    System.out.println("Switch took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                }

            }
        };

        sim.setTicker(ticker);


//      the rest of the code is the drawer (makes simulation images) and data logger (makes csv files)

        /*********************************************************
         * Set up the drawer
         */
        BSimDrawer drawer = new BSimP3DDrawer(sim, 600, 350) {
            /**
             * Draw the default cuboid boundary of the simulation as a partially transparent box
             * with a wireframe outline surrounding it.
             */
            @Override
            public void boundaries() {
                p3d.noFill();
                p3d.stroke(128, 128, 255);
                p3d.pushMatrix();
                p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
                p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
                p3d.popMatrix();
                p3d.noStroke();
            }

            @Override
            public void draw(Graphics2D g) {
                p3d.beginDraw();

                if(!cameraIsInitialised){
                    // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
                    p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f,
                            // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                            simX > simY ? (float)simX : (float)simY,
//                            10,
                            (float)bound.x*0.5f, (float)bound.y*0.5f, 0,
                            0,1,0);
                    cameraIsInitialised = true;
                }

                p3d.textFont(font);
                p3d.textMode(PConstants.SCREEN);

                p3d.sphereDetail(10);
                p3d.noStroke();
                p3d.background(255, 255,255);

                scene(p3d);
                boundaries();
                time();

                p3d.endDraw();
                g.drawImage(p3d.image, 0,0, null);
            }

            /**
             * Draw the formatted simulation time to screen.
             */
            @Override
            public void time() {
                p3d.fill(0);
//                p3d.text(sim.getFormattedTimeHours(), 50, 50);
                p3d.text(sim.getFormattedTime(), 50, 50);
            }

            @Override
            public void scene(PGraphics3D p3d) {
                p3d.ambientLight(128, 128, 128);
                p3d.directionalLight(128, 128, 128, 1, 1, -1);

//                draw(bac_act, new Color(55, 126, 184));
//                draw(bac_rep, new Color(228, 26, 28));

                for(Bacterium b : bac) {
                    draw(b, Color.blue);
                }

            }
        };
        sim.setDrawer(drawer);

        export=true;
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


            // Export a csv file that matches CellProfiler's output
            CellProfilerLogger cp_logger = new CellProfilerLogger(sim, filePath + "MyExpt_EditedObjects8_simulation.csv", bacteriaAll);
            cp_logger.setDt(10);			// Set export time step
            sim.addExporter(cp_logger);

            sim.export();

            /**
             * Drawing a java plot once we're done?
             * See TwoCellsSplitGRNTest
             */

        } else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();

        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");
    }
}
