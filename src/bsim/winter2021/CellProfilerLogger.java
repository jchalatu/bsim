package bsim.winter2021;

import bsim.BSim;
import bsim.export.BSimLogger;
import bsim.winter2021.Bacterium;

import java.io.IOException;
import java.util.ArrayList;

public class CellProfilerLogger extends BSimLogger {

    ArrayList<? extends Bacterium> bac;
    double pixel_to_um_ratio;

    public CellProfilerLogger(BSim sim, String filename, ArrayList<? extends Bacterium> bac, double pixel_to_um_ratio) {
        super(sim, filename);
        this.bac = bac;
        this.pixel_to_um_ratio = pixel_to_um_ratio;
    }

    @Override
    public void before() {
        super.before();
        String Age_fields = "CellAge_Generation,Node_x1_x,Node_x1_y,Node_x2_x,Node_x2_y";
        String AreaShape_fields = "ImageNumber,ObjectNumber,AreaShape_Area,AreaShape_BoundingBoxArea," +
                "AreaShape_BoundingBoxMaximum_X,AreaShape_BoundingBoxMaximum_Y,AreaShape_BoundingBoxMinimum_X," +
                "AreaShape_BoundingBoxMinimum_Y,AreaShape_Center_X,AreaShape_Center_Y,AreaShape_Compactness," +
                "AreaShape_Eccentricity,AreaShape_EquivalentDiameter,AreaShape_EulerNumber,AreaShape_Extent," +
                "AreaShape_FormFactor,AreaShape_MajorAxisLength,AreaShape_MaxFeretDiameter,AreaShape_MaximumRadius," +
                "AreaShape_MeanRadius,AreaShape_MedianRadius,AreaShape_MinFeretDiameter,AreaShape_MinorAxisLength," +
                "AreaShape_Orientation,AreaShape_Perimeter,AreaShape_Solidity";
        String TrackObjects_fields = "TrackObjects_Displacement_50,TrackObjects_DistanceTraveled_50,TrackObjects_FinalAge_50," +
                "TrackObjects_IntegratedDistance_50,TrackObjects_Label_50,TrackObjects_Lifetime_50," +
                "TrackObjects_Linearity_50,TrackObjects_ParentImageNumber_50,TrackObjects_ParentObjectNumber_50," +
                "TrackObjects_TrajectoryX_50,TrackObjects_TrajectoryY_50";
        String other_fields = "Location_Center_X,Location_Center_Y,Location_Center_Z,Number_Object_Number,Parent_filtered_objects_1_labeled";
        String other_fields2 = "Location_Center_X,Location_Center_Y,Number_Object_Number,Parent_EditedObjects8,Parent_IdentifyPrimaryObjects7,Division_Angle";
        //write(CellProfilerFields + "," + other2 + "\n");
        write(AreaShape_fields + "," + other_fields + "," + TrackObjects_fields + "," + Age_fields + "\n");
    }

    @Override
    public void during() {
        String str = "-,";

        String buffer = "";
        for(Bacterium b : bac) {
            double area = (Math.PI * b.radius * b.radius + 2 * b.radius * b.L) * 13.89 * 13.89;
            int image_number = (int) (sim.getTimestep() * sim.getDt() / dt + 0.5) + 1;

            String AreaShape_fields =  image_number + "," + (b.id + 1) + "," + area + "," + str.repeat(5) +
                    (b.position.x * pixel_to_um_ratio) + "," +
                    (b.position.y  * pixel_to_um_ratio) + "," + str.repeat(6) +
                    (b.L + 2 * b.radius) * pixel_to_um_ratio + "," + str.repeat(5) +
                    (2 * b.radius) * pixel_to_um_ratio + "," +
                    b.cell_profiler_angle() + ",-,-";

            /*
            String other_fields = (b.position.x * pixel_to_um_ratio) + "," +
                    (b.position.y * pixel_to_um_ratio) + "," +
                    (b.id + 1) + "," + (b.id + 1) + "," + "-," +
                    b.angle_initial;
             */
            String other_fields2 = (b.position.x * pixel_to_um_ratio) + "," +
                    (b.position.y * pixel_to_um_ratio) + "," + "-," +
                    (b.id + 1) + ",-";

            String TrackObjects_fields = str.repeat(4) + (b.origin_id + 1) + "," + str.repeat(2) +
                    (image_number - 1) + "," + (b.lifetime == 0 ? b.parent_id + 1 : b.id + 1) + ",-,-";

            String Age_fields = String.valueOf(b.generation) + "," + (b.x1.x * pixel_to_um_ratio) + "," + (b.x1.y * pixel_to_um_ratio) + "," + (b.x2.x * pixel_to_um_ratio) + "," + (b.x2.y * pixel_to_um_ratio);
            buffer += AreaShape_fields + "," + other_fields2 + "," + TrackObjects_fields + "," + Age_fields + "\n";
        }
        write(buffer);
    }

    @Override
    public void write(String text) {
        try {
            bufferedWriter.write(text);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
