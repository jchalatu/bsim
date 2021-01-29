package PhysModBsim;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.export.BSimLogger;

import java.io.IOException;
import java.util.ArrayList;

public class CellProfilerLogger extends BSimLogger {

    ArrayList<BSimCapsuleBacterium> bacteriaAll;

    public CellProfilerLogger(BSim sim, String filename, ArrayList<BSimCapsuleBacterium> bacteriaAll) {
        super(sim, filename);
        this.bacteriaAll = bacteriaAll;
    }

    @Override
    public void before() {
        super.before();
        String CellProfilerFields = "ImageNumber,ObjectNumber,AreaShape_Area,AreaShape_BoundingBoxArea," +
                "AreaShape_BoundingBoxMaximum_X,AreaShape_BoundingBoxMaximum_Y,AreaShape_BoundingBoxMinimum_X," +
                "AreaShape_BoundingBoxMinimum_Y,AreaShape_Center_X,AreaShape_Center_Y,AreaShape_Compactness," +
                "AreaShape_Eccentricity,AreaShape_EquivalentDiameter,AreaShape_EulerNumber,AreaShape_Extent," +
                "AreaShape_FormFactor,AreaShape_MajorAxisLength,AreaShape_MaxFeretDiameter,AreaShape_MaximumRadius," +
                "AreaShape_MeanRadius,AreaShape_MedianRadius,AreaShape_MinFeretDiameter,AreaShape_MinorAxisLength," +
                "AreaShape_Orientation,AreaShape_Perimeter,AreaShape_Solidity,Children_EditedObjects8_Count," +
                "Children_RelateObjects9_Count,Location_Center_X,Location_Center_Y,Number_Object_Number," +
                "Parent_EditedObjects8,Parent_IdentifyPrimaryObjects7\n";
        /*
        String CellProfilerFields = "ImageNumber,ObjectNumber,Intensity_IntegratedIntensityEdge_DilateImage6," +
                "Intensity_IntegratedIntensity_DilateImage6,Intensity_LowerQuartileIntensity_DilateImage6," +
                "Intensity_MADIntensity_DilateImage6,Intensity_MassDisplacement_DilateImage6," +
                "Intensity_MaxIntensityEdge_DilateImage6,Intensity_MaxIntensity_DilateImage6," +
                "Intensity_MeanIntensityEdge_DilateImage6,Intensity_MeanIntensity_DilateImage6," +
                "Intensity_MedianIntensity_DilateImage6,Intensity_MinIntensityEdge_DilateImage6," +
                "Intensity_MinIntensity_DilateImage6,Intensity_StdIntensityEdge_DilateImage6," +
                "Intensity_StdIntensity_DilateImage6,Intensity_UpperQuartileIntensity_DilateImage6," +
                "Location_CenterMassIntensity_X_DilateImage6,Location_CenterMassIntensity_Y_DilateImage6," +
                "Location_CenterMassIntensity_Z_DilateImage6,Location_Center_X,Location_Center_Y," +
                "Location_MaxIntensity_X_DilateImage6,Location_MaxIntensity_Y_DilateImage6," +
                "Location_MaxIntensity_Z_DilateImage6,Number_Object_Number,Parent_IdentifyPrimaryObjects7," +
                "RadialDistribution_FracAtD_DilateImage6_1of4,RadialDistribution_FracAtD_DilateImage6_2of4," +
                "RadialDistribution_FracAtD_DilateImage6_3of4,RadialDistribution_FracAtD_DilateImage6_4of4," +
                "RadialDistribution_MeanFrac_DilateImage6_1of4,RadialDistribution_MeanFrac_DilateImage6_2of4," +
                "RadialDistribution_MeanFrac_DilateImage6_3of4,RadialDistribution_MeanFrac_DilateImage6_4of4," +
                "RadialDistribution_RadialCV_DilateImage6_1of4,RadialDistribution_RadialCV_DilateImage6_2of4," +
                "RadialDistribution_RadialCV_DilateImage6_3of4,RadialDistribution_RadialCV_DilateImage6_4of4," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_0_0," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_1_1," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_2_0," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_2_2," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_3_1," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_3_3," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_4_0," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_4_2," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_4_4," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_5_1," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_5_3," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_5_5," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_6_0," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_6_2," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_6_4," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_6_6," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_7_1," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_7_3," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_7_5," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_7_7," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_8_0," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_8_2," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_8_4," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_8_6," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_8_8," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_9_1," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_9_3," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_9_5," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_9_7," +
                "RadialDistribution_ZernikeMagnitude_DilateImage6_9_9," +
                "RadialDistribution_ZernikePhase_DilateImage6_0_0,RadialDistribution_ZernikePhase_DilateImage6_1_1," +
                "RadialDistribution_ZernikePhase_DilateImage6_2_0,RadialDistribution_ZernikePhase_DilateImage6_2_2," +
                "RadialDistribution_ZernikePhase_DilateImage6_3_1,RadialDistribution_ZernikePhase_DilateImage6_3_3," +
                "RadialDistribution_ZernikePhase_DilateImage6_4_0,RadialDistribution_ZernikePhase_DilateImage6_4_2," +
                "RadialDistribution_ZernikePhase_DilateImage6_4_4,RadialDistribution_ZernikePhase_DilateImage6_5_1," +
                "RadialDistribution_ZernikePhase_DilateImage6_5_3,RadialDistribution_ZernikePhase_DilateImage6_5_5," +
                "RadialDistribution_ZernikePhase_DilateImage6_6_0,RadialDistribution_ZernikePhase_DilateImage6_6_2," +
                "RadialDistribution_ZernikePhase_DilateImage6_6_4,RadialDistribution_ZernikePhase_DilateImage6_6_6," +
                "RadialDistribution_ZernikePhase_DilateImage6_7_1,RadialDistribution_ZernikePhase_DilateImage6_7_3," +
                "RadialDistribution_ZernikePhase_DilateImage6_7_5,RadialDistribution_ZernikePhase_DilateImage6_7_7," +
                "RadialDistribution_ZernikePhase_DilateImage6_8_0,RadialDistribution_ZernikePhase_DilateImage6_8_2," +
                "RadialDistribution_ZernikePhase_DilateImage6_8_4,RadialDistribution_ZernikePhase_DilateImage6_8_6," +
                "RadialDistribution_ZernikePhase_DilateImage6_8_8,RadialDistribution_ZernikePhase_DilateImage6_9_1," +
                "RadialDistribution_ZernikePhase_DilateImage6_9_3,RadialDistribution_ZernikePhase_DilateImage6_9_5," +
                "RadialDistribution_ZernikePhase_DilateImage6_9_7,RadialDistribution_ZernikePhase_DilateImage6_9_9," +
                "TrackObjects_Displacement_50,TrackObjects_DistanceTraveled_50,TrackObjects_FinalAge_50," +
                "TrackObjects_IntegratedDistance_50,TrackObjects_Label_50,TrackObjects_Lifetime_50," +
                "TrackObjects_Linearity_50,TrackObjects_ParentImageNumber_50,TrackObjects_ParentObjectNumber_50," +
                "TrackObjects_TrajectoryX_50,TrackObjects_TrajectoryY_50,Length,Orientation\n";
         */
        write(CellProfilerFields);
    }

    @Override
    public void during() {
        /*
        String temp = "";
        for(int i = 0; i < 18; i++) {
            temp += "-,";
        }
        final String Fields_C_to_T = new String(temp);

        temp = "";
        for(int i = 0; i < 3; i++) {
            temp += "-,";
        }
        final String Fields_W_to_Y = new String(temp);

        temp = "";
        for(int i = 0; i < 72; i++) {
            temp += "-,";
        }
        final String Fields_AB_to_CU = new String(temp);
         */
        String temp = "";
        for(int i = 0; i < 6; i++) {
            temp += "-,";
        }
        final String Fields_C_to_H = new String(temp);

        temp = "";
        for(int i = 0; i < 6; i++) {
            temp += "-,";
        }
        final String Fields_K_to_P = new String(temp);

        temp = "";
        for(int i = 0; i < 5; i++) {
            temp += "-,";
        }
        final String Fields_R_to_V = new String(temp);

        temp = "";
        for(int i = 0; i < 4; i++) {
            temp += "-,";
        }
        final String Fields_Y_to_AB = new String(temp);

        String buffer = "";
        for(BSimCapsuleBacterium b : bacteriaAll) {
            /*
            buffer += (sim.getTimestep() * sim.getDt() / dt + 1) + "," + (b.id + 1) + "," + Fields_C_to_T +
                    b.position.x + "," + b.position.y + "," + Fields_W_to_Y +
                    (b.id + 1) + "," + "-" + "," + Fields_AB_to_CU +
                    "-,-,-,-," + (b.id + 1) + "," + "-,-,-,-," +
                    "-,-," + b.L + "," + b.direction() + "\n";
            */
            buffer += ((int) (sim.getTimestep() * sim.getDt() / dt) + 1) + "," + (b.id + 1) + "," + Fields_C_to_H +
                    (b.position.x * NeighbourInteractions.pixel_to_um_ratio) + "," +
                    (b.position.y  * NeighbourInteractions.pixel_to_um_ratio) + "," + Fields_K_to_P +
                    (b.L + 2 * b.radius) * NeighbourInteractions.pixel_to_um_ratio + "," + Fields_R_to_V +
                    (2 * b.radius) * NeighbourInteractions.pixel_to_um_ratio + "," +
                    (b.direction() - Math.PI / 2) + "," + Fields_Y_to_AB +
                    (b.position.x * NeighbourInteractions.pixel_to_um_ratio) + "," +
                    (b.position.y * NeighbourInteractions.pixel_to_um_ratio) + "," +
                    (b.id + 1) + "," + (b.id + 1) + "," + "-\n";
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
