package bsim.winter2021;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class RawReader {
    BufferedReader csvReader;
    final double pixel_to_um_ratio;

    public RawReader(double pixel_to_um_ratio) {
        this.pixel_to_um_ratio = pixel_to_um_ratio;
    }

    // TODO: improve
    // The csv file is closed after reading
    public ArrayList<double[]> readcsv(String filepath) {
        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader(filepath));
        } catch (FileNotFoundException e) {
            // if that doesn't work, print out an error
            e.printStackTrace();
        }

        double[][] initEndpoints = new double[4][];
        ArrayList<double[]> initEndpoints_arrlist = new ArrayList();

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

        for(int i = 0; i < initEndpoints[i].length; i++) {
            initEndpoints_arrlist.add(new double[] {initEndpoints[0][i] / pixel_to_um_ratio,
                    initEndpoints[1][i] / pixel_to_um_ratio,
                    initEndpoints[2][i] / pixel_to_um_ratio,
                    initEndpoints[3][i] / pixel_to_um_ratio});
        }
        return initEndpoints_arrlist;
    }
}
