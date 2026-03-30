/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package bio.coil.CoilEdinburgh;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Roi;
import ij.plugin.ChannelSplitter;
import ij.plugin.ZProjector;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import io.scif.*;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imagej.ops.special.computer.Computers;
import net.imagej.ops.special.computer.UnaryComputerOp;
import net.imagej.roi.ROIService;
import net.imglib2.FinalDimensions;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.MaskInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import org.apache.commons.io.FilenameUtils;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;
import java.io.IOException;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * This example illustrates how to create an ImageJ {@link Command} plugin.
 * <p>
 * </p>
 * <p>
 * You should replace the parameter fields with your own inputs and outputs,
 * and replace the {@link run} method implementation with your own logic.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins>Users Plugins>HweiLing_Centromere_Quantification")
public class HweiLing_Centromere_Quantification<T extends RealType<T>> implements Command {
    //
    // Feel free to add more parameters here...
    //
    @Parameter
    private FormatService formatService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private UIService uiService;

    @Parameter
    private OpService ops;

    @Parameter
    private ROIService roiService;

    @Parameter(label = "Open Folder: ", style="directory")
    public File filePath;

    @Parameter(label = "Threshold 1 (partial): ")
    public double T1 = 1.5;

    @Parameter(label = "Threshold 2 (full): ")
    public double T2 = 2;

    @Parameter(label = "Centromere channel (0-3): ")
    public int centromereCh = 3;

    @Parameter(label = "DNA channel (0-3): ")
    public int DNACh = 0;

    @Parameter(label = "Endogenous channel (0-3): ")
    public int endoCh = 1;

    @Parameter(label = "TransGene channel (0-3): ")
    public int transCh = 2;



    @Parameter(label = "File Extension: ", choices = {".ims", ".nd2", ".czi", ".tif"})
    public String extension;

    RoiManager roiManager;
    double pixelSize;

    String filename;
    double pixel_size;
    @Override
    public void run() {

        //IJ.run("Cellpose setup...", "cellposeenvdirectory=C:\\Users\\tmchugh.ED\\.conda\\envs\\cellpose " +
        //        "envtype=conda usegpu=true usemxnet=false usefastmode=false useresample=false version=2.0");

            File[] files = filePath.listFiles();
            roiManager = new RoiManager();
            boolean first = true;
            for (File file : files) {
                if (file.toString().contains(extension) && !file.toString().contains(".xml")) {
                    //Open file and get filename and filepath
                    ImagePlus imp = IJ.openImage(file.getPath());
                    ImagePlus impMaxProj = ZProjector.run(imp,"Max");
                    filename = FilenameUtils.removeExtension(file.getName());
                    pixel_size = imp.getCalibration().pixelWidth;

//                    if(first){
//                    //IJ.run( ImageJFunctions.wrap(img, "test"), "Cellpose setup...", "");
//                    first=false;}

                    //Get cell outlines

                    ImagePlus sliceCent = ChannelSplitter.split(impMaxProj)[centromereCh];
                    ImagePlus outputImp = ChannelSplitter.split(impMaxProj)[centromereCh];
                    ImagePlus overlay = IJ.createImage("Overlay", "16-bit", outputImp.getWidth(), outputImp.getHeight(), outputImp.getNChannels(), outputImp.getNSlices(), outputImp.getNFrames());
                    overlay.show();

                    //Use the centromere channel and cellpose to get cell outlines
                    outputImp.show();
                    IJ.run("Cellpose Advanced", "diameter=60 cellproba_threshold=0.0 flow_threshold=0.4" +
                            " anisotropy=1.0 diam_threshold=12.0 model=cyto2 nuclei_channel=0 cyto_channel=1 " +
                            "dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
                    getROIsfromMask();
                    Roi[] cells = roiManager.getRoisAsArray();
                    roiManager.reset();

                    //Split in to single channels
                    ImagePlus sliceDNA = ChannelSplitter.split(impMaxProj)[DNACh];
                    //removed gauss filter - not necessary if image not deconvolved
                    ImagePlus sliceEndo = ChannelSplitter.split(impMaxProj)[endoCh];
                    ImagePlus sliceTrans = ChannelSplitter.split(impMaxProj)[transCh];

                    ArrayList<double[]> output = new ArrayList<>();
                    for (Roi cell : cells) {
                        //For each cell find the nucleus in the DAPI channel
                        ArrayList<Roi> nucleii = findNucleus(sliceDNA, cell);
                        roiManager.reset();

                        for (Roi nucleus:nucleii) {
                            //for each nucleus get the total intensity in the endo, trans and cent channels and the nuclear area
                            double[] intensities = new double[8];
                            intensities[0] = getRoiIntensityIMP(sliceEndo,nucleus)*nucleus.getStatistics().area;
                            intensities[2] = getRoiIntensityIMP(sliceTrans,nucleus)*nucleus.getStatistics().area;
                            intensities[4] = nucleus.getStatistics().area;
                            intensities[6] = getRoiIntensityIMP(sliceCent,nucleus)*nucleus.getStatistics().area;

                            //For each nucleus find centromeres in the centromere channel, threshold per nucleus
                            ArrayList<Roi> centromeres = findCentromeres(sliceCent, nucleus);
                            roiManager.reset();
                            double area = 0;
                            for (Roi centromere: centromeres){
                                area += centromere.getStatistics().area;
                                    intensities[1] += getRoiIntensityIMP(sliceEndo,centromere)*centromere.getStatistics().area;
                                    intensities[3] += getRoiIntensityIMP(sliceTrans,centromere)*centromere.getStatistics().area;
                                    intensities[7] += getRoiIntensityIMP(sliceCent,centromere)*centromere.getStatistics().area;
                            }
                            intensities[0]= (intensities[0]-intensities[1])/(intensities[4]-area);
                            intensities[2]= (intensities[2]-intensities[3])/(intensities[4]-area);
                            intensities[6]= (intensities[6]-intensities[7])/(intensities[4]-area);
                            intensities[1] = intensities[1]/area;
                            intensities[3] = intensities[3]/area;
                            intensities[7] = intensities[7]/area;
                            intensities[5] = area;
                            //Draw centromeres and label cell
                            output.add(intensities);
                            roiManager.reset();
                            drawNumbers(output.size(),overlay,nucleus,centromeres);
                        }
                    }

                    //Create output image
                    sliceEndo.setTitle("Endo");
                    sliceEndo.show();
                    sliceTrans.setTitle("Trans");
                    sliceTrans.show();
                    sliceCent.setTitle("Centro");
                    sliceCent.show();
                    overlay.show();
                    //c1, c2 etc here correspond to the colours in the imageJ Image>Channels>Merge Channels... tool
                    //c1=red, c2=green, c3=blue, c4=grey, c5=cyan, c6=magenta, c7=yellow
                    IJ.run("Merge Channels...", "c1=[Centro] c2=[Endo] c3=[Trans] c4=[Overlay] create");
                    IJ.save(WindowManager.getCurrentImage(), Paths.get(String.valueOf(filePath), filename + "_Overlay.tif").toString());

                    //create or update results table
                    try {
                        MakeResults(output, T1, T2);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    IJ.run("Close All", "");
                }
            }
        }

    private double getRoiIntensityIMP(ImagePlus image, Roi roi){
        image.setRoi(roi);
        ImageProcessor ip = image.getProcessor();
        return ip.getStatistics().mean;
    }

    //Finds Kinetochore by applying a Yen threshold within an ROI
    private ArrayList<Roi> findNucleus(ImagePlus image, Roi roi){
        IJ.run(image, "Gaussian Blur...", "sigma=0.3 scaled");
        image.setRoi(roi);
        IJ.setAutoThreshold(image, "Default dark");

        IJ.run(image, "Analyze Particles...", "size=10-35 circularity=0.30-1.00 include add");
        Roi[] allROIs = roiManager.getRoisAsArray();
        ArrayList<Roi> finalROIs = new ArrayList<>();
        for (Roi roIs : allROIs) {
            int x = (int) roIs.getContourCentroid()[0];
            int y = (int) roIs.getContourCentroid()[1];
            if (roi.contains(x, y)) {
                finalROIs.add(roIs);
            }
        }
        roiManager.reset();
        return finalROIs;
    }

    private ArrayList<Roi> findCentromeres(ImagePlus image, Roi roi){
        roiManager.reset();
        //ImagePlus imp = image.duplicate();
        //IJ.run(imp, "Subtract Background...", "rolling=50");
        image.setRoi(roi);
        IJ.setAutoThreshold(image, "Yen dark");
        double min = 5*pixel_size*pixel_size;
        double max = 500*pixel_size*pixel_size;
        IJ.run( image, "Analyze Particles...", "size=size="+min+"-"+max+" exclude include add");
        Roi[] allROIs = roiManager.getRoisAsArray();
        ArrayList<Roi> finalROIs = new ArrayList<>();
        for (Roi roIs : allROIs) {
            int x = (int) roIs.getContourCentroid()[0];
            int y = (int) roIs.getContourCentroid()[1];
            if (roi.contains(x, y)) {
                finalROIs.add(roIs);
            }
        }
        roiManager.reset();
        return finalROIs;
    }

    private void drawNumbers(int Counter, ImagePlus ProjectedWindow, Roi roi, ArrayList<Roi> centromeres) {
        ProjectedWindow.show();
        ImageProcessor ip = ProjectedWindow.getProcessor();
        Font font = new Font("SansSerif", Font.PLAIN, 12);
        ip.setFont(font);

        String cellnumber = String.valueOf(Counter);
        ip.setColor(Color.yellow);
        ip.draw(roi);
        ip.setColor(Color.blue);
        for(Roi cent:centromeres){
            ip.draw(cent);
        }
        ip.setColor(Color.white);
        ip.drawString(cellnumber, (int) roi.getContourCentroid()[0], (int) roi.getContourCentroid()[1]);
        ProjectedWindow.updateAndDraw();
    }


    public void MakeResults(ArrayList<double[]> intensities, double T1, double T2) throws IOException {
        Date date = new Date(); // This object contains the current date value
        SimpleDateFormat formatter = new SimpleDateFormat("dd-MM-yyyy, hh:mm:ss");
        String CreateName = Paths.get(String.valueOf(filePath), "_Results.csv").toString();
        try {
            FileWriter fileWriter = new FileWriter(CreateName, true);
            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
            bufferedWriter.newLine();
            bufferedWriter.write(formatter.format(date));
            bufferedWriter.newLine();
            bufferedWriter.write("Pixel size:,"+pixel_size);
            bufferedWriter.newLine();
            bufferedWriter.write("DAPI CH:,"+DNACh+",Centromere CH:,"+centromereCh+",Endogenous CH:,"+endoCh+",Transfected CH:,"+transCh);
            bufferedWriter.newLine();
            bufferedWriter.write("Threshold 1 = ,"+ T1 + ",Threshold 2 =," + T2);
            bufferedWriter.newLine();
            bufferedWriter.write("File, Number,Nuclear_Area, Cent_Area, 488_total, 488_centromere, 555_total, 555_centromere,647_total, 647_centromere, T2>448>T1, T2>555>T1, 448>T2, 555>T2");//write header 1
            bufferedWriter.newLine();
            int partial_488 = 0;
            int partial_555 = 0;
            int full_488 = 0;
            int full_555 = 0;

            double average_488_area = 0;
            double average_555_area = 0;
            double average_488_cell = 0;
            double average_555_cell = 0;
            double average_647_cell = 0;
            double average_488_cent = 0;
            double average_555_cent = 0;
            double average_647_cent = 0;
            int count = 0;

            for (int i =0; i < intensities.size(); i++){//for each slice create and write the output string
                double[] values = intensities.get(i);

                if(values[1]>0) {
                    count++;
                    average_488_area += values[4];
                    average_555_area += values[5];
                    average_488_cell += values[0];
                    average_555_cell += values[2];
                    average_647_cell += values[6];
                    average_488_cent += values[1];
                    average_555_cent += values[3];
                    average_647_cent += values[7];

                    if(aboveThreshold(values[1], values[0],T2)){
                        full_488++;
                    }else if (aboveThreshold(values[1], values[0],T1)){
                        partial_488++;
                    }
                    if(aboveThreshold(values[3], values[2],T2)){
                        full_555++;
                    }else if (aboveThreshold(values[3], values[2],T1)){
                        partial_555++;
                    }

                    bufferedWriter.newLine();
                    bufferedWriter.write(filename + "," + (i + 1) + "," + values[4] + "," + values[5] + "," + values[0] + "," + values[1] + "," + values[2] + "," + values[3] + ","+values[6] + "," + values[7] +
                            "," + (aboveThreshold(values[1], values[0],T1)&!aboveThreshold(values[1] , values[0], T2))  + "," + (aboveThreshold(values[3] , values[2], T1)&!aboveThreshold(values[3] , values[2], T2))
                            + "," + aboveThreshold(values[1] , values[0], T2)+ "," + aboveThreshold(values[3] , values[2], T2));
                }
            }
            bufferedWriter.newLine();
            bufferedWriter.write(filename + ", Total/Average,"+average_488_area/count + "," +average_555_area/count + "," +average_488_cell/count + "," +average_488_cent/count +
                    "," +average_555_cell/count + "," +average_555_cent/count + ","+average_647_cell/count + "," +average_647_cent/count + "," +partial_488 + "," +partial_555 + "," +full_488 + "," +full_555);
            bufferedWriter.close();
        } catch (IOException ex) {
            System.out.println(
                    "Error writing to file '" + CreateName + "'");
        }
    }

    private boolean aboveThreshold(double cent, double cell, double threshold){
        if(cent/cell>=threshold){
            return true;
        }
        return false;
    }

    //Adds the Masks created by cellpose to the ROI manager
    public void getROIsfromMask() {

        //Gets the current image (the mask output from cellpose)
        ImagePlus mask = WindowManager.getCurrentImage();
        ImageStatistics stats = mask.getStatistics();
        //For each ROI (intensity per cell mask is +1 to intensity
        for (int i = 1; i < stats.max + 1; i++) {
            //Set the threshold for the cell and use analyse particles to add to ROI manager
            IJ.setThreshold(mask, i, i);
            IJ.run(mask, "Analyze Particles...", "add");
        }
    }




    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();
        ij.command().run(HweiLing_Centromere_Quantification.class, true);
    }

}
