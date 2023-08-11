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
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import io.scif.*;
import io.scif.bf.BioFormatsFormat;
import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import io.scif.services.FormatService;
import io.scif.config.SCIFIOConfig.ImgMode;
import io.scif.ome.OMEMetadata;
import loci.common.Location;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.axis.CalibratedAxis;
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
import net.imglib2.type.logic.BitType;
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
import java.io.IOException;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.List;

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

    RoiManager roiManager;
    double pixelSize;

    String filename;
    @Override
    public void run() {

        //IJ.run("Cellpose setup...", "cellposeenvdirectory=C:\\Users\\tmchugh.ED\\.conda\\envs\\cellpose " +
        //        "envtype=conda usegpu=true usemxnet=false usefastmode=false useresample=false version=2.0");

            File[] files = filePath.listFiles();
            roiManager = new RoiManager();
            boolean first = true;
            for (File file : files) {
                if (file.toString().contains(".ims") && !file.toString().contains(".xml")) {
                    //Open file and get filename and filepath

                    ImagePlus imp = IJ.openImage(file.getPath());
                    Img<T> img = ImageJFunctions.wrapReal(imp);
                    //Img<T> img = openDataset(file);
                    uiService.show(img);
                    filename = FilenameUtils.removeExtension(file.getName());

                    //Maximum Z-project
                    RandomAccessibleInterval<T> zProj = (RandomAccessibleInterval) zProject(img);

//                    if(first){
//                    //IJ.run( ImageJFunctions.wrap(img, "test"), "Cellpose setup...", "");
//                    first=false;}

                    //Get cell outlines
                    RandomAccessibleInterval<T> slice647 = getChannel((RandomAccessibleInterval<T>) zProj, 3);
                    ImagePlus outputImp = ImageJFunctions.wrap(slice647, "Slice_647");
                    ImagePlus overlay = IJ.createImage("Overlay", "RGB", outputImp.getWidth(), outputImp.getHeight(), outputImp.getNChannels(), outputImp.getNSlices(), outputImp.getNFrames());
                    overlay.show();

                    outputImp.show();
                    IJ.run("Cellpose Advanced", "diameter=60 cellproba_threshold=0.0 flow_threshold=0.4" +
                            " anisotropy=1.0 diam_threshold=12.0 model=cyto2 nuclei_channel=0 cyto_channel=1 " +
                            "dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
                    getROIsfromMask();
                    Roi[] cells = roiManager.getRoisAsArray();
                    roiManager.reset();
                    outputImp.show();
                    IJ.run("RGB Color");

                    //Create output image
                    RandomAccessibleInterval<T> slice405 = getChannel(zProj, 0);
                    slice405 = ops.filter().gauss(slice405, 5);
                    RandomAccessibleInterval<T> slice488 = getChannel(zProj, 1);
                    RandomAccessibleInterval<T> slice555 = getChannel(zProj, 2);

                    ArrayList<double[]> output = new ArrayList<>();
                    for (Roi cell : cells) {
                        //Find nuclei ch0
                        //thresholding (Default) particles 10-Inf. Circ. 0.3-1 or cellpose?
                        //uiService.show(slice405);
                        ArrayList<Roi> nucleii = findNucleus(ImageJFunctions.wrap(slice405, "kSlice"), cell);
                        //uiService.show(slice647);
                        roiManager.reset();
                        //For each nucleus find centromeres in ch3 threshold per nucleus
                        for (Roi nucleus:nucleii) {
                            double[] intensities = new double[8];
                            intensities[0] = getRoiIntensity(slice488,nucleus)*nucleus.getStatistics().area;
                            intensities[2] = getRoiIntensity(slice555,nucleus)*nucleus.getStatistics().area;
                            intensities[4] = nucleus.getStatistics().area;
                            intensities[6] = getRoiIntensity(slice647,nucleus)*nucleus.getStatistics().area;
                            ArrayList<Roi> centromeres = findCentromeres(ImageJFunctions.wrap(slice647,"slice"), nucleus);
                            roiManager.reset();
                            double area = 0;
                            for (Roi centromere: centromeres){
                                area += centromere.getStatistics().area;
                                    intensities[1] += getRoiIntensity(slice488,centromere)*centromere.getStatistics().area;
                                    intensities[3] += getRoiIntensity(slice555,centromere)*centromere.getStatistics().area;
                                    intensities[7] += getRoiIntensity(slice647,centromere)*centromere.getStatistics().area;
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
                    ImagePlus imp488 = ImageJFunctions.wrap(slice488,"488");
                    imp488.show();
                    ImagePlus imp555 = ImageJFunctions.wrap(slice555,"555");
                    imp555.show();
                    ImagePlus imp647 = ImageJFunctions.wrap(slice647,"647");
                    imp647.show();
                    overlay.show();
                    IJ.run(overlay,"32-bit","");

                    IJ.run("Merge Channels...", "c1=[555] c2=[488] c3=[647] c4=[Overlay] create");
                    IJ.save(WindowManager.getCurrentImage(), Paths.get(String.valueOf(filePath), filename + "_Overlay.tif").toString());
                    try {
                        MakeResults(output, T1, T2);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    IJ.run("Close All", "");
                }
            }
        }

    private double getRoiIntensity(RandomAccessibleInterval<T> image, Roi roi){
        MaskInterval maskInterval = roiService.toMaskInterval(roi);
        IterableInterval<T> mask = Views.interval(image, maskInterval);
        int count = 0;
        double  intensity = 0;
        net.imglib2.Cursor<T> cursor = mask.localizingCursor();

        //Loop through the pixels in the mask
        for(int k =0; k< mask.size(); k++){
            //RealType<T> value = cursor.get();
            int x = (int) cursor.positionAsDoubleArray()[0];
            int y = (int) cursor.positionAsDoubleArray()[1];

            //If the pixel is in the bounding ROI
            if(roi.contains(x,y)) {
                count++;
                intensity = intensity + cursor.get().getRealDouble();
            }
            //Move the cursors forwards
            cursor.fwd();
        }
        return intensity/count;
    }

    //Finds Kinetochore by applying a Yen threshold within an ROI
    private ArrayList<Roi> findNucleus(ImagePlus image, Roi roi){
        image.setRoi(roi);
        IJ.setAutoThreshold(image, "Default dark");
        IJ.run(image, "Analyze Particles...", "size=250-Infinity circularity=0.30-1.00 include add");
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
        image.setRoi(roi);
        IJ.setAutoThreshold(image, "Triangle dark");
        IJ.run( image, "Analyze Particles...", "size=size=5-500 pixel exclude include add");
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

    private RandomAccessibleInterval<T> getChannel(RandomAccessibleInterval<T> stack, int channel){
        return ops.transform().hyperSliceView(stack, 2, channel);
    }

    //Takes a nD image dims[x,y,c,z,t] and performs a mean z-projection
    private IterableInterval<T> zProject(RandomAccessibleInterval<T> stack){
        int dim = 3;
        long[] dims = stack.dimensionsAsLongArray();
        long[] newDims = new long[dims.length-1];
        int count = 0;
        for(int i = 0; i < dims.length; i++) {
            if(i!=dim){
                newDims[count]=dims[i];
                count++;
            }
        }
        FinalDimensions dimensions = new FinalDimensions(newDims);
        Img<FloatType> proj = ops.create().img(dimensions, new FloatType());
        UnaryComputerOp maxOp = Computers.unary(ops, Ops.Stats.Max.class,RealType.class, Iterable.class);
        return ops.transform().project(proj, stack, maxOp, dim);
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
            bufferedWriter.write("Threshold 1 = "+ T1 + " Threshold 2 = " + T2);
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
