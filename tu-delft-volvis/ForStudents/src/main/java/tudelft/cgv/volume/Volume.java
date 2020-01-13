/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tudelft.cgv.volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel modified by Anna
 */

//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////

public class Volume {
    

    //Do NOT modify these attributes
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;

    // Do NOT modify this function
    // This function returns the nearest neighbour given a position in the volume given by coord.
    public short getVoxelNN(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.round(coord[0]); 
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
    
        return getVoxel(x, y, z);
    }
        

    //Do NOT modify this function
    //This function linearly interpolates the value g0 and g1 given the factor (t) 
    //the result is returned. It is used for the tri-linearly interpolation the values 
    private float interpolate(float g0, float g1, float factor) {
        float result = (1 - factor)*g0 + factor*g1;
        return result; 
    }
             
    //Do NOT modify this function
    // This function returns the trilinear interpolated value of the position given by  position coord.
    public float getVoxelLinearInterpolate(double[] coord) {
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return 0;
        }
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.floor(coord[0]); 
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        float fac_x = (float) coord[0] - x;
        float fac_y = (float) coord[1] - y;
        float fac_z = (float) coord[2] - z;

        float t0 = interpolate(getVoxel(x, y, z), getVoxel(x+1, y, z), fac_x);
        float t1 = interpolate(getVoxel(x, y+1, z), getVoxel(x+1, y+1, z), fac_x);
        float t2 = interpolate(getVoxel(x, y, z+1), getVoxel(x+1, y, z+1), fac_x);
        float t3 = interpolate(getVoxel(x, y+1, z+1), getVoxel(x+1, y+1, z+1), fac_x);
        float t4 = interpolate(t0, t1, fac_y);
        float t5 = interpolate(t2, t3, fac_y);
        float t6 = interpolate(t4, t5, fac_z);
        
        return t6; 
    }


    float a = -0.75f; // global variable that defines the value of a used in cubic interpolation.
    // you need to chose the right value
        
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
        
    // Function that computes the weights for one of the 4 samples involved in the 1D interpolation 
    private float weight (float x, Boolean one_two_sample)
     {
         //Calculates the absolute value of X, X^2 and X^3.
         float X = Math.abs(x);
         float X2 = (float) Math.pow(X,2);
         float X3 = (float) Math.pow(X,3);

          // Compute the cubic interpolation kernel
         if(0<=X && X<1){
             return (float)(a+2)*X3-(a+3)*X2+1;
         } else if(1<=X && X<2) {
             return (float)a*X3-5*a*X2+8*a*X-4*a;
         } else {
             return (float)0;
         }
    } 
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Function that computes the 1D cubic interpolation. g0,g1,g2,g3 contain the values of the voxels that we want to interpolate
    // factor contains the distance from the value g1 to the position we want to interpolate to.
    // We assume the out of bounce checks have been done earlier
    
    public float cubicinterpolate(float g0, float g1, float g2, float g3, float factor) {
        // Compute the weights for the samples of the cubic interpolation kernel
        float h0 = weight(factor+1,true);
        float h1 = weight(factor,true);
        float h2 = weight(factor-1,true);
        float h3 = weight(factor-2,true);

        float fx = g0*h0+g1*h1+g2*h2+g3*h3;
        //take into account the exceptions
        if (fx < 0) {
            return 0;
        }
        if (fx > 255) {
            return 255;
        }
        float result = fx;                   
        return result; 
    }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 2D cubic interpolation implemented here. We do it for plane XY. Coord contains the position.
    // We assume the out of bounce checks have been done earlier
    public float bicubicinterpolateXY(double[] coord, int z) {
        // Interpolation can be seperated per axis 2D
        float Gx0, Gx1, Gx2, Gx3;
        /* notice that in this framework we assume that the distance between neighbouring voxels is 1 in all directions*/
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        
        //These factors can be computed in the same way as was done with the linearinterpolation
        float fac_X = (float) (coord[0] - x);
        float fac_Y = (float) (coord[1] - y);

        //First we interpolate points with the same y-coord on the sides along the x-axis
        Gx0 = cubicinterpolate(
                (float) getVoxel(x - 1, y - 1, z),
                (float) getVoxel(x, y - 1, z),
                (float) getVoxel(x + 1, y - 1, z),
                (float) getVoxel(x + 2, y - 1, z),
                fac_X);
        Gx1 = cubicinterpolate(
                (float) getVoxel(x - 1, y, z),
                (float) getVoxel(x, y, z),
                (float) getVoxel(x + 1, y, z),
                (float) getVoxel(x + 2, y, z),
                fac_X);
        Gx2 = cubicinterpolate(
                (float) getVoxel(x - 1, y + 1, z),
                (float) getVoxel(x, y + 1, z),
                (float) getVoxel(x + 1, y + 1, z),
                (float) getVoxel(x + 2, y + 1, z),
                fac_X);
        Gx3 = cubicinterpolate(
                (float) getVoxel(x - 1, y + 2, z),
                (float) getVoxel(x, y + 2, z),
                (float) getVoxel(x + 1, y + 2, z),
                (float) getVoxel(x + 2, y + 2, z),
                fac_X);
        //Secondly, we interpolate the points along the y-axis with the found interpolations from the x-axis
        float result = cubicinterpolate(Gx0, Gx1, Gx2, Gx3, fac_Y);
        return result;
    }
            
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // 3D cubic interpolation implemented here given a position in the volume given by coord.
    
    public float getVoxelTriCubicInterpolate(double[] coord) {
        if (coord[0] < 1 || coord[0] > (dimX-3) || coord[1] < 1 || coord[1] > (dimY-3)
                || coord[2] < 1 || coord[2] > (dimZ-3)) {
            return 0;
        }
       /*Again interpolation can be seperated per axis. Therefore, below we only have to do a cubicinterpolation
       with the values that can be found by doing the bicubicinterpolation for the 4 different XY planes.
        */
        int z = (int) Math.floor(coord[2]);
        float fac_Z = (float) (coord[2] - z);
        
        float y0, y1, y2, y3;
        
        y0 = bicubicinterpolateXY(coord, z - 1);
        y1 = bicubicinterpolateXY(coord, z);
        y2 = bicubicinterpolateXY(coord, z + 1);
        y3 = bicubicinterpolateXY(coord, z + 2);

        float result = cubicinterpolate(y0, y1, y2, y3, fac_Z);
        //Take into account the out of bounds coordinates
        result = Math.max(0, result);
        result = Math.min(255, result);
        
        return result; 
    }
    
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

	
    //Do NOT modify this function
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    //Do NOT modify this function
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    //Do NOT modify this function
    public short getVoxel(int x, int y, int z) {
    	int i = x + dimX*(y + dimY * z);
        return data[i];
    }
    
    //Do NOT modify this function
    public void setVoxel(int x, int y, int z, short value) {
    	int i = x + dimX*(y + dimY * z);
        data[i] = value;
    }
    
	//Do NOT modify this function
    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    //Do NOT modify this function
    public short getVoxel(int i) {
        return data[i];
    }
    
	//Do NOT modify this function
    public int getDimX() {
        return dimX;
    }
    
    //Do NOT modify this function
    public int getDimY() {
        return dimY;
    }
    
    //Do NOT modify this function
    public int getDimZ() {
        return dimZ;
    }

    //Do NOT modify this function
    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }
    
    //Do NOT modify this function
    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    //Do NOT modify this function
    public int[] getHistogram() {
        return histogram;
    }
    
    //Do NOT modify this function
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
}
