import java.io.*;
//import java.util.*;

// for image output
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import java.awt.Color;

public class Wave2D{
	static double c, c2, time;
	static int snap1, snap2, snapX, snapY;
	double[][][] u, uV, uA;
	double[][] temp;
	double[] sample;
	int x0, y0, steps, lx, ly;
	double u0, h;

	public Wave2D(int lx, int ly, int stepCount, double stepDuration) {
		// lx and ly are the number of grid points in those directions. These are the sizes of the 1st and 2nd array dimensions
		this.lx = lx;
		this.ly = ly;

		// steps is the number of time steps to perform. This is the size of the third array dimension
		// h is the duration of each step (seconds)
		this.steps = stepCount;
		this.h = stepDuration;

		// c is the wave speed c2 is c squared
		c2 = 15;	c = Math.sqrt(c2);
		
		time = 0;
		
		// snap are locations to print tracing/debugging
		snap1 = 0;	snap2 = snap1+1/*(int)(snap1*Math.sqrt(2))*/;
		snapX = 4;
		snapY = 4;

		// (x0,y0) is the grid point of the initial peak disturbance (used to generate the initial conditions)
		x0 = (int)(0.5*lx)/*(int)(Math.random()*(lx-1)+1)*/;
		y0 = (int)(0.5*ly)/*(int)(Math.random()*(ly-1)+1)*/;
		u0 = 1;

		// allocate the arrays given the lx,ly, and steps values
		u = new double[lx][ly][steps];
		uV = new double[lx][ly][steps];
		uA = new double[lx][ly][steps];
		temp = new double[lx][ly];
		sample = new double[steps];

		// Init all the allocated data structures to all 0 (java's runtime probably does this so we could skip if it takes significant time)
		for (int x = 0; x < lx; x++) {
			for (int y = 0; y < ly; y++) {
				temp[x][y] = 0;
				for (int t = 0; t < steps; t++) {
					u[x][y][t] = 0;
					uV[x][y][t] = 0;
					uA[x][y][t] = 0;
				}
			}
		}
		for (int t = 0; t < steps; t++) {
			sample[t] = 0;
		}
	}

	public void initGaussian(int x0, int y0, double u0) {
		this.x0 = x0;
		this.y0 = y0;
		this.u0 = u0;
		for(int x=0; x<lx; x++)
			for(int y=0; y<ly; y++)
				u[x][y][0] = 0;
		u[x0][y0][0] = u0;
	}

	public void initPyramid(int x0, int y0, double u0) {
		this.x0 = x0;
		this.y0 = y0;
		this.u0 = u0;
		for (int x = 1; x < lx - 1; x++) {
			for (int y = 1; y < ly - 1; y++) {
				// lower quad
				if (y < (y0 + 0.0) / x0 * (x) && y < (y0 + 0.0) / (x0 - lx) * (x - lx)) {
					u[x][y][0] = 1.0 * u0 * (y + 0.0) / y0;
				}
				// left quad
				else if (y >= (y0 + 0.0) / x0 * (x) && y < ly - 1.0 * (ly - y0) / (x0 - 0) * x) {
					u[x][y][0] = 1.0 * u0 * (x + 0.0) / x0;
				}
				// upper quad
				else if (y >= (ly - y0 + 0.0) / (lx - x0) * (x - lx) + ly) {
					u[x][y][0] = 1.0 * u0 * (ly - y - 1.0) / (ly - y0 - 1);
				}
				// right quad
				else {
					u[x][y][0] = 1.0 * u0 * (lx - x - 1.0) / (lx - x0 - 1);
				}
			}
		}
	}

	public void writeSurface(int timeStep) throws IOException {
		int height = 300;
		int width = (int)(height * 1.0*lx/ly);
		double scaleX = 1.0*lx/width;
		double scaleY = 1.0 * ly / height;
		double scaleU = 1.0 * 255 / u0;

		// Create a new image
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

		// Draw the image
		for (int px = 0; px < width; px++) {
			for (int py = 0; py < height; py++) {
				// Map pixel coords to surface coords and scale the u value to a grayscale color
				int x = (int)(px*scaleX);
				int y = ly -1 - (int) (py * scaleY);
				int grayValue = (int)(u[x][y][timeStep] * scaleU);

				// clip grayValue to 255 in case there is an error it wont crash the run 
				grayValue = (grayValue>255) ? 255 : grayValue;

				// this makes 1 pixel wide X and Y axis so we know where the origin is.
				if (px==0 || py==height-1)
					grayValue = 255;

				Color color = new Color(grayValue, grayValue, grayValue);

				// Set the pixel color
				image.setRGB(px, py, color.getRGB());
			}
		}

		// Write the image to a file
		String baseFilename = new String("out/2DSurf-") + timeStep;
		File outputfile = new File(baseFilename + ".png");
		ImageIO.write(image, "png", outputfile);

		// write it as a txt file with a pont (x,y,u) on each line
		PrintStream txtOut = new PrintStream(new File(baseFilename + ".txt"));
		for (int y = 0; y < ly; y++) {
			for (int x = 0; x < lx; x++) {
				txtOut.printf("%d,%d,%.4f\n", x, y, u[x][y][timeStep]);
			}
		}
	}

	double u2(int x, int y){								//2D
		return c2*(temp[x-1][y]+temp[x+1][y]+temp[x][y-1]+temp[x][y+1])/4-temp[x][y];
	}

	double u2(int x, int y, int t){							//2D
		return c2*(u[x-1][y][t]+u[x+1][y][t]+u[x][y-1][t]+u[x][y+1][t])/4-u[x][y][t];
	}

	double u2(double x, int y, double t){					//2D
		int xL = (int)Math.floor(x);
		int tL = (int)Math.floor(t);
		if(xL <= 0 || xL >= lx-1){
			return 0;
		}
		/*double uLow = ((x-xL)*u[xL+1][y][tL]+(1-x+xL)*u[xL][y][tL]);
		double uHigh = ((x-xL)*temp[xL+1][y]+(1-x+xL)*temp[xL][y]);
		double uC = (t-tL)*uHigh+(1-t+tL)*uLow;
		double uLeftLow = ((x-xL)*u[xL][y][tL]+(1-x+xL)*u[xL-1][y][tL]);
		double uLeftHigh = ((x-xL)*temp[xL][y]+(1-x+xL)*temp[xL-1][y]);
		double uL = (t-tL)*uLeftHigh+(1-t+tL)*uLeftLow;
		double uRightLow = ((x-xL)*u[xL+2][y][tL]+(1-x+xL)*u[xL+1][y][tL]);
		double uRightHigh = ((x-xL)*temp[xL+2][y]+(1-x+xL)*temp[xL+1][y]);
		double uR = (t-tL)*uRightHigh+(1-t+tL)*uRightLow;
		
		double uLowUp = ((x-xL)*u[xL+1][y+1][tL]+(1-x+xL)*u[xL][y+1][tL]);
		double uHighUp = ((x-xL)*temp[xL+1][y+1]+(1-x+xL)*temp[xL][y+1]);
		double uCUp = (uHighUp+uLowUp)/4;
		double uLowDn = ((x-xL)*u[xL+1][y-1][tL]+(1-x+xL)*u[xL][y-1][tL]);
		double uHighDn = ((x-xL)*temp[xL+1][y-1]+(1-x+xL)*temp[xL][y-1]);
		double uCDn = (uHighDn+uLowDn)/4;
		
		return c2*((uL+uR+uCUp+uCDn)/3-uC);*/
		
		double uLLL = u[xL][y-1][tL];		double z1 = 1/Math.sqrt(Math.pow(x-xL,2)+1+Math.pow(t-tL,2));			//(unitLength)/r
		double uLLH = temp[xL][y-1];		double z2 = 1/Math.sqrt(Math.pow(x-xL,2)+1+Math.pow(tL+1-t,2));
		double uLHL = u[xL][y+1][tL];		double z3 = 1/Math.sqrt(Math.pow(x-xL,2)+1+Math.pow(t-tL,2));
		double uLHH = temp[xL][y+1];		double z4 = 1/Math.sqrt(Math.pow(x-xL,2)+1+Math.pow(tL+1-t,2));
		double uHLL = u[xL+1][y-1][tL];		double z5 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1+Math.pow(t-tL,2));
		double uHLH = temp[xL+1][y-1];		double z6 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1+Math.pow(tL+1-t,2));
		double uHHL = u[xL+1][y+1][tL];		double z7 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1+Math.pow(t-tL,2));
		double uHHH = temp[xL+1][y+1];		double z8 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1+Math.pow(tL+1-t,2));
		
		double z = z1+z2+z3+z4+z5+z6+z7+z8;
		double uC = (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4+uHLL*z5+uHLH*z6+uHHL*z7+uHHH*z8)/z;
		
		double uLLLL = u[xL-1][y-1][tL];
		double uLLLH = temp[xL-1][y-1];
		double uLLHL = u[xL-1][y+1][tL];
		double uLLHH = temp[xL-1][y+1];
		double uLHLL = u[xL][y-1][tL];
		double uLHLH = temp[xL][y-1];
		double uLHHL = u[xL][y+1][tL];
		double uLHHH = temp[xL][y+1];
		double uL = (uLLLL*z1+uLLLH*z2+uLLHL*z3+uLLHH*z4+uLHLL*z5+uLHLH*z6+uLHHL*z7+uLHHH*z8)/z;
		
		double uH;
		if(xL >= lx-2){
			uH = 0;
		}
		else{
			double uHLLL = u[xL+1][y-1][tL];
			double uHLLH = temp[xL+1][y-1];
			double uHLHL = u[xL+1][y+1][tL];
			double uHLHH = temp[xL+1][y+1];
			double uHHLL = u[xL+2][y-1][tL];
			double uHHLH = temp[xL+2][y-1];
			double uHHHL = u[xL+2][y+1][tL];
			double uHHHH = temp[xL+2][y+1];
			uH = (uHLLL*z1+uHLLH*z2+uHLHL*z3+uHLHH*z4+uHHLL*z5+uHHLH*z6+uHHHL*z7+uHHHH*z8)/z;
		}
		return c2*((uH+uL)/2-uC);
	}
	
	double u2(int x, double y, double t){					//2D
		int yL = (int)Math.floor(y);
		int tL = (int)Math.floor(t);
		if(yL <= 0 || yL >= ly-1){
			return 0;
		}
		
		double uLLL = u[x-1][yL][tL];		double z1 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(t-tL,2));			//(unitLength)/r
		double uLLH = temp[x-1][yL];		double z2 = 1/Math.sqrt(1+Math.pow(y-yL,2)+1+Math.pow(tL+1-t,2));
		double uLHL = u[x+1][yL][tL];		double z3 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(t-tL,2));
		double uLHH = temp[x+1][yL];		double z4 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(tL+1-t,2));
		double uHLL = u[x-1][yL+1][tL];		double z5 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(t-tL,2));
		double uHLH = temp[x-1][yL+1];		double z6 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(tL+1-t,2));
		double uHHL = u[x+1][yL+1][tL];		double z7 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(t-tL,2));
		double uHHH = temp[x+1][yL+1];		double z8 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(tL+1-t,2));
		
		double z = z1+z2+z3+z4+z5+z6+z7+z8;
		double uC = (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4+uHLL*z5+uHLH*z6+uHHL*z7+uHHH*z8)/z;
		
		double uLLLL = u[x-1][yL-1][tL];
		double uLLLH = temp[x-1][yL-1];
		double uLLHL = u[x+1][yL-1][tL];
		double uLLHH = temp[x+1][yL-1];
		double uLHLL = u[x-1][yL][tL];
		double uLHLH = temp[x-1][yL];
		double uLHHL = u[x+1][yL][tL];
		double uLHHH = temp[x+1][yL];
		double uL = (uLLLL*z1+uLLLH*z2+uLLHL*z3+uLLHH*z4+uLHLL*z5+uLHLH*z6+uLHHL*z7+uLHHH*z8)/z;
		
		double uH;
		if(yL >= ly-2){
			uH = 0;
		}
		else{
			double uHLLL = u[x-1][yL+1][tL];
			double uHLLH = temp[x-1][yL+1];
			double uHLHL = u[x+1][yL+1][tL];
			double uHLHH = temp[x+1][yL+1];
			double uHHLL = u[x-1][yL+2][tL];
			double uHHLH = temp[x-1][yL+2];
			double uHHHL = u[x+1][yL+2][tL];
			double uHHHH = temp[x+1][yL+2];
			uH = (uHLLL*z1+uHLLH*z2+uHLHL*z3+uHLHH*z4+uHHLL*z5+uHHLH*z6+uHHHL*z7+uHHHH*z8)/z;
		}
		return c2*((uH+uL)/2-uC);
		
		/*int yL = (int)Math.floor(y);
		int tL = (int)Math.floor(t);
		if(yL <= 0 || yL >= ly-2){
			return 0;
		}
		double uLow = ((y-yL)*u[x][yL+1][tL]+(1-y+yL)*u[x][yL][tL]);
		double uHigh = ((y-yL)*temp[x][yL+1]+(1-y-yL)*temp[x][yL]);
		double uC = (t-tL)*uHigh+(1-t+tL)*uLow;
		double uLeftLow = ((y-yL)*u[x][yL][tL]+(1-y-yL)*u[x][yL-1][tL]);
		double uLeftHigh = ((y-yL)*temp[x][yL]+(1-y-yL)*temp[x][yL-1]);
		double uL = (t-tL)*uLeftHigh+(1-t+tL)*uLeftLow;
		double uRightLow = ((y-yL)*u[x][yL+2][tL]+(1-y-yL)*u[x][yL+1][tL]);
		double uRightHigh = ((y-yL)*temp[x][yL+2]+(1-y-yL)*temp[x][yL+1]);
		double uR = (t-tL)*uRightHigh+(1-t+tL)*uRightLow;
		
		double uLowUp = ((y-yL)*u[x+1][yL+1][tL]+(1-y-yL)*u[x+1][yL][tL]);
		double uHighUp = ((y-yL)*temp[x+1][yL+1]+(1-y-yL)*temp[x+1][yL]);
		double uCUp = (uHighUp+uLowUp)/4;
		double uLowDn = ((y-yL)*u[x-1][yL+1][tL]+(1-y-yL)*u[x-1][yL][tL]);
		double uHighDn = ((y-yL)*temp[x-1][yL+1]+(1-y-yL)*temp[x-1][yL]);
		double uCDn = (uHighDn+uLowDn)/4;
		
		return c2*((uL+uR+uCUp+uCDn)/3-uC);*/
	}

	double u2(double x, int y){								//2D
		int xL = (int)Math.floor(x);
		if(xL <= 0 || xL >= lx-1){
			return 0;
		}
		
		double uLL = temp[xL][y-1];			double z1 = 1/Math.sqrt(Math.pow(x-xL,2)+1);
		double uLH = temp[xL][y+1];			double z2 = 1/Math.sqrt(Math.pow(x-xL,2)+1);
		double uHL = temp[xL+1][y-1];		double z3 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1);
		double uHH = temp[xL+1][y+1];		double z4 = 1/Math.sqrt(Math.pow(xL+1-x,2)+1);
		
		double z = z1+z2+z3+z4;
		double uC = (uLL*z1+uLH*z2+uHL*z3+uHH*z4)/z;
		
		double uLLL = temp[xL-1][y-1];
		double uLLH = temp[xL-1][y+1];
		double uLHL = temp[xL][y-1];
		double uLHH = temp[xL][y+1];
		double uL = (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4)/z;
		
		double uH;
		if(xL >= lx-2){
			uH = 0;
		}
		else{
			double uHLL = temp[xL+1][y-1];
			double uHLH = temp[xL+1][y+1];
			double uHHL = temp[xL+2][y-1];
			double uHHH = temp[xL+2][y+1];
			uH = (uHLL*z1+uHLH*z2+uHHL*z3+uHHH*z4)/z;
		}
		return c2*((uH+uL)/2-uC);
		
		/*int xL = (int)Math.floor(x);
		if(xL <= 0 || xL >= lx-2){
			return 0;
		}
		double uC = (x-xL)*temp[xL+1][y]+(1-x+xL)*temp[xL][y];
		double uL = (x-xL)*temp[xL][y]+(1-x+xL)*temp[xL-1][y];
		double uR = (x-xL)*temp[xL+2][y]+(1-x+xL)*temp[xL+1][y];
		double uCUp = ((x-xL)*temp[xL+1][y+1]+(1-x+xL)*temp[xL][y+1])/2;
		double uCDn = ((x-xL)*temp[xL+1][y-1]+(1-x+xL)*temp[xL][y-1])/2;
		return c2*((uL+uR+uCUp+uCDn)/3-uC);*/
	}
	
	double u2(int x, double y){								//2D
		int yL = (int)Math.floor(y);
		if(yL <= 0 || yL >= ly-1){
			return 0;
		}
		
		double uLL = temp[x-1][yL];			double z1 = 1/Math.sqrt(Math.pow(y-yL,2)+1);
		double uLH = temp[x+1][yL];			double z2 = 1/Math.sqrt(Math.pow(y-yL,2)+1);
		double uHL = temp[x-1][yL+1];		double z3 = 1/Math.sqrt(Math.pow(yL+1-y,2)+1);
		double uHH = temp[x+1][yL+1];		double z4 = 1/Math.sqrt(Math.pow(yL+1-y,2)+1);
		
		double z = z1+z2+z3+z4;
		double uC = (uLL*z1+uLH*z2+uHL*z3+uHH*z4)/z;
		
		double uLLL = temp[x-1][yL-1];
		double uLLH = temp[x+1][yL-1];
		double uLHL = temp[x-1][yL];
		double uLHH = temp[x+1][yL];
		double uL = (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4)/z;
		
		double uH;
		if(yL >= ly-2){
			uH = 0;
		}
		else{
			double uHLL = temp[x-1][yL+1];
			double uHLH = temp[x+1][yL+1];
			double uHHL = temp[x-1][yL+2];
			double uHHH = temp[x+1][yL+2];
			uH = (uHLL*z1+uHLH*z2+uHHL*z3+uHHH*z4)/z;
		}
		return c2*((uH+uL)/2-uC);
		
		/*int yL = (int)Math.floor(y);
		if(yL <= 0 || yL >= ly-2){
			return 0;
		}
		double uC = (y-yL)*temp[x][yL+1]+(1-y-yL)*temp[x][yL];
		double uL = (y-yL)*temp[x][yL]+(1-y-yL)*temp[x][yL-1];
		double uR = (y-yL)*temp[x][yL+2]+(1-y-yL)*temp[x][yL+1];
		double uCUp = ((y-yL)*temp[x+1][yL+1]+(1-y-yL)*temp[x+1][yL])/2;
		double uCDn = ((y-yL)*temp[x-1][yL+1]+(1-y-yL)*temp[x-1][yL])/2;
		return c2*((uL+uR+uCUp+uCDn)/3-uC);*/
	}
	
	void rkA(int x, int y, int t){							//2D
		double k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1, k2, k3, k4;
		k1x = h*this.u2(x,y,t-1);				k1y = h*this.u2(x,y,t-1);				k1 = (k1x+k1y)/2;
		k2x = h*this.u2(x+k1x/2,y,t-0.5);		k2y = h*this.u2(x,y+k1y/2,t-0.5);		k2 = (k2x+k2y)/2;
		k3x = h*this.u2(x+k2x/2,y,t-0.5);		k3y = h*this.u2(x,y+k2y/2,t-0.5);		k3 = (k3x+k3y)/2;
		k4x = h*this.u2(x+k3x,y);				k4y = h*this.u2(x,y+k3y);				k4 = (k4x+k4y)/2;
		uV[x][y][t] = uV[x][y][t-1]+(k1+2*k2+2*k3+k4)/6;
	}

	double u1(double x, int y, double t){					//2D
	
		int xL = (int)Math.floor(x);
		int tL = (int)Math.floor(t);
		if(xL <= 0 || xL >= lx){
			return 0;
		}
		
		double uLLL = uV[xL][y-1][tL];		double z1 = 1/Math.sqrt(1+Math.pow(x-xL,2)+Math.pow(t-tL,2));			//(unitLength)/r
		double uLLH = uV[xL][y-1][tL+1];	double z2 = 1/Math.sqrt(1+Math.pow(x-xL,2)+1+Math.pow(tL+1-t,2));
		double uLHL = uV[xL][y+1][tL];		double z3 = 1/Math.sqrt(1+Math.pow(x-xL,2)+Math.pow(t-tL,2));
		double uLHH = uV[xL][y+1][tL+1];	double z4 = 1/Math.sqrt(1+Math.pow(x-xL,2)+Math.pow(tL+1-t,2));
		double z5 = 1/Math.sqrt(1+Math.pow(xL+1-x,2)+Math.pow(t-tL,2));
		double z6 = 1/Math.sqrt(1+Math.pow(xL+1-x,2)+Math.pow(tL+1-t,2));
		double z7 = 1/Math.sqrt(1+Math.pow(xL+1-x,2)+Math.pow(t-tL,2));
		double z8 = 1/Math.sqrt(1+Math.pow(xL+1-x,2)+Math.pow(tL+1-t,2));
		double z = z1+z2+z3+z4+z5+z6+z7+z8;
		
		if(xL >= lx-1){
			return (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4)/z;
		}
		
		double uHLL = uV[xL+1][y-1][tL];	
		double uHLH = uV[xL+1][y-1][tL+1];	
		double uHHL = uV[xL+1][y+1][tL];	
		double uHHH = uV[xL+1][y+1][tL+1];	
		
		return (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4+uHLL*z5+uHLH*z6+uHHL*z7+uHHH*z8)/z;
		
		/*int xL = (int)Math.floor(x);
		int tL = (int)Math.floor(t);
		if(xL <= 0 || xL >= lx-1){
			return 0;
		}
		double uVLowUp = ((x-xL)*uV[xL+1][y+1][tL]+(1-x+xL)*uV[xL][y+1][tL]);
		double uVHighUp = ((x-xL)*uV[xL+1][y+1][tL+1]+(1-x+xL)*uV[xL][y+1][tL+1]);
		double uVLowDn = ((x-xL)*uV[xL+1][y-1][tL]+(1-x+xL)*uV[xL][y-1][tL]);
		double uVHighDn = ((x-xL)*uV[xL+1][y-1][tL+1]+(1-x+xL)*uV[xL][y-1][tL+1]);
		return (t-tL)*(uVHighUp+uVHighDn)/2+(1-t+tL)*(uVLowUp+uVLowDn)/2;			//Potential inaccuracy*/
	}
	
	double u1(int x, double y, double t){					//2D
		
		int yL = (int)Math.floor(y);
		int tL = (int)Math.floor(t);
		if(yL <= 0 || yL >= ly){
			return 0;
		}
		
		double uLLL = uV[x-1][yL][tL];		double z1 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(t-tL,2));			//(unitLength)/r
		double uLLH = uV[x-1][yL][tL+1];	double z2 = 1/Math.sqrt(1+Math.pow(y-yL,2)+1+Math.pow(tL+1-t,2));
		double uLHL = uV[x+1][yL][tL];		double z3 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(t-tL,2));
		double uLHH = uV[x+1][yL][tL+1];	double z4 = 1/Math.sqrt(1+Math.pow(y-yL,2)+Math.pow(tL+1-t,2));
		double z5 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(t-tL,2));
		double z6 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(tL+1-t,2));
		double z7 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(t-tL,2));
		double z8 = 1/Math.sqrt(1+Math.pow(yL+1-y,2)+Math.pow(tL+1-t,2));
		double z = z1+z2+z3+z4+z5+z6+z7+z8;
		
		if(yL >= ly-1){
			return (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4)/z;
		}
		
		double uHLL = uV[x-1][yL+1][tL];
		double uHLH = uV[x-1][yL+1][tL+1];
		double uHHL = uV[x+1][yL+1][tL];
		double uHHH = uV[x+1][yL+1][tL+1];
		
		return (uLLL*z1+uLLH*z2+uLHL*z3+uLHH*z4+uHLL*z5+uHLH*z6+uHHL*z7+uHHH*z8)/z;
		
		/*int yL = (int)Math.floor(y);
		int tL = (int)Math.floor(t);
		if(yL <= 0 || yL >= ly-1){
			return 0;
		}
		double uVLowUp = ((y-yL)*uV[x+1][yL+1][tL]+(1-y-yL)*uV[x+1][yL][tL]);
		double uVHighUp = ((y-yL)*uV[x+1][yL+1][tL+1]+(1-y-yL)*uV[x+1][yL][tL+1]);
		double uVLowDn = ((y-yL)*uV[x-1][yL+1][tL]+(1-y-yL)*uV[x-1][yL][tL]);
		double uVHighDn = ((y-yL)*uV[x-1][yL+1][tL+1]+(1-y-yL)*uV[x-1][yL][tL+1]);
		return (t-tL)*(uVHighUp+uVHighDn)/2+(1-t+tL)*(uVLowUp+uVLowDn)/2;			//Potential inaccuracy*/
	}

	void rkV(int x, int y, int t){
		double k1x, k2x, k3x, k4x, k1y, k2y, k3y, k4y, k1, k2, k3, k4;
		k1x = h*uV[x][y][t-1];					k1y = h*uV[x][y][t-1];					k1 = (k1x+k1y)/2;
		k2x = h*this.u1(x+k1x/2,y,t-0.5);		k2y = h*this.u1(x,y+k1y/2,t-0.5);		k2 = (k2x+k2y)/2;
		k3x = h*this.u1(x+k2x/2,y,t-0.5);		k3y = h*this.u1(x,y+k2y/2,t-0.5);		k3 = (k3x+k3y)/2;
		k4x = h*this.u1(x+k3x,y,t);				k4y = h*this.u1(x,y+k3y,t);				k4 = (k4x+k4y)/2;
		u[x][y][t] = u[x][y][t-1]+(k1+2*k2+2*k3+k4)/6;
	}

	//Sampler Theoretically 2D
	void samp(int t){
		double x = lx*Math.sqrt(2)/2;
		double y = ly*Math.sqrt(2)/2;
		int xL = (int)Math.floor(x);	int yL = (int)Math.floor(y);
		sample[t] = (x-xL)*(y-yL)*u[xL+1][yL+1][t] + (1-x+xL)*(y-yL)*u[xL][yL+1][t] + (x-xL)*(1-y+yL)*u[xL+1][yL][t] + (1-x+xL)*(1-y+yL)*u[xL][yL][t];
	}

	private void runSimulation() throws IOException {
		PrintStream spacialTemp1 = new PrintStream(new File("out/2DSpacialTemp1.txt"));
		PrintStream spacialTemp2 = new PrintStream(new File("out/2DSpacialTemp2.txt"));
		for (int t = 1; t < steps - 1; t++) {
			for (int x = 1; x < lx - 1; x++) {
				for (int y = 1; y < ly - 1; y++) {
					temp[x][y] = u[x][y][t - 1] + uV[x][y][t - 1] * h + uA[x][y][t - 1] * Math.pow(h, 2) / 2; // temp
																														// u(x,y,t)
																														// for
																														// calculating
																														// uA(x,y,t)
					if (t == snap1) {
						spacialTemp1.printf("%.4f ", temp[x][y]);
					} else if (t == snap2) {
						spacialTemp2.printf("%.4f ", temp[x][y]);
					}
				}
				if (t == snap1) {
					spacialTemp1.printf("%n");
				} else if (t == snap2) {
					spacialTemp2.printf("%n");
				}
			}
			for (int x = 1; x < lx - 2; x++) {
				for (int y = 1; y < ly - 2; y++) {
					uA[x][y][t] = u2(x, y); /* In terms of temp */
					rkA(x, y, t); /* Defines uV(x,t) */
					rkV(x, y, t); /* Defines u(x,t) */
				}
			}
			for (int x = 1; x < lx - 2; x++) {
				for (int y = 1; y < ly - 2; y++) {
					uA[x][y][t] = u2(x, y, t); /* In terms of u */
				}
			}
			samp(t); /* Defines u(lsqrt(2)/2,t); */
			time += h;
			// output.printf("%f %f %n", time, sample[t]);

			writeSurface(t);
		}
	}

	private void writeResults() throws FileNotFoundException {
		PrintStream spacial1 = new PrintStream(new File("out/2DSpacial1.txt"));
		PrintStream spacialV1 = new PrintStream(new File("out/2DSpacialV1.txt"));
		PrintStream spacialA1 = new PrintStream(new File("out/2DSpacialA1.txt"));
		PrintStream spacial2 = new PrintStream(new File("out/2DSpacial2.txt"));
		PrintStream spacialV2 = new PrintStream(new File("out/2DSpacialV2.txt"));
		PrintStream spacialA2 = new PrintStream(new File("out/2DSpacialA2.txt"));
		PrintStream outputT = new PrintStream(new File("out/2DTime.txt"));

		for (int y = 0; y < ly; y++) { // plots u(x,y) at const t
			for (int x = 0; x < lx; x++) {
				spacial1.printf("%.4f ", u[x][y][snap1]);
				spacialV1.printf("%.4f ", uV[x][y][snap1]);
				spacialA1.printf("%.4f ", uA[x][y][snap1]);
				spacial2.printf("%.4f ", u[x][y][snap2]);
				spacialV2.printf("%.4f ", uV[x][y][snap2]);
				spacialA2.printf("%.4f ", uA[x][y][snap2]);
			}
			spacial1.printf("%n");
			spacialV1.printf("%n");
			spacialA1.printf("%n");
			spacial2.printf("%n");
			spacialV2.printf("%n");
			spacialA2.printf("%n");
		}
		time = 0;
		for (int t = 0; t < steps; t++) { // plots u(t) at const x,y
			outputT.printf("%.3f %.4f %n", time, u[snapX][snapY][t]);
			time += h;
		}
	}

	public static void main(String[] args) throws IOException
	{
		
		Wave2D a = new Wave2D(21, 21, 10, 0.01);


		// set the initial t=0 surface shape
		a.initGaussian(10,10,1.0);
		//a.initPyramid(10,10,1.0);
		
		a.writeSurface(0);

		a.runSimulation();

		a.writeResults();
		
		// for(int x=0; x<l; x++){
		// 	output.printf("%d ", x);
		// 	for(int t=1; t<20; t++){
		// 		output.printf("%f ", a.u[x][t]);
		// 	}
		// 	output.println();
		// }

		/*int tOutStart = steps - (100/tOutZoom)*800;			Output for zoomed region 1D
		for(int t=tOutStart; t<steps; t+=100/tOutZoom){
			for(int x=0; x<l; x++){
				output.printf("%f ", a.u[x][t]);
			}
			output.println();
		}*/


		/*Complex[] complexSamples = new Complex[steps];
		for(int t=0; t<steps; t++){
			complexSamples[t] = new Complex(a.sample[t],0);
		}
		Complex[] comp = FFT.fft(complexSamples);
		PrintStream fft2 = new PrintStream(new File("out/fft2.txt"));
		PrintStream sampFile2 = new PrintStream(new File("out/SampFile2.txt"));
		time = 0;
		for(int t=500000; t<steps; t++){
			fft2.println(time+" "+comp[t].abs());
			sampFile2.println(time+" "+a.sample[t]);
			time+=h;
		}*/
	}

}
