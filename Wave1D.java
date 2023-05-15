/*
Professor,
Presentation
https://docs.google.com/presentation/d/1VOkbA3EL5aSKoPNnr-o9xNYslc8rDHAxGbpeDjBm6vs/edit?usp=sharing

My father and I have discovered since class that the FFT was in fact working, but was obscured by extra data.
Attached to the Moodle are jpgs of the FFT zoomed into the meaningful section and the plot of the array that was input into the FFT method. 
The point is at l*sqrt(2)/2 and is plotted vs. time.

Also attached are the Complex and FFT classes that I got from the Prinston CS website.


*/

import java.io.*;
import java.util.*;
public class Wave1D{
	static double c, c2, h, time;
	static int steps, l, x0;
	double[][] u, uV, uA;
	double[] temp, sample;

	public Wave1D(){
		c2 = 15;	c = Math.sqrt(c2);	h = 0.01;	time = 0;
		steps = 524288;	l = 100;	x0 = (int)(0.15*l)/*(int)(Math.random()*(l-1)+1)*/;
		u = new double[l][steps];	u[x0][0] = 1;
		uV = new double[l][steps];
		uA = new double[l][steps];
		temp = new double[l];	temp[0] = 0;	temp[l-1] = 0;
		sample = new double[steps];
	}

	void squareWaveInit(int nodeCount) {
		double[] u0 = new double[l];
		double v = 0;
		for(int x=1; x<l; x++){
			if ( (x% (l/nodeCount)) == 0 ) {
				if (v>0.5) {
					v = 0;
				} else {
					v = 1;
				}
			}
			u[x][0] = v;
		}
	}

	void triangleWaveInit(int x1, int x2, int x3) {
		double[] u0 = new double[l];
		double v = 0;
		for(int x=1; x<l; x++){
			if ( x > x3) {
				u[x][0] = 0;
			} else if (x > x2) {
				u[x][0] = 1.0 - (x-x2) * u[x2][0]/(x3-x2);
			} else if ( x> x1) {
				u[x][0] =  (x-x1) * 1.0/(x2-x1) ;
			} else {
				u[x][0] = 0;
			}
		}
	}

	double u2(int x){
		return c2*(temp[x-1]+temp[x+1])/2-temp[x];
	}

	double u2(int x, int t){
		return c2*((u[x-1][t]+u[x+1][t])/2-u[x][t]);
	}

	double u2(double x, double t, boolean test){
		if(Math.floor(x-1) < 0){
			return 0;
		}
		if(Math.ceil(x+1) > l-1){
			return 0;
		}
		int xL = (int)Math.floor(x);	int xH = xL+1;
		int tL = (int)Math.floor(t);
		double uLow = ((x-xL)*u[xH][tL]+(1-x+xL)*u[xL][tL]);
		double uHigh = ((x-xL)*temp[xH]+(1-x+xL)*temp[xL]);
		double uC = (t-tL)*uHigh+(1-t+tL)*uLow;
		double uLeftLow = ((x-xL)*u[xL][tL]+(1-x+xL)*u[xL-1][tL]);
		double uLeftHigh = ((x-xL)*temp[xL]+(1-x+xL)*temp[xL-1]);
		double uL = (t-tL)*uLeftHigh+(1-t+tL)*uLeftLow;
		double uRightLow = ((x-xL)*u[xH+1][tL]+(1-x+xL)*u[xH][tL]);
		double uRightHigh = ((x-xL)*temp[xH+1]+(1-x+xL)*temp[xH]);
		double uR = (t-tL)*uRightHigh+(1-t+tL)*uRightLow;
		if(test) System.out.printf("uLow = %f %n", uLow);
		if(test) System.out.printf("uHigh = %f %n", uHigh);
		if(test) System.out.printf("uC = %f %n", uC);
		if(test) System.out.printf("uLeftLow = %f %n", uLeftLow);
		if(test) System.out.printf("uLeftHigh = %f %n", uLeftHigh);
		if(test) System.out.printf("uL = %f %n", uL);
		if(test) System.out.printf("uRightLow = %f %n", uRightLow);
		if(test) System.out.printf("uRightHigh = %f %n", uRightHigh);
		if(test) System.out.printf("uR = %f %n", uR);
		if(test) System.out.printf("xLtL = %f %n", u[xL][tL]);
		if(test) System.out.printf("xLtH = %f %n", temp[xL]);
		if(test) System.out.printf("xHtL = %f %n", u[xH][tL]);
		if(test) System.out.printf("xHtH = %f %n", temp[xH]);
		return c2*((uL+uR)/2-uC);
	}

	double u2(double x, int t, boolean test){
		if(Math.floor(x-1) < 0 || Math.ceil(x+1) > l-1){
			return 0;
		}
		int xL = (int)Math.floor(x);	int xH = xL+1;
		double uC = (x-xL)*temp[xH]+(1-x+xL)*temp[xL];
		double uL = (x-xL)*temp[xL]+(1-x+xL)*temp[xL-1];
		double uR = (x-xL)*temp[xH+1]+(1-x+xL)*temp[xH];
		if(test) System.out.printf("uC = %f %n", uC);
		if(test) System.out.printf("uL = %f %n", uL);
		if(test) System.out.printf("uR = %f %n", uR);
		if(test) System.out.printf("xL-1 = %f %n", temp[xL-1]);
		if(test) System.out.printf("xL = %f %n", temp[xL]);
		if(test) System.out.printf("xH = %f %n", temp[xH]);
		if(test) System.out.printf("xH = %f %n", temp[xH]);
		if(test) System.out.printf("xH+1 = %f %n", temp[xH+1]);
		return c2*((uL+uR)/2-uC);
	}

	void rkA(int x, int t){
		double k1, k2, k3, k4;
		boolean test = false;
		if(x == 96 && t == 1){
			test=false;
		}
		k1 = h*this.u2(x,t-1);
		k2 = h*this.u2(x+k1/2,t-0.5, false);
		k3 = h*this.u2(x+k2/2,t-0.5, false);
		k4 = h*this.u2(x+k3,t, test);
		uV[x][t] = uV[x][t-1]+(k1+2*k2+2*k3+k4)/6;
		/*if(test) System.out.printf("k1 = %f %n", k1);
		if(test) System.out.printf("k2 = %f %n", k2);
		if(test) System.out.printf("k3 = %f %n", k3);
		if(test) System.out.printf("k4 = %f %n", k4);
		if(test){
			System.exit(0);
		}*/
	}

	double u1(double x, double t){
		int xL = (int)Math.floor(x);	int xH = xL+1;
		int tL = (int)Math.floor(t);	int tH = tL+1;
		if(xL < 1 || xL > l-1 || xH > l-1 || xH <1){
			return 0;
		}
		double uVLow = ((x-xL)*uV[xH][tL]+(1-x+xL)*uV[xL][tL]);
		double uVHigh = ((x-xL)*uV[xH][tH]+(1-x+xL)*uV[xL][tH]);
		return (t-tL)*uVHigh+(1-t+tL)*uVLow;
	}

	void rkV(int x, int t){
		double k1, k2, k3, k4;
		boolean test = false;
		if(x == 96 && t == 2){
			test=false;
		}
		k1 = h*uV[x][t-1];
		k2 = h*this.u1(x+k1/2,t-0.5);
		k3 = h*this.u1(x+k2/2,t-0.5);
		k4 = h*this.u1(x+k3,t);
		u[x][t] = u[x][t-1]+(k1+2*k2+2*k3+k4)/6;
		if(test) System.out.printf("k1 = %f %n", k1);
		if(test) System.out.printf("k2 = %f %n", k2);
		if(test) System.out.printf("k3 = %f %n", k3);
		if(test) System.out.printf("k4 = %f %n", k4);
		if(test) System.out.printf("u[x][t-1] = %f %n", u[x][t-1]);
		if(test) System.out.printf("u[x][t] = %f %n", u[x][t]);
		if(test){
			System.exit(0);
		}
	}

	void samp(int t){
		double x = l*Math.sqrt(2)/2;
		int xL = (int)Math.floor(x);	int xH = xL+1;
		sample[t] = ((x-xL)*u[xH][t]+(1-x+xL)*u[xL][t]);
	}

	public static void main(String[] args) throws FileNotFoundException{
		PrintStream output = new PrintStream(new File("Wave1D.txt"));
		Wave1D a = new Wave1D();
		for(int t=0; t<steps; t++){						/*boundary conditions*/
			a.u[0][t] = 0;	a.u[l-1][t] = 0;
			a.uV[0][t] = 0;	a.uV[l-1][t] = 0;
			a.uA[0][t] = 0;	a.uA[l-1][t] = 0;
		}
		//output.printf("%d %f %n %d %f %n", x0, a.u[x0][0], 0, a.u[0][0]);
		for(int x=1; x<l; x++){						/*initial conditions*/
			if(x<=x0){
				a.u[x][0] = a.u[x0][0]/x0*x;
			}
			else{
				a.u[x][0] = a.u[x0][0]/(x0-l)*(x-l);
			}
			a.uV[x][0] = 0;
			a.uA[x][0] = 0;
			//output.printf("%d %f %n", x, a.u[x][0]);
		}

		// uncomment one of these at a time
		a.triangleWaveInit(0, 50, 100);
		//a.squareWaveInit(10);
		//a.triangleWaveInit(30, 50, 70);
		int tOutZoom = 10; // can be set to '1' or '10'

		//output.printf("%f %f %n", time, a.sample[0]);
		for(int t=1; t<steps-1; t++){
			a.temp[0] = 0; a.temp[l-1] = 0;
			for(int x=1; x<l-1; x++){
				a.temp[x] = a.u[x][t-1]+a.uV[x][t-1]*h+a.uA[x][t-1]*Math.pow(a.h,2)/2;			/*temp u(x,t) for calculating uA(x,t)*/
			}
			for(int x=1; x<l-2; x++){
				a.uA[x][t] = a.u2(x);								/*In terms of temp*/
				a.rkA(x,t);											/*Defines uV(x,t)*/
				a.rkV(x,t);											/*Defines u(x,t)*/
				a.uA[x][t] = a.u2(x,t);								/*In terms of u*/
			}
			a.samp(t);										/*Defines u(lsqrt(2)/2,t);*/
			time+=h;
			//output.printf("%f %f %n", time, a.sample[t]);
		}
		// for(int x=0; x<l; x++){
		// 	output.printf("%d ", x);
		// 	for(int t=1; t<20; t++){
		// 		output.printf("%f ", a.u[x][t]);
		// 	}
		// 	output.println();
		// }

		int tOutStart = steps - (100/tOutZoom)*800;
		for(int t=tOutStart; t<steps; t+=100/tOutZoom){
			for(int x=0; x<l; x++){
				output.printf("%f ", a.u[x][t]);
			}
			output.println();
		}


		Complex[] complexSamples = new Complex[steps];
		for(int t=0; t<steps; t++){
			complexSamples[t] = new Complex(a.sample[t],0);
		}
		Complex[] comp = FFT.fft(complexSamples);
		PrintStream fft = new PrintStream(new File("fft.txt"));
		PrintStream sampFile = new PrintStream(new File("SampFile.txt"));
		time = 0;
		for(int t=500000; t<steps; t++){
			fft.println(time+" "+comp[t].abs());
			sampFile.println(time+" "+a.sample[t]);
			time+=h;
		}
	}
}
