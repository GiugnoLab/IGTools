package igtools.common.util;

public class AvgSD {

	
	public static double avg(double a[]){
		double avg = 0.0;
		for(int i=0; i<a.length; i++){
			avg += a[i];
		}
		return avg / ((double)a.length);
	}
	
	public static double avg(double a[], double thr){
		double avg = 0.0;
		double count = 0.0;
		for(int i=0; i<a.length; i++){
			if(a[i] > thr){
				avg += a[i];
				count++;
			}
		}
		return avg / count;
	}
	
	
	public static double sd(double[] a, double avg, double thr){
		double sd = 0.0;
		for(int i=0; i<a.length; i++){
			if(a[i] > thr){
				sd += (a[i] - avg) * (a[i] - avg);
			}
		}
		return Math.sqrt(sd / ((double) a.length));
	}
	public static double sd(double[] a, double avg){
		double sd = 0.0;
		for(int i=0; i<a.length; i++){
			sd += (a[i] - avg) * (a[i] - avg);
		}
		return Math.sqrt(sd / ((double) a.length));
	}
	public static double sd(double[] a){
		double avg = avg(a);
		return sd(a,avg);
	}
	
	
	public static void avg_sd(double[] a, double[] o_avg_sd){
		o_avg_sd[0] = 0.0;
		o_avg_sd[1] = 0.0;
		
		double e_x = 0.0;
		double e_x2 = 0.0;
		
		for(int i=0; i<a.length; i++){
			e_x += a[i];
			e_x2 += a[i] * a[i];
		}
		
		o_avg_sd[0] = e_x / ((double)(a.length));
		
		o_avg_sd[1] = Math.sqrt(( e_x2 / ((double)(a.length)) )  - (o_avg_sd[0] * o_avg_sd[0]));
	}
	
	public static void avg_sd(double[] a, double thr, double[] o_avg_sd){
		o_avg_sd[0] = 0.0;
		o_avg_sd[1] = 0.0;
		
		double count = 0.0;
		double e_x = 0.0;
		double e_x2 = 0.0;
		
		for(int i=0; i<a.length; i++){
			if(a[i] >= thr){
				count += 1.0;
				e_x += a[i];
				e_x2 += a[i] * a[i];
			}
		}
		
		o_avg_sd[0] = e_x / count;
		
		o_avg_sd[1] = Math.sqrt(( e_x2 / count )  - (o_avg_sd[0] * o_avg_sd[0]));
	}
	
	//TODO others
}
