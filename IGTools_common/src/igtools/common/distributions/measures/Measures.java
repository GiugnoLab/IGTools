package igtools.common.distributions.measures;

import java.util.Map;

import igtools.common.distributions.DistributionUtils;

public class Measures {

	public static class DMeasures{
		public double mean = 0.0;
		public double sd = 0.0;
		public double skew = 0.0;
		public double kurtosis = 0.0;
		public void print(){
			System.out.println("mean "+mean);
			System.out.println("sd "+sd);
			System.out.println("skew "+skew);
			System.out.println("kurtosis "+kurtosis);
		}
	}
	
	
	
	public static void measure(Map<Double,Double> distr, DMeasures mes){
		double[][] aa = DistributionUtils.toArray(distr);
		measure(aa,mes);
	}
	public static void measureRaw(Map<Double,Double> distr, DMeasures mes){
		double[][] aa = DistributionUtils.toArray(distr);
		measureRaw(aa,mes);
	}
	
	
	public static void measure(double[][] distr, DMeasures mes){
		
		double p_x = distr[0][0];
		double count = 1.0;
		mes.mean = distr[0][1];
		for(int i=1; i<distr.length; i++){
			count++;
			for(double j=p_x+1.0; j<distr[i][0]; j++){
				mes.mean += 0.0;
				count++;
			}
			mes.mean += distr[i][1];
			p_x = distr[i][0];
		}
		mes.mean /= count;
		
		
		p_x = distr[0][0];
		mes.sd = (distr[0][1] - mes.mean)*(distr[0][1] - mes.mean);
		for(int i=1; i<distr.length; i++){
			for(double j=p_x+1.0; j<distr[i][0]; j++){
				mes.sd += (0.0-mes.mean)*(0.0-mes.mean);
			}
			mes.sd += (distr[i][1] - mes.mean)*(distr[i][1] - mes.mean);
			p_x = distr[i][0];
		}
		mes.sd = Math.sqrt(mes.sd/count);
		
		
		
		
		p_x = distr[0][0];
		mes.skew = Math.pow((distr[0][1] - mes.mean)/mes.sd,3);
		mes.kurtosis = Math.pow((distr[0][1] - mes.mean)/mes.sd,4);
		for(int i=1; i<distr.length; i++){
			for(double j=p_x+1.0; j<distr[i][0]; j++){
				mes.skew += Math.pow((0.0 - mes.mean)/mes.sd,3);
				mes.kurtosis += Math.pow((0.0 - mes.mean)/mes.sd,4);
			}
			mes.skew += Math.pow((distr[i][1] - mes.mean)/mes.sd,3);
			mes.kurtosis += Math.pow((distr[i][1] - mes.mean)/mes.sd,4);
			p_x = distr[i][0];
		}
		mes.sd = Math.sqrt(mes.sd/count);
	}
	
		
	public static void measureRaw(double[][] distr, DMeasures mes){
		mes.mean = 0.0;
		mes.sd = 0.0;
		mes.skew = 0.0;
		mes.kurtosis = 0.0;
		
		double count = 1.0;
		for(int i=0; i<distr.length; i++){
			count++;
			mes.mean += distr[i][1];
		}
		mes.mean /= count;
		
		
		for(int i=0; i<distr.length; i++){
			mes.sd += (distr[i][1] - mes.mean)*(distr[i][1] - mes.mean);
		}
		mes.sd = Math.sqrt(mes.sd/count);
		
		for(int i=0; i<distr.length; i++){
			mes.skew += Math.pow((distr[i][1] - mes.mean)/mes.sd,3);
			mes.kurtosis += Math.pow((distr[i][1] - mes.mean)/mes.sd,4);
		}
		mes.sd = Math.sqrt(mes.sd/count);
	}
}
