package igtools.analyses.distribution.peaks;

import igtools.common.distributions.DistributionUtils;
import igtools.common.distributions.EstimatedDistribution;

import java.util.Arrays;
import java.util.TreeMap;

public class ExpDistributionPeakFinder {
	
	public static double expGDistrFirst(double[][] idistr){
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		//double[][] odistr = new double[idistr.length][2];
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(idistr);
			return estimator.getValue(idistr[0][0]);
			
		}catch(Exception e){
			gestimator.estimateDistrParameter(idistr);
			return gestimator.getValue(idistr[0][0]);
		}
	}
	
	public static double[][] expGDistr(double[][] idistr){
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		double[][] odistr = new double[idistr.length][2];
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(idistr);
			for(int i=0; i<idistr.length; i++){
				odistr[i][0] = idistr[i][0];
				odistr[i][1] =estimator.getValue(idistr[i][0]);
			}
		}catch(Exception e){
			gestimator.estimateDistrParameter(idistr);
			for(int i=0; i<idistr.length; i++){
				odistr[i][0] = idistr[i][0];
				odistr[i][1] =gestimator.getValue(idistr[i][0]);
			}
		}
		return odistr;
	}
	
	
	
	
	
	
	public static double getAreaValue(double perc, double[] ovalues){
		double sum = 0.0;
		for(int i=0; i<ovalues.length; i++)
			sum += ovalues[i];
		
		double psum = 0.0;
		for(int i=ovalues.length-1; i>=0; i--){
			psum += ovalues[i];
			if(psum/sum > perc){
				return ovalues[i];
			}
		}
		return ovalues[0];
	}
	public static double getPercentileValue(double percentile, double[] ovalues){
		int p = 0;
		for(int i=ovalues.length-1; i>=0; i--){
			p++;
			if(p/(double)ovalues.length >= percentile){
				return ovalues[i];
			}
		}
		return ovalues[0];
	}
	
	
	
	
	public static TreeMap<Double,Double> findPeakPositions2(double[][] distr, double resampleP){
//		double[][] rdistr = new double[distr.length][2];
//		for(int i=0; i<distr.length; i++){
//			rdistr[i][0] = distr[i][0];
//			rdistr[i][1] = distr[i][1];
//		}
//		DistributionUtils.resample(rdistr, resampleP);
//		rdistr = expGDistr(rdistr);
//		
//		for(int i=0; i<distr.length; i++){
//			//rdistr[i][0] = distr[i][0];
//			rdistr[i][1] = distr[i][1] - rdistr[i][1];  
//		}
//		distr = rdistr;
		
		TreeMap<Double,Double> ret = new TreeMap<Double,Double>();
		
		int wlength = (int)Math.ceil(((double)distr.length)*resampleP/2.0);
		
		
		double avg, count, sd;
		
		for(int i=0; i<distr.length; i++){
			avg = 0.0;
			count = 0.0;
			sd = 0.0;
			
			for(int j=i-wlength; j<i+wlength; j++){
				if(j>=0 && j<distr.length){
					count++;
					avg += distr[j][1];
				}
			}
			avg = avg / count;
			
//			for(int j=i-wlength; j<i+wlength; j++){
//				if(j>=0 && j<distr.length){
//					sd += (distr[i][1]-avg) * (distr[i][1]-avg);
//				}
//			}
//			sd = Math.sqrt(sd / count);
			
			if(distr[i][1] > avg){
				//ret.put(distr[i][0], Math.log(distr[i][0])* distr[i][1]);
				ret.put(distr[i][0],  distr[i][1] - avg);
			}
			else{
				ret.put(distr[i][0], distr[i][1] * (-1.0));
			}
		}
		
		
		return ret;
	}
	
	
	
	
	public static TreeMap<Double,Double> findPeakPositions(double[][] distr, double resampleP){
		
//		for(int i=0; i<distr.length; i++){
//			System.out.println(i+" "+distr[i][0]+" "+distr[i][1]);
//		}
		
		double[][] rdistr = new double[distr.length][2];
		for(int i=0; i<distr.length; i++){
			rdistr[i][0] = distr[i][0];
			rdistr[i][1] = distr[i][1];
		}
		DistributionUtils.resample(rdistr, resampleP);
		//double[][] rdistr = expGDistr(distr);
		
//		for(int i=0; i<rdistr.length; i++){
//			System.out.println(i+" "+rdistr[i][0]+" "+rdistr[i][1]);
//			if(rdistr[i][0] >=20)
//				break;
//		}
    	
		int c1 = 0;
		double c2 = 0.0;
		
		double[] evalues = new double[distr.length];
		for(int i=0; i<distr.length; i++){
			evalues[i] = rdistr[i][1];
		}
		
		for(int i=0; i<distr.length; i++){
			if(evalues[i] < distr[i][1]){
				//evalues[i] = (distr[i][1] - evalues[i]) / distr[i][1];
				if(distr[i][0]<=20){
					System.out.println(distr[i][1] +" - "+ evalues[i] +" "+(distr[i][1] - evalues[i]));
				}
				evalues[i] = (distr[i][1] - evalues[i]);
				c1++;
			}
			else{
				evalues[i] = 0.0;
			}
		}
		
		
		TreeMap<Double,Double> ret = new TreeMap<Double,Double>();
		
		double[] oevalues = new double[c1];
		for(int i=0,j=0; i<evalues.length; i++){
			if(evalues[i] > 0.0){
				oevalues[j] =evalues[i];
				j++;
			}
		}
		Arrays.sort(oevalues);
		//double thr = getAreaValue(0.025, oevalues);
		//double thr = getPercentileValue(((double)oevalues.length)/((double)evalues.length), oevalues);
		//double thr = oevalues[(int)Math.floor(((double)oevalues.length)/2.0)];
		//System.out.println(((double)oevalues.length)/((double)evalues.length));
		double thr = igtools.common.util.AvgSD.avg(oevalues) + igtools.common.util.AvgSD.sd(oevalues);
		
		for(int i=0; i<distr.length; i++){
			if(evalues[i] >= thr){
				ret.put(distr[i][0], evalues[i]);
			}
			else{
				ret.put(distr[i][0], evalues[i] * (-1.0));
			}
		}
		
		return ret;
	}

}
