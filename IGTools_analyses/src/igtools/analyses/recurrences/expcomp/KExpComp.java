package igtools.analyses.recurrences.expcomp;

import java.util.Map;
import java.util.TreeMap;

import igtools.analyses.recurrences.distances.ProperMinimalRecurrenceDistancesExtractor;
import igtools.common.distributions.DistributionUtils;
import igtools.common.distributions.EstimatedDistribution;
import igtools.common.distributions.distance.DistributionDistance;
import igtools.dictionaries.elsa.IELSAIterator;

public class KExpComp {
	
	
	public static EstimatedDistribution getExponentialEstimator(Map<Double,Double> it_distr, double resampleP){
		
		double[][] aa = DistributionUtils.toArray(it_distr);
		
		if(resampleP>0.0){
			DistributionUtils.resample(aa, resampleP);
		}
		
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		//EstimatedDistribution gestimator = new EstimatedDistribution.ExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		
		//Map<Double,Double> e_distr = new TreeMap<Double,Double>();
		//double[][] aa = DistributionUtils.toArray(it_distr);
		
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(aa);
		}catch(Exception e){
			gestimator.estimateDistrParameter(aa);
		}
		
		return estimator;
	}
	
	public static Map<Double,Double> estimateExponential(Map<Double,Double> it_distr, double resampleP){
		
		double[][] aa = DistributionUtils.toArray(it_distr);
		
		if(resampleP>0.0){
			DistributionUtils.resample(aa, resampleP);
		}
		
		
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		//EstimatedDistribution gestimator = new EstimatedDistribution.ExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		
		
		Map<Double,Double> e_distr = new TreeMap<Double,Double>();
		//double[][] aa = DistributionUtils.toArray(it_distr);
		
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(aa);
			for(Map.Entry<Double, Double> en : it_distr.entrySet()){
				e_distr.put(en.getKey(), estimator.getValue(en.getKey()));
			}
		}catch(Exception e){
			gestimator.estimateDistrParameter(aa);
			for(Map.Entry<Double, Double> en : it_distr.entrySet()){
				e_distr.put(en.getKey(), gestimator.getValue(en.getKey()));
			}
		}
		
		return e_distr;
	}
	
	
	/**
	 * Return the distance to the estimated exponential distribution
	 * 
	 * @param it
	 * @param dd
	 * @return
	 */
	public static Double distanceToExponential(IELSAIterator it, DistributionDistance dd, double resampleP){
		if(it.multiplicity() > 1){
			Map<Double,Double> it_distr = ProperMinimalRecurrenceDistancesExtractor.factory(false, true).recurrenceDistanceDistributionMap(it);
			//DistributionUtils.normalize(it_distr);
			
			double[][] aa = DistributionUtils.toArray(it_distr);
			
			if(resampleP>0.0){
				DistributionUtils.resample(aa, resampleP);
			}
			
			
			EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
			//EstimatedDistribution gestimator = new EstimatedDistribution.ExponentialBased();
			EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
			
			
			Map<Double,Double> e_distr = new TreeMap<Double,Double>();
			//double[][] aa = DistributionUtils.toArray(it_distr);
			
			
			try{//sometimes it sucks
				estimator.estimateDistrParameter(aa);
				for(Map.Entry<Double, Double> en : it_distr.entrySet()){
					e_distr.put(en.getKey(), estimator.getValue(en.getKey()));
				}
			}catch(Exception e){
				gestimator.estimateDistrParameter(aa);
				for(Map.Entry<Double, Double> en : it_distr.entrySet()){
					e_distr.put(en.getKey(), gestimator.getValue(en.getKey()));
				}
			}
			
			
			DistributionUtils.normalize(it_distr);
			DistributionUtils.normalize(e_distr);
			
			return dd.distance(it_distr, e_distr);
		}
		return 0.0;
	}
	
	
	
	
	
	public static Double exponentialDistance(IELSAIterator it, DistributionDistance dd, double resampleP){
		if(it.multiplicity() > 1){
			Map<Double,Double> it_distr = ProperMinimalRecurrenceDistancesExtractor.factory(false, true).recurrenceDistanceDistributionMap(it);
			DistributionUtils.normalize(it_distr);
			
			double[][] aa = DistributionUtils.toArray(it_distr);
			
			if(resampleP>0.0){
				DistributionUtils.resample(aa, resampleP);
			}
			
			
			EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
			//EstimatedDistribution gestimator = new EstimatedDistribution.ExponentialBased();
			EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
			
			
			Map<Double,Double> e_distr = new TreeMap<Double,Double>();
			//double[][] aa = DistributionUtils.toArray(it_distr);
			
			
			try{//sometimes it sucks
				estimator.estimateDistrParameter(aa);
				for(Map.Entry<Double, Double> en : it_distr.entrySet()){
					e_distr.put(en.getKey(), estimator.getValue(en.getKey()));
				}
			}catch(Exception e){
				gestimator.estimateDistrParameter(aa);
				for(Map.Entry<Double, Double> en : it_distr.entrySet()){
					e_distr.put(en.getKey(), gestimator.getValue(en.getKey()));
				}
			}
			
			//DistributionUtils.normalize(e_distr);
			
			return dd.distance(e_distr, it_distr);
		}
		return 0.0;
	}
	
	
	
	
//	public static Double distanceToExponential(IELSAIterator it, DistributionDistance dd, double resampleP, double minCoRecc){
//		if(it.multiplicity() > 1){
//			double[][] it_distr = ProperMinimalRecurrenceDistancesExtractor.factory(false, true).recurrenceDistanceDistribution(it);
//			DistributionUtils.normalize(it_distr);
//			DistributionUtils.resample(it_distr, resampleP);
//			
//			double[][] e_distr = expGDistr(it_distr);
//			
//			
//			
//			
//			return dd.distance(it_distr, e_distr);
//		}
//		return 0.0;
//	}
//	
//	
//	public static boolean checkMinCoRecc(double[][] distr){
//		for(int i=0; i<distr.length; i++){
//			
//		}
//		return false;
//	}
	
	
	
	public static double[][] expGDistr(double[][] idistr){
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		double[][] odistr = new double[idistr.length][2];
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(idistr);
			for(int i=0; i<idistr.length; i++){
				odistr[i][1] =gestimator.getValue(odistr[i][0]);
			}
		}catch(Exception e){
			gestimator.estimateDistrParameter(idistr);
			for(int i=0; i<idistr.length; i++){
				odistr[i][1] =gestimator.getValue(odistr[i][0]);
			}
		}
		return odistr;
	}
	
	public static double expGDistrFirst(double[][] idistr){
		EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
		EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
		double[][] odistr = new double[idistr.length][2];
		
		try{//sometimes it sucks
			estimator.estimateDistrParameter(idistr);
			return estimator.getValue(odistr[0][0]);
			
		}catch(Exception e){
			gestimator.estimateDistrParameter(idistr);
			return gestimator.getValue(odistr[0][0]);
		}
	}
	
}
