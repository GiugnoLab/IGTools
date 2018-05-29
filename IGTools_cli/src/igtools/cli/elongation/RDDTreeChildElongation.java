package igtools.cli.elongation;

import igtools.analyses.elongation.ElongationAlgorithm;
import igtools.analyses.elongation.L2RElongation;
import igtools.analyses.elongation.ElongationAlgorithm.ElongationEvent;
import igtools.analyses.elongation.ElongationAlgorithm.ElongationListener;
import igtools.analyses.elongation.ElongationAlgorithm.IteratorMeasure;
import igtools.analyses.recurrences.distances.ProperMinimalRecurrenceDistancesExtractor;
import igtools.analyses.recurrences.expcomp.KExpComp;
import igtools.common.distributions.DistributionUtils;
import igtools.common.distributions.EstimatedDistribution;
import igtools.common.distributions.distance.DistributionDistance;
import igtools.common.distributions.distance.KLDistance;
import igtools.common.distributions.distance.KSDistance;
import igtools.common.distributions.distance.PointDiffSUMDistance;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.AvgSD;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.util.Map;
import java.util.TreeMap;

public class RDDTreeChildElongation {

	public static void usage(){
		System.out.println("Usage: cmd b3seq nelsa <kl|ks|sum> startK");
	}
	
	
	enum DistType {KL,KS,SUM};
	
	public static void main(String[] args){
		String a_b3seq = "";
		String a_nelsa = "";
		DistType a_dtype = DistType.KL; 
		int k = 0;
		
		try{
			a_b3seq = args[0];
			a_nelsa = args[1];
			
			
			if(args[2].compareTo("kl") == 0){
				a_dtype = DistType.KL;
			}
			else if(args[2].compareTo("ks") == 0){
				a_dtype = DistType.KS;
			} 
			else if(args[2].compareTo("sum") == 0){
				a_dtype = DistType.SUM;
			}
			else
				throw new Exception();
			
			k = Integer.parseInt(args[3]);
			if(k<1)
				throw new Exception();
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_b3seq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_nelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			DistributionDistance dd;
			switch(a_dtype){
			case KL:
				dd = new KLDistance.MaxKLDistance();
				break;
			case KS:
				dd = new KSDistance();
				break;
			case SUM:
				dd = new PointDiffSUMDistance();
				break;
			default:
				dd = null;
				break;
			}
			
			final DistributionDistance fdd = dd;
			final int fstartK = k;
			
			
			L2RElongation ealgo = new L2RElongation(nelsa);
			ealgo.elongate(k, 
					new ElongationAlgorithm.IteratorMeasure() {
						@Override
						public Double value(IELSAIterator it) {
							return KExpComp.distanceToExponential(it, fdd, 0.01);
						}
					}, 
					new ElongationAlgorithm.ElongationListener() {
						
						@Override
						public void event(IELSAIterator it, ElongationEvent ev, Double parentValue, Double itValue) {
							if(ev == ElongationEvent.EXTENDED){
								System.out.println("E: "+B3Nucleotide.toString(it.kmer())+" "+parentValue+" "+itValue);
							}
							else{
								if(it.k() != fstartK){
									System.out.println("S: "+B3Nucleotide.toString(it.kmer())+" "+itValue);
								}
							}
						}
							
					});
			
			
//			int i = 0;
//			it = nelsa.begin(k);
//			//double value;
//			while(it.next()){
//				if(it.multiplicity() > 1){
//					Map<Double,Double> it_distr = ProperMinimalRecurrenceDistancesExtractor.factory(false, true).recurrenceDistanceDistributionMap(it);
//					DistributionUtils.normalize(it_distr);
//					
//					EstimatedDistribution estimator = new EstimatedDistribution.FExponentialBased();
//					EstimatedDistribution gestimator = new EstimatedDistribution.GeometricBased();
//					
//					Map<Double,Double> e_distr = new TreeMap<Double,Double>();
//					double[][] aa = DistributionUtils.toArray(it_distr);
//					
//					
//					try{//sometimes it sucks
//						estimator.estimateDistrParameter(aa);
//						for(Map.Entry<Double, Double> en : it_distr.entrySet()){
//							e_distr.put(en.getKey(), estimator.getValue(en.getKey()));
//						}
//					}catch(Exception e){
//						gestimator.estimateDistrParameter(aa);
//						for(Map.Entry<Double, Double> en : it_distr.entrySet()){
//							e_distr.put(en.getKey(), gestimator.getValue(en.getKey()));
//						}
//					}	
//				}
//			}
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
