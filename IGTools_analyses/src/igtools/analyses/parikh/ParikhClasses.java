package igtools.analyses.parikh;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.util.AvgSD;
import igtools.common.util.Maths;
import igtools.common.util.Vector4Int;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;

public class ParikhClasses {

	
	
	public static void getClassesNMultiplicities(
			NELSA nelsa,
			int k,
			Map<Vector4Int, Integer> parikhClasses,	
			Map<Vector4Int, Integer> parikhMults)
	{
		Vector4Int s_vect = new Vector4Int();
		Vector4Int c_vect;
		B3Nucleotide[] kmer = new B3Nucleotide[k];
		Integer g, m;
		
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			it.kmer(kmer);
			igtools.common.util.Parikh.getParikh(kmer, s_vect);
			
//			System.out.println(B3Nucleotide.toString(kmer)+ " " +s_vect);
			
			g = parikhClasses.get(s_vect);
			if(g == null){
//				System.out.println("new ");
				
				c_vect = s_vect.clone();
				
				parikhClasses.put(c_vect, 1);
				parikhMults.put(c_vect, it.multiplicity());
			}
			else{
				m = parikhMults.get(s_vect);
				
				parikhClasses.put(s_vect, g+1);
				parikhMults.put(s_vect, m + it.multiplicity());
			}
		}
	}
	
	
	
	
	
	public static void getParikhClasses(NELSA nelsa, int k, Set<Vector4Int> classes) {
		
		Vector4Int s_vect = new Vector4Int();
		B3Nucleotide[] kmer = new B3Nucleotide[k];
		
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			it.kmer(kmer);
			igtools.common.util.Parikh.getParikh(kmer, s_vect);
			classes.add(s_vect.clone());
		}
		
	}
	
	
	
	public static class PClassStats{
		public int nofKmers = 0;
		public int totMultiplicity = 0;
		public double tRatio = 0.0;
		public double seqCoverage = 0.0;
		public double punctualCoverageAvg = 0.0;
		public double punctualCoverageSD = 0.0;
		
		public PClassStats(){};
	}
	
	
	public static void parikhAnalysis(
			NELSA nelsa,
			int k,
			Map<Vector4Int, PClassStats> res,
			boolean getCoverage
			){
		
		Vector4Int s_vect = new Vector4Int();
		B3Nucleotide[] kmer = new B3Nucleotide[k];
		PClassStats cs;
		
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			it.kmer(kmer);
			igtools.common.util.Parikh.getParikh(kmer, s_vect);
			
			cs = res.get(s_vect);
			if(cs == null){
				cs = new PClassStats();
				res.put(s_vect.clone(), cs);
			}
			
			cs.nofKmers++;
			cs.totMultiplicity += it.multiplicity();
		}
		
		
		for(Map.Entry<Vector4Int, PClassStats> entry : res.entrySet()){
			entry.getValue().tRatio = ((double)entry.getValue().nofKmers) / ((double)Maths.nofAnagrams(entry.getKey().values, k)) ;  
		}
		
		
		if(getCoverage){
			
			double seqLength = nelsa.b3seq().length() - nelsa.b3seq().countBads();
			
			Vector4Int c;
			double[] cov = new double[nelsa.length()];
			int[] sa = nelsa.sa();
			
			for(Map.Entry<Vector4Int, PClassStats> entry : res.entrySet()){
				c = entry.getKey();
				cs = entry.getValue();
				Arrays.fill(cov, 0.0);
				
				it = nelsa.begin(k);
				while(it.next()){
					it.kmer(kmer);
					igtools.common.util.Parikh.getParikh(kmer, s_vect);
					
					if(s_vect.compareTo(c) == 0){
						for(int i=it.istart(); i<it.iend(); i++){
							for(int j=sa[i]; j<sa[i]+k; j++){
								cov[j]++;
							}
						}
					}
				}
				
				
				double scov = 0.0;
				for(int i=0; i<cov.length; i++){
					if(cov[i] > 0){
						scov++;
					}
				}
				cs.seqCoverage = scov / seqLength;
				
				cs.punctualCoverageAvg = AvgSD.avg(cov, 0.0);
				cs.punctualCoverageSD = AvgSD.sd(cov, cs.punctualCoverageAvg, 0.0);
			}
			
		}
	}
}
