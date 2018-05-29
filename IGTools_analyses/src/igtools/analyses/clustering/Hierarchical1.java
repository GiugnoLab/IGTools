package igtools.analyses.clustering;

import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.dictionaries.elsa.DLNELSA;

import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

public class Hierarchical1 {

	
	
	
	
	public static void run(
		B3LLSequence seq,
		int[] start,
		int[] end,
		double thr, 
		int min_k,
		int max_k,
		DkDistances dist, 
		Vector<Cluster> clusters){
		
		
		int nofSeqs = start.length;
		
		B3Sequence seqs = DLNELSA.extract(seq, start, end);
		
		int[] lengths = DLNELSA.getLengths(start, end);
		
		System.out.println("Building DLNELSA...");
		DLNELSA dlnelsa = new DLNELSA(seqs,lengths);
		
		System.out.println("Building DkVectors...");
		DkVectors dkv = new DkVectors(dlnelsa);
		if(min_k>0)
			dkv.setMinK(min_k);
		if(max_k>0)
			dkv.setMaxK(max_k);
		
		
		dlnelsa.print(dkv.getMrl()+1);
		
		
		//double[] weights = dkv.getKLogWeights();
		double[] weights = dkv.getKWeights();
		
		double[] s_vect_i = new double[weights.length];
		double[] s_vect_j = new double[weights.length];
		
		SimilarityInfDiagonalTable sims = new SimilarityInfDiagonalTable(nofSeqs, Double.NaN);
		
		double d;
		
		for(int i=0; i<nofSeqs; i++){
			dkv.getVector(i, s_vect_i);
			
			for(int j=0; j<i; j++){
				System.out.println("sim "+i+"/"+nofSeqs+" "+j+"/"+i);
				dkv.getVector(j, s_vect_j);
				d = dist.distance(s_vect_i, s_vect_j, weights);
				sims.set(i, j, d);
				System.out.println("sim "+i+"/"+nofSeqs+" "+j+"/"+i+" : "+d);
			}
		}
		
		
		
		
		int[] iclu = new int[nofSeqs];
		for(int i=0; i<iclu.length; i++){
			iclu[i] = i;
		}
		
		for(int i=nofSeqs-1; i>=0; i--){
			
			for(int j=0; j<i; j++){
				
				if(sims.get(i, j) >= thr){
					if(iclu[j] == j){
						iclu[j] = i;
					}
					else{
						
						for(int k=0; k<iclu.length; k++){
							if(iclu[k] == j)
								iclu[k] = i;
						}
						
					}
				}
			}
		}
		
		Set<Integer> tclusters = new TreeSet<Integer>();
		for(int i=0; i<iclu.length; i++){
			tclusters.add(iclu[i]);
		}
		
		for(Integer t : tclusters){
			Cluster c = new Cluster();
			for(int i=0; i<iclu.length; i++){
				if(iclu[i] == t){
					c.add(i);
				}
			}
			clusters.add(c);
		}
	}
	
}
