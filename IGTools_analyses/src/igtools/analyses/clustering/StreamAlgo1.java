package igtools.analyses.clustering;

import java.util.Vector;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.dictionaries.elsa.DLNELSA;

public class StreamAlgo1 {

	
	public static void run(
			B3LLSequence seq,
			int[] start,
			int[] end,
			double thr, 
			int min_k,
			int max_k,
			DkDistances dist, 
			Vector<Cluster> clusters){
		
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
		
		
		//dkv.setMaxK(dkv.getMrl());
		//dkv.setMinK(dkv.getMrl()-3);
		
		//clusters = new Vector<Cluster>();
		
		int min_c;
		double min_d, d;
//		, d_thr;
		double[] weights = dkv.getKLogWeights();
		double[] s_vect = new double[weights.length];
		
		System.out.println("Running...");
		for(int s=0; s<lengths.length; s++){
			
			B3Nucleotide[] kmer = new B3Nucleotide[end[s]-start[s]+1];
			seq.getB3(start[s], kmer);
			System.out.print(B3Nucleotide.toString(kmer)+" ");
			
			//double[] s_vect = dkv.getVector(s);
			dkv.getVector(s, s_vect);
			
			if(clusters.size() == 0){
				System.out.println(s+" new cluster 0");
				
				clusters.add(new Cluster(s, s_vect));
			}
			else{
				min_c = Integer.MAX_VALUE;
				min_d = Double.MAX_VALUE;
				
				for(int c =0; c<clusters.size(); c++){
					d = dist.distance(s_vect, clusters.get(c).centroid, weights);
					//d_thr = dkv.nofFeatures(s_vect, clusters.get(c).centroid) * thr;
					
					//System.out.println(c+" "+d+" "+thr);
					
					if(d >= thr && d < min_d){
						min_d = 1;
						min_c = c;
					}
				}
				
				if(min_d != Double.MAX_VALUE){
					System.out.println(s+" to cluster "+min_c);
					clusters.get(min_c).add(s, s_vect);
				}
				else{
					System.out.println(s+" new cluster "+clusters.size());
					clusters.add(new Cluster(s, s_vect));
				}
			}
			
		}
		System.out.println("done.");
		
	}
}
