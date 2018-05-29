package igtools.dictionaries.elsa;

import igtools.common.nucleotide.B3Nucleotide;


/**
 * It stops at suffix-tree nodes with different elongability or with elongability equal to 0.
 * 
 * @author vbonnici
 *
 */
public class NELSAMaximalPrefixIntersection {

	
	public interface Listener{
		public void maximal(int m1, int m2);
	}
	
	public interface ExtListener{
		public void maximal(B3Nucleotide[] kmer, int m1,int m2);
	}
	
	
	private NELSA nelsa1;
	private NELSA nelsa2;
	
	public NELSAMaximalPrefixIntersection(NELSA nelsa1, NELSA nelsa2){
		this.nelsa1 = nelsa1;
		this.nelsa2 = nelsa2;
	}
	
	
	//end is not inclusive
	private int elo(NELSA nelsa, int start, int end, int[][] elos, int k){
		int nof_elos = 0;
		
		
		while(start < end &&   ( (nelsa.sa()[start]+k)  >= nelsa.b3seq().length())  ){
			start++;
		}
		
		int c = 0;
		for(c=0; c<4 && start < end; c++){
			if(nelsa.b3seq().getB3(nelsa.sa()[start] + k) == c){
				nof_elos++;
				elos[c][0] = start;
				while(start < end && nelsa.b3seq().getB3(nelsa.sa()[start] + k) == c){
					start++;
				}
				elos[c][1] = start;
			}
			else{
				elos[c][0] = -1;
				elos[c][1] = -1;
			}
		}
		for(;c<4;c++){
			elos[c][0] = -1;
			elos[c][1] = -1;
		}
		
		return nof_elos;
	}
	
	
	
	public void ext_run(ExtListener l, int max_k){
		ext_rec_run(0, nelsa1.length(), 0, nelsa2.length(), 0, l, max_k);
	}
	
	private void ext_rec_run(int s1, int e1, int s2, int e2, int k, ExtListener l, int max_k){
		if(k == max_k){
			B3Nucleotide[] kmer = new B3Nucleotide[k];
			nelsa1.b3seq().getB3(nelsa1.sa()[s1], kmer);
			l.maximal(kmer, e1 - s1, e2 - s2);
			return;
		}
//		System.out.println("@ext_rec_run(s1 "+s1+", e1 "+e1+", s2 "+s2+", e2 "+e2+", k "+k+")");
		
		int[][] elos1 = new int[4][2];
		int nof_elos1 = elo(nelsa1, s1, e1, elos1, k);
		
		int[][] elos2 = new int[4][2];
		int nof_elos2 = elo(nelsa2, s2, e2, elos2, k);
		
//		for(int i=0; i<4; i++){
//			System.out.println("\t"+elos1[i][0]+" "+elos1[i][1]+" | "+elos2[i][0]+" "+elos2[i][1]);
//		}
		
		if(nof_elos1==0 || nof_elos2==0){
			B3Nucleotide[] kmer = new B3Nucleotide[k];
			nelsa1.b3seq().getB3(nelsa1.sa()[s1], kmer);
			l.maximal(kmer, e1 - s1, e2 - s2);
		}
		else{
			
			int n_count = 0;
			int u_count = 0;
			int d_count = 0;
			for(int i=0; i<4; i++){
				if(elos1[i][0]==-1 && elos2[i][0]==-1)
					n_count++;
				else if(elos1[i][0]!=-1 && elos2[i][0]!=-1)
					u_count++;
				else
					d_count++;
				
			}
			
			if(n_count == 4 || d_count > 0){
				B3Nucleotide[] kmer = new B3Nucleotide[k];
				nelsa1.b3seq().getB3(nelsa1.sa()[s1], kmer);
				l.maximal(kmer, e1 - s1, e2 - s2);
			}
			else{
				for(int i=0; i<4; i++){
					//if(elos1[i][0]!=-1 && elos2[i][0]!=-1){ //because d_count == 0
					if(elos1[i][0]!=-1 && elos2[i][0]!=-1){
						ext_rec_run(elos1[i][0], elos1[i][1], elos2[i][0], elos2[i][1], k+1, l, max_k);
					}
				}
			}
			
		}
	}
	
}
