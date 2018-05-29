package igtools.cli.entropy;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.NELSAMaximalPrefixIntersection;

public class JKL {

	public static void usage(){
		System.out.println("Usage: cmd seq1 nelsa1 seq2 nelsa2");
	}
	
	
	public static void main(String[] args){
		String a_seq1 = "";
		String a_nelsa1 = "";
		String a_seq2 = "";
		String a_nelsa2 = "";
		
		try{
			a_seq1 = args[0];
			a_nelsa1 = args[1];
			a_seq2 = args[2];
			a_nelsa2 = args[3];
		}catch(Exception e){
			usage();
			System.exit(0);
		}
		
		try{
			
			B3LLSequence seq1 = B3LLSequence.load(a_seq1);
			B3LLSequence seq2 = B3LLSequence.load(a_seq2);
			NELSA nelsa1 = new NELSA();
			nelsa1.load(a_nelsa1);
			nelsa1.setSequence(seq1);
			NELSA nelsa2 = new NELSA();
			nelsa2.load(a_nelsa2);
			nelsa2.setSequence(seq2);
			
			
			int m1 =  (int)Math.ceil(Math.log(seq1.length()) / Math.log(2));
			int m2 =  (int)Math.ceil(Math.log(seq2.length()) / Math.log(2));
			int m = 0;
			if(m1 <= m2) m = m1;
			else m = m2;
			
			if(m1 <= m2) m = m2;
			else m = m1;
			
			//do not limit extension depth
			//m = Integer.MAX_VALUE;
			
			System.out.println("ceil(log2("+seq1.length()+")) = "+m1+", ceil(log2("+seq2.length()+")) = "+m2+"");
			System.out.println("m "+ m);
			
			
			final double[] nofs = {0.0};
			final double[] sums = {0.0,0.0};
			
			NELSAMaximalPrefixIntersection mi = new NELSAMaximalPrefixIntersection(nelsa1, nelsa2);
			NELSAMaximalPrefixIntersection.ExtListener li = new NELSAMaximalPrefixIntersection.ExtListener() {
				@Override
				public void maximal(B3Nucleotide[] kmer, int m1, int m2) {
					System.out.println(B3Nucleotide.toString(kmer)+" "+m1+" "+m2);
					nofs[0]++;
					sums[0] += m1;
					sums[1] += m2;
				}
			};
			mi.ext_run(li, m);
			System.out.println(nofs[0]+" "+sums[0]+" "+sums[1]);
			
			
			final double[] jkl = {0.0};
			final double[] p = {0.0,0.0};
			final double log2 = Math.log(2.0);
			li = new NELSAMaximalPrefixIntersection.ExtListener() {
				@Override
				public void maximal(B3Nucleotide[] kmer, int m1, int m2) {
					p[0] = ((double)m1) / sums[0];
					p[1] = ((double)m2) / sums[1];
					jkl[0] += (p[0] - p[1]) * Math.log(p[0] / p[1]) * log2; 
				}
			};
			mi.ext_run(li, m);
			
			System.out.println("jkl "+jkl[0]);
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
}
