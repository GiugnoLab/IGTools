package igtools.cli.teaching;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.NELSAMaximalPrefixIntersection;
import igtools.dictionaries.elsa.NELSAMaximalSubPrefixIntersection;

public class EsJP {

	
	public static void main(String[] args){
		
		try{
			
			
			String  g1 = "ACGTTACGATTGGCCAGTTCAACGTTCCAGTTTGGGAAAACACGTGCAAGTTCAG";
			String  g2 = "ACGTTACGATTGGCCAGAACAACGTTCCAGTTTCCCAAAACACGTGCTTGTTCAG";
			//String  g2 = "ACGT";
			//String  g2 = "AACGTAGAT";
			
			B3LLSequence s1 = new B3LLSequence(g1);
			B3LLSequence s2 = new B3LLSequence(g2);
			
			NELSA n1 = new NELSA(s1);
			NELSA n2 = new NELSA(s2);
			
			System.out.println("================================================================================");
			int k = 20;
			B3Nucleotide[] kmer = new B3Nucleotide[k];
			
			System.out.println(s1);
			for(int i=0; i<n1.length(); i++){
				//s1.getB3(n1.sa()[i], kmer);
				//System.out.println(i+"\t"+n1.sa()[i]+"\t"+ B3Nucleotide.toString(kmer));
				//System.out.println(i+"\t"+n1.sa()[i]+"\t"+ B3Nucleotide.toString(kmer, 0,    k < s1.length() - n1.sa()[i] ?  k  : s1.length() - n1.sa()[i] ));
				System.out.println(i+"\t"+n1.sa()[i]+"\t"+ s1.subSequence(n1.sa()[i], s1.length));
			}
			
			
			System.out.println("================================================================================");
			
			System.out.println(s2);
			for(int i=0; i<n2.length(); i++){
				System.out.println(i+"\t"+n2.sa()[i]+"\t"+ s2.subSequence(n2.sa()[i], s2.length));
			}
			
			
			System.out.println("================================================================================");
			//stops at suffix-tree nodes with different elongability or with elongability euqal to 0
			NELSAMaximalPrefixIntersection mi = new NELSAMaximalPrefixIntersection(n1, n2);
			NELSAMaximalPrefixIntersection.ExtListener li = new NELSAMaximalPrefixIntersection.ExtListener() {
				@Override
				public void maximal(B3Nucleotide[] kmer, int m1, int m2) {
					System.out.println(B3Nucleotide.toString(kmer)+" "+m1+" "+m2);
				}
			};
			mi.ext_run(li,  Integer.MAX_VALUE);
			
			System.out.println("================================================================================");
			//visit all common tree paths 
			NELSAMaximalSubPrefixIntersection msi = new NELSAMaximalSubPrefixIntersection(n1, n2);
			NELSAMaximalSubPrefixIntersection.ExtListener lsi = new NELSAMaximalSubPrefixIntersection.ExtListener() {
				@Override
				public void maximal(B3Nucleotide[] kmer, int m1, int m2) {
					System.out.println(B3Nucleotide.toString(kmer)+" "+m1+" "+m2);
				}
			};
			msi.ext_run(lsi,  Integer.MAX_VALUE);
			
			
			
		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
		}
		
	}
}
