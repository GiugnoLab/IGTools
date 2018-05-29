package igtools.cli.teaching;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.intersection.BigNumericUnion;
import igtools.dictionaries.intersection.BigNumericUnionIntersection;
import igtools.dictionaries.intersection.NumericUnion;
import igtools.dictionaries.intersection.NumericUnionIntersection;

public class EsSets {

	
	public static void main(String[] args){
		String  g1 = "ACGTTACGATTGGCCAGTTCAACGTTCCAGTTTGGGAAAACACGTGCAAGTTCAG";
		String  g2 = "ACGTTACGATTGGCCAGAACAACGTTCCAGTTTCCCAAAACACGTGCTTGTTCAG";
		//String  g2 = "ACGT";
		//String  g2 = "AACGTAGAT";
		
		B3LLSequence s1 = new B3LLSequence(g1);
		B3LLSequence s2 = new B3LLSequence(g2);
		
		NELSA n1 = new NELSA(s1);
		NELSA n2 = new NELSA(s2);
		
		System.out.println("================================================================================");
		
		System.out.println(s1);
		for(int i=0; i<n1.length(); i++){
			System.out.println(i+"\t"+n1.sa()[i]+"\t"+ s1.subSequence(n1.sa()[i], s1.length));
		}
		
		System.out.println("================================================================================");
		
		System.out.println(s2);
		for(int i=0; i<n2.length(); i++){
			System.out.println(i+"\t"+n2.sa()[i]+"\t"+ s2.subSequence(n2.sa()[i], s2.length));
		}
		
		
		System.out.println("================================================================================");
		
		
		int k = 5;
		System.out.println("================================================================================");
		
		
		IELSAIterator it1 = n1.begin(k);
		IELSAIterator it2;
		while(it1.next()){
			it2 = n2.find(it1.kmer());
			if(it2 != null){
				System.out.println("I\t"+B3Nucleotide.toString(it1.kmer()));
			}
			else{
				System.out.println("U\t"+B3Nucleotide.toString(it1.kmer()));
			}
		}
		it2 = n2.begin(k);
		while(it2.next()){
			it1 = n1.find(it2.kmer());
			if(it1 == null){
				System.out.println("U\t"+B3Nucleotide.toString(it2.kmer()));
			}
		}
		
		
		
		System.out.println("================================================================================");
		B3LLSequence[] seqs = {s1,s2};
		IELSAIterator[] its = new IELSAIterator[2];
		its[0] = n1.begin(k);
		its[1] = n2.begin(k);
		//its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		
		
		if (k < 32){
			NumericUnionIntersection.UnionListerner listener = new NumericUnionIntersection.UnionListerner() {
				@Override
				public void intersection(B3Nucleotide[] kmer, int mult, int nof_seqs) {
					System.out.println("I\t"+B3Nucleotide.toString(kmer) +"\t"+ mult+"\t"+nof_seqs);
				}
				@Override
				public void union(B3Nucleotide[] kmer, int mult, int nof_seqs) {
					System.out.println("U\t"+B3Nucleotide.toString(kmer) +"\t"+ mult+"\t"+nof_seqs);
				}
			};
			NumericUnionIntersection union =  new NumericUnionIntersection( seqs , its, k, listener);
			union.union();
		}
		else{
			BigNumericUnionIntersection.UnionListerner listener = new BigNumericUnionIntersection.UnionListerner() {
				@Override
				public void intersection(B3Nucleotide[] kmer, int mult, int nof_seqs) {
					System.out.println("I\t"+B3Nucleotide.toString(kmer) +"\t"+ mult+"\t"+nof_seqs);
				}
				@Override
				public void union(B3Nucleotide[] kmer, int mult, int nof_seqs) {
					System.out.println("U\t"+B3Nucleotide.toString(kmer) +"\t"+ mult+"\t"+nof_seqs);
				}
			};
			BigNumericUnionIntersection union =  new BigNumericUnionIntersection( seqs , its, k, listener);
			union.union();
		}
		
		
		System.out.println("================================================================================");
		
		its[0] = n1.begin(k);
		its[1] = n2.begin(k);
		//its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		
		NumericUnion.UnionListerner listener = new NumericUnion.UnionListerner() {
			@Override
			public void intersection(B3Nucleotide[] kmer) {
			}
			@Override
			public void union(B3Nucleotide[] kmer, int mult) {
				System.out.println("U\t"+B3Nucleotide.toString(kmer));
			}
		};
		NumericUnion union = null;
		if(k < 32)
			union = new NumericUnion( seqs , its, k, listener);
		else 
			union = new BigNumericUnion(seqs, its, k, listener);
		union.union();
		
	}
}
