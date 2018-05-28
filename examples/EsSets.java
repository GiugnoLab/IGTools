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
		
		//finding the union and the interseciotn of D_5(s1)  and  D_5(s2)
		IELSAIterator it1 = n1.begin(k);
		IELSAIterator it2;
		while(it1.next()){
			//list kmers in D_5(s1)
			it2 = n2.find(it1.kmer());
			if(it2 != null){
				//if they exist in s2, then they are in the intersection only
				System.out.println("I\t"+B3Nucleotide.toString(it1.kmer()));
				System.out.println("U\t"+B3Nucleotide.toString(it1.kmer()));
			}
			else{
				//otherwise they are also included in the union of the two dictionaries
				System.out.println("U\t"+B3Nucleotide.toString(it1.kmer()));
			}
		}
		//to find the complete union we need to search for kmers in s2 that are not in s1
		it2 = n2.begin(k);
		while(it2.next()){
			it1 = n1.find(it2.kmer());
			if(it1 == null){
				System.out.println("U\t"+B3Nucleotide.toString(it2.kmer()));
			}
		}
		
		
		
		System.out.println("================================================================================");
		//IGTools API provide a more efficient way to calculate intersection and union
		//The methods are also able to perform set operations without loading the NELSAs in main memory
		B3LLSequence[] seqs = {s1,s2};
		IELSAIterator[] its = new IELSAIterator[2];
		its[0] = n1.begin(k);
		its[1] = n2.begin(k);
		//the OnDiskNELSAIteratorV2 iterates over kmers without loading the NELSA in main memory, but jut the 3-bit sequence
		//its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		
		//set-theoretic operations take the advantage of suffix order over the suffix arrays
		//moreover, they compare kmers by means of their position in the lexicographic order, rather than using string comparison
		//in this way, kmers are coverted into integers, and integers are compared several times
		//time is saved w.r.t. string comparison, just because objects are compared several times
		if (k < 32){
			//if kmers are smaller than 32 nucleotides we can use the implicit integer representation given by the 3-bit encoding
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
			//otherwise we need to calculate the lexicographic order
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
		//the previous methods give both intersection and union
		//however, if onyl union is needed, this twos perform better
		its[0] = n1.begin(k);
		its[1] = n2.begin(k);
		//its[i] = new OnDiskNELSAIteratorV2(seqs[i], inelsas.get(i), k);
		
		NumericUnion.UnionListerner listener = new NumericUnion.UnionListerner() {
			@Override
			public void intersection(B3Nucleotide[] kmer) {
				//included in the interface just for compatibility with the other operators 
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

