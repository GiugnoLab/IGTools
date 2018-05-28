import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.CompleteIterator;



public class EsSegments{
	public static void main(String[] args){
		try{

			//loading a pre-built 3-bit sequence from a file
			String sequence_file_path = "nanoa.3bit";
			B3LLSequence b3seq = B3LLSequence.load(sequence_file_path);

			
			//loading a pre-built NELSA from a file
			String a_inelsa =  "nanoa.nelsa";
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			//remember to reassing the 3-bitsequence to the NELSA
			nelsa.setSequence(b3seq);

			System.out.println("sequence length: "+ b3seq.length());
			int seg_length = 100000;
			int nof_segments = (int)Math.ceil((double)b3seq.length() / (double)seg_length);
			System.out.println("number of segments: "+ nof_segments);

			//the goal of this example is to calculate k-mers multiplicities in segments of the original strings
			//thus, a table with #segments rows and 4^k columns is retrieved
		
			int k = 2;
			int table[][] = new int[nof_segments][ (int)Math.pow(4,k) ];
			B3Nucleotide kmer[] = new B3Nucleotide[k];
			
			
			for(int seg_index = 0; seg_index < nof_segments; seg_index++){
				//for each segments
				//we exctract it from the original sequence
				//no copies are made because the B3LLSequence.subSequence method creates an object that refers to the original sequence data
				//then a segment-specific NELSA is created and an iterator over it is performed
				//The lexicographic order is used to specific the position/index of the kmer within the table
				
				int start = seg_index * seg_length;
				int end = start + seg_length > b3seq.length() ? b3seq.length() : start + seg_length;

				B3LLSequence segment = new B3LLSequence(b3seq.subSequence(start,end));

				NELSA seg_nelsa = new NELSA(segment);
				IELSAIterator it = seg_nelsa.begin(k);
				while(it.next()){
					it.kmer(kmer);
					table[seg_index] [(int) B3Nucleotide.toLexicoOrder(kmer)] = it.multiplicity();
				}

			}

			System.out.print("#\t");
			//The CompleteIterator is a way to list all the theoretical k-mers in T^k
			CompleteIterator cit = new CompleteIterator(k);
			while(cit.next()){
				cit.kmer(kmer);
				System.out.print(B3Nucleotide.toString(kmer) +"\t");
			}
			System.out.println("SUM");


			for(int i=0; i<nof_segments; i++){
				System.out.print(i+"\t");
				int sum = 0;
				for(int j=0; j< (int)Math.pow(4,k); j++){
					System.out.print(table[i][j]+"\t");
					sum += table[i][j];
				}
				System.out.println(sum);
			}


		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
		}
	}
}



/*public static long toLexicoOrder(B3Nucleotide[] ns){
	long ret = 0;
	
	long p4 = 1;
	for(int i=0; i<ns.length - 1; i++){
		p4 *= 4;
	}
	
	for(int i=0; i<ns.length; i++){
		ret += ns[i].code *  p4;
		p4 = p4 / 4;
	}
	
	return ret;
}*/
