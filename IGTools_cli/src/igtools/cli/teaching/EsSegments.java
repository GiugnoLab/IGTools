package igtools.cli.teaching;


import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.dictionaries.elsa.NELSA;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.CompleteIterator;



public class EsSegments{
	public static void main(String[] args){
		try{


			String sequence_file_path = "nanoa.3bit";
			B3LLSequence b3seq = B3LLSequence.load(sequence_file_path);

			String a_inelsa =  "nanoa.nelsa";
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			nelsa.setSequence(b3seq);

			
			System.out.println("sequence length: "+ b3seq.length());
			int seg_length = 100000;
			int nof_segments = (int)Math.ceil((double)b3seq.length() / (double)seg_length);
			System.out.println("number of segments: "+ nof_segments);

		
			int k = 2;
			int table[][] = new int[nof_segments][ (int)Math.pow(4,3) ];
			B3Nucleotide kmer[] = new B3Nucleotide[k];
			
			
			for(int seg_index = 0; seg_index < nof_segments; seg_index++){
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


			/*for(int seg_index = 0; seg_index < nof_segments; seg_index++){
				int start = seg_index * seg_length;
				int end = start + seg_length > b3seq.length() ? b3seq.length() : start + seg_length;
				B3Sequence segment = b3seq.subSequence(start,end);
				NELSA seg_nelsa = new NELSA(segment);


				IELSAIterator fit;
				int i=0;
				CompleteIterator it = new CompleteIterator(k);
				while(it.next()){
					it.kmer(kmer);
					fit = nelsa.find(kmer);
					if(fit != null){
						table[seg_index][i] += fit.multiplicity();
					}
					i++;
				}
			}*/


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
