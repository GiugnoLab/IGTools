package igtools.dictionaries.elsa;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.common.versioning.DataVersioning;


/**
 * Multi-sequence NELSA.
 * (Sequences are merged together to form a single indexed sequence, by concatenating them with a N at the end of each sequence) 
 * 
 * (Not implemented) The array cumLengths stores the ID of the original sequence at position i of the merged sequence.
 * (implemented) The array ssIDs (ss_ids) stores the ID of the original sequence at position SA[i].
 * 
 * @author vbonnici
 *
 */
public class DLNELSA extends NELSA{
	public static final long serialVersionUID = DataVersioning.DLNELSA |	0x0000000000000001L;
	
	
	/**
	 * Extract sub-sequences in [start,end[ and concat them in a single sequence by putting N in the middle.
	 * 
	 * @param iseq
	 * @param start
	 * @param end
	 * @return
	 */
	public static B3LLSequence extract(B3Sequence iseq, int[] start, int[] end){
		int t_length = 0;
		for(int i=0; i<start.length; i++){
			t_length += end[i] - start[i] + 1;
		}
		B3LLSequence oseq = new B3LLSequence(t_length);
		
		//System.out.println("##"+start.length+" "+end.length);
		
		int p = 0;
		for(int i=0; i<start.length; i++){
			for(int j=start[i]; j<end[i]; j++){
				oseq.setB3(p, iseq.getB3(j));
				p++;
			}
			oseq.setB3(p, B3Nucleotide.N_CODE);
			p++;
		}
		return oseq;
	}
	
	public static int[] getLengths(int[] start, int[] end){
		int[] lengths = new int[start.length];
		for(int i=0; i<start.length; i++){
			lengths[i] += end[i] - start[i] + 1;
		}
		return lengths;
	}
	
	
	
	
	
	/**
	 * number of sequences indexed by the DLNELSA
	 */
	protected int nof_ss;
	/**
	 * length of each indexed sequence
	 */
	protected int[] ss_lengths;
	/**
	 * ss_ids[i] = id of the original sequence 
	 */
	protected int[] ss_ids; 
	
	
	public DLNELSA(){
		super();
	}
	
	public DLNELSA(B3Sequence b3seq, int[] ss_lengths){
		this.nof_ss = ss_lengths.length;
		this.ss_lengths = ss_lengths;
		
		//int t_length = 0;
		
//		ss_lengths = new int[nof_ss];
//		for(int i=0; i<ss_lengths.length; i++){
//			t_length += ss_lengths[i];
//		}
		
//		int p = 0;
//		for(int i=0; i<sequences.length; i++){
//			for(int j=0; j<sequences[i].length; j++){
//				b3seq.setB3(p, sequences[i][j].code());
//				p++;
//			}
//			b3seq.setB3(p, B3Nucleotide.N_CODE);
//			p++;
//		}
		
		super.init(b3seq);
		
		
		computeSSIDS();
	}
	
	
	public int nofSequences(){
		return this.nof_ss;
	}
	public int[] sequenceLengths(){
		return this.ss_lengths;
	}
	public int[] ssIDs(){
		return this.ss_ids;
	}
	
	
	protected void computeSSIDS(){
		int[] sa = this.sa;
		this.ss_ids = new int[sa.length];
		
//		int[] cum_lengths = new int[nof_ss];
//		
//		cum_lengths[0] = ss_lengths[0] ;
//		for(int i=1; i<ss_lengths.length; i++){
//			cum_lengths[i] = cum_lengths[i-1] + ss_lengths[i];
//		}
//		
//		int sid;
//		for(int i=0; i<sa.length; i++){
//			sid = 0;
//			while(cum_lengths[sid] < sa[i])
//				sid++;
//			ss_ids[i] = sid-1;
//		}
		int[] sid = new int[sa.length];
		int ii = 0;
		for(int i=0; i<ss_lengths.length; i++){
			for(int j=0; j<ss_lengths[i]; j++){
				sid[ii] = i;
				ii++;
			}
		}
		
		for(int i=0; i<sa.length; i++){
			ss_ids[i] = sid[sa[i]];
		}
	}
	
	
	@Override
	public void print(int k){
	  System.out.println("i\tSA\tLCP\tFN\tS\tsuff");
	  for(int i=0; i< sa.length; i++){
		  System.out.print(i +"\t"+ sa[i] +"\t"+ lcp[i] +"\t"+ ns[i] +"\t"+ss_ids[i]+"\t");
		  
		  int nn = s.length() - sa[i] < k ? s.length() - sa[i] : k;
		  B3Nucleotide[] suff = new B3Nucleotide[nn];
		  s.getB3(sa[i], suff);
		  for(int j=0; j<suff.length; j++)
			  System.out.print(suff[j]);
		  
		  System.out.println();
	  }
    }
}
