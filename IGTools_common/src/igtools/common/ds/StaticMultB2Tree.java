package igtools.common.ds;

import igtools.common.nucleotide.B2Nucleotide;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;

public class StaticMultB2Tree {
	public class StaticB2TreeNode{
		protected StaticB2TreeNode[] childs = null;
		//= new StaticB2TreeNode[4];
		public long multiplicity = 0;
		
		protected StaticB2TreeNode(){
		}
		public StaticB2TreeNode firstChild(){
			if(this.childs != null){
				for(int i=0; i<childs.length; i++){
					if(childs[i] != null)
						return childs[i];
				}
			}
			return null;
		}
		public int firstChildCode(){
			if(this.childs != null){
				for(int i=0; i<childs.length; i++){
					if(childs[i] != null)
						return i;
				}
			}
			return -1;
		}
	}
	
	
	
	
	
	
	private StaticB2TreeNode root = new StaticB2TreeNode();
	
	public StaticMultB2Tree(){
	}
	
	
	public StaticB2TreeNode root(){
		return this.root;
	}
	
	public StaticB2TreeNode add(int n_code, StaticB2TreeNode tn, int multiplicity){
		if(tn.childs == null){
			tn.childs = new StaticB2TreeNode[4];
		}
		if(tn.childs[n_code] == null){
			tn.childs[n_code] = new StaticB2TreeNode();
		}
		tn.childs[n_code].multiplicity += multiplicity;
		return tn.childs[n_code];
	}
	public StaticB2TreeNode add(B2Nucleotide n, StaticB2TreeNode tn, int multiplicity){
		return add(n.code(), tn, multiplicity);
	}
	public StaticB2TreeNode add(B2Nucleotide n, int multiplicity){
		return add(n.code(), root, multiplicity);
	}
	public StaticB2TreeNode add(B3Sequence seq, int start, int end, int multiplicity){
		StaticB2TreeNode tn = add(B2Nucleotide.by(seq.getB2(start)), multiplicity);
		for(int i=start+1; i<end; i++){
			tn = add(seq.getB2(i), tn, multiplicity);
		}
		return tn;
	}
	public StaticB2TreeNode add(B3Nucleotide[] seq, int multiplicity){
		StaticB2TreeNode tn = add(seq[0].code(), root, multiplicity);
		for(int i=1; i< seq.length; i++){
			tn = add(seq[i].code(), tn, multiplicity);
		}
		return tn;
	}
	
	
	
	public StaticB2TreeNode get(StaticB2TreeNode tn, int n_code){
		if(tn.childs != null)
			return tn.childs[n_code];
		return null;
	}
	public StaticB2TreeNode get(StaticB2TreeNode tn, B2Nucleotide n){
		return get(tn, n.code());
	}
	public StaticB2TreeNode get(B3LLSequence seq){
		StaticB2TreeNode n = get(root, seq.getB2(0));
		if(n != null){
			for(int i=1; i<seq.length(); i++){
				n = get(n, seq.getB2(i));
				if(n == null)
					return null;
			}
		}
		return n;
	}
	
	
}
