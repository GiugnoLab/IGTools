package igtools.common.ds;

import java.util.TreeMap;

import igtools.common.nucleotide.B2Nucleotide;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;

public class StaticB2Tree {
	public class StaticB2TreeNode{
		protected StaticB2TreeNode[] childs = new StaticB2TreeNode[4];
		public int occurrences = 0;
		
		protected StaticB2TreeNode(){
		}
		public StaticB2TreeNode firstChild(){
			for(int i=0; i<childs.length; i++){
				if(childs[i] != null)
					return childs[i];
			}
			return null;
		}
		public int firstChildCode(){
			for(int i=0; i<childs.length; i++){
				if(childs[i] != null)
					return i;
			}
			return -1;
		}
	}
	
	
	
	
	
	
	private StaticB2TreeNode root = new StaticB2TreeNode();
	
	public StaticB2Tree(){
	}
	
	
	public StaticB2TreeNode root(){
		return this.root;
	}
	
	public StaticB2TreeNode add(int n_code, StaticB2TreeNode tn){
		if(tn.childs[n_code] == null){
			tn.childs[n_code] = new StaticB2TreeNode();
		}
		//else{
		//	tn.childs[n_code].occurrences++;
		//}
		return tn.childs[n_code];
	}
	public StaticB2TreeNode add(B2Nucleotide n, StaticB2TreeNode tn){
		return add(n.code(), tn);
	}
	public StaticB2TreeNode add(B2Nucleotide n){
		return add(
				n.code(),
				root);
	}
	public StaticB2TreeNode add(B3Sequence seq, int start, int end){
		StaticB2TreeNode tn = add(B2Nucleotide.by(seq.getB2(start)));
		for(int i=start+1; i<end; i++){
			tn = add(seq.getB2(i), tn);
		}
		tn.occurrences++;
		return tn;
	}
	public StaticB2TreeNode addReverse(B3Sequence seq, int start, int end){
		StaticB2TreeNode tn = add(B2Nucleotide.by(seq.getB2(end-1)));
		for(int i=end-2; i>=start; i--){
			tn = add(seq.getB2(i), tn);
		}
		tn.occurrences++;
		return tn;
	}
	public StaticB2TreeNode addAll(B3Sequence seq, int start, int end){
		StaticB2TreeNode tn = add(B2Nucleotide.by(seq.getB2(start)));
		tn.occurrences++;
		for(int i=start+1; i<end; i++){
			tn = add(seq.getB2(i), tn);
			tn.occurrences++;
		}
		return tn;
	}
	
	
	
	public StaticB2TreeNode get(StaticB2TreeNode tn, int n_code){
		return tn.childs[n_code];
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
	
	
	
	public interface ExtensionListener{
		public void extension(String extension);
	}
	public void extend(StaticB2TreeNode tn, String extension, ExtensionListener listener){
		if(tn != null){
			if(tn.occurrences > 0)
				listener.extension(extension);
			for(int i=0; i<4; i++){
				extend(tn.childs[i], extension + B3Nucleotide.by(i), listener);
			}
		}
	}
	
	
	
	
	public void printWords(String prefix){
		printWords(root, prefix);
	}
	public void printWords(StaticB2TreeNode tn, String s){
		for(int i=0; i<4; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					System.out.println(s+B2Nucleotide.by(i));
				}
				else{
					printWords(tn.childs[i], s+B3Nucleotide.by(i));
				}
			}
		}
	}
	
	
	
	public void stats(TreeMap<Integer,Integer> comults){
		stats(root, 0, comults);
	}
	public void stats(StaticB2TreeNode tn, int level, TreeMap<Integer,Integer> comults){
		for(int i=0; i<4; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					Integer value = comults.get(level+1);
					if(value == null)
						value= 0;
					comults.put(level+1,value+1);
				}
				else{
					stats(tn.childs[i], level+1, comults);
				}
			}
		}
	}
	
	
	
	public interface WordListener{
		public void word(String s);
	}
	
	
	public void listWords(WordListener listener){
		listWords(root, "", listener);
	}
	public void listWords(StaticB2TreeNode tn, String s, WordListener listener){
		for(int i=0; i<4; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					//System.out.println(s+B2Nucleotide.by(i));
					listener.word(s+B2Nucleotide.by(i));
				}
				else{
					listWords(tn.childs[i], s+B3Nucleotide.by(i), listener);
				}
			}
		}
	}
	
	
	public void listOriginalWords(WordListener listener){
		listOriginalWords(root, "", listener);
	}
	public void listOriginalWords(StaticB2TreeNode tn, String s, WordListener listener){
		if(tn != null){
			if(tn.occurrences > 0)
				listener.word(s);
			
			for(int i=0; i<4; i++){
				listOriginalWords(tn.childs[i], s+B3Nucleotide.by(i), listener);
			}
		}
	}
	
	
	public void listSharedPrefixes(WordListener listener){
		listSharedPrefixes(root, "", listener);
	}
	public void listSharedPrefixes(StaticB2TreeNode tn, String s, WordListener listener){
		if(tn != null){
			boolean noGoodChilds = true;
			
			
			for(int i=0; i<4; i++){
				if(tn.childs[i] != null && tn.childs[i].occurrences > 1)
					noGoodChilds = false;
				//listOriginalWords(tn.childs[i], s+B3Nucleotide.by(i), listener);
			}
			
			if(noGoodChilds)
				listener.word(s);
			else{
				for(int i=0; i<4; i++){
					if(tn.childs[i] != null && tn.childs[i].occurrences > 1)
						listSharedPrefixes(tn.childs[i], s+B3Nucleotide.by(i), listener);
				}
			}
		}
	}
	
	
	public void listMaximalSharedPrefixes(WordListener listener){
		listSharedPrefixes(root, "", listener);
	}
	public void listMaximalSharedPrefixes(StaticB2TreeNode tn, String s, WordListener listener){
		if(tn != null){
			int noGoodChilds = 0;
			
			
			for(int i=0; i<4; i++){
				if(tn.childs[i] != null && tn.childs[i].occurrences > 1)
					noGoodChilds++;
			}
			
			if(noGoodChilds == 0)
				listener.word(s);
			else{
				for(int i=0; i<4; i++){
					if(tn.childs[i] != null && tn.childs[i].occurrences > 1)
						listMaximalSharedPrefixes(tn.childs[i], s+B3Nucleotide.by(i), listener);
				}
			}
		}
	}
	
	
	public boolean containsOriginal(B3Sequence seq){
		StaticB2TreeNode n = root;
		int code;
		for(int i=0; i<seq.length(); i++){
			code = (B2Nucleotide.by(seq.getB2(i))).code();
			if(n.childs[code] == null)
			    return false;
			else
				n = n.childs[code];
		}
		return n.occurrences > 0;
	}
	
	public int nof(B3Sequence seq){
		StaticB2TreeNode n = root;
		int code;
		for(int i=0; i<seq.length(); i++){
			code = (B2Nucleotide.by(seq.getB2(i))).code();
			if(n.childs[code] == null)
			    return 0;
			else
				n = n.childs[code];
		}
		return n.occurrences;
	}
	
	public boolean contains(B3Sequence seq){
		StaticB2TreeNode n = root;
		int code;
		for(int i=0; i<seq.length(); i++){
			code = (B2Nucleotide.by(seq.getB2(i))).code();
			if(n.childs[code] == null)
			    return false;
			else
				n = n.childs[code];
		}
		return true;
	}
	
	
}
