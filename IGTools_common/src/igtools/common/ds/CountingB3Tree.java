package igtools.common.ds;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;

public class CountingB3Tree {

	public class Node{
		public B3Nucleotide c;
		public int count = 0;
		public boolean word = false;
		
		Node[] childs = new Node[5];
		
		public Node addChild(B3Nucleotide n){
			return addChild(n.code());
		}
		public Node addChild(int code){
			if(childs[code] == null){
				childs[code] = new Node();
				childs[code].c = B3Nucleotide.by(code);
			}
			return childs[code];
		}
	}
	
	
	
	private Node root = new Node();
	
	public CountingB3Tree(){
		
	}
	
	public void add(B3Nucleotide[] ns){
		Node n = root;
		for(int i=0; i<ns.length; i++){
			n = n.addChild(ns[i]);
			n.count++;
		}
		n.word = true;
	}
	
	public void add(B3Sequence seq){
		Node n = root;
		for(int i=0; i<seq.length(); i++){
			n = n.addChild(seq.getB3(i));
			n.count++;
		}
		n.word = true;
	}
	
	public void addRC(B3Sequence seq){
		Node n = root;
		for(int i=seq.length()-1; i>=0; i--){
			n = n.addChild(B3Nucleotide.rcs_code[seq.getB3(i)]);
			n.count++;
		}
		n.word = true;
	}
	
	
	public boolean contains(B3Nucleotide[] ns){
		Node n = root;
		for(int i=0; i<ns.length; i++){
			n = n.childs[ns[i].code()];
			if(n == null) return false;
		}
		return n.word;
	}
	
	public boolean containsRC(B3Nucleotide[] ns){
		Node n = root;
		for(int i=ns.length-1; i>=0; i--){
			n = n.childs[B3Nucleotide.rcs_code[ns[i].code()]];
			if(n == null) return false;
		}
		return n.word;
	}
	
	
	public boolean contains(B3Sequence ns){
		Node n = root;
		for(int i=0; i<ns.length(); i++){
			n = n.childs[ns.getB3(i)];
			if(n == null) return false;
		}
		return n.word;
	}
	
	public boolean containsRC(B3Sequence ns){
		Node n = root;
		for(int i=ns.length()-1; i>=0; i--){
			n = n.childs[B3Nucleotide.rcs_code[ns.getB3(i)]];
			if(n == null) return false;
		}
		return n.word;
	}
	
	public long nofWords(){
		return nofWords(root);
	}
	private long nofWords(Node n){
		long ret = n.word ? 1 : 0;
		for(int i=0; i<n.childs.length; i++){
			if(n.childs[i] != null)
				ret += nofWords(n.childs[i]); 
		}
		return ret;
	}
	
	
	public static interface WordListener{
		public void word(String w, int w_count, int c_ount);
	}
	
	public void listWords(WordListener l){
		for(int i=0; i<4; i++){
			listWords(root.childs[i], "", l);
		}
	}
	public void listWords(Node n, String w, WordListener l){
		if(n != null){
			if(n.word){
				int w_count = n.count;
				int c_count = 0;
				for(int i=0; i<n.childs.length; i++){
					if(n.childs[i]!=null){
						w_count -= n.childs[i].count;
						c_count += n.childs[i].count;
					}
				}
				l.word(w + n.c, w_count, c_count);
			}
			for(int i=0; i<4; i++){
				listWords(n.childs[i], w + n.c, l);
			}
		}
	}
	
	public void listMaximalWords(WordListener l){
		for(int i=0; i<4; i++){
			listMaximalWords(root.childs[i], "", l);
		}
	}
	public void listMaximalWords(Node n, String w, WordListener l){
		if(n != null){
			boolean maximal = true;
			for(int i=0; i<4; i++){
				listWords(n.childs[i], w + n.c, l);
				if(n.childs[i] != null)
					maximal = false;
			}
			if(maximal){
				l.word(w + n.c, n.count, 0);
			}
		}
	}
}
