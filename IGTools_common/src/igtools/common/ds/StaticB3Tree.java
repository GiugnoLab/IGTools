package igtools.common.ds;

import java.util.Map;


import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;

public class StaticB3Tree {

	public class StaticB3TreeNode{
		protected StaticB3TreeNode[] childs = new StaticB3TreeNode[5];
		public int occurrences = 1;
		
		protected StaticB3TreeNode(){
		}
		public StaticB3TreeNode firstChild(){
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
	
	private StaticB3TreeNode root = new StaticB3TreeNode();
	
	public StaticB3Tree(){
	}
	
	
	public StaticB3TreeNode root(){
		return this.root;
	}
	
	
	public StaticB3TreeNode add(int n_code, StaticB3TreeNode tn){
//		if(tn==null)
//			tn = root;
		if(tn.childs[n_code] == null){
			tn.childs[n_code] = new StaticB3TreeNode();
		}
		else{
			tn.childs[n_code].occurrences++;
		}
		return tn.childs[n_code];
	}
	public StaticB3TreeNode add(B3Nucleotide n, StaticB3TreeNode tn){
		return add(n.code(), tn);
	}
	public StaticB3TreeNode add(B3Nucleotide n){
		return add(n.code(), root);
	}
	public StaticB3TreeNode add(B3LLSequence seq, int start, int end){
//		System.out.print("add["+start+","+end+"]:"+B3Nucleotide.by(seq.getB3(start)));
		StaticB3TreeNode tn = add(B3Nucleotide.by(seq.getB3(start)));
		for(int i=start+1; i<end; i++){
//			System.out.print(B3Nucleotide.by(seq.getB3(i)));
			tn = add(seq.getB3(i), tn);
		}
//		System.out.println();
		return tn;
	}
	
	
	public StaticB3TreeNode get(StaticB3TreeNode tn, int n_code){
		return tn.childs[n_code];
	}
	public StaticB3TreeNode get(StaticB3TreeNode tn, B3Nucleotide n){
		return get(tn, n.code());
	}
	
	public StaticB3TreeNode get(B3LLSequence seq, int start, int end){
//		System.out.println("S:"+start+":"+seq.getB3(start));
		
		StaticB3TreeNode tn =  root.childs[seq.getB3(start)];
		
		int i = start +1;
		while(tn!=null && i<end){
			tn = tn.childs[seq.getB3(i)];
			i++;
		}
		
		return tn;
	}
	
	public void print(int[] stats){
		print(root, "",stats);
	}
	public void print(StaticB3TreeNode tn, String s, int[] stats){
		for(int i=0; i<5; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					stats[0]++;
					System.out.println(s+B3Nucleotide.by(i));
				}
				else{
					stats[1]++;
					print(tn.childs[i], s+B3Nucleotide.by(i), stats);
				}
			}
		}
	}
	
	public void statsP(long[] stats, Map<Integer, Integer> comul){
		statsP(root, 1,stats, comul, "");
	}
	public void statsP(StaticB3TreeNode tn, int depth, long[] stats, Map<Integer, Integer> comul, String s){
		for(int i=0; i<5; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					stats[0]++;
					stats[2]+=depth;
					
					System.out.println(s+B3Nucleotide.by(i) +":"+depth);
					
					if(comul.containsKey(depth)){
						comul.put(depth, comul.get(depth)+1);
					}
					else{
						comul.put(depth, 1);
					}
				}
				else{
					stats[1]++;
					statsP(tn.childs[i], depth+1, stats, comul, s+B3Nucleotide.by(i));
				}
			}
		}
	}
	
	public void stats(long[] stats, Map<Integer, Integer> comul){
		stats(root, 1,stats, comul);
	}
	public void stats(StaticB3TreeNode tn, int depth, long[] stats, Map<Integer, Integer> comul){
		for(int i=0; i<5; i++){
			if(tn.childs[i] != null){
				if(tn.childs[i].firstChild()==null){
					stats[0]++;
					stats[2]+=depth;
					
					if(comul.containsKey(depth)){
						comul.put(depth, comul.get(depth)+1);
					}
					else{
						comul.put(depth, 1);
					}
				}
				else{
					stats[1]++;
					stats(tn.childs[i], depth+1, stats, comul);
				}
			}
		}
	}
}
