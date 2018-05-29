package igtools.analyses.clustering;

import java.util.Arrays;

import igtools.analyses.Inform;
import igtools.dictionaries.elsa.DLNELSA;
import igtools.dictionaries.elsa.IELSAIterator;

public class DkVectors {

	
	private DLNELSA dlnelsa;
	private int mrl;
	
	private int min_k = 0;
	private int max_k = 0;
	
	private int v_length;
	
	private double[] weigths = null;
	
	public DkVectors(DLNELSA dlnelsa){
		this.dlnelsa = dlnelsa;
		
		this.mrl = Inform.Mrl(dlnelsa);
		this.min_k = 1;
		//this.min_k = getCurveFlex();
		//this.min_k = getCompleteDLK();
		this.max_k = mrl;
		
		computeVectorLength();
		
		System.out.println("mrl "+mrl);
		System.out.println("vLength "+v_length);
	}
	
	public int getCurveFlex(){
		int p_nof = 0;
		int k_nof;
		for(int k=1; k<mrl; k++){
			k_nof = dlnelsa.nof(k);
			if(k_nof < p_nof){
				return k_nof-1;
			}
		}
		return mrl;
	}
	
	public int getCompleteDLK(){
		int[] sids = dlnelsa.ssIDs();
		
		int[] count = new int[dlnelsa.nofSequences()];  
//		int nof_k;
		
		
//		int p_id;
		int is, ie;
		
		int k = 1;
		for(; k<= mrl; k++){
			Arrays.fill(count, 0);
//			nof_k=0;
			IELSAIterator it = dlnelsa.begin(k);
			while(it.next()){
				is = it.istart();
				ie = it.iend();
				
//				p_id = -1;
				for(int i=is; i<ie; i++){
					
					System.out.println(i+" "+sids[i]);
					
					count[sids[i]]++;
				}
			}
			
			for(int i=0; i<count.length; i++){
				if(count[i] == 0)
					return k;
			}
		}
		return 1;
	}
	
	public int getMrl(){
		return this.mrl;
	}
	public int minK(){
		return this.min_k;
	}
	public int maxK(){
		return this.max_k;
	}
	public void setMinK(int min_k){
		this.min_k = min_k;
		computeVectorLength();
	}
	public void setMaxK(int max_k){
		this.max_k = max_k;
		computeVectorLength();
	}
	
	
	
//	public void maxKByAll(){
//		boolean[] pre = new boolean[dlnelsa.nofSequences()];
//		Arrays.fill(pre, false);
//		
//		for(int k=0; k<=mrl; k++){
//			Arrays.fill(pre, false);
//			IELSAIterator it = dlnelsa.begin(k);
//			while(it.next()){
//				
//			}
//		}
//	}
	
	public int vectorLength(){
		return v_length;
	}
	
	private void computeVectorLength(){
		v_length = 0;
		for(int k=min_k; k<=max_k; k++){
			System.out.println("k "+k+" ;nof "+dlnelsa.nof(k));
			v_length += dlnelsa.nof(k);
		}
	}
	
	
	
	public double[] getKLogWeights(){
		if(this.weigths == null){
			this.weigths = new double[v_length];
			int v=0;
			
			for(int k=min_k; k<=max_k; k++){
				int nof =  dlnelsa.nof(k);
				for(int j=0; j<nof; j++){
					this.weigths[v] = Math.log(k);
					v++;
				}
			}	
		}
		return this.weigths;
	}
	public double[] getKWeights(){
		if(this.weigths == null){
			this.weigths = new double[v_length];
			int v=0;
			
			for(int k=min_k; k<=max_k; k++){
				int nof =  dlnelsa.nof(k);
				for(int j=0; j<nof; j++){
					this.weigths[v] = k;
				}
			}	
		}
		return this.weigths;
	}
	
	
	public double[] getVector(int sID){
		double[] vect = new double[v_length];
		getVector(sID, vect);
		return vect;
	}
	public void getVector(int sID, double[] vect){
		//double[] vect = new double[v_length];
		int v = 0;
		
		int[] ids = this.dlnelsa.ssIDs();
		
		double count;
		
		for(int k=min_k; k<=max_k; k++){
			IELSAIterator it = dlnelsa.begin(k);
			while(it.next()){
				count = 0;
				for(int i=it.istart(); i<it.iend(); i++){
					if(ids[i] == sID){
						count++;
					}
				}
				vect[v] = count;
				v++;
			}
		}
	}
	
	public double nofFeatures(double[] v1, double[]v2){
		double ret = 0.0;
		for(int i=0; i<v1.length; i++){
			if(v1[i] > v2[i]) 
				ret += v1[i];
			else
				ret += v2[i];
		}
		return ret;
	}
}
