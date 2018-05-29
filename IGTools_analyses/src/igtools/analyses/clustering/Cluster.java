package igtools.analyses.clustering;

import java.util.LinkedList;

public class Cluster {
	public double nofElements;
	public LinkedList<Integer> elementIds;
	public double[] centroid;
	
	
	public Cluster(){
		nofElements = 0.0;
		elementIds = new LinkedList<Integer>();
	}
	
	public void add(int i){
		elementIds.add(i);
		nofElements++;
	}
	
	
	public Cluster(int sID, double[] vector){
		nofElements = 1.0;
		elementIds = new LinkedList<Integer>();
		elementIds.push(sID);
		centroid = new double[vector.length];
		for(int i=0; i<vector.length; i++)
			centroid[i] = vector[i];
	}
	
	public void add(int sID, double[] vector){
		elementIds.push(sID);
		for(int i=0; i<centroid.length; i++){
			centroid[i] = ((centroid[i] * nofElements) + vector[i]) / (nofElements+1); 
		}
		nofElements++;
	}
	
}
