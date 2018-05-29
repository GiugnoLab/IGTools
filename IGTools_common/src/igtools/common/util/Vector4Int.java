package igtools.common.util;

public class Vector4Int implements Comparable<Vector4Int>{

	public int[] values = new int[4]; 
	
	public Vector4Int(){
		
	}
	public Vector4Int(int[] values){
		setValues(values);
	}
	
	public void setValues(int[] values){
		for(int i=0; i<4; i++)
			this.values[i] = values[i];
		
	}
	
	public Vector4Int clone(){
		return new Vector4Int(this.values);
	}
	
	@Override
	public int compareTo(Vector4Int o) {
		int res = 0;
		for(int i=0; i<4; i++){
			res = this.values[i] - o.values[i];
			if(res != 0)
				break;
		}
		return res;
	}
	
	@Override
	public boolean equals(Object o){
		if(o instanceof Vector4Int){
			return this.compareTo((Vector4Int)o) == 0;
		}
		return false;
	}
	
	@Override
	public String toString(){
		return "("+values[0]+","+values[1]+","+values[2]+","+values[3]+")";
	}
	
}
