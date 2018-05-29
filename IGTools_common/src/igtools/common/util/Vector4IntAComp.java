package igtools.common.util;

import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

public class Vector4IntAComp implements Comparator<Vector4Int>{
	
	
	private Map<Vector4Int, Long> ta = new TreeMap<Vector4Int, Long>();
	
	public Vector4IntAComp(){
		
	}
	
	private Long nofAnagrams(Vector4Int v){
		Long a = ta.get(v);
		if(a == null){
			a = Maths.nofAnagrams(v.values, v.values[0] + v.values[1] + v.values[2] + v.values[3]);
			ta.put(v.clone(), a);
		}
		return a;
	}

	@Override
	public int compare(Vector4Int o1, Vector4Int o2) {
		Long res = (nofAnagrams(o1) - nofAnagrams(o2));
		if(res == 0){
			return o1.compareTo(o2);
		}
		return res.intValue();
	}

}
