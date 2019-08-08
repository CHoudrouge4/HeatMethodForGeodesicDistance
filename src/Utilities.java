/**************************************************************************
 * @author Hussein Houdrouge
 * Master m1 
 * Ecole Ploytechnique de Paris
 * 
 *************************************************************************/

import java.util.ArrayList;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;

public class Utilities {
	
	public static double dist(Point_3 p, Point_3 q) {
		double d = 0.0;
		
		double px = p.getX().doubleValue();
		double py = p.getY().doubleValue();
		double pz = p.getZ().doubleValue();
		
		double qx = q.getX().doubleValue();
		double qy = q.getY().doubleValue();
		double qz = q.getZ().doubleValue();
		
		d = (px - qx) * (px - qx) + (py - qy) * (py - qy) + (pz -qz) * (pz - qz);
		return Math.sqrt(d);
	}
	
	public static double mean(ArrayList<Double> array) {
		int n = array.size();
		double sum = 0.0;
		for(int i = 0; i < n; ++i) 
			sum += array.get(i);
		return sum / n;
	}
	
	public static double meanSqaure(ArrayList<Double> array) { 
		double m = mean(array);
		return m * m;
	}
	
	/**
	 * compute the area of a triangle given its three 
	 * sides a, b, and c.
	 * s = (a + b + c)/2;
	 * A = sqrt(s * (s - a) * (s - b) * (s - c));
	 * @param a
	 * @param b
	 * @param c
	 * @return A
	 */
	public static double heron(double a, double b, double c) {
		double s  = (a + b + c)/2.0;
		double a2 = s * (s - a) * (s - b) * (s - c);
		return Math.sqrt(a2);
	}
	/**
	 * Given a face as a triangle compute
	 * its area.
	 */
	public static double triangleArea(Face<Point_3> f) {
		Halfedge<Point_3> h = f.getEdge();
		double s1 = dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
		h = h.next;
		double s2 = dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
		h = h.next;
		double s3 = dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
		return heron(s1, s2, s3);
	}
	/**
	 * compute the vector norm 
	 * as |v| = sqrt(<v, v>)
	 * @param v
	 * @return
	 */
	public static double vectorNorm(Vector_3 v) {
		double norm = v.innerProduct(v).doubleValue();
		return Math.sqrt(norm);
	}

	/**
	 * normalize a vector v  by devising each component by vector norm. 
	 * @param v
	 * @return
	 */
	public static Vector_3 normalize(Vector_3 v) {
		double norm = vectorNorm(v);
		double nx = v.getX().doubleValue()/norm;
		double ny = v.getY().doubleValue()/norm;
		double nz = v.getZ().doubleValue()/norm;
		v.setX(nx);
		v.setY(ny);
		v.setZ(nz);
		return v;
	}
	
	/**
	 * 
	 * @param e
	 * @return
	 */
	public static Vector_3 getVector(Halfedge<Point_3> e) {
		double v[] = new double[3];
		v[0] = e.vertex.getPoint().getX().doubleValue() - e.opposite.vertex.getPoint().getX().doubleValue();
		v[1] = e.vertex.getPoint().getY().doubleValue() - e.opposite.vertex.getPoint().getY().doubleValue();
		v[2] = e.vertex.getPoint().getZ().doubleValue() - e.opposite.vertex.getPoint().getZ().doubleValue();
		return new Vector_3(v[0], v[1], v[2]);
	}
	
	public static double vectorLength(Vector_3 v) {
		double l = v.innerProduct(v).doubleValue();
		return Math.sqrt(l);
	}
	
	/**
	 * compute the angle between to vector 
	 * @param u
	 * @param v
	 * @return
	 */
	public static double computeAngele(Vector_3 u, Vector_3 v) {
		double inerProduct = u.innerProduct(v).doubleValue();
		double lu = vectorLength(u);
		double lv = vectorLength(v); 
		double cos = inerProduct / (lu * lv);
		return Math.acos(cos);
	}
}
