import java.util.HashMap;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;

public class Color_Map {
	
	public static float[] linear_interpolation(float x1, float y1, float x2, float y2) {
		float [] ab = new float[2];
		ab [0] = (y2 - y1) / (x2 - x1);
		ab [1] = y1 - (ab[0]*x1);
		return ab;
	}
	
	public static float fl(float a, float b, float x) {
		return a * x + b;
	}
	
	// colors from low to high
	public static float[][] create_color_map(float [][] colors, int map_size) { 
		int number_of_colors = colors.length; 
		float [][] map = new float[map_size][3];
		map[0] = colors[0];
		map[map_size - 1] = colors[number_of_colors - 1];
		int step = map_size/(number_of_colors - 1);
		int appand = map_size % step;
		for(int i = map_size - 1; i >= map_size - 1 - appand; i--) {
			System.out.println("Yes I am running!!");
			map[i] = colors[number_of_colors - 1];
		}
		
		int j = 1;
		for(int i = step; i <= map_size - step; i += step) {
			System.out.println(i);
			map[i - 1] = colors[j];
			j++;			
		} 
		
		for(int k = 0; k < number_of_colors - 1; ++k) { 
			int i = Math.max(k * step - 1, 0);
			int ii = (k + 1) * step - 1;
			float[] ab0 = linear_interpolation(i, map[i][0], ii, map[ii][0]);
			float[] ab1 = linear_interpolation(i, map[i][1], ii, map[ii][1]);
			float[] ab2 = linear_interpolation(i, map[i][2], ii, map[ii][2]);
			for(j = i + 1; j < i + step; ++j) {
				map[j][0] = fl(ab0[0], ab0[1], j);
				map[j][1] = fl(ab1[0], ab1[1], j);
				map[j][2] = fl(ab2[0], ab2[1], j);
			}
		}
		return map;
	}
	
	/** 
	 *  compute the barycentric coordinates given the centroid p and the a, b, c the vertices of a triangle
	 */
	public static float[] compute_barycentric_coordinate(Point_3 p, Point_3 a, Point_3 b, Point_3 c) {
		float b_coord [] = new float[3];
		Vector_3 v0 = (Vector_3) b.minus(a);
	    System.out.println(v0);               	
				
	    Vector_3 v1 = (Vector_3) c.minus(a);
	    Vector_3 v2 = (Vector_3) p.minus(a);
		float d00 =  v0.innerProduct(v0).floatValue();
		float d01 =  v0.innerProduct(v1).floatValue();
		float d11 =  v1.innerProduct(v1).floatValue();
		float d20 =  v2.innerProduct(v0).floatValue();
		float d21 =  v2.innerProduct(v1).floatValue();
		float denom = d00 * d11 - d01 * d01;
		b_coord[1] = (d11 * d20 - d01 * d21) / denom;
		b_coord[2] = (d00 * d21 - d01 * d20) / denom;
		b_coord[0] = 1.0f - b_coord[1] - b_coord[2];
	
		return b_coord;
	}
	
	public static float[] color_scheme(float [][] color_map, double color, double maxColor) {
		int i = (int)Math.floor((color / maxColor) * (color_map.length - 1));
		System.out.println(i);
		return color_map[i];
	}
	
	public static double computeColorValue_Average(Face<Point_3> f, HashMap<Vertex<Point_3>, Double> lv) {
		Halfedge<Point_3> h = f.getEdge();
		double sum = lv.get(h.getVertex());
		h = h.next;
		sum += lv.get(h.getVertex());
		h = h.next;
		sum += lv.get(h.getVertex());
		sum /= 3;
		return sum;
	}
	
	public static double computeColorValue_Baricentric(Face<Point_3> f, HashMap<Vertex<Point_3>, Double> lv) {
		Halfedge<Point_3> h = f.getEdge();
		Vertex<Point_3> p1 = h.getVertex();
		h = h.next;
		Vertex<Point_3> p2 = h.getVertex();
		h = h.next;
		Vertex<Point_3> p3 = h.getVertex();
		
		float x = (p1.getPoint().getX().floatValue() + p2.getPoint().getX().floatValue() + p3.getPoint().getX().floatValue())/3; 
		float y = (p1.getPoint().getY().floatValue() + p2.getPoint().getY().floatValue() + p3.getPoint().getY().floatValue())/3;
		float z = (p1.getPoint().getZ().floatValue() + p2.getPoint().getZ().floatValue() + p3.getPoint().getZ().floatValue())/3;
		
		Point_3 p = new Point_3(x, y, z);
		float [] barycenters = compute_barycentric_coordinate(p, p1.getPoint(), p2.getPoint(), p3.getPoint());
		
		return barycenters[0] * lv.get(p1) + barycenters[1] * lv.get(p2) + barycenters[2] * lv.get(p3);
	}
	
	public static float[][] create_color_map_level_set(float [][] colors, int map_size, int period) {
		float color_map[][] = create_color_map(colors, map_size);
		float color [] = {100, 100, 100};
		for(int i = 0; i < color_map.length; i += period)
			for(int j = i; j < i + 3 && j < color_map.length; ++j) {
				color_map[j] = color;
			}
		
		return color_map;
	}
}
