/**************************************************************************
 * @author Hussein Houdrouge
 * Master m1 
 * Ecole Ploytechnique de Paris
 * 
 * This class implements the heat method.
 * 
 *************************************************************************/

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleLUDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class HeatMethod {
	
	public Polyhedron_3<Point_3> polyhedron3D;
	public double [] colors; 
	public double MaxColor;
	public MeshViewer view;
	int s = 0;
	double t;
	private DoubleMatrix2D A;
	private SparseDoubleMatrix2D Lc;
	private SparseDoubleMatrix2D M;
	private DenseDoubleLUDecomposition lu;
	private DoubleMatrix1D L;
	private DoubleMatrix1D phi;
	private DenseDoubleLUDecomposition lulc; 
	public HeatMethod(Polyhedron_3<Point_3> polyhedron3D, int s) {
		long pb = System.currentTimeMillis();
		this.s = s;
		this.polyhedron3D = polyhedron3D;
		t = 100 * computeTime();
		colors = new double[this.polyhedron3D.facets.size()];
		A  = computeA(); //compute the vertices area matrix
		Lc = computeLc();
		M  = computeM();
		long pa = System.currentTimeMillis();
		System.out.println("elapsed time for pre computation phase one = " + (pa - pb)/1000.0);
		
		long b = System.currentTimeMillis();
		lu = new DenseDoubleLUDecomposition(M);
		long a = System.currentTimeMillis();
		System.out.println("elapsed time to compute LU of m = " + (a - b)/1000.0);
		
		long lcb = System.currentTimeMillis();
		lulc = new DenseDoubleLUDecomposition(Lc);
		long lca = System.currentTimeMillis();
		System.out.println("elapsed time to compute LU of Lc = " + (lca - lcb)/1000.0 );
	}
	
	/**
	 * 
	 * t = h^2   
	 * h is the mean spacing between adjacent nodes 
	 * @return t 
	 */
	private double computeTime() {
		HashSet<Halfedge<Point_3>>    edges = new HashSet<Halfedge<Point_3>>();
		ArrayList<Halfedge<Point_3>> hedges = this.polyhedron3D.halfedges;
		double sum = 0.0;
		int n = 0;
		for(Halfedge<Point_3> h : hedges) { 
			if(edges.contains(h)) continue;
			if(edges.contains(h.opposite)) continue;
			edges.add(h);
			double d = Utilities.dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
			sum += d;
			n++;
		}
		return (sum /n ) * (sum /n);
	}
	
	public double getTime() {
		return t;
	}
	
	/**
	 * compute the area of each face in the mesh
	 * and store it in a hash map.
	 * @return hash map from face to its area.
	 */
	private HashMap<Face<Point_3>, Double> computeFacesArea() {
		HashMap<Face<Point_3>, Double> mp = new HashMap<Face<Point_3>, Double>(); // map from Face -> area
		ArrayList<Face<Point_3>> faces = this.polyhedron3D.facets;
		for(Face<Point_3> f: faces) {
			mp.put(f, Utilities.triangleArea(f)); // compute the area of the face which is a triangle 
		}
		return mp;
	}
	
	/**
	 * 
	 * vertex area for a vertex v is one third the area of all the triangles incident on vertex v.
	 * given a map that store the area of each face.
	 */
	private double computeVertexArea(Vertex<Point_3> v, HashMap<Face<Point_3>, Double> mp) {
		double va = 0.0;
		Halfedge<Point_3> h = v.getHalfedge();
		Halfedge<Point_3> current = h.next.opposite;
		while(current != h) {
			if(mp.get(current.getFace()) != null)
			va += mp.get(current.getFace());
			current = current.next.opposite;
		}
		va += mp.get(current.getFace());
		return va / 3.0;
	}
		
	/**
	 * if i = j then
	 *  	Aij is one third of the area of all triangles incident on vertex i
	 * Aij = 0
	 * A is |V| X |V|
	 * 
	 */
	public DoubleMatrix2D computeA() {
		ArrayList<Vertex<Point_3>> vertices = this.polyhedron3D.vertices;
		DoubleMatrix2D M = new SparseDoubleMatrix2D(vertices.size(), vertices.size()); //create M |V| X |V| matrix
		HashMap<Face<Point_3>, Double> mp = computeFacesArea();                       // mp is map to store the faces area.
		int i = 0;
		for(Vertex<Point_3> v: vertices) {
			M.set(i, i, computeVertexArea(v, mp));
			i++;
		}
		return M;
	}
	
	/**
	 * Given a vertex v return all of its neighbor in a hashSet
	 * @param v
	 * @return HashSet
	 */
	private HashSet<Vertex<Point_3>> getNeighbor(Vertex<Point_3> v) {
		HashSet<Vertex<Point_3>> s = new HashSet<Vertex<Point_3>>();
		Halfedge<Point_3> h = v.getHalfedge().opposite;
		Halfedge<Point_3> current = h.opposite.next;
		while(current != h) {
			s.add(current.vertex);
			current = current.opposite.next;
		}
		s.add(current.vertex);
		return s;
	}
	
	/**
	 * compute the cotan formula without computing the angles
	 * 1/(8 * A1) (d0^2 - d1^2 - d2^2) + 1/(8 * A2)(d0^2 - d3^2 - d4^2) 
	 * @param u
	 * @param v
	 * @return
	 */
	private double cotan1(Vertex<Point_3> u, Vertex<Point_3> v) {
		//direct h to from u to v 
		Halfedge<Point_3>    h = u.getHalfedge().opposite;
		while(h.vertex != v) h = h.opposite.next;
		Halfedge<Point_3> o    = h.opposite;		
		double d0 = Utilities.dist(h.vertex.getPoint(), o.vertex.getPoint());
		
		h = h.next;
		o = o.next;
		double d2 = Utilities.dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
		double d3 = Utilities.dist(o.vertex.getPoint(), o.opposite.vertex.getPoint());
		
		h = h.next;
		o = o.next;
		double d1 = Utilities.dist(h.vertex.getPoint(), h.opposite.vertex.getPoint());
		double d4 = Utilities.dist(o.vertex.getPoint(), o.opposite.vertex.getPoint());
		
		double A1 = Utilities.heron(d0, d1, d2); 
		double A2 = Utilities.heron(d0, d3, d4); 
		A1 = 1/(8 * A1);  
		A2 = 1/(8 * A2);
		double d02 = d0 * d0;
		return A1 * (d02 - d1 * d1 - d2 * d2) + A2 * (d02 - d3 * d3 - d4 * d4);
	}
	
	/**
	 * given two vertices vi vj compute the contan formula
	 * 1/2 (cot(ai) + cot(bi). 
	 * @param u
	 * @param v
	 * @return
	 */
	private double computeCotan(Vertex<Point_3> u, Vertex<Point_3> v) {
		return cotan1(u, v);
	}
	
	/**
	 * Lc (i, j) = - 1/2 cot(aij) + 1/2 cot(bij) if i and j are adjacent
 	 * Lc (i, i) sum over j for L(i, j).
 	 * Lc is |V| X |V| matrix
	 */
	public SparseDoubleMatrix2D computeLc() {
		ArrayList<Vertex<Point_3>> vertices = this.polyhedron3D.vertices;
		SparseDoubleMatrix2D M = new SparseDoubleMatrix2D(vertices.size(), vertices.size()); //create M |V|x|V| matrix
		for(int i = 0; i < vertices.size(); ++i) {
			HashSet<Vertex<Point_3>> s = getNeighbor(vertices.get(i)); // get the neighbor of the vertex vi
			double sum = 0.0;
			for(Vertex<Point_3> v : s) { 
				double cotan = computeCotan(vertices.get(i), v); // compute cotan formula
				M.set(i, v.index, -1 * cotan);
				sum += cotan;
			}
			M.set(i, i, sum);
		}
		return M;
	}
	
	/**
	 * Compute the matrix M = A - t*Lc
	 * @return
	 */
	public SparseDoubleMatrix2D computeM() {
		SparseDoubleMatrix2D M = new SparseDoubleMatrix2D(Lc.rows(), Lc.columns());
		for(int i = 0; i < A.rows(); ++i) {
			HashSet<Vertex<Point_3>> s = getNeighbor(this.polyhedron3D.vertices.get(i)); // get the neighbor of the vertex vi
			double tmp = A.get(i, i) -  t * Lc.get(i, i);
			M.set(i, i, tmp);
			for(Vertex<Point_3> v: s) {
				tmp = - t * Lc.get(i, v.index);
				if (Math.abs(tmp) != 0)
					M.set(i, v.index, tmp);
			}
		}
		return M;
	}
	
	/**
	 * 
	 * 
	 * @param s the index of the heat source
	 */
	public void heatDiffusion(int s) {
		this.s = s;
		DoubleMatrix1D K = new SparseDoubleMatrix1D(A.rows());
		K.set(s, 1);
		L = lu.solve(K);
	}
	
	/**
	 * compute the color value for a face the average value of its vertices 
	 * @param f
	 * @param lv map that store the value of each face
	 * @return 
	 */
	private double computeColorValue(Face<Point_3> f, HashMap<Vertex<Point_3>, Double> lv) {
		Halfedge<Point_3> h = f.getEdge();
		double sum = lv.get(h.getVertex());
		h = h.next;
		sum += lv.get(h.getVertex());
		h = h.next;
		sum += lv.get(h.getVertex());
		sum /= 3;
		return sum;
	}

	/**
	 * given a choice c that can be either h for heat or d for distance define a heat map
	 * for each face.
	 * @param c
	 */
	public void colorScheme(char c) {
		DoubleMatrix1D cs;
		if(c == 'h') cs = L;
		else if(c == 'd') cs = phi;
		else return;
			
	
		HashMap<Vertex<Point_3>, Double> lv = new HashMap<Vertex<Point_3>,Double>(); //hash map from each vertex to its value 
		ArrayList<Vertex<Point_3>> vertices = this.polyhedron3D.vertices;
		for(int i = 0; i < vertices.size(); ++i) lv.put(vertices.get(i), cs.get(i));
		
		ArrayList<Face<Point_3>> faces = this.polyhedron3D.facets;
		
		//compute the color of each face and set the max color value
		for(int i = 0; i < faces.size(); ++i) {
//			colors[i] = Color_Map.computeColorValue_Baricentric(faces.get(i), lv); 
			colors[i] = computeColorValue(faces.get(i), lv);
			MaxColor = Math.max(MaxColor, colors[i]);
		}
	}
	
	/**
	 * compute the normal vector of a given face
	 * @param f
	 * @return
	 */
	private Vector_3 normalVector(Face<Point_3> f) {
		Halfedge<Point_3> e1 = f.getEdge();
		Halfedge<Point_3> e2 = e1.next.next.opposite;
		//compute the vector corresponded to e1.
		double v1[] = new double[3];
		v1[0] = e1.vertex.getPoint().getX().doubleValue() - e1.opposite.vertex.getPoint().getX().doubleValue();
		v1[1] = e1.vertex.getPoint().getY().doubleValue() - e1.opposite.vertex.getPoint().getY().doubleValue();
		v1[2] = e1.vertex.getPoint().getZ().doubleValue() - e1.opposite.vertex.getPoint().getZ().doubleValue();
		
		//compute the vector corresponded to e2.
		double v2[] = new double[3];
		v2[0] = e2.vertex.getPoint().getX().doubleValue() - e2.opposite.vertex.getPoint().getX().doubleValue();
		v2[1] = e2.vertex.getPoint().getY().doubleValue() - e2.opposite.vertex.getPoint().getY().doubleValue();
		v2[2] = e2.vertex.getPoint().getZ().doubleValue() - e2.opposite.vertex.getPoint().getZ().doubleValue();
		
		Vector_3 v_1 = new Vector_3(v1[0], v1[1], v1[2]);
		Vector_3 v_2 = new Vector_3(v2[0], v2[1], v2[2]);
		
		Vector_3 n = v_1.crossProduct(v_2);
		return Utilities.normalize(n);
	}
	
	/**
	 * compute the gradient of a face as Delta u = 1/(2 * Af) sum_i (N ^ ei)
	 * where ei is the edge oppose the vertex i and N is the unit normal;
	 * @param f
	 * @param vu
	 * @return
	 */
	private Vector_3 gradiantFace(Face<Point_3> f, HashMap<Vertex<Point_3>, Double> vu) {
		double Af = Utilities.triangleArea(f); // compute the area of a triangle
		Vector_3 N = normalVector(f);		   // compute the normal of f.
		
		
		Halfedge<Point_3> e1 = f.getEdge();
		Halfedge<Point_3> e2 = e1.next;
		Halfedge<Point_3> e3 = e2.next;		
		
		Vector_3 v1 = Utilities.getVector(e1);
		Vector_3 v2 = Utilities.getVector(e2);
		Vector_3 v3 = Utilities.getVector(e3);
		
	//	System.out.println("is it zero = " + ((Math.abs(in - 0) < 0.000000000001)? 0 : -1));
		
		Vertex<Point_3> u1 = e2.vertex;
		Vertex<Point_3> u2 = e3.vertex;
		Vertex<Point_3> u3 = e1.vertex;
		
		double u = vu.get(u1);
		Vector_3 c1 = N.crossProduct(v1);
		c1 = c1.multiplyByScalar(u);
		
		u = vu.get(u2);
		Vector_3 c2 = N.crossProduct(v2);
		c2 = c2.multiplyByScalar(u);
		
		u = vu.get(u3);
		Vector_3 c3 = N.crossProduct(v3);
		c3 = c3.multiplyByScalar(u);
		
		double delta_x = c1.getX().doubleValue() + c2.getX().doubleValue() + c3.getX().doubleValue();
		double delta_y = c1.getY().doubleValue() + c2.getY().doubleValue() + c3.getY().doubleValue();
		double delta_z = c1.getZ().doubleValue() + c2.getZ().doubleValue() + c3.getZ().doubleValue(); 
		
		double area = 1/(2 * Af);
		delta_x = delta_x * area;
		delta_y = delta_y * area;
		delta_z = delta_z * area;
	//	System.out.println("c = quiver3(" + N.getX().doubleValue() + "," + N.getY().doubleValue() + "," + N.getZ().doubleValue() +")");
	//	System.out.println("hold on");
		return new Vector_3(delta_x, delta_y, delta_z);
	}
	
	
	/**
	 * given a half edge of a face compute the delta X on that face
	 * @param h
	 * @param Xj
	 * @return
	 */
	private double computeDeltaXuFace(Halfedge<Point_3> h, Vector_3 Xj) {
		Halfedge<Point_3> e1 = h;
		Halfedge<Point_3> e2 = e1.next;
		Halfedge<Point_3> e3 = e2.next;
		
		Vector_3 v1  = Utilities.getVector(e1);
		Vector_3 v2  = Utilities.getVector(e2.opposite);
		Vector_3 v3  = Utilities.getVector(e3);
		
		double theta2 = Utilities.computeAngele(v3, v2);	
		v3 = Utilities.getVector(e3.opposite);
		double theta1 = Utilities.computeAngele(v1, v3);	

		Vector_3 E2 = Utilities.getVector(e1.opposite);
		Vector_3 E1 = Utilities.getVector(e2);
		double inner1 = E1.innerProduct(Xj).doubleValue();
		double inner2 = E2.innerProduct(Xj).doubleValue();
		
		double cot1 = 1/Math.tan(theta1);
		double cot2 = 1/Math.tan(theta2);
		
		return (cot1 * inner1 + cot2 * inner2);
	}
	
	/**
	 * compute the delta xu for each vertex
	 * @param u
	 * @param fv is the gradient in every face 
	 * @return
	 */
	private double computeDeltaXu(Vertex<Point_3> u, HashMap<Face<Point_3>, Vector_3> fv) {
		Halfedge<Point_3> h = u.getHalfedge();
		double delta = computeDeltaXuFace(h, fv.get(h.getFace())); 
		Halfedge<Point_3> current = h.next.opposite;
		while(current != h) {
			delta += computeDeltaXuFace(current, fv.get(current.getFace()));
			current = current.next.opposite;
		}
		
		return 0.5 * delta;
	}
	
	/**
	 * Compute the integrated divergence
	 * @param s is the source index 
	 * @return
	 */
	private DoubleMatrix1D computeDeltaX(int s) {
		HashMap<Vertex<Point_3>, Double> vu = new HashMap<Vertex<Point_3>, Double>();
		ArrayList<Vertex<Point_3>> vertices = this.polyhedron3D.vertices;
		heatDiffusion(s); // compute the heat diffusion 
		// store in a map the heat associated with each vector
		for(int i = 0; i < vertices.size(); ++i) vu.put(vertices.get(i), L.get(i));
		
		//compute the gradient vector for each face and store it in a map.
		ArrayList<Face<Point_3>> faces = this.polyhedron3D.facets;
		HashMap<Face<Point_3>, Vector_3> fv = new HashMap<Face<Point_3>, Vector_3>();
		for(int i = 0; i < faces.size(); ++i) {
			Vector_3 x = gradiantFace(faces.get(i), vu);
			//normalize x
			x = Utilities.normalize(x);
			fv.put(faces.get(i), x);
		}
		//compute delta X
		DoubleMatrix1D X = new DenseDoubleMatrix1D(vertices.size());
		for(int i = 0; i < vertices.size(); ++i) {
			X.set(i, -1 * computeDeltaXu(vertices.get(i), fv));
		}
		return X;
	}
	
	private DoubleMatrix1D computeB(int s) {
		return computeDeltaX(s);
	}
	
	public void computePhi(int s) {
		this.s = s;
		long before = System.currentTimeMillis();
		DoubleMatrix1D b  = computeB(s);
		long after = System.currentTimeMillis();
		System.out.println("the elapsed time to compute X = " + (after - before)/1000.0);
	
		before = System.currentTimeMillis();
		phi = lulc.solve(b);
		after = System.currentTimeMillis();
		System.out.println("the elapsed time to solve for phi = " + (after - before)/1000.0);
		
		double mini = 0;
		for(int i = 0; i < phi.size(); ++i) {
			mini = Math.min(phi.get(i), mini);
		}
		if(mini < 0)  {
			for(int i = 0; i < phi.size(); i++) {
				phi.set(i, phi.get(i) + Math.abs(mini));
			}
		} 
		System.out.println("\n phi  = \n " +  phi);
	}
}
