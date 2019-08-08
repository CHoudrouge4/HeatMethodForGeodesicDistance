import Jcg.geometry.*;
import Jcg.mesh.SharedVertexRepresentation;
import Jcg.polyhedron.*;

/**
 * Class for rendering a surface triangle mesh (using Processing)
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 * edited By Hussein Houdrouge for INF555 project.
 *
 */
public class SurfaceMesh {
	
	double scaleFactor = 60; // scaling factor: useful for 3d rendering
	MeshViewer view; // Processing 3d frame (where meshes are rendered)
	public Polyhedron_3<Point_3> polyhedron3D; // triangle mesh
	float[][] color_map; //
	float colors [][] = {{0, 0, 255}, {0, 255, 0}, {255, 140, 0}, {255, 0, 0}, {0, 0, 0}};  //

	/**
	 *
	 * Create a surface mesh from an OFF file
	 */	
	public SurfaceMesh(MeshViewer view, String filename) {
		this.view = view;
		
		//color_map = Color_Map.create_color_map_level_set(this.colors, 100, 10);

//		shared vertex representation of the mesh
    	SharedVertexRepresentation sharedVertex=new SharedVertexRepresentation(filename);
    	LoadMesh<Point_3> load3D=new LoadMesh<Point_3>();
    	
    	polyhedron3D=load3D.createTriangleMesh(sharedVertex.points,sharedVertex.faceDegrees,
				sharedVertex.faces,sharedVertex.sizeHalfedges);

    	//System.out.println(polyhedron3D.verticesToString());   	
    	//System.out.println(polyhedron3D.facesToString());
    	polyhedron3D.isValid(false);
    	    	
    	this.scaleFactor=this.computeScaleFactor();
	}
	
	/**
	 * Draw a segment between two points
	 */	
	public void drawSegment(Point_3 p, Point_3 q) {
		float s=(float)this.scaleFactor;
		float x1=(float)p.getX().doubleValue()*s;
		float y1=(float)p.getY().doubleValue()*s;
		float z1=(float)p.getZ().doubleValue()*s;
		float x2=(float)q.getX().doubleValue()*s;
		float y2=(float)q.getY().doubleValue()*s;
		float z2=(float)q.getZ().doubleValue()*s;
		this.view.line(	x1, y1, z1, x2, y2, z2 );		
	}
	
	/**
	 * Draw a vertex (as a small sphere)
	 */	
	public void drawVertex(Point_3 p, Point_3 src) {
		if(p.equals(src)) view.fill(255, 0, 0);
		else view.fill(0, 0, 255);
		float s=(float)this.scaleFactor;
		float x1=(float)p.getX().doubleValue()*s;
		float y1=(float)p.getY().doubleValue()*s;
		float z1=(float)p.getZ().doubleValue()*s;
		
		view.translate(x1, y1, z1);
		view.sphere(s/25f);
		view.translate(-x1, -y1, -z1);
	}

	/**
	 * Draw a triangle
	 */	
	public void drawTriangle(Point_3 p, Point_3 q, Point_3 r) {
		float s= (float)this.scaleFactor;
		view.vertex( (float)(p.getX().doubleValue()*s), (float)(p.getY().doubleValue()*s), (float)(p.getZ().doubleValue()*s));
		view.vertex( (float)(q.getX().doubleValue()*s), (float)(q.getY().doubleValue()*s), (float)(q.getZ().doubleValue()*s));
		view.vertex( (float)(r.getX().doubleValue()*s), (float)(r.getY().doubleValue()*s), (float)(r.getZ().doubleValue()*s));
	}

	/**
	 * Draw a (triangle or polygonal) face
	 */	
	public void drawFace(Face<Point_3> f, double color, double maxColor) {
		Halfedge<Point_3> h=f.getEdge();
		Halfedge<Point_3> pEdge=h.getNext();
		
		Point_3 u=h.getOpposite().getVertex().getPoint();
		view.noStroke();
		//view.fill(200,200,200,255); // color of the triangle
		
		float[] rgb = fiveColorsHeatMap(color, maxColor); // ccall code 
		System.out.println(color);
		//float[] my_rgb = Color_Map.color_scheme(color_map, color, maxColor);
		
		view.fill(rgb[0], rgb[1], rgb[2], 255);

//		view.fill(my_rgb[0], my_rgb[1], my_rgb[2], 255);
		
		while(pEdge.getVertex()!=h.getOpposite().getVertex()) {
			Point_3 v=pEdge.getOpposite().getVertex().getPoint();
			Point_3 w=pEdge.getVertex().getPoint();
			this.drawTriangle(u, v, w); // draw a triangle face
			
			pEdge=pEdge.getNext();
		}
	}
	
	/**
	 * Draw the entire mesh
	 */
	public void draw(int type, double color[], double maxColor, Point_3 s) {
		//this.drawAxis();
		
		// draw all faces
		view.beginShape(view.TRIANGLES);
		int i = 0;
		for(Face<Point_3> f: this.polyhedron3D.facets) {
				this.drawFace(f, color[i], maxColor);
				i++;
		}
		view.endShape();
		
		if(type==1) return; // no rendering of edges
		
		// draw all edges
		view.strokeWeight(2); // line width (for edges)
		view.stroke(20);
		for(Halfedge<Point_3> e: this.polyhedron3D.halfedges) {
			Point_3 p=e.vertex.getPoint();
			Point_3 q=e.opposite.vertex.getPoint();
			
			this.drawSegment(p, q); // draw edge (p,q)
		}
		//view.strokeWeight(1);
		
		if(type==0) return; // no rendering for vertices
		
		view.noStroke();
		view.fill(0f, 0f, 250f);
		for(Vertex<Point_3> v: this.polyhedron3D.vertices) {
			this.drawVertex(v.getPoint(), s);
		}
		view.strokeWeight(1);
	}
	
	/**
	 * Draw the X, Y and Z axis
	 */
	public void drawAxis() {
		double s = 1;
		Point_3 p000=new Point_3(0., 0., 0.);
		Point_3 p100=new Point_3(s, 0., 0.);
		Point_3 p010=new Point_3(0.,s, 0.);
		Point_3 p011=new Point_3(0., 0., s);
		
		drawSegment(p000, p100);
		drawSegment(p000, p010);
		drawSegment(p000, p011);
	}

	/**
	 * Return the value after truncation
	 */
	public static double round(double x, int precision) {
		return ((int)(x*precision)/(double)precision);
	}
	
	/**
	 * Compute the scale factor (depending on the max distance of the point set)
	 */
	public double computeScaleFactor() {
		if(this.polyhedron3D==null || this.polyhedron3D.vertices.size()<1)
			return 1;
		double maxDistance=0.;
		Point_3 origin=new Point_3(0., 0., 0.);
		for(Vertex<Point_3> v: this.polyhedron3D.vertices) {
			double distance=Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
			maxDistance=Math.max(maxDistance, distance);
		}
		return Math.sqrt(3)/maxDistance*150;
	}
	
	/**
	 * Update the scale factor
	 */
	public void updateScaleFactor() {
		this.scaleFactor=this.computeScaleFactor();
	}
	
	/**
     * Heatmap function based on 5 colors
     *
     * @param v
     *         the input value, which must be in [0,max]
     * @return
     *         an array containing the three RGB components (between 0..255)
     */
    public static float[] fiveColorsHeatMap(double v, double max) {
            double value;
            if(max>0.)
            	value=v/max; // normalize value, between 0..1
            else value=0;

            int NUM_COLORS = 4;
            float[][] color = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
            // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
            int idx1;         // |-- Our desired color will be between these two indexes in "color".
            int idx2;         // |
            float fractBetween = 0;  // Fraction between "idx1" and "idx2" whereour value is.
            if(value <= 0) {  idx1 = idx2 = 0;            }    // accountsfor an input <=0
            else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
            else
            {
            	value = value * (NUM_COLORS-1);        // Will multiply value by 3.
            	idx1  = (int)(value);                  // Our desired color will be after this index.
            	idx2  = idx1+1;                        // ... and before this index (inclusive).
            	fractBetween = (float)value - (idx1);    // Distance between the two indexes(0-1).
            }

            float red   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
            float green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
            float blue  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
            float[] result=new float[]{red*255, green*255, blue*255};
            return result;
    }

    

    
}
