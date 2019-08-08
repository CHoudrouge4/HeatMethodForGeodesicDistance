import processing.core.*;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 * and edited by Hussein Houdrouge for the purpose of the project.
 * 
 */
public class MeshViewer extends PApplet {

	SurfaceMesh mesh; // 3d surface mesh
	HeatMethod hm;
	int renderType  = 0; // choice of type of rendering
	int renderModes = 3; // number of rendering modes
	
	//String filename="high_genus.off";
	String filename="sphere.off";
	//String filename="cube.off";
	//String filename="torus_33.off";
	//String filename="tore.off";
	//String filename="tri_hedra.off";
	//String filename="letter_a.off";
	//String filename="OFF/star.off";
	//String filename = "bague.off";
	//String filename = "out.off";
	//String filename = "bunny.off";
	//String filename = "tetrahedron.off";
	//String filename = "tri_cow.off";	
	//String filename = "camel.off";
	public void setup() {		
		  size(800,600,P3D);
		  ArcBall arcball = new ArcBall(this);
		  
		  this.mesh = new SurfaceMesh(this, filename);
		  
		  hm = new HeatMethod(this.mesh.polyhedron3D, 0);
		 /// System.out.println("The time t = " + hm.getTime());
		  hm.heatDiffusion(0);
		  hm.computePhi(0);
		  hm.colorScheme('d');

		  // System.out.println(" " + hm.polyhedron3D.facesToString());
		  //ms.polyhedron3D.isValid(false);
		 
	}
		 
		public void draw() {
		  background(0);
		  //this.lights();
		  directionalLight(101, 204, 255, -1,  0, 0);
		  directionalLight(51,  102, 126,  0, -1, 0);
		  directionalLight(51,  102, 126,  0,  0, -1);
		  directionalLight(102,  50, 126,  1,  0, 0);
		  directionalLight(51,   50, 102,  0,  1, 0);
		  directionalLight(51,   50, 102,  0,  0, 1);
		 
		  translate(width/2.f, height/2.f, -1 * height/2.f);
		  this.strokeWeight(1);
		  stroke(150, 150, 150); 
		  this.mesh.draw(renderType, hm.colors, hm.MaxColor, hm.polyhedron3D.vertices.get(hm.s).getPoint());
		}
		
		public void keyPressed(){
			  switch(key) {
			    case('r'):this.renderType = (this.renderType + 1) % this.renderModes;
			    break;
			  }
		}
		
		/**
		 * For running the PApplet as Java application
		 */
		public static void main(String args[]) {
			//PApplet pa=new MeshViewer();
			//pa.setSize(400, 400);
			PApplet.main(new String[] { "MeshViewer" });
		}
}