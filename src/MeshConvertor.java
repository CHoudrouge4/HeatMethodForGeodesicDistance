import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

public class MeshConvertor {
	
	public static void convertToOFFfromMesh(String fileName) throws FileNotFoundException {
		Scanner scan = new Scanner(new File(fileName));
		PrintWriter print = new PrintWriter("out.off");
		int vert = 0, elem = 0;
		scan.nextLine();
		String str = "";
		while(scan.hasNextLine()) {
			String s = scan.nextLine();
			if(s == "elements") {
				int size = scan.nextInt();
				for(int i = 0; i < size; ++i) {
					scan.nextInt(); 
					int type = scan.nextInt();
					if(type == 2) {
						str += scan.nextInt() + " " + scan.nextInt() + " " + scan.nextInt() + '\n'; 		
						elem++;
					} else {
						scan.nextLine();
					}
				} 
			}
			if(s == "vertices") {
				vert = scan.nextInt();
				while(!scan.hasNextDouble()) scan.nextLine();
				for(int i = 0; i < vert; ++i) {
					s += scan.nextDouble() + " " + scan.nextDouble() + " " + scan.nextDouble() + '\n';
				}
			}
		}
		
		
		print.close();
		scan.close();
	}
}
