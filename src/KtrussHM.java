/*
 * k-truss decomposition
 * Serial algorithm without optimization (HashMap version)
 * @author Alison
 */
import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.lang.Math;
import it.unimi.dsi.webgraph.ImmutableGraph;

public class KtrussHM {
	ImmutableGraph G;
	int n;
	int m;
	int ms; //max support
	int[] sup;
	int[] A; //stores first edge point
    int[] B; //stores second edge point
	//List<Integer>[] keys;
	HashMap<List<Integer>, Integer> edge_ids;
	public KtrussHM(String basename) throws Exception {
		G = ImmutableGraph.load(basename);

		n = G.numNodes();
		m = (int) (G.numArcs() / 2); //sym graph includes edge for each direction

		ms = 0;
		sup = new int[m];
		A = new int[m];
		B = new int[m];
		//keys = new List<Integer>[m];

		edge_ids = new HashMap<List<Integer>, Integer>(m);
		//List<Integer> key;
		//int[] vN, uN; //successor arrays of adjacent vertices
		//Set<Integer>  v_set;
		int vNLen, uNLen;
		int edge_num = 0;
		int j = 0;
		int k = 0;
		int edge_uv_sup = 0;
		for(int v=0; v<n; v++) {
			int[] vN = G.successorArray(v);
			vNLen = vN.length;
			for(int i=0; i < vNLen; i++){
				int u = vN[i];
				if(u > v){
					/*
					v_set = new HashSet<Integer>(Arrays.stream(vN).boxed().collect( Collectors.toList() ));
					v_set.retainAll(new HashSet<Integer>(Arrays.stream(G.successorArray(u)).boxed().collect( Collectors.toList() )));
					edge_uv_sup = v_set.size() + 2; //adding two to each support value to get truss numbers in final result
					*/
					int[] uN = G.successorArray(u);
					uNLen = uN.length;

					//count the number of shared neighbours using two-pointer method
					while(j < vNLen && k < uNLen){
						if(vN[j] < uN[k]){
							j++;
						}else if(uN[k] < vN[j]){
							k++;
						}else{
							edge_uv_sup++;
							j++;
							k++;
						}
					}
					sup[edge_num] = edge_uv_sup;
					//initialize edge arrays in order of edge numbers
					A[edge_num] = v;
					B[edge_num] = u;

					//store edge_num in hash map
					List<Integer> key = new ArrayList<Integer>(2);
					key.add(v);
					key.add(u);
					edge_ids.put(key, edge_num);
					//keys[edge_num] = key;

					edge_num++;

					if(ms < edge_uv_sup)
						ms = edge_uv_sup;

					//for next intersection
					edge_uv_sup = 0;
					j = 0;
					k = 0;
				}
			}
		}
	}
    public int[] KTrussCompute () {
		long t0 = System.currentTimeMillis();
		System.out.println("Starting step 2: sorting edges by support");

    	int[] pos = new int[m]; //position of edge numbers in sortedEdges
		int[] sortedEdges = new int[m];
    	int[] bin = new int[ms+1];

    	for(int s=0; s<=ms; s++)
    		bin[s] = 0;
    	for(int e=0; e<m; e++)
    		bin[ sup[e] ]++;

    	int start = 0;
    	for(int s=0; s<=ms; s++) {
    		int num = bin[s];
    		bin[s] = start;
    		start += num;
    	}

    	//bin-sort edges by support
    	for(int e=0; e<m; e++) {
    		pos[e] = bin[ sup[e] ];

			//store edge numbers in order in separate array
			sortedEdges[ pos[e] ] = e;

			bin[ sup[e] ]++;
    	}
    	//recover bin[]
    	for(int s=ms; s>=1; s--)
    		bin[s] = bin[s-1];
		bin[0] = 0;

		long t1 = System.currentTimeMillis();
		System.out.println("Step 2: Time elapsed (sec) = " + (t1 - t0)/1000.0);

		System.out.println("Starting step 3: computing the truss numbers of edges");
		//step 3: "pealing" edges
		Integer edge_id;
		List<Integer> key;
		for(int e=0; e<m; e++){
			int u = A[ sortedEdges[e] ];
			int v = B[ sortedEdges[e] ];

			for(int w : G.successorArray(u)){
				key = new ArrayList<Integer>();
				key.add(Math.min(v,w));
				key.add(Math.max(v,w));
				edge_id = edge_ids.get(key);
				if(edge_id != null){
					//edge (v,w) exists
					int edge_vw = edge_id.intValue();
					int sup_vw = sup[edge_vw]; //support of (v,w)
					if(sup_vw > sup[e]){
						int pos_vw = pos[edge_vw]; //position of (v,w) in sortedEdges
						int pos_firstEdge = bin[sup_vw]; //position of first edge with same support
						int firstEdge = sortedEdges[pos_firstEdge];
						if(edge_vw != firstEdge){
							//swap
							pos[edge_vw] = pos_firstEdge;
							sortedEdges[pos_vw] = firstEdge;
							pos[firstEdge] = pos_vw;
							sortedEdges[pos_firstEdge] = edge_vw;
						}
						bin[sup_vw]++;
						sup[edge_vw]--;
					}

					//repeat for edge (u,w)
					key = new ArrayList<Integer>();
					key.add(Math.min(u,w));
					key.add(Math.max(u,w));
					edge_id = edge_ids.get(key);
					int edge_uw = edge_id.intValue();
					int sup_uw = sup[edge_uw]; //support of (u,w)
					if(sup_uw > sup[e]){
						int pos_uw = pos[edge_uw]; //position of (u,w) in sortedEdges
						int pos_firstEdge = bin[sup_uw]; //position of first edge with same support
						int firstEdge = sortedEdges[pos_firstEdge];
						if(edge_uw != firstEdge){
							//swap
							pos[edge_uw] = pos_firstEdge;
							sortedEdges[pos_uw] = firstEdge;
							pos[firstEdge] = pos_uw;
							sortedEdges[pos_firstEdge] = edge_uw;
						}
						bin[sup_uw]++;
						sup[edge_uw]--;
					}
				}
			}
		}
		System.out.println("Step 3: Time elapsed (sec) = " + (System.currentTimeMillis() - t1)/1000.0);
		return sup;
	}
	public static void main(String[] args) throws Exception {
		System.out.println("Please enter graph's base name");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

		System.out.println("Starting " + basename);

		long t0 = System.currentTimeMillis();
		System.out.println("Starting step 1: computing support of edges");
		KtrussHM kt = new KtrussHM(basename+"-noself-sym");
		System.out.println("Step 1: Time elapsed (sec) = " + (System.currentTimeMillis() - t0)/1000.0);

		int[] res = kt.KTrussCompute();

		// //storing the truss number of each edge in a file
		// PrintStream ps = new PrintStream(new File(basename + "_BZresults"));
		// ps.println("Truss number for each edge in " + basename);
		//
		// for(int i=0; i<res.length; i++) {
		// 	ps.println("edge " + i + ":" + (res[i] + 2) + " "); //adding two to each support value to get truss numbers in final result
		// }
		System.out.println("n	m	smax");
		System.out.println(kt.n + "\t" + (kt.m) + "\t" + kt.ms);

	}
}
