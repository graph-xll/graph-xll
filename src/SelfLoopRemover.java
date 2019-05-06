public class SelfLoopRemover {
        it.unimi.dsi.webgraph.ImmutableGraph G;
        String basename1, basename2;

        public SelfLoopRemover(String basename1, String basename2) throws Exception {
                this.basename1 = basename1;
                this.basename2 = basename2;

                G = it.unimi.dsi.webgraph.Transform.filterArcs(
                                it.unimi.dsi.webgraph.ImmutableGraph.load(basename1),
                                it.unimi.dsi.webgraph.Transform.NO_LOOPS);
        }

        public void store() throws Exception {
                System.out.println("Storing graph without selfloops...");
                it.unimi.dsi.webgraph.BVGraph.store(G, basename2);
                System.out.println("Graph without selfloops stored!");
        }


        public static void main(String[] args) throws Exception {
                //args = new String[] {"edges02", "edges02_noself"};
                //args = new String[] {"data1_astrocnet", "data1_astrocnet_noself"};
        		//args = new String[] {"ljournal-2008-noself-t", "ljournal-2008-noself-t"};
        	
    			if(args.length != 2 ) {
    				System.out.println("Specify: basename basename_output");
    				System.exit(0);
    			}
    			
                SelfLoopRemover kc = new SelfLoopRemover(args[0], args[1]);
                kc.store();
        }
}
