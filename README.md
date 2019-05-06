# Graph-XLL
> A Graph Library for Extra Large Graph Analytics on a Single Machine
***

## Input

The input graphs are in WebGraph format.

There are three files in this format:

*basename.graph* 
*basename.properties* 
*basename.offsets*

Available datasets in this format can be found in: <http://law.di.unimi.it/datasets.php>

The first two files are for the forward (regular) graph. The other two are for the transposed graph. What's missing is the "offsets" file. This can be easily created by running: 
```
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000
```

### *Converting Edgelist Format to WebGraph Format*

This section is for the case when your graph is given a text file of edges (known as edgelist). *If your graph is already in WebGraph format, skip to the next section.*

It is very easy to convert an edgelist file into WebGraph format. I am making the folloiwng assumptions:

1. The graph is unlabeled and the vertices are given by consecutive numbers, 0,1,2,... 
   (If there are some vertices "missing", e.g. you don't have a vertex 0 in your file, it's not a problem. WebGraph will create dummy vertices, e.g. 0, that does not have any neighbor.)
2. The edgelist file is TAB separated (not comma separated).

Now, to convert the edgelist file to WebGraph format execute the following steps:

Sort the file, then remove any duplicate edges:

```
sort -nk 1 edgelistfile | uniq > edgelistsortedfile
```

(If you are on Windows, download *sort.exe* and *uniq.exe* from <http://gnuwin32.sourceforge.net/packages/coreutils.htm>)

Run:

```
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -g ArcListASCIIGraph edgelistsortedfile  dummyBasename
```

For example:

```
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -g ArcListASCIIGraph car-2000.txt  cnr-2000
```

### *Undirected Graphs in WebGraph*

Undirected graphs are used in some algorithms, e.g., k-truss decomposition. After creating the offset files, remove self-loops in the graphs by running:

```
java -cp "lib/*:bin" SelfLoopRemover cnr-2000 cnr-2000-noself
java -cp "lib/*:bin" SelfLoopRemover cnr-2000-t cnr-2000-noself-t
```

Symmetrize by taking union:

```
java -cp "lib/*" it.unimi.dsi.webgraph.Transform union cnr-2000-noself cnr-2000-noself-t cnr-2000-noself-sym
```

Now we obtain the undirected simple graph in WebGraph as the input:

*cnr-2000-noself-sym.graph* 
*cnr-2000-noself-sym.properties*
*cnr-2000-noself-sym.offsets*

## Compiling

The programs are already compiled with javac 1.8 (Java 8). Please use Java 8 or newer version.

```
javac -cp "lib/*" -d bin src/*
```

## Running

```
java -cp "lib/*:bin" PageRankWG
```
For larger graphs, allocate more memory to avoid OutOfMemoryError by appending -Xmx4g or -Xmx8g or -Xmx16g etc. to the above line.

## Acknowledgements

We would like to thank:
* Wissam Khaouid and Marina Barsky for developing the [k-core decomposition](https://github.com/athomo/kcore) program.
* Diana Popova for developing the [influence-maximization](https://github.com/dianapopova/InfluenceMax) program.
* Yudi Santoso for developing the [triad enumeration](https://github.com/ySant/triads) program.
* Michael Simpson for developing the [feedback-arc-set](https://github.com/stamps/FAS) program.
