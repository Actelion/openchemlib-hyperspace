Hyperspace CLI 0.3


To run any of the Hyperspace CLI services, you have to run java programs from within the hyperspace-cli.jar file, see below for example calls.

To specify the input of the algorithms you have to provide usually two files:
1. a config file (a json file), see sss_conf_a.json for a template
2. a file containing input structures. This can be either a dwar file, or a text file containing idcodes that is referenced from the config file, e.g. see test_idcodes.txt or test_structures.dwar.


SUBSTRUCTURE SEARCH

To run substructure search, use the following call, where "sss_conf_a.json" is your configuration file for the substructure search. The provided example file references the input file "test_idcodes.txt" containing the query idcodes:

java -Xmx16G -cp hyperspace-cli.jar com.idorsia.research.chem.hyperspace.services.HyperspaceSubstructureSearchService sss_conf_a.json

Note: if things work correclty, the results file should contain hits for the first structure, the second structure cannot be found in REAL space.



SAR-DIRECTED HIT EXPANSION

To run sar-directed hit expansion, exactly the same kind of config file can be used. Please use a call similar to:

java -Xmx16G -cp hyperspace-cli.jar com.idorsia.research.chem.hyperspace.services.HyperspaceSARDirectedHitExpansionService sss_conf_a.json


Note: here, the output contains one additional column, indicating the specific substructure search that was performed as part of the sar-directed hit expansion.


PAIRS OF CONNECTED SUBGRAPHS SEARCH

This is a method that can be employed in case that SAR-directed hit expansion does not find any hits.
The method searches for similar molecules that share certains substructures. The method works in multiple steps:
1. an algorithm creates a specific number of representative substructures of the query molecule
2. for pairs of representative structures that are within a specific distance, a Hyperspace substructure search is performed, searching for the two fragments and some in-between part.


To this search method, there are 5 additional integer parameters that have to be specified. Please use a call similar to:

java -Xmx16G -cp hyperspace-cli.jar com.idorsia.research.chem.hyperspace.services.HyperspaceConnectedFragmentsSearchService sss_conf_connectedSubgraphs.json

NOTE: This search may generate very large number of results, depending on the query.



RELEVANT SUBGRAPHS SEARCH

This is a method that can be employed in case that SAR-directed hit expansion does not find any hits.
The method uses a method similar to the first step in the "pairs of connected subgraphs search", i.e. it computes subgraphs of the query molecule.
These subgraphs are then used for substructure search in the specified combinatorial library.

To this search method, there are 5 additional integer parameters that have to be specified. Please use a call similar to:

java -Xmx16G -cp hyperspace-cli.jar com.idorsia.research.chem.hyperspace.services.HyperspaceRelevantSubgraphSearchService sss_conf_relevantSubgraphs.json

NOTE: This search may generate very large number of results, depending on the query.




Options in the configuration file for substructure search and sardhx :

"input": path to input file containing structures
"output": path to output file (file format is dwar)
"space": specifies the space to load ( currently for real space, please use "real_space_2020_04_FragFp.data"
"maxCombiHits": maximum number of combinatorial hits, 32 will be reached only rarely
"maxExpandedHits": maximum number of total enumerated hits for a query, or -1 for no limit
"maxExpandedHitsPerCombiHit": maximum number of enumerated hits for each combinatorial hit, or -1 for no limit
"minFFPQueryBits": this is a simple feature to avoid "too simple / general" queries, by counting the number of bits set in the ffp descriptor to estimate complexity of the query. You can set this to zero if you are not interested in this. The current value (24) should be exceeded by I think all reasonable queries, and only sort out very general queries (e.g. only a chain of c atoms, or only an unsubstituted single ring system etc.)
"maxQueryTimeMillis": maximum time for a query in millisceonds until the query is aborted and skipped
"threads": number of threads to use for computation, or -1 to ask the system for the number of cores
"maxSplits": maximum number of splits used in hyperspace. Default value is 3, and most likely you should use this. However, valid options are 1,2 or 3 (because in REAL Space there are synthon reactions with up to 3 splits). Decreasing the value from 3 to 2 or 1 will exclude search in synthon reactions with higher number of splits and therefore miss results, but max work much faster (but this should not really be needed).


Additional parameters for pairs of connected subgraphs search :
"minFragSize": minimum size (atoms) of considered fragments
"maxFragSize": maximum size (atoms) of considered fragments
"minDist": min dist in between fragments to consider the pair
"maxDist": max dist in between fragments to consider the pair
"numClusters": number of representative fragments that are computed per query structure


Additional parameters for relevant subgraphs search :

"minFragSize": minimum size (atoms) of considered fragments
"maxFragSize": maximum size (atoms) of considered fragments
"minFFPBits": min dist in between fragments to consider the pair
"maxFFPBits": max dist in between fragments to consider the pair
"numClusters": number of representative fragments that are computed per query structure



IMPORTANT REMARKS:

There is one thing to keep in mind when running large number of queries in a single run (large meaning hundreds or thousands of substructure queries): currently, all enumerated structures for all hits are stored in memory, and only at the very end they will be written into the output file. This means, if you have thousands of queries with each > thousands of results, then the java program may run out of memory.

