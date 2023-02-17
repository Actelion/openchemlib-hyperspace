# openchemlib-hyperspace
Cheminformatics tools, workflows and pipelines for substructure search, virtual screening and data analysis

# Developers
- Thomas Liphardt

# Build
## General

Maven is used as build tool.
The project is organized as a Maven project "openchemlib-hyperspace" containing a number of different submodules.
The parent project has a fixed version "DEV-snapshot" / "dev".
The actual version of the different submodules is specified via the propery ${revision}.

# Submodules

**hyperspace-core**
Provides datastructures and algorithms to represent and search in combinatorial libraries (aka. synthon spaces).

**hyperspace-core-gui**
Provides GUI elements and server/client tools that enable search / virtual screening.  

**hyperspace-tools**
Contains tools to parse combinatorial libraries in synthon representation and create the hyperspace data files.

**hyperspace-sar**
Coming soon.

# hyperspace-core

Provides algorithms and datastructures to represent combinatorial libraries in synthon representation. Provides algorithms for fast fingerprint filtering, enabling fast substructure searching in enumerated libraries. Provides full algorithm for substructure searching in enumerated libraries. Provides implementations of similarity search algorithms for combinatorial libraries.

Most important classes:
com.idorsia.research.chem.hyperspace.SynthonSpace
Serializable class that represents a combinatorial library in synthon representation. Provides methods for fast substructure search.

com.idorisa.research.chem.hyperspace.SynthonSimilaritySpace3
Serializable class that extends a SynthonSpace by methods for similarity searching.

# hyperspace-core-gui
Provides a basic extendable GUI to provide search functionality. The GUI can be extended by implementing "search providers" that provide GUI and logic for searching. Contains default implementations for hyperspace similarity and substructure search. The GUI also provides tools to help with the implementation of server / client based searching. It further contains a substructure search service that provides substructure search capability for the datawarrior hyperspace substructure plugin.

## Hyperspace Server Manual

Example launch of server

`java -Xmx16G -jar <path_to_hs_server_jar> com.idorsia.research.chem.hyperspace.gui.HyperspaceServer -t <num_threads_per_query> -s <max_num_simultaneous_queries> -p <port> --init <path_to_json_config_file>`

`java -Xmx16G -jar <path_to_hs_server_jar> com.idorsia.research.chem.hyperspace.gui.HyperspaceServer -t 12 -s 4 -p 8090 --init server_config.json`

The file server_config.json just contains the different available search services. A search service has a specific specific search service name, the serviceProvider must be specified and the configuration for the service provider must be provided.

```
{
   "ServiceProviders":[
      {
         "ServiceName":"REAL Space SSS",
         "ServiceProvider":"HyperspaceSSS",
         "Config":{
            "SpaceName":"REAL Space",
            "File":"/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data",
            "MaxNumberOfThreads":-1
         }
      },
      {
         "ServiceName":"Idorsia VCS 3.2 (2S) SSS",
         "ServiceProvider":"HyperspaceSSS",
         "Config":{
            "SpaceName":"IVCS 3 (2S)",
            "File":"/home/liphath1/hyperspace_base_2/data_hyperspace/ivcs3_2s_FragFp.data",
            "MaxNumberOfThreads":-1
         }
      },
      {
         "ServiceName":"REAL Space SSS DWService",
         "ServiceProvider":"HyperspaceSSSForDW",
         "Config":{
            "SpaceName":"REAL Space",
            "File":"/home/liphath1/hyperspace_base_2/data_hyperspace/REAL_Space_latest_FragFp.data"
         }
      }
   ]
}
```


## Hyperspace GUI Manual
Example launch of GUI:

`java -Xmx16G -jar <path_to_hyperspace_server_jar> com.idorsia.research.chem.hyperspace.gui.HyperspaceSearchGUI <path_to_json_config_file>`

The structure of the config file for the GUI is the same as for the server, i.e. it is a list of service provider configurations. There are remote service providers available, that allow to send queries to a hyperspace server. The "SearchServiceName" has to agree with the configured service name of the service provider running on the server.

Example Standard (local) configuration

Exactly the same as the configuration of the server (however, only search providers that provide a GUI can be used, but most search providers do).

Example Remote configuration

```
{
   "ServiceProviders":[
      {
         "ServiceName":"REAL Space SSS",
         "ServiceProvider":"HyperspaceSSSRemote",
         "Config":{
            "Server":"hyperspaceserver",
            "Port":"8090",
            "SearchServiceName":"REAL Space SSS",
            "GUIServiceName":"REAL Space SSS"
         }
      },
      {
         "ServiceName":"IVCS3 2S SSS",
         "ServiceProvider":"HyperspaceSSSRemote",
         "Config":{
            "Server":"hyperspaceserver",
            "Port":"8090",
            "SearchServiceName":"Idorsia VCS 3.2 (2S) SSS",
            "GUIServiceName":"IVCS 3 (2s) SSS"
         }
      },
      {
         "ServiceName":"REAL Space Sim",
         "ServiceProvider":"HyperspaceSimilarityRemote",
         "Config":{
            "Server":"hyperspaceserver",
            "Port":"8090",
            "SearchServiceName":"REAL Space Similarity",
            "GUIServiceName":"REAL Space Similarity"
         }
      },
      {
         "ServiceName":"IVCS3 2S Sim",
         "ServiceProvider":"HyperspaceSimilarityRemote",
         "Config":{
            "Server":"hyperspaceserver",
            "Port":"8090",
            "SearchServiceName":"Idorsia VCS 3.2 (2S) Similarity",
            "GUIServiceName":"IVCS 3 (2s) Similarity"
         }
      }
   ]
}
```

# hyperspace-tools

Provides tools to create the hyperspace data files that can be loaded by the search providers of the tools in hyperspace-server.

Example for parsing / creating hyperspace data files for substructure search:

`java -Xmx16G  com.idorsia.research.chem.hyperspace.io.SynthonSpaceParser2 /home/liphath1/hyperspace_base_2/data/REAL_Space_latest/2021-02_REAL_synthons_SMALL.txt <name_of_space_file> FragFp <number_of_threads>`


This creates a synthon space with "FragFp" descriptors. This is the descriptor / space that is necessary for the substructure search services.

The input file format is a tsv (tab separated values) file with the following format:
first line is header line (and is not parsed)
all other lines are:
1. synthon in smiles, 2. a synthon id, 3. an integer specifying the synthon set of the synthon reaction, 4. the synthon reaction

i.e. the following file would create two synthon reactions, each with two synthon sets (synthon reaction r1 has 2x2 synthons, synthon reaction r2 has 1x3 synthons to assemble). Of course the synthons must have matching connectors (i.e. matching transuranium elements (atomic mass 92-99), that uniquely identify each possible connection point). All synthons of a synthon set of a reaction must have the same connectors (i.e. a single U, or Np and Pu, etc.). The correctness of synthons is checked during the creation of the datastructures and in case of problems, the software will output error messages for the problematic synthons / synthon reactions. 

```
smiles   synthon  synton#   reaction_id
smiles1	r1_1_s1	1	r1
smiles2	r1_1_s2	1	r1
smiles3	r1_2_s1	2	r1
smiles4	r1_2_s2	2	r1
smiles5	r2_1_s1	1	r2
smiles6	r2_2_s1	2	r2
smiles7	r2_2_s2	2	r2
smiles8	r2_2_s3	2	r2
