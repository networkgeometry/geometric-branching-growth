### Codes and data sets for "Scaling up real networks by geometric branching growth".

#### Citation
[_Scaling up real networks by geometric branching growth_](https://doi.org/10.1073/pnas.2018994118)<br>
Muhua Zheng, Guillermo García-Pérez, Marián Boguñá, and M. Ángeles Serrano<br>
Proceedings of the National Academy of Sciences USA 118, e2018994118 (2021)

#### Abstract
Real networks often grow through the sequential addition of new nodes that connect to older ones in the graph. However, many real systems evolve through the branching of fundamental units, whether those be scientific fields, countries, or species. Here, we provide empirical evidence for self-similar growth of network structure in the evolution of real systems---the journal citation network and the world trade web---and present the Geometric Branching Growth model, which predicts this evolution and explains the symmetries observed. The model produces multiscale unfolding of a network in a sequence of scaled-up replicas preserving network features, including clustering and community structure, at all scales. Practical applications in real instances include the tuning of network size for best response to external influence and finite-size scaling to assess critical behavior under random link failures.

#### Dependencies
libstable: http://www.lpi.tel.uva.es/stable

 
#### 1. Generate_kappas
	In this folder, we will illustrate how to obtain the kappas of descendants from libstable software. 
	(1) Compiling the libstable software in current folder with the make command:
			$ make
		After compilation, both shared libstable.so and static (libstable.a) versions of the library are produced.
	(2) To illustrate the usage better, we give an example with citation network. There are 3 files about the original snapshot in folder "data":
		"**_edgelist.txt"-------edgelist file for each line inditates an edge node_i---node_j
		"**_coordinates.txt"----coordinates for each node. Each line inditates node kappa theta. Noted that node index has been reasigned by ascending theta. The corresponding edgelist are in consistent.  
		"**_parameters.txt"-----embedding parameters:  N  beta  mu

	(3) run function "Find_kappa.m" to generate kappa in each descendants' layer with command:
			[pars_fit]=Find_kappa(filepath, layers, ps) 
		where filepath is the path and rootname of files about the original snapshot, layers is the total descendants' layer, ps~[0,1] is the blanching probability so that b=1+ps.
		for example:
		[pars_fit]=Find_kappa('./data/TW_65', 10, 0.5) will load files from filepath './data/TW_65' and generate 10 layer descendants with blanching probability 0.5.
		The function returns some related files in 'data' folder:
		"TW_65_kappa_l_**.txt"----------------is the kappas in the descendants layer** with node_descendants node_ancestor kappa_descendants. These files are very important in the GBG.
		"TW_65_fitting_performance.txt"-------is the fitting performance of the complementary cumulative distribution with [z  Pc(z)_empirical Pc(z)_fitting] in each colomn
		"TW_65_fitting_stable_parameter.txt"--is the stalbe distribution fitting parameters [alpha, eta, c, d]	
		"TW_65_z_kappa_l_**.txt"--------------is the complementary cumulative distribution of z and kappa in each layer. Each colomn indicates [z  Pc(z) kappa Pc(kappa)].

#### 2. Estimate_a
	In this folder, we give the code "Estimate_a.cpp" to estimate parameter $a$. 
	(1)To run this code, you should obtain the kappas from stable distribution with Find_kappa.m in step 1.
    (2)You also need to obtain the slope of k as a function of N in log-log scale as well as the standard deviation. 
	(3)Change the input parameters and file name in main function in Estimate_a.cpp. 
	(4)Combile the code with: g++ Estimate_a.cpp -o Estimate and run with: ./Estimate

#### 3. Evolution	
	In this folder, we give the code "Growth_p_split_evolution.cpp" to Scale up a network. 
	(1) given a the blanching probability and get kappas in step 1. 
		For example: we have some related files in 'data' folder about snapshot TW_65.
	(2) combile the code with: g++ Growth_p_split_evolution.cpp -o evolution
	(3) run the code with: ./evolution <filepath> <a> <total_layers>. filepath is the path and rootname of files about the original snapshot, it should be same as in step 1. a is the parameter estimated in
	step 2. total_layers is the total layers you want to scale up. 
	For example: run the code: ./evolution ./data/TW_65 1.415 8. It will evolute the JCN network to layer 8 from time window TW_65 with a=1.415.

#### 4. Empirical_Data
	In this folder, we give the empirical data and the hyperbolic maps used in this paper.
		"**_edgelist.txt"-------edgelist file for each line inditates an edge node_i---node_j
		"**_coordinates.txt"----coordinates for each node. Each line inditates node kappa theta. Noted that node index has been reasigned by ascending theta. The corresponding edgelist are in consistent.  
		"**_parameters.txt"-----embedding parameters:  N  beta  mu
