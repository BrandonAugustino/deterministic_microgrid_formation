# deterministic_microgrid_formation
This repository contains Python scripts for the deterministic microgrid formation optimization problem presented in "Allocation of Resources Using a Microgrid Formation Approach for Resilient Electric Grids" by Kwami S. A. Sedzro, Alberto J. Lamadrid and Luis F. Zuluaga. Currently there are 3 files in the repository. 

1) "model_total_new.py" A python script which solves the MILP as presented in the paper. The model is constructed and solved using the gurobipy module, thus in order to execute the script one must (i) Ensure that gurobi is installed on the computer they are using to run the code (ii) Install gurobipy: 

https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html

2) "ieee30_s5_python.dat" A data file corresponding to the IEEE 30-Bus sytem. This system has a mesh topology and to solve this instance, one should ensure that the data is loaded (around line 8) as follows:

file = open("../resil/ieee30_s5_python.dat","r")

3) "ieee118_s5_python.dat" A data file corresponding to the IEEE 118-Bus sytem. This system has a mesh topology and to solve this instance, one should ensure that the data is loaded (around line 8) as follows:

file = open("../resil/ieee118_s5_python.dat","r")

If you have any questions regarding the code please feel free to contact me:

Brandon Augustino
bra216@lehigh.edu
