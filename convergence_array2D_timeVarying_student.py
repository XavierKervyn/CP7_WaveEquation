# Ce code est une code simplifier pour l'etude de convergence d'un problem 2D temps variantes
import os
import numpy as np
import scipy.stats.mstats as mstats
import matplotlib.pyplot as plt
import AnalyticalSolutions as ansol

# Data
do_simulations = True
exec_name = "Exercice7_solution"
delimiter = " "
font_size = 32
marker_size = 18
line_width = 3
legend_location = 1
filename_grid = "output_mesh_"
filename_solution = "output_f_"
output_grid = "output_mesh"
output_file = "output_file"
filename_config = 'configuration.in'
convergence_study = 'CFL'; # values: CFL, Nx,Ny
mesh_values = [-2.,0.,10]#[1.,2.,10]
mesh_type = " "#"int"
time = 0.353553390593274
# input array
inputs = {'tfin':time,'xL':0.,'xR':10.,'yL':0.,'yU':10.,\
'pert_amplitude':0.,'pert_velocity':0.,'u':4.,'g':9.81,'h0':2.,\
'h':1.,'Nx':300,'Ny':300,'CFL':0.95,'type_u2':'const',\
'outputDataProfiles':'outputDataProfiles.out','ecrire_f':'true',\
'mode_num_x':1,'mode_num_y':1,'bc_left':'dirichlet','bc_right':'dirichlet',\
'bc_lower':'dirichlet','bc_upper':'dirichlet','impulsion':'false',\
'type_init':'harmonic','F0':1.,'A':0.,'omega':0.e0,'write_mesh':'true',\
'write_f':'true','n_stride':0,'output_mesh':'output_mesh.out',\
'output_file':'output_f.out','output_energy':'output_E.out',\
'output_velocity':'output_u.out'}

# initialisation
x_mesh = []
y_mesh = []
error_list = []
solution = []

# TODO: Generate Mesh
convergence_mesh = [0.]
if(mesh_type == "int"):
  convergence_mesh = np.ceil(convergence_mesh)

# Execute simulation and read data
for element in convergence_mesh:
  print(element)
  output_mesh_name = "".join([filename_grid,convergence_study,str(element),".out"])
  output_f_name = "".join([filename_solution,convergence_study,str(element),".out"])
  inputs[output_grid] = output_mesh_name
  inputs[output_file] = output_f_name 
  inputs[convergence_study] = element
  with open(filename_config,'w') as fid_input:
    for key,value in inputs.items():
      fid_input.write(''.join([key,' = ',str(value),'\n']))
  # check whether new simulations should be performed
  if(do_simulations):
      print("".join(["simulating element: ",str(element)]))
      os.system(''.join(['./',exec_name,' ',filename_config]))  
  # read the results
  print("".join(["Reading solution element: ",str(element)]))
  with open(output_mesh_name,'r') as mesh_file:
    x_mesh.append(np.fromstring(mesh_file.readline(),dtype=np.float64,sep=delimiter))
    y_mesh.append(np.fromstring(mesh_file.readline(),dtype=np.float64,sep=delimiter))
  results = np.loadtxt(output_f_name,dtype=np.float64,delimiter=delimiter)     
  solution.append(results[results[:,0]==time,1:])
    
#TODO Compute error
error_list = [-999]
print(error_list)

# TODO Compute linear regression
slope_value = [-999]
print(slope_value)

# Plot convergence test
plt.rcParams.update({'font.size': font_size})
plt.figure(1)
plt.loglog(convergence_mesh,error_list,'o-',linewidth=line_width,markersize=marker_size)
plt.legend([slope_value],loc=legend_location)
plt.title("".join(["Convergence study w.r.t ",convergence_study]))
plt.xlabel(convergence_study)
plt.ylabel("error")
plt.grid()
plt.show()
