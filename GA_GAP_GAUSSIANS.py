#correct  binding energies with areas.

#load areas 
import pickle
import os
vasp_or_pnn='pnn'
main_dir = os.getcwd()
method='gdac'

#############################load binding energies###########################
def configurations(system,sub_system,suby_system):
    standardize_with_vdw = 'yes'
    if sub_system=='GAP-GAP':
        name_of_system=system
    else:
        name_of_system=suby_system+'_'+system
    vdw_cutoff=2
    if system == 'AA':
        max_num_configs = 576
        num_atoms=26 + 26
    if system == 'PP':
        max_num_configs = 576
        num_atoms= 35 + 35
    if system == 'GAP':
        max_num_configs = 576
        num_atoms= 140 + 26 + 35
    if system == 'GA':
        max_num_configs = 24
        num_atoms = 140 +26
    if system == 'GP':
        max_num_configs = 24
        num_atoms= 140 + 35
    if system == 'AP':
        max_num_configs = 576
        num_atoms= 26+35
    import os


    system_path = '/nas/longleaf/home/jarkeith/PNN_EDA/'+str(method)+'/'+str(system)
    os.chdir('/nas/longleaf/home/jarkeith/calculate_area_of_molecule/PNN/'+str(name_of_system))
   

    # Load the data from the file
    file_path = str(name_of_system)+'_shadow_area.pickle'
    with open(file_path, 'rb') as file:
        list1_loaded, list2_loaded = pickle.load(file)
    good_geometries=list1_loaded[0]
    shadow_area=list2_loaded[0]
    #print(good_configurations_list)
    #print(vdw_overlap_areas)


    os.chdir('/nas/longleaf/home/jarkeith/PNN_EDA')
    os.chdir( os.getcwd()+'/'+str(method)+'/'+str(system)+'/results')

    filename = "good_geometries.txt"
    good_geometries =[]
    with open(system_path+'/'+filename, "r") as f:
        for line in f:
            good_geometries.append(line.strip())

    from numpy import loadtxt
    import numpy as np
    from ase.visualize import view
    import matplotlib.pyplot as plt
    from ase.visualize.plot import plot_atoms
    #import text file into NumPy array
    if sub_system == 'none':
        tot = loadtxt('total_energy.txt')
        disp = loadtxt('dispersion.txt')
        repul = loadtxt('repulsion.txt')
        es = loadtxt('electrostatics.txt')
    if sub_system == 'GA-GAP':
        num_atoms= 140+26
        tot = loadtxt('total_energy_12.txt')
        disp = loadtxt('dispersion_12.txt')
        repul = loadtxt('repulsion_12.txt')
        es = loadtxt('electrostatics_12.txt')
    if sub_system == 'GP-GAP':
        num_atoms= 140+35
        tot = loadtxt('total_energy_13.txt')
        disp = loadtxt('dispersion_13.txt')
        repul = loadtxt('repulsion_13.txt')
        es = loadtxt('electrostatics_13.txt')
    if sub_system == 'AP-GAP':
        num_atoms= 35+26
        tot = loadtxt('total_energy_23.txt')
        disp = loadtxt('dispersion_23.txt')
        repul = loadtxt('repulsion_23.txt')
        es = loadtxt('electrostatics_23.txt')
    if sub_system == 'GAP-GAP':
        num_atoms= 140+26 +35
        tot_1 = loadtxt('total_energy_12.txt')
        disp_1 = loadtxt('dispersion_12.txt')
        repul_1 = loadtxt('repulsion_12.txt')
        es_1 = loadtxt('electrostatics_12.txt')
        tot_2 = loadtxt('total_energy_13.txt')
        disp_2 = loadtxt('dispersion_13.txt')
        repul_2 = loadtxt('repulsion_13.txt')
        es_2 = loadtxt('electrostatics_13.txt')
        tot=np.add(tot_1,tot_2)
        es=np.add(es_1,es_2)
        disp=np.add(disp_1,disp_2)
        repul=np.add(repul_1,repul_2)

    os.chdir(main_dir)
    disordered_list_of_indexes = good_geometries
    my_list = disordered_list_of_indexes
    disordered_list_of_indexes = [int(x) for x in my_list]
    ordered_list_of_indexes = np.arange(max_num_configs)

    disordered_tot_values = tot
    ordered_tot_values = np.zeros(max_num_configs)
    ordered_tot_values[disordered_list_of_indexes] = disordered_tot_values

    disordered_es_values = es
    ordered_es_values = np.zeros(max_num_configs)
    ordered_es_values[disordered_list_of_indexes] = disordered_es_values

    disordered_disp_values = disp
    ordered_disp_values = np.zeros(max_num_configs)
    ordered_disp_values[disordered_list_of_indexes] = disordered_disp_values

    disordered_repul_values = repul
    ordered_repul_values = np.zeros(max_num_configs)
    ordered_repul_values[disordered_list_of_indexes] = disordered_repul_values

    tot=ordered_tot_values
    es=ordered_es_values
    disp=ordered_disp_values
    repul=ordered_repul_values

    good_disp = []
    good_repul = []
    good_es = []
    #good_affinity=[]
    import numpy as np
    negative_tot=[]
    index = []
    positive_tot=[]

    if standardize_with_vdw != 'no':
        area_corrected_negative_tot=[]
        area_corrected_good_disp=[]
        area_corrected_good_repul=[]
        area_corrected_good_es=[]
        #area_corrected_good_affinity=[]

    num_atoms=1
    for i in range(0,max_num_configs): 
        #if tot[i] <= -0.01:
        #temp_tot=tot[i]*0.0104/num_atoms
        #temp_disp=disp[i]*0.0104/num_atoms
        #temp_repul=repul[i]*0.0104/num_atoms
        #temp_es=es[i]*0.0104/num_atoms
        negative_tot.append(tot[i]*0.0104/shadow_area[i]*num_atoms)
        good_disp.append(disp[i]*0.0104/shadow_area[i]*num_atoms)
        good_repul.append(repul[i]*0.0104/shadow_area[i]*num_atoms)
        good_es.append(es[i]*0.0104/shadow_area[i]*num_atoms)
        index.append(i)
        #good_affinity.append((disp[i]*0.0104+es[i]*0.0104)/repul[i]*0.0104)
        #good_affinity.append((temp_es+temp_disp/temp_repul))

    final_good_geometries=[]

    print(str(vasp_or_pnn)+' - '+str(sub_system))
    print("number of good geometries",len(index))
    print("vdw cutoff",vdw_cutoff)
    print('#########################################')
    print('############ Direct Energies ############')
    print('#########################################')


    std_negative_tot=np.average(np.std(negative_tot, axis=None, dtype=None, out=None, ddof=0))
    std_good_es=np.average(np.std(good_es, axis=None, dtype=None, out=None, ddof=0))
    std_good_disp=np.average(np.std(good_disp, axis=None, dtype=None, out=None, ddof=0))
    std_good_repul=np.average(np.std(good_repul, axis=None, dtype=None, out=None, ddof=0))
    #std_non_local_affinity=np.average(np.std(good_affinity, axis=None, dtype=None, out=None, ddof=0))
    #std_non_local_affinity=(std_good_es+std_good_disp)+std_good_repul

    print("E_int : ",np.average(negative_tot) ,std_negative_tot)
    print("E_es : ",np.average(good_es), std_good_es)
    print("E_disp : ",np.average(good_disp), std_good_disp)
    print("E_rep : ",np.average(good_repul), std_good_repul)
    #print((np.average(good_disp)+np.average(good_es))/np.average(good_repul),std_non_local_affinity)
    #print("Affinity : ",(np.average(good_affinity),std_non_local_affinity))

    print('#########################################')
    print("Min E_int : " , np.min(negative_tot))
    print("Min E_es : " , np.min(good_es))
    print("Min E_disp : " , np.min(good_disp))
    print("Min E_rep : " , np.min(good_repul))
    #print("Min Affinity : " , np.min(good_affinity))

    with open("/nas/longleaf/home/jarkeith/VASP/GAP/results/"+str(sub_system)+'_'+str(system)+".txt", "w") as f:
        f.write(str(vasp_or_pnn)+'/'+str(system)+' '+str(sub_system)+'\n')
        f.write("number of good geometries: " + str(len(index))+'\n')
        f.write('#########################################\n')
        f.write('############ Direct Energies ############\n')
        f.write('#########################################\n')
        f.write("Avg E_int : " + str(np.average(negative_tot)) + " " + str(std_negative_tot)+'\n')
        f.write("Avg E_es : " + str(np.average(good_es)) + " " + str(std_good_es)+'\n')
        f.write("Avg E_disp : " + str(np.average(good_disp)) + " " + str(std_good_disp)+'\n')
        f.write("Avg E_rep : " + str(np.average(good_repul)) + " " + str(std_good_repul)+'\n')
        #f.write("Avg Affinity : " + str(np.average(good_affinity)) + " " + str(std_non_local_affinity)+'\n')
        f.write('#########################################\n')
        f.write("Min E_int : " + str(np.min(negative_tot)) +'\n')
        f.write("Min E_es : " + str(np.min(good_es)) +'\n')
        f.write("Min E_disp : " + str(np.min(good_disp)) +'\n')
        f.write("Min E_rep : " + str(np.min(good_repul)) +'\n')
       # f.write("Min Affinity : " + str(np.min(good_affinity)) +'\n')


    return negative_tot , good_es , good_repul, good_disp 

vasp_or_pnn='VASP'
standardize_with_vdw='no'
#method= 'gdac' # gdac

system ='GA'
sub_system = 'none'
suby_system=''
negative_tot , good_es , good_repul, good_disp = configurations(system,sub_system,suby_system)
GA_negative_tot=negative_tot
GA_good_es=good_es
GA_good_repul=good_repul
GA_good_disp=good_disp

system ='GP'
sub_system = 'none'
suby_system=''
negative_tot , good_es , good_repul, good_disp =configurations(system,sub_system,suby_system)
GP_negative_tot=negative_tot
GP_good_es=good_es
GP_good_repul=good_repul
GP_good_disp=good_disp


#system ='GAP'
#sub_system = 'GA-GAP'
#suby_system='GA'
#negative_tot , good_es , good_repul, good_disp=configurations(system,sub_system,suby_system)
#GA_GAP_negative_tot=negative_tot
#GA_GAP_good_es=good_es
#GA_GAP_good_repul=good_repul
#GA_GAP_good_disp=good_disp

#system ='GAP'
#sub_system = 'GP-GAP'
#suby_system='GP'
#negative_tot , good_es , good_repul, good_disp =configurations(system,sub_system,suby_system)
#GP_GAP_negative_tot=negative_tot
#GP_GAP_good_es=good_es
#GP_GAP_good_repul=good_repul
#GP_GAP_good_disp=good_disp

system ='GAP'
sub_system = 'GAP-GAP'
suby_system=''
negative_tot , good_es , good_repul, good_disp =configurations(system,sub_system,suby_system)
GAP_GAP_negative_tot=negative_tot
GAP_GAP_good_es=good_es
GAP_GAP_good_repul=good_repul
GAP_GAP_good_disp=good_disp

#system ='GAP'
#sub_system = 'AP-GAP'
#negative_tot , good_es , good_repul, good_disp =configurations(system,sub_system)
#AP_GAP_negative_tot=negative_tot
#AP_GAP_good_es=good_es
#AP_GAP_good_repul=good_repul
#AP_GAP_good_disp=good_disp


########################################## export corrected values #################################################

num_bins = 10

# Set the dpi parameter to a high value for high resolution
dpi = 300

import numpy as np
import matplotlib.pyplot as plt

# Assuming you have three lists, each containing the binding energy values for a system
system1_binding_energies = np.array(GP_negative_tot)
system2_binding_energies = np.array(GA_negative_tot) 
system3_binding_energies = np.array(GAP_GAP_negative_tot)

# Determine the range of binding energies
min_energy = min(min(system1_binding_energies),
                 min(system2_binding_energies),
                 min(system3_binding_energies))
max_energy = max(max(system1_binding_energies),
                 max(system2_binding_energies),
                 max(system3_binding_energies))

# Define the number of bins and calculate the bin width

bin_width = (max_energy - min_energy) / num_bins

# Bin the data and count the occurrences for each system
system1_counts, _ = np.histogram(system1_binding_energies, bins=num_bins, range=(min_energy, max_energy))
system2_counts, _ = np.histogram(system2_binding_energies, bins=num_bins, range=(min_energy, max_energy))
system3_counts, _ = np.histogram(system3_binding_energies, bins=num_bins, range=(min_energy, max_energy))

# Normalize the counts by the total number of states in each system
total_states_system1 = len(system1_binding_energies)
total_states_system2 = len(system2_binding_energies)
total_states_system3 = len(system3_binding_energies)
system1_pdf = system1_counts / total_states_system1
system2_pdf = system2_counts / total_states_system2
system3_pdf = system3_counts / total_states_system3

# Calculate the x-axis values (center of each bin)
x = np.linspace(min_energy + bin_width/2, max_energy - bin_width/2, num_bins)

# Plot the PDFs
plt.plot(x, system1_pdf, label='GP')
#plt.plot(x, system2_pdf, label='GA')
plt.plot(x, system3_pdf, label='GAP')
plt.xlabel('Binding Energy [eV/Å²]')
plt.ylabel('Probability Density')
plt.title('Comparison of Probability Density Functions for total binding Energy VASP')
plt.legend()
# Save the figure with high resolution
plt.savefig('/nas/longleaf/home/jarkeith/calculate_area_of_molecule/GAP_total_prob_density.png', dpi=dpi)
plt.show()

import numpy as np
from itertools import combinations
from scipy.stats import ks_2samp

# Example data for three systems
system1_data = system1_binding_energies
system2_data = system2_binding_energies
system3_data = system3_binding_energies

# Combine the systems for pairwise comparisons
systems = [system1_data, system2_data, system3_data]
system_pairs = combinations(range(len(systems)), 2)

# Perform the Kolmogorov-Smirnov test for each pair
for pair in system_pairs:
    system1_index, system2_index = pair
    system1_data, system2_data = systems[system1_index], systems[system2_index]
    statistic, p_value = ks_2samp(system1_data, system2_data)
    print(f"Comparing System {system1_index + 1} and System {system2_index + 1}")
    print("Test Statistic:", statistic)
    print("P-value:", p_value)
    print()

# ################### ################### ################### ################### ################### #################

import numpy as np

# Combine the arrays column-wise
data = np.column_stack((x, system1_pdf, system3_pdf))

# Save the data to a text file
np.savetxt('data.txt', data, delimiter=' ')

print("File 'data.txt' created successfully.")

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Let's create a function to model and create data
def func(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))  # Gaussian function
y1 = system1_pdf
y2 = system3_pdf
# Plot out the current state of the data and model
fig, ax = plt.subplots()
#ax.plot(x, y1, c='k', label='GP')
#ax.scatter(x, y1)
#ax.plot(x, y2, c='b', label='GAP')
#ax.scatter(x, y2)

# Executing curve_fit on data
popt1, pcov1 = curve_fit(func, x, y1)
popt2, pcov2 = curve_fit(func, x, y2)

# Calculate the best fit lines for system1_pdf and system2_pdf
ym1 = func(x, popt1[0], popt1[1], popt1[2])
ym2 = func(x, popt2[0], popt2[1], popt2[2])

# Plot the best fit lines
ax.plot(x, ym1, c='r', label='GP')
ax.plot(x, ym2, c='g', label='GAP')
ax.legend()

# Set labels and title
ax.set_xlabel('Binding Energy')
ax.set_ylabel('prob density')
#ax.set_title('Best Fit Lines')

# Save the plot
fig.savefig('model_fit.png')

# Show the plot
plt.show()
