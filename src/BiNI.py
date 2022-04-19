import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from os import listdir
from os.path import join, isfile


# Assert file loaded
error1 = 'error_1: ' + ' <missing one or more input files>'
error2 = 'error_2: ' + sys.argv[1] + ' <not a csv file>'
error3 = 'error_3: ' + sys.argv[2] + ' <not a csv file>'

if len(sys.argv) < 3:
    print(error1)
    exit()
if sys.argv[1].endswith('.csv') == False:
    print(error2)
    exit()
if sys.argv[2].endswith('.csv') == False:
    print(error3)
    exit()


# Load BIGFAM file
BIGFAM_file = pd.read_csv(sys.argv[1])
length_file = pd.read_csv(sys.argv[2])

# Complete Genomes Filter
BIGFAM_file = BIGFAM_file[BIGFAM_file['completeness'] == 'complete']
length_file = length_file.astype({"length": float})

# Dictionary Creation 
strains_BF = set(BIGFAM_file['strain'])
strains_LF = set(length_file['strain'])
strains = strains_BF.intersection(strains_LF)

def get_distance(strain, BIGFAM_file = BIGFAM_file):
  """
    This function returns the total sum of all the corresponding complete genomes strain's distances.

    Args:
        strain (str): This parameters corresponds to the strain we are interested in getting its total distance.
        BIGFAM_file (df): This parameters corresponds to the DataFrame created from the BIGFAM_file.

    Returns:
        float: The return value (sum_distances) is a float number that corresponds to the sumation of every distance linked to the selected strain.

    """
  sum_distances = sum(BIGFAM_file[BIGFAM_file['strain'] == strain]['distance'])
  return sum_distances

def get_number_of_clusters(strain, BIGFAM_file = BIGFAM_file):
  """
    This function returns the number of clusters from a specific strain.

    Args:
        strain (str): This parameters corresponds to the strain we are interested in getting its total distance.
        BIGFAM_file (df): This parameters corresponds to the DataFrame created from the BIGFAM_file.

    Returns:
        int: The return value (number_of_clusters) is a integer number that corresponds to the total number of clusters linked to the selected strain.

  """
  number_of_clusters = BIGFAM_file[BIGFAM_file['strain'] == strain].shape[0]
  return number_of_clusters

def get_distance_mean(sum_distances, number_of_clusters):
  """
    This function returns the mean of the distance from a specific strain based on the sumber of clusters linked to that strain and the total distance.

    Args:
        sum_distances (float): This parameters corresponds to the strain we are interested in getting its total distance.
        number_of_clusters (int): This parameters corresponds to the DataFrame created from the BIGFAM_file.

    Returns:
        float: The return value (distance_mean) is a float number that corresponds to reason between number of clusters and the sumation fo the distances from a specific strain.

  """
  distance_mean = round(sum_distances/number_of_clusters, 4)
  return distance_mean

def get_strain_lenght(strain, length_file = length_file):
  """
    This function returns the length of the assembly from a specific strain.

    Args:
        strain (str): This parameters corresponds to the strain we are interested in getting its total distance.
        length_file (df): This parameters corresponds to the DataFrame created from the length_file.
        
    Returns:
        float: The return value (length) is a float number that corresponds to length of specific strain assembly.
  """
  if len(length_file[length_file['strain'] == strain]['length']) >= 1:
    length = length_file[length_file['strain'] == strain]['length'].values[0]
    return float(length)

def get_BiNI(distance_mean, length):
  """
    This function returns the length of the assembly from a specific strain.

    Args:
        distance_mean (float): This parameters corresponds to the mean of the distance from a specific strain based on the sumber of clusters linked to that strain and the total distance.
        length (float): This parameters corresponds to the length of the assembly from a specific strain.
        
    Returns:
        float: The return value (BiNI) is a float number that corresponds to BiNI score of specific strain set.
  """
  BiNI = round(distance_mean/length, 4)
  return BiNI

def get_DataFrame(strains):
  strains_list = []
  distances_list = []
  clusters_list = []
  lengths_list = []
  distances_mean_list = []
  BiNIs_list = []

  for strain in strains:
    strains_list.append(strain)

    sum_distances = get_distance(strain)
    distances_list.append(sum_distances)

    number_of_clusters = get_number_of_clusters(strain)
    clusters_list.append(number_of_clusters)

    distance_mean = get_distance_mean(sum_distances, number_of_clusters)
    distances_mean_list.append(distance_mean)

    length = get_strain_lenght(strain)
    lengths_list.append(length)

    BiNI = get_BiNI(distance_mean, length)
    BiNIs_list.append(BiNI)

    BiNI_df = pd.DataFrame()
    BiNI_df['strain'] = strains_list
    BiNI_df['sum_distances'] = distances_list
    BiNI_df['#_Clusters'] = clusters_list
    BiNI_df['Length'] = lengths_list
    BiNI_df['Mean_distance'] = distances_mean_list
    BiNI_df['BiNI_Complete_BCGS'] = BiNIs_list

  return BiNI_df

def plot_BiNI(BiNI_df, file_name):
  with plt.style.context('seaborn-paper'):
    x = BiNI_df['strain']
    y = BiNI_df['BiNI_Complete_BCGS']/max(BiNI_df['BiNI_Complete_BCGS'])
    BiNI_mean = y.mean()
    colors = y
    plt.plot(x, [BiNI_mean]*len(x), label='BiNI score mean', linestyle='--', zorder=1, linewidth=0.5, color='g')
    plt.scatter(x, y, sizes = y*100, c=colors ,alpha=0.5, cmap='RdPu', zorder=2)
    plt.xticks(rotation = 90)
    plt.xlabel('Strains', fontsize = 10, fontweight='bold')
    plt.ylabel('BiNI Complete BGC', fontsize = 10, fontweight='bold')
    plt.colorbar()
    plt.legend()
    plt.tight_layout()
    plt.savefig(file_name + '_BiNI_fig.png', dpi=300)
    #plt.show()

def main():
  file_name  = input('Provide output file name: ')

  directory = os.getcwd()
  
  folders = [x[0] for x in os.walk(directory)]

  if file_name in directory:
    file_name  = input('name already exist, please provide another output file name: ')

  path = os.path.join(directory, file_name)
  
  os.mkdir(path)

  df = get_DataFrame(strains)
  df.to_csv(path + '/' + file_name + '_BiNI_table_parameters.csv', index=False)
  plot_BiNI(df, path + '/' + file_name)

if __name__ == "__main__":
  main()