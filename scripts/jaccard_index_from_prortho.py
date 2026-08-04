import pandas as pd
import numpy as np
from itertools import combinations
import os



##before starting you need to modify the proteinortho .tsv file for this script to work

#remove .faa
#change #Species in 'Cluster', in that column change with cluster 1, cluster 2 etc
#remove column 2 and 3
#The file should look like this


### Cluster	sequenceA	sequenceB	sequenceC	sequenceD
### cluster1	IOKCPI_025	*	*	LPIKCG_041
### cluster2	*	*	*	LPIKCG_003
### cluster3	IOKCPI_045	*	*	LPIKCG_040
 

# Change this to your folder path
working_dir = "/path_to_proteinortho_folder/"

os.chdir(working_dir)
print("Working directory set to: {working_dir}")

# Load ProteinOrtho output
po = pd.read_csv("your.proteinortho.file.modified.tsv", sep="\t", index_col=0)


# modify the file again: '*' -> 0, everything else -> 1
pa_matrix = po.replace('*', 0)   # '*' becomes 0
pa_matrix = pa_matrix.map(lambda x: 1 if x != 0 else 0)  # everything else -> 1

# Display first few rows to check
pa_matrix.head()


#traspose table
# sequence names as rows, clusters as columns
pa_matrix = pa_matrix.T
pa_matrix.head()


#define jaccard function
import numpy as np

def jaccard_similarity(vec1, vec2):
    intersection = np.logical_and(vec1, vec2).sum()
    union = np.logical_or(vec1, vec2).sum()
    if union == 0:
        return 0
    else:
        return intersection / union
        
        
#Step 4: Compute pairwise Jaccard index

from itertools import combinations

edges = []

for sequence, sequence2 in combinations(pa_matrix.index, 2):
    jaccard = jaccard_similarity(pa_matrix.loc[sequence1], pa_matrix.loc[sequence2])
    edges.append([sequence1, sequence2, jaccard])

edge_df = pd.DataFrame(edges, columns=['sequence1', 'sequence2', 'Jaccard_index'])
edge_df.head()


#save the file 
edge_df.to_csv('MGEs_jaccard_edges.csv', index=False)

#filter edges for plotting
threshold = 0.3
edge_df_filtered = edge_df[edge_df['Jaccard'] > threshold]
edge_df_filtered.to_csv('MGE_jaccard_edges_filtered0.3.csv', index=False)
