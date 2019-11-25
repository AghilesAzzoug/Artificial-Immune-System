# Artificial Immune Recognition System V2
# AIRS2 Algorithm

# What does it do?
* Classifying irises (Flowers)

# How?
* Using the artificial immune recognition system algorithm v2
    * Initialization
    * Memory cell (MC) identification
    * ARB generation 
    * Competition for resources
    * (Potential) Candidate memory introduction to the set of Cells 
    * K-Nearest Neighbors approach for classification

# Algorithm parameters
    * Hyper Clonal Rate : Define the number of clones an ARB is allowed to produce
    * Clonal Rate : Define the number of ressources an ARB can obtain and used also to determine the number of clones allowed to produce
    * Class Number : The number of classes (3 in this case)
    * Mc Init Rate : Define the number of training data to be copied in memory cells 
    * Total Num Resources : The total number of resources to share between ARBs
    * Affinity Threshold Scalar : Give a cut-off value for cell replacement
    * K : The number of memory cells to use for classification
    * Test Size : The percentage of global data to take as test data
    * Mutation rate : The probability for a feature (in any given antigene) to mutate
    
# Accuracy
Mean accuracy after 1000 executions is : ``` 93,78 %```

In ``` 44.06 seconds```

Using ```Intel i7-6700HQ CPU```

With this parameters :
```
    * Hyper Clonal Rate : 20
    * Clonal Rate : 0.8
    * Mc Init Rate : 0.4 
    * Total Num Resources : 10
    * Affinity Threshold Scalar : 0.8
    * K : 3
    * Test Size : 0.4
    * Mutation rate : 0.2
```

# Output Example
```
Execution time : 0.0389 seconds
Accuracy : 98.33 %
``` 

# Some notes 
The algorithm is kind-of a mix between KNearestNeighbors and Genetic Algorithms. The Memory Cells (MC) that are used for final classification tend better represent the data, by choosing data points that are closer to their true classes (the true distribution).
The picture below shows the separation and the shape of data point using T-SNE representation. The algorithm tends to 
generate new points to better represent the data and remove outliers (points from a certain class located in other ones).

![alt text](https://github.com/AghilesAzzoug/Artificial-Immune-System/blob/master/tsne_img.png)
 

# Requirements

* Python >= 3


