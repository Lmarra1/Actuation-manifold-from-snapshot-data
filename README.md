## Overview

This repository accompanies the article titled "Actuation Manifold from Snapshots Data" by Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.

The paper proposes a data-driven methodology to learn a low-dimensional actuation manifold of controlled flows. The methodology starts by resolving snapshot flow data for a representative ensemble of actuations. Key enablers for the actuation manifold are isometric mapping as an encoder and a combination of a neural network and a $k$-nearest neighbour interpolation as a decoder. This methodology is tested for the fluidic pinball at Re = 30. The proposed methodology yields a five-dimensional manifold describing a wide range of dynamics with small representation error. The manifold is shown to be a key enabler for control-oriented flow estimation.



## Paper Information
- Title: Actuation Manifold from Snapshots Data
- Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro
- Journal: Journal of Fluid Mechanics
- Year: 2024
- DOI: https://doi.org/10.1017/jfm.2024.593

## Key Features
- Encoding using ISOMAP procedure for manifold learning from snapshots data of fluidic pinball with actuation.
- Decoder that uses neural networks and KNN for snapshots reconstruction
- Functions to plot the manifold and reconstruction results.


### Downloading the Dataset
Remember to download the dataset from Zenodo using the following link: https://zenodo.org/doi/10.5281/zenodo.12802191

### Dependencies
Please download and install the following functions and toolboxes:

1. **Toolbox Image**
   - Download from: Gabriel Peyre (2024). Toolbox image (https://www.mathworks.com/matlabcentral/fileexchange/16201-toolbox-image), MATLAB Central File Exchange. Recuperato luglio 23, 2024. Only consider the folder "toolbox_image"


2. **L2_distance Function**
   - Download from: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fww2.mathworks.cn%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2F6ff224c7-d91d-460e-b695-086355a59386%2F0963789e-d6f8-46ce-ad06-ec3a9f5c8a91%2Ffiles%2FFSLib_v6.1_2018%2Flib%2Fdrtoolbox%2Ftechniques%2FL2_distance.m&embed=web. 
(C) Laurens van der Maaten, Delft University of Technology

----

-- **IsoMap Function** --
IsoMap Function has been downloaded from: https://www.mathworks.com/matlabcentral/fileexchange/62449-isomap-d-n_fcn-n_size-options
(J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global geometric framework for nonlinear dimensionality reduction. Science 290 (5500): 2319-2323, 22 December 2000.)




### Functions Overview

1. **Encoding.m**
   - Uses the matrix of snapshots of the fluidic pinball for the given explored actuations and applies ISOMAP for manifold learning.

2. **isomap_decoder.m**
   - Performs a Knn interpolation to transition from manifold coordinates to snapshots as proposed by Farzamnik E, Ianiro A, Discetti S, et al. in "From snapshots to manifolds – a tale of shear flows."

3. **PLOT_encoding**
   - Plots related to the manifold of the test case in question (Figures 1-5 of the paper).

4. **PLOT_SnapReconstruction**
   - Plots the results of the reconstruction using the two sensing methods (forces and velocity) and two reconstruction methods (MLP+KNN interpolation and KNN+KNN regression).

5. **MLP_1 and MLP_2**
   - The folders are included in a NN.zip file
   - Folders containing codes for training the neural networks (TrainNN.py) and making predictions (EvalTest.py) on the test dataset.
   - Required packages to install:
     - numpy; pandas; tensorflow; scikit-learn; matplotlib.

7. **ReadH5.m**
   - Code to read the data from the paper downloadable from Zenodo.

8. **Reconstruction_KNN_KNN.m**
   - Performs the reconstruction of snapshots using two KNN regressions (actuations+sensors -> manifold coordinates and subsequently manifold coordinates -> snapshots).

9. **Reconstruction_MLP_KNN.m**
   - Performs the reconstruction of snapshots using MLP (actuations+sensors -> manifold coordinates) and subsequently KNN interpolation (manifold coordinates -> snapshots).

## Contact Information
- Email: luigi.marra@uc3m.es


## Issues and Feedback
If you encounter any issues or have feedback regarding this code, please open an issue on our GitHub repository. Your insights and suggestions are valuable and appreciated.
