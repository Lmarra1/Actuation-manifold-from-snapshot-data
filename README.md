## Overview

This repository accompanies the article titled "Actuation Manifold from Snapshots Data" by Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro.

The paper proposes a data-driven methodology to learn a low-dimensional actuation manifold of controlled flows. The methodology starts by resolving snapshot flow data for a representative ensemble of actuations. Key enablers for the actuation manifold are isometric mapping as an encoder and a combination of a neural network and a $k$-nearest neighbour interpolation as a decoder. This methodology is tested for the fluidic pinball at Re = 30. The proposed methodology yields a five-dimensional manifold describing a wide range of dynamics with small representation error. The manifold is shown to be a key enabler for control-oriented flow estimation.



## Paper Information
- Title: Actuation Manifold from Snapshots Data
- Authors: Luigi Marra, Guy Y. Cornejo Maceda, Andrea Meilán-Vila, Vanesa Guerrero, Salma Rashwan, Bernd R. Noack, Stefano Discetti, and Andrea Ianiro
- Journal: Journal of Fluid Mechanics
- Year: 2024
- DOI: To be 

## Key Features

- Encoding using ISOMAP procedure for manifold learning from snapshots data of fluidic pinball with actuation.
- Decoder that uses neural networks and KNN for snapshots reconstruction
- Functions to plot the manifold and reconstruction results.


### Downloading the Dataset

Remember to download the dataset from Zenodo using the following DOI: [Zenodo DOI]

### Dependencies

Please download and install the following functions and toolboxes:

1. **ISOMAP Function**
   - Download from: https://www.mathworks.com/matlabcentral/fileexchange/62449-isomap-d-n_fcn-n_size-options
(J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global geometric framework for nonlinear dimensionality reduction. Science 290 (5500): 2319-2323, 22 December 2000.)

2. **Toolbox Image**
   - Download from: Gabriel Peyre (2024). Toolbox image (https://www.mathworks.com/matlabcentral/fileexchange/16201-toolbox-image), MATLAB Central File Exchange. Recuperato luglio 23, 2024. Only consider the folder "toolbox_image"


3. **L2_distance Function**
   - Download from: https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fww2.mathworks.cn%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2F6ff224c7-d91d-460e-b695-086355a59386%2F0963789e-d6f8-46ce-ad06-ec3a9f5c8a91%2Ffiles%2FFSLib_v6.1_2018%2Flib%2Fdrtoolbox%2Ftechniques%2FL2_distance.m&embed=web. 
(C) Laurens van der Maaten, Delft University of Technology





### Functions Overview

1. **Main Encoding**
   - Uses the matrix of snapshots of the fluidic pinball for the given explored actuations and applies ISOMAP for manifold learning.

2. **ISOMAP Decoder**
   - Performs a Knn interpolation to transition from manifold coordinates to snapshots as proposed by Farzamnik E, Ianiro A, Discetti S, et al. in "From snapshots to manifolds – a tale of shear flows."

3. **Plot Encoding**
   - Plots related to the manifold of the test case in question (Figures 1-5 of the paper).

4. **Plot Snapshots Reconstruction**
   - Plots the results of the reconstruction using the two sensing methods (forces and velocity) and two reconstruction methods (MLP+KNN interpolation and KNN+KNN regression).

5. **MLP_1 and MLP_2**
   - Folders containing codes for training the neural networks (TrainNN.py) and making predictions (EvalTest.py) on the test dataset.
   - Required packages to install:
     - numpy
     - pandas
     - tensorflow
     - scikit-learn
     - matplotlib

6. **ReadH5**
   - Code to read the data from the paper downloadable from Zenodo.

7. **Reconstruction_KNN_KNN**
   - Performs the reconstruction of snapshots using two KNN regressions (actuations+sensors -> manifold coordinates and subsequently manifold coordinates -> snapshots).

8. **Reconstruction_MLP_KNN**
   - Performs the reconstruction of snapshots using MLP (actuations+sensors -> manifold coordinates) and subsequently KNN interpolation (manifold coordinates -> snapshots).

## Contact Information
- Email: luigi.marra@uc3m.es


## Issues and Feedback
If you encounter any issues or have feedback regarding this code, please open an issue on our GitHub repository. Your insights and suggestions are valuable and appreciated.

## Funding
This work is supported by the National Science Foundation of China (NSFC) through grants 12172109 and 12302293, by the Guangdong Basic and Applied Basic Research Foundation under grant 2022A1515011492, and by the Shenzhen Science and Technology Program under grants JCYJ20220531095605012, KJZD20230923115210021 and 29853MKCJ202300205, by the projects PITUFLOW-CM-UC3M and PREDATOR-CM-UC3M, funded by the call ``Estímulo a la Investigación de Jóvenes Doctores/as'' within the frame of the Convenio Plurianual CM-UC3M and the V PRICIT (V Regional Plan for Scientific Research and Technological Innovation), by the project ARTURO, ref. PID2019-109717RB-I00/AEI/10.13039/501100011033, funded by the Spanish State Research Agency, by the project ACCREDITATION (Grant No TED2021-131453B-I00), funded by MCIN/AEI/ 10.13039/501100011033 and by the “European Union NextGenerationEU/PRTR”, by the funding under "Orden 3789/2022, del Vicepresidente, Consejero de Educación y Universidades, por la que se convocan ayudas para la contratación de personal investigador predoctoral en formación para el año 2022", by the European Research Council (ERC) under the European Union’s Horizon H2020 research and innovation programme (grant agreement No 949085) and by the project EXCALIBUR (Grant No PID2022-138314NB-I00), funded by MCIU/AEI/ 10.13039/501100011033 and by“ERDF A way of making Europe”.
