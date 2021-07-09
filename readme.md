# Analysis and modelling of Locomotion Homework 3 : PCA analysis of motion-capture and EMG gait parameters

## Data processing pipeline

To run the data processing pipeline proceed as follow : 

1. The submitted code does not contain the datasets, place the src file in the same folder as the dataset in order to correctly import the data
2. There are three different scripts for the 3 different species, each script can be made to process all of the datasets from a given species. To select a dataset, set the *name* variable to the dataset to analyze, the list of available names is presented as a comment in the header of each script. The script can be made to plot results by setting the show_plot boolean to **true**. The resulting feature vector is automatically saved in the **features** folder and can in turn be loaded by the PCA script. The scripts files for the different species are the following :
   1. *Human_Feature_Extraction.m* processes the human datasets
   2. *NHP_Feature_Extraction.m* processes the non-human primates datasets
   3. *Rodent_Feature_Extraction.m* processes the rodent datasets
3. To run the PCA analysis run the *pca_analysis.m* script, it will automatically load the results from the feature extraction scripts and present biplot and eigenvector plots of the relevant datasets. For the PCA to run correctly, all feature vectors must have been computed, else it will raise an error.

## Investigated time series

| **Dataset name**                  | Species           | Condition                               | Context                    |
| --------------------------------- | ----------------- | --------------------------------------- | -------------------------- |
| H01_TDM_2kmh                      | Human             | Healthy                                 | 2kmh walk                  |
| H01_TDM_35kmh                     | Human             | Healthy                                 | 3.5kmh walk                |
| H01_TDM_2kmh_20_incl              | Human             | Healthy                                 | 2kmh walk on a slope (20°) |
| DM002_TDM_08_2kmh                 | Human             | Spinal Cord Injury, EES treatement      | 2kmh walk                  |
| DM002_TDM_1kmh_NoEES              | Human             | Spinal Cord Injury, no EES treatement   | 1kmh walk                  |
| DM002_TDM_08_1kmh                 | Human             | Spinal Cord Injury, EES treatement      | 1khm walk                  |
| Elektra_20190425_TM20_004         | Non-Human Primate | Healthy                                 | 2kmh walk                  |
| Elektra_20190425_TM30_002         | Non-Human Primate | Healthy                                 | 3kmh walk                  |
| Elektra_20190425_TM40_005         | Non-Human Primate | Healthy                                 | 4kmh walk                  |
| Healthy_330_BIP_RW_06             | Rodent            | Healthy                                 | bipedal locomotion         |
| Healthy_330_BIP_RW_07             | Rodent            | Healthy                                 | bipedal locomotion         |
| Healthy_332_BIP_RW_05             | Rodent            | Healthy                                 | bipedal locomotion         |
| Healthy_332_BIP_RW_11             | Rodent            | Healthy                                 | bipedal locomotion         |
| SCI_Trained_207_RW_STIM_25_04     | Rodent            | Spinal Cord Injury, Trained with EES    | bipedal locomotion         |
| SCI_Trained_207_RW_STIM_25_07     | Rodent            | Spinal Cord Injury, Trained with EES    | bipedal locomotion         |
| SCI_Trained_207_RW_STIM_35_02     | Rodent            | Spinal Cord Injury, Trained with EES    | bipedal locomotion         |
| SCI_Trained_207_RW_STIM_BWS40_10  | Rodent            | Spinal Cord Injury, Trained with EES    | bipedal locomotion         |
| SCI_trained_207_RW_SPONT_30_08    | Rodent            | Spinal Cord Injury, Trained without EES | bipedal locomotion         |
| SCI_trained_207_RW_SPONT_30_10    | Rodent            | Spinal Cord Injury, Trained without EES | bipedal locomotion         |
| SCI_trained_207_RW_SPONT_BWS40_03 | Rodent            | Spinal Cord Injury, Trained without EES | bipedal locomotion         |
| SCI_trained_207_RW_SPONT_BWS45_05 | Rodent            | Spinal Cord Injury, Trained without EES | bipedal locomotion         |