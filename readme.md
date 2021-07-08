# Analysis and modelling of Locomotion Homework 3 : PCA analysis of motion-capture and EMG gait parameters

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



## Data processing pipeline

*List of parameters obtained from EMG and mo-cap :*

| parameter                               | implemented                           |
| --------------------------------------- | ------------------------------------- |
| average cycle time                      | computed using the get_cycle function |
| variance cycle time                     | computed using the get_cycle function |
| walking velocity                        | fixed for each dataset                |
| stance duration in cycle proportion     | computed using swing_stance           |
| step height normalized to stance height |                                       |
| stance width                            |                                       |

*Motion capture data from the H01_TDM_2kmh timeseries :*
![mocap_cloud_test](./figures/mocap_cloud_test.jpg)