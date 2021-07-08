# Analysis and modelling of Locomotion Homework 3 : PCA analysis of motion-capture and EMG gait parameters

## Investigated time series

| **Dataset name**         | Species | Condition                             | Context                    |
| ------------------------ | ------- | ------------------------------------- | -------------------------- |
| H01_TDM_2kmh.mat         | Human   | Healthy                               | 2kmh walk                  |
| H01_TDM_35kmh.mat        | Human   | Healthy                               | 3.5kmh walk                |
| H01_TDM_2kmh_20_incl.mat | Human   | Healthy                               | 2kmh walk on a slope (20°) |
| DM002_TDM_08_2kmh.mat    | Human   | Spinal Cord Injury, EES treatement    | 2kmh walk                  |
| DM002_TDM_1kmh_NoEES.mat | Human   | Spinal Cord Injury, no EES treatement | 1kmh walk                  |
| DM002_TDM_08_1kmh.mat    | Human   | Spinal Cord Injury, treatement        | 1khm walk                  |
|                          |         |                                       |                            |
|                          |         |                                       |                            |



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