# Introduction

miRPara_model_trainer.pl was designed to train the miRPara models based on newly released miRBase data. It supports miRBase 13.0 or higher. Basically, it will: * extract the experimentally verified pri-miRNAs * calculate the parameters for each pri-miRNA * group the parameters into overall, sub-groups and species name * create 20 levels negative parameters for each group * train the models

The selection of miRPara_model_trainer.pl have to match the miRPara.pl version: 
|Model Trainer Version|miRPara Version| 
|:--------------------|:--------------| 
|3.0 or above |6.0 or above | 
|2.3 or below |5.3 or below |

# Requirement

**Version 3.0 or above** 
* File::chdir 
* Cwd 
* threads 
* threads::shared 
* miRPara.pl (any version) 
* libsvm ([http://www.csie.ntu.edu.tw/~cjlin/libsvm/](http://www.csie.ntu.edu.tw/~cjlin/libsvm/))

**Version 2.3 or below** 
* perl 
* File::chdir 
* Algorithm::SVM 
* Algorithm::SVM::DataSet 
* Cwd 
* miRPara.pl

# Usage

**Version 3.0 or above**

The script is written in perl and no need to install.

**_Usage:_** 

```
Usage: perl miRPara_model_trainer_v3.0.pl

Example: perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 1,2,3,4 12 all

Steps(more than one steps can be given at one time, separated by ","): 1 extract experimentally varified miRNA sequences 2 calculate the miRPara parameters 3 generate the parameters for training 4 model training

Species (one of following): - all species - eg. hsa ``` **_Example:_**

1.  Only extract the experimentally verified miRNA from miRBase version 19, the number of cores and the species does not matter in this step, any number and any species can be given: `perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 1 1 all perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 1 2 hsa`

2.  Only calculate the miRPara parameters, the species does not matter in this step, any species can be given: `perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 2 1 all`

3.  Only generate the parameters for training, here 12 cores and hsa parameters was generated: `perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 3 12 hsa`

4.  Only train the models for cel with 12 cores: `perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 4 12 cel`

5.  Completely train the models for the miRBase 18 with 12 cores: `perl miRPara_model_trainer_v3.0.pl 18 /usr/bin/miRPara.pl 1,2,3,4 12 all`

6.  Only generate parameters and train the models for the human from miRBase 19 with 12 cores: `perl miRPara_model_trainer_v3.0.pl 19 /usr/bin/miRPara.pl 3,4 12 hsa`
```
**Version 2.3 or below**

The script is written in perl and no need to install.

**_Usage:_** 

```perl miRPara_model_trainer.pl <miRBase Version> <Full miRPara address>```

**_Example:_** 

```perl miRPara_model_trainer.pl 19 /home/xxx/xxx/miRPara.pl```

# Notes

**Version 3.0 or above** 

The training is a slow process, it take around 72 hours in: Ubuntu 12.04 CPU 2.40GB X 12 Memory: 24 GB

**Version 2.3 or below**

The training is a slow process, it might take as long as one week to finish all trains.

# Release

*   version 3.3 (2012-09-25): 1) added a function to sort parameter file; 2) fixed a bug in reading organisms.txt in mirbase20 and other bugs in creating the positive and negative sequences
*   version 3.2 (2012-11-24 15:32:17 ): fix a bug that terminal the program after step 3
*   version 3.1 (2012-11-23 20:21:12 ): enable multicore processing for the parameter calculation
*   version 3.0 (2012-11-17 09:40:06 ): 1) use libsvm instead of Algorithm::SVM; 2) able to run the whole process separatly for different steps; 3) able to run for one specific species; 4) multiple core support
*   version 2.3 (2012-10-12 18:30:07 ): fix a bug that not able to produce parameters for some species
*   version 2.2 (2012-10-09 12:51:39 ): ignore empty data set while training, fix a bug to include species group
*   version 2.1 (2012-10-03 17:52:47 ): fix the bug to download miRBase17 data
*   version 2.0 (2012-10-02 15:44:02 )
