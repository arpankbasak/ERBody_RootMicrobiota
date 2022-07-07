#!/bin/bash

# Originally by @Rui Guan @Ruben Garrido Otter modified by @Arpan Kumar Basak
# Wrapper to fish out the SPIKE-in from the dataset
# scripts for demultiplexing raw sequencing data for 16S
# exits whenever a function returns 1
# set -e

# # # load functions
# scripts_dir=`dirname $0`

# log() {
#     echo $(date -u)": "$1 >> $logfile
# }

# if ($1=='Bacteria')
# then
# ./Bacteria/scripts/CAS/6_GLM_Statistics_combat.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_family_all.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_family_all.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_family_CAS.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_family_CAS.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_family_CAS13.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_family_CAS11.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_family_CAS11.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_family_CAS13.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_asv_all.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_asv_all.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_asv_CAS.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_asv_CAS.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_asv_CAS11.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_asv_CAS11.R
./Bacteria/scripts/all/GLM/RP_GLM_Statistics_asv_CAS13.R
./Bacteria/scripts/all/GLM/E_GLM_Statistics_asv_CAS13.R

# ./Bacteria/scripts/STREX/6_GLM_Statistics_combat.R
# ./Bacteria/scripts/STREX/7_GLM_Kmeans.R
# ./Bacteria/scripts/STREX/8_GLM_Summary.R

# fi

# if ($1=='Fungi')
# then
# ./Fungi/scripts/CAS/6_GLM_Statistics_combat.R
# ./Fungi/scripts/CAS/6_GLM_Statistics_combat.R
# ./Fungi/scripts/CAS/7_GLM_Kmeans.R
# ./Fungi/scripts/CAS/8_GLM_Summary.R

# ./Fungi/scripts/STREX/6_GLM_Statistics_combat.R
# ./Fungi/scripts/STREX/7_GLM_Kmeans.R
# ./Fungi/scripts/STREX/8_GLM_Summary.R

# fi
