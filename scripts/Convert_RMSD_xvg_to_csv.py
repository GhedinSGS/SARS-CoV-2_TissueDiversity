"""
Written by: Rita Afriyie Baoteng

# input: xvg file with root mean square deviations of WT and mutant simulations

#To run: python3 convert_RMSD_xvg_to_csv.py

#This code extract the frames in all .xvg to csv files (Root Mean Square Deviation)
"""

import os
xvg=os.listdir('https://github.com/GhedinSGS/SARS-CoV-2_TissueDiversity/tree/main/data/RMSD_xvg')
for i in xvg:
    if ".xvg" in i:
        xvg_name=i[:-4]
        os.system("cat "+i+" | tail -n+19 | awk '{print $2}' > "+xvg_name+".csv")
    else:
        pass
