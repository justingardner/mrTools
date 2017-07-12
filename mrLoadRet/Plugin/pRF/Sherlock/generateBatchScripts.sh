#!/bin/bash

saveName=$1
sessionpath=$2
suid=$3
whichSplit=$4

name="${saveName}_split${whichSplit}"
email="$suid@stanford.edu"

script1="#!/bin/bash\n#\n#SBATCH --job-name=split${whichSplit}_${saveName}\n#SBATCH --output=/share/PI/jlg/log/$name.%j.out\n#SBATCH --error=/share/PI/jlg/log/$name.%j.err\n#SBATCH --time=90\n#SBATCH --qos=long\n#SBATCH -p normal\n#SBATCH --nodes=3\n#SBATCH --mem=8G\n#SBATCH -c 1\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=$email\nmodule load matlab\nmatlab -nodesktop <<EOF\naddpath(genpath('/share/PI/jlg/mrTools'));\ncd $sessionpath\npRFRunSplits('$1', $4);\nEOF"

script2="sbatch $sessionpath/Splits/Scripts/$name.sbatch"

echo $script1 > Splits/Scripts/$name.sbatch

echo $script2 >> Splits/Scripts/runAll.sh
