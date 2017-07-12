#!/bin/bash

name=$1
sessionpath=$2
suid=$3

email="$suid@stanford.edu"

script1="#!/bin/bash\n#\n#SBATCH --job-name=$name\n#SBATCH --output=/share/PI/jlg/log/$name.%j.out\n#SBATCH --error=/share/PI/jlg/log/$name.%j.err\n#SBATCH --time=30\n#SBATCH --qos=long\n#SBATCH -p normal\n#SBATCH --nodes=2\n#SBATCH --mem=8G\n#SBATCH -c 1\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=$email\nmodule load matlab\nmatlab -nodesktop <<EOF\nsplitfile = '$name.mat';\naddpath(genpath('/share/PI/jlg/mrTools'));\ncd $sessionpath\npRFRunSplits(splitfile);\nEOF"


script2="sbatch $sessionpath/Splits/Scripts/$name.sbatch"

echo $script1 > Splits/Scripts/$name.sbatch

echo $script2 >> Splits/Scripts/runAll.sh
