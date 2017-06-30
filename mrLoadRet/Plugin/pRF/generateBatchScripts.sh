#!/bin/bash

name=$1
echo $scriptsDir

script1="#!/bin/bash\n#\n#SBATCH --job-name=$name\n#SBATCH --output=out/$name.%j.out\n#SBATCH --error=out/$name.%j.err\n#SBATCH --time=30\n#SBATCH --qos=long\n#SBATCH -p normal\n#SBATCH --nodes=2\n#SBATCH --mem=8G\n#SBATCH -c 1\n#SBATCH --mail-type=END,FAIL # notifications for job done & fail\n#SBATCH --mail-user=akshayj@stanford.edu\nmodule load matlab\nmatlab -nodesktop <<EOF\nsplitfile = '$name.mat';\naddpath(genpath('/share/PI/jlg/mrTools'));\ncd /share/PI/jlg/data/mgldoublebars/s036020170331/\npRFRunSplits(splitfile);\nEOF"

echo $script1 > Splits/Scripts/$name.sbatch
