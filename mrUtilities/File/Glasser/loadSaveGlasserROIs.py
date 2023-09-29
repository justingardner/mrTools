import numpy as np
import nibabel as nib
import scipy.io as sio
import pandas as pd
import sys

FREESURFER_DIR = '/Users/gru/data/freesurfer'
if len(sys.argv) < 2:
    SUBJECT = 's0424'
else:
    SUBJECT = sys.argv[1]

print(f'Saving glasser ROIs as .mat files for subject {SUBJECT}')
#######
img_lh = nib.load('{}/{}/atlas/lh.Glasser2016.nii'.format(FREESURFER_DIR, SUBJECT))
img_rh = nib.load('{}/{}/atlas/rh.Glasser2016.nii'.format(FREESURFER_DIR, SUBJECT))

df = pd.read_csv('/Users/gru/proj/mrTools/mrUtilities/File/Glasser/glasser_atlas/Glasser2016.txt', sep=' ', names=['index', 'roi'])

sio.savemat('{}/{}/atlas/glasser.mat'.format(FREESURFER_DIR, SUBJECT), 
        {'rh': np.array(img_rh.dataobj), 'lh': np.array(img_lh.dataobj), 'labels': list(df['roi'])})
