import numpy as np
import nibabel as nib
import scipy.io as sio
import pandas as pd

FREESURFER_DIR = '/Users/gru/data/freesurfer'
SUBJECT = 's0423'


#######
img_lh = nib.load('{}/{}/atlas/lh.Glasser2016.nii'.format(FREESURFER_DIR, SUBJECT))
img_rh = nib.load('{}/{}/atlas/rh.Glasser2016.nii'.format(FREESURFER_DIR, SUBJECT))

df = pd.read_csv('/Users/gru/Downloads/atlasmgz/Glasser2016.txt', sep=' ', names=['index', 'roi'])

sio.savemat('{}/{}/atlas/glasser.mat'.format(FREESURFER_DIR, SUBJECT), 
        {'rh': np.array(img_rh.dataobj), 'lh': np.array(img_lh.dataobj), 'labels': list(df['roi'])})
