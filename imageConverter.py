### ----------------------------------- ONE PATH IMAGE CONVERTER TO 32,16,8 BIT -------------------------------------- ###
# For one SR-micro-CT stack a time
# ---------------------------------------------------------------------------------------------------------------------- #

import os
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from skimage import io

# ---------------------------------------------------------------------------------------------------------------------- #
# Input parameters
# ---------------------------------------------------------------------------------------------------------------------- #

# input directory/ies
indir = r"G:\Giulia\skin_tissues_agar\neo_corr_phrt"
#"E:\\Melanoma_may23\\B21_21086_01B_2um_horiz\\"
indir1 = r'G:\Giulia\skin_tissues_agar\no_neo_corr_phrt'
# indir1 = 'E:\\Extra\\newTanninFoams\\recons32bit\\sulfuric_acid_2_Z0.0mm_corr_phrt\\'
# indir2 = 'E:\\Extra\\newTanninFoams\\recons32bit\\sulfuric_acid_2_Z4.5mm_corr_phrt\\'
# indir3 = 'E:\\Extra\\newTanninFoams\\recons32bit\\sulfuric_acid_3_Z0.0mm_corr_phrt\\'

# data type
export_data_type = 0 # 0=8bit,1=16bit,2=float

# single stack
single_stack = False # False
check_limits = False

# ---------------------------------------------------------------------------------------------------------------------- #
# Main script
# ---------------------------------------------------------------------------------------------------------------------- #

def create_input_dir(input_dir):
    """
    Return the path of the input directory
    Par:
        - input_dir: a generic path to a folder with images inside
    """
    print(input_dir)
    return os.path.join(input_dir,'slices\\')

def create_export_dir(input_dir,export_data_type):
    """
    Return the path of the output directory while creating the folder
    Par:
        - input_dir: a generic path to a folder with images inside
        - export_data_type: an integer value which stands for the image conversion type (0:8bit, 1:16bit, 2:32bit)
    """
    if export_data_type == 0:
        export_dir = os.path.join(input_dir,'slices8bit\\')
    elif export_data_type == 1:
        export_dir = os.path.join(input_dir,'slices16bit\\')
    else:
        export_dir = os.path.join(input_dir,'slices32bit\\')
    if not os.path.exists(export_dir):
        os.makedirs(export_dir)
    return export_dir

def get_filelist_in_dir(input_dir):
    """
    Return a list of files with the name of all slices inside the input_dir parameter while sorting in numerical order
    Par:
        - input_dir: a generic path to a folder with images inside
    """
    return np.sort(os.listdir(input_dir))

def compute_histogram_limits(input_dir):
    """
    Return the miminum and the maximum histogram values of a number of slices decimated by a file_steps,
    Return the list of files in the directory from get_filelist_in_dir function,
    Return the composed input path from create_input_dir function
    Par:
        - input_dir: a generic path to a folder with images inside
    """
    max_gval = -1000000
    min_gval = 1000000
    indir_slices = create_input_dir(input_dir)
    dirnames = get_filelist_in_dir(indir_slices)
    file_steps = int(np.floor(len(dirnames)/20))
    # check all files in the stack decimated by a file step starting from first files
    for i in range(0,len(dirnames), file_steps):
        if '' in dirnames[i]:
            full_path = os.path.join(indir_slices, dirnames[i])
            im = np.array(io.imread(full_path),dtype='float')
            max_gval = np.max([max_gval, np.max(im)])
            min_gval = np.min([min_gval, np.min(im)])
    # check all files in the stack decimated by a file step starting from last files
    for i in range(len(dirnames)-1,0,-file_steps):
        if '' in dirnames[i]:
            full_path = os.path.join(indir_slices, dirnames[i])
            im = np.array(io.imread(full_path),dtype='float')
            max_gval = np.max([max_gval, np.max(im)])
            min_gval = np.min([min_gval, np.min(im)])
    print('done: max_gval = ',max_gval,' , min_gval = ',min_gval)
    return max_gval, min_gval, dirnames, indir_slices

def export_converted_images(input_dir,export_data_type,check_limits=False):
    """
    Return THE END when it lasts running
    Par:
        - input_dir: a generic path to a folder with images inside
        - export_data_type: an integer value which stands for the image conversion type (0:8bit, 1:16bit, 2:32bit)
    """
    outdir = create_export_dir(input_dir,export_data_type)
    max_hist, min_hist, dirnames, indir_slices = compute_histogram_limits(input_dir)
    count = 0
    for entry in dirnames:
        if check_limits:
            if count == 0:
                break
        im = np.array(io.imread(os.path.join(indir_slices,entry)),dtype='float')
        if export_data_type==0:
            im = np.asarray((im-min_hist)/(max_hist-min_hist)*255,dtype='ubyte')
        if export_data_type==1:
            im = np.asarray((im-min_hist)/(max_hist-min_hist)*65535,dtype='ushort')
        if export_data_type==2:
            im = np.asarray((im-min_hist)/(max_hist-min_hist),dtype='single')*(max_hist-min_hist)+min_hist
        print('saving converted ',entry)
        io.imsave(outdir +'/'+entry, im) #'/slice_%04d' % num+'.tif',I1)
        count = count + 1
        if count == len(dirnames):
            print('THE END')
            #os.rmdir(indir_slices) # to be tested
    return max_hist, min_hist

def export_converted_images2(input_dir,export_data_type,max_hist,min_hist):
    """
    Return THE END when it lasts running
    Par:
        - input_dir: a generic path to a folder with images inside
        - export_data_type: an integer value which stands for the image conversion type (0:8bit, 1:16bit, 2:32bit)
        - max_hist: maximum histogram limit
        - min_hist: minimum histogram limit
    """
    outdir = create_export_dir(input_dir,export_data_type)
    indir_slices = create_input_dir(input_dir)
    dirnames = get_filelist_in_dir(indir_slices)
    count = 0
    for entry in dirnames:
        im = np.array(io.imread(os.path.join(indir_slices,entry)),dtype='float')
        if export_data_type==0:
            im = np.asarray((im-min_hist)/(max_hist-min_hist)*255,dtype='ubyte')
        if export_data_type==1:
            im = np.asarray((im-min_hist)/(max_hist-min_hist)*65535,dtype='ushort')
        if export_data_type==2:
            im = np.asarray((im-min_hist)/(max_hist-min_hist),dtype='single')*(max_hist-min_hist)+min_hist
        print('saving converted ',entry)
        io.imsave(outdir +'/'+entry, im) #'/slice_%04d' % num+'.tif',I1)
        count = count + 1
        if count == len(dirnames):
            print('THE END')
            #os.rmdir(indir_slices) # to be tested
    return max_hist, min_hist
    
if single_stack == True:
    export_converted_images(indir,export_data_type,check_limits)
if single_stack == False:
    max_hist1,min_hist1 = export_converted_images(indir,export_data_type,check_limits)
    if indir1:
        export_converted_images2(indir1,export_data_type,max_hist1,min_hist1)
    if indir2:
        export_converted_images2(indir2,export_data_type,max_hist1,min_hist1)
    if indir3:
        export_converted_images2(indir3,export_data_type,max_hist1,min_hist1)

