### ----------------------------------- STITCHING -------------------------------------- ###
# For one or more stacks of synchrotron x-rays micro-CT slices

"""
B22_11609_01C_norm_Z0.0mm
B22_11609_01C_norm_Z3.6mm
B22_11609_01C_norm_Z7.2mm
B22_11609_01C_patho_Z0.0mm_corr_phrt
B22_11609_01C_patho_Z3.6mm_corr_phrt
B22_11609_01C_patho_Z7.2mm_corr_phrt
Min and Max of all stacks: max_gval =  0.12637895345687866  , min_gval =  -0.1357395350933075
"""

# ----------------------------------------------------------------------------------------
# Input parameters

# source dir, all stacks need to have separate folder (can contain additional subfolder)
input_dir = r'X:\reconstructions\20240604_InHouse_GS\2035-2G\2035_2G_right' 
#'E:\\Melanoma_may23\\recons\\B21-20177-01C\\B21_20177_01C_2um_dx\\' 

# export dir 
export_dir = r'Z:\giulia.saccomano\20240604_InHouse_GS\2035-2G\2035_2G_right-stitched'
#r'X:\reconstructions\20240604_InHouse_GS\12634_5E\12634_5E\12634-5E-singlecells_stitched' 
#'E:\\Melanoma_may23\\mergedStacks\\B21-20177-01C\\B21_20177_01C_2um_dx\\'

# True: if projections are upside-down (such as OrcaFlash v.4.0), scans have been done in reversed order. 
reverse_stacks = True
# if also the slices are reconstructed in reversed order
reverse_files = False

# True: let it compute min and max in all stacks and then break the program
compute_minmax = True
start_stitching = True

# consider the medium stack (for example if z0.0, z4.6, z9.2, z13.8: take z9.2)
n_slice = 169

# reference slice number in stack (starting with 0)
# consider the last stack (for example if z0.0, z4.6, z9.2, z13.8: take z13.8)
end_slice = 1437 # 1881


# data type
export_data_type = 1 # 0=8bit,1=16bit,2=float

# range of slices to check for optimal overlap

range_slice = n_slice - 50 # generally a 20-30 range
overlap_range = [n_slice - range_slice, n_slice + range_slice] 

if not compute_minmax:
    max_provided = 0.12637895345687866
    min_provided = -0.1357395350933075

# ---------------------------------------------------------------------------------------- #
# Functions

import os
import numpy as np
from skimage import io
#from skimage import data
from skimage.registration import phase_cross_correlation
from skimage.transform import warp_polar, rotate, rescale, SimilarityTransform, warp
#from skimage.util import img_as_float
import warnings
warnings.filterwarnings('ignore')

if not os.path.exists(export_dir):
    os.makedirs(export_dir)

# def write_log(lock, filenames, export_dir, logfilename):    	      
#     lock.acquire()
#     try: 
# 		# Print out execution time:
#         log = open(logfilename,"a")
#         log.write(os.linesep + 'saving ' + filenames[i] + ' -> ' + export_dir + 'slice_%04d' % num + '.tif')
# 		#log.write(os.linesep + "\t%s converted in %0.3f sec." % (os.path.basename(fname), iotime))
#         #('save ',filenames[i],' -> ',export_dir+'slice_%04d' % num+'.tif')
#         log.close()	
#     finally:
#         lock.release()	

def cropped(I1,cropped_rows_and_cols):
    ''' 
    Return a the image cropped from the center
    '''
    rpos = np.max([0,int((I1.shape[0]-cropped_rows_and_cols)/2)])
    cpos = np.max([0,int((I1.shape[1]-cropped_rows_and_cols)/2)])
    I1=I1[rpos:rpos+cropped_rows_and_cols+1,cpos:cpos+cropped_rows_and_cols+1]
    return I1

def corr(I1,I2):
    ''' 
    Return the correlation index as optimization criteria for the image correlation of two stacks 
    Par:
        - im1, reference image
        - im2, 
    '''
    mean1=np.average(I1)
    mean2=np.average(I2)
    I1=I1-mean1
    I2=I2-mean2
    _corr = np.sum(I1*I2)/np.sqrt(np.sum(I1*I1)*np.sum(I2*I2))
    print('correlation = ',mean1,' ',mean2,' ',_corr)
    return _corr

def detect_rotation(I1,I2,display):
    ''' 
    Return the index of correlation and the shift coefficients, detecting the rotation angle between two images using fourier shift in polar coordinates.
    The angle is detected with a precision of 0.01 degree. ensure that the important image content is brighter than the background displaying results is optional 
    also the correlation coefficient is calculated and reported
    '''
    I1_polar = warp_polar(I1, radius=I1.shape[0]*.5)
    I2_polar = warp_polar(I2, radius=I1.shape[0]*.5)

    shifts, error, phasediff = phase_cross_correlation(I1_polar, I2_polar,upsample_factor=100)

    I2_rot = rotate(I2,-shifts[0]);
    Diff = I1-I2_rot;      
    _corr = corr(I1,I2_rot)

    print('corr = ',f'{_corr}',' angle =',f'{shifts[0]}')
    return [_corr, shifts[0]]

def getListOfFiles(dirName):
    ''' 
    Return a list all filenames in a directory including the provided path results may not be sorted
    from having a list of file and sub directories names in the given directory
    '''
    #listOfFile = os.listdir(dirName)
    listOfFile = [file for file in os.listdir(dirName) if file.endswith('.tif')]
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)      
    return allFiles


def export_files(filenames,start_slice,end_slice,angle,dx,dy,num,export_dir,cropped_rows_and_cols):
    ''' 
    Return the number of exported files. 
    It copies all images and applies rotation and cropping a global filenumber is used as counter 
    '''
    for i in range(start_slice,end_slice+1):
        I1 = cropped(np.array(io.imread(filenames[i]),dtype='float'),cropped_rows_and_cols)
        #I1 = np.array(io.imread(filenames[i]),dtype='float')
        if (angle!=0):
            I1=rotate(I1,-angle)
        if dx!=0 or dy!=0:
            tform = SimilarityTransform(translation=(-dx, -dy))
            I1 = warp(I1,tform)
        
        if export_data_type==0:
            I1 = np.asarray((I1-min_gval)/(max_gval-min_gval)*255,dtype='uchar')
        if export_data_type==1:
            I1 = np.asarray((I1-min_gval)/(max_gval-min_gval)*65535,dtype='ushort')
        if export_data_type==2:
            I1 = np.asarray((I1-min_gval)/(max_gval-min_gval),dtype='single')*(max_gval-min_gval)+min_gval

        print('save ',filenames[i],' -> ',export_dir+'/slice_%04d' % num+'.tif')
        #write_log(lock, filenames, export_dir, logfilename)
        io.imsave(export_dir+'/slice_%04d' % num+'.tif',I1)
        num+=1
    return num

def detect_overlap_and_rotation(I1,filenames_moving,overlap_range,cropped_rows_and_cols):
    ''' 
    Return the number and the correlation shift. 
    It detects the rotation and correlation coefficient comparing one fixed image with a range of moving images
    the image number with the highest correlation coefficient and the resulting angle is reported 
    '''
    res=np.zeros((overlap_range[1]-overlap_range[0]+1, 2))
    for i in range(overlap_range[0],overlap_range[1]+1):
        I2 = cropped(np.array(io.imread(filenames_moving[i]),dtype='float'),cropped_rows_and_cols)
        res[i-overlap_range[0],:]=detect_rotation(I1,I2,False)
    pos = np.argmax(res[:,0])
    rot = res[pos,1]
    slice_num = pos+overlap_range[0]    
    print('best guess ',slice_num,' with a rotation of ',rot)   
    
    #validate
    I2 = cropped(np.array(io.imread(filenames_moving[pos]),dtype='float'),cropped_rows_and_cols)
    I2 = rotate(I2,-rot)
    return [slice_num,rot]

def detect_corr(img1,filenames_moving,overlap_range,cropped_rows_and_cols):
    ''' 
    Return a graph of the behaviour of the correlation coefficient through the images.
    It also detects the rotation and correlation coefficient comparing one fixed image with a range of moving images
    the image number with the highest correlation coefficient and the resulting angle is reported
    Par:
        - img
    '''
    res=np.zeros((overlap_range[1]-overlap_range[0]+1, 2))
    for i in range(overlap_range[0],overlap_range[1]+1):
        #print('analysing ',filenames_moving[i])
        img2 = cropped(np.array(io.imread(filenames_moving[i]),dtype='float'),cropped_rows_and_cols)
        #print('mean image n.2 ',np.average(img2))
        res[i-overlap_range[0],:]=corr(img1,img2)
    pos = np.argmax(res[:,0])
    slice_num = pos+overlap_range[0]    
    print('best guess ',slice_num,' corr = ',res[pos])   
    
    #validate
    img2 = cropped(np.array(io.imread(filenames_moving[pos]),dtype='float'),cropped_rows_and_cols)
    return [slice_num,0]

#def showImage(img, dpi):
#    """
#    Return the preview of the image in input
#    Par: 
#        - img , image as an array 2D
#        - dpi , screen resolution as a integer
#    """
#    dim1, dim2 = img.shape
#    fig = plt.figure(frameon = False)
#    fig.set_size_inches(dim1/dpi, dim2/dpi)
#    ax = plt.Axes(fig, [0., 0., 1., 1.])
#    ax.set_axis_off()
#    fig.add_axes(ax)
#    return plt.imshow(img, cmap = 'gray')

# ----------------------------------------------------------------------------------------
#try:
if compute_minmax == True:
    print("Start computing Min and Max . . . ")
    ### --- check histo limits
    # steps
    file_steps = 100
    # limits
    max_gval = -1000000;
    min_gval = 1000000;

    dirnames=np.sort(os.listdir(input_dir))  

    # check all stacks
    for i in range(0,len(dirnames)):
        if '' in dirnames[i]:
            print(dirnames[i])

            filenames=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[i])))
            for j in range(0,len(filenames),file_steps):
                I = np.array(io.imread(filenames[j]),dtype='float')
                max_gval = np.max([max_gval, np.max(I)])
                min_gval = np.min([min_gval, np.min(I)])
else:
    max_gval = max_provided
    min_gval = min_provided

print('Min and Max of all stacks: max_gval = ',max_gval,' , min_gval = ',min_gval)
#except:
#    print('Something went wrong: it could be directories, or missing files')


# Runs only if sticthing procedure i
#try:
if not compute_minmax or (start_stitching and compute_minmax):
    print("Start testing overlap . . .")
    # test overlap
    if reverse_stacks :
        dirnames=os.listdir(input_dir)[::-1]  #np.sort([file for file in os.listdir(input_dir) if file.endswith('.tif')])[::-1]  
    else:
        dirnames=os.listdir(input_dir) #np.sort([file for file in os.listdir(input_dir) if file.endswith('.tif')])  

    if reverse_files:
        filenames_fixed=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[0])))[::-1]
    else:
        filenames_fixed=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[0])))
    
    if reverse_files:
        filenames_float=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[1])))[::-1]
    else:
        filenames_float=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[1])))
    
    I1 = np.array(io.imread(filenames_fixed[end_slice]),dtype='float')
    I2 = np.array(io.imread(filenames_float[overlap_range[0]]),dtype='float')
    I3 = np.array(io.imread(filenames_float[overlap_range[1]]),dtype='float')

    # fig,ax=plt.subplots(figsize=(19,10),ncols=3),
    # ax[0].imshow(I1,cmap='Greys_r')
    # ax[0].set_title('reference')
    # ax[1].imshow(I2,cmap='Greys_r')
    # ax[1].set_title('start overlap region')
    # ax[2].imshow(I3,cmap='Greys_r')
    # ax[2].set_title('end overlap region')
    # plt.show()

    print('reference slice: ', overlap_range[0], ', start overlap region slice: ', overlap_range[0], ', end overlap region slice: ', overlap_range[1])
#except:
#    print('Something went wrong: check if directory paths are ok')
        
try:
    if not compute_minmax or (start_stitching and compute_minmax):
        print('Start stitching procedure . . .')
        # initializing 
        # ----------------------------------------------------------------------------------------
        start_slice = 0
        angle = 0
        num = 0
        if reverse_stacks :
            #dirnames = np.sort([file for file in os.listdir(input_dir) if file.endswith('.tif')])[::-1] 
            dirnames=np.sort(os.listdir(input_dir))[::-1] 
        else:
            dirnames=np.sort(os.listdir(input_dir)) 
            #dirnames = np.sort([file for file in os.listdir(input_dir) if file.endswith('.tif')])
        
        if reverse_files:
            filenames_fixed=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[0])))[::-1]
        else:
            filenames_fixed=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[0])))

        # we crop everything smaller to be on the safe side
        I1 = np.array(io.imread(filenames_fixed[end_slice]),dtype='float')
        cropped_rows_and_cols = I1.shape[0]-10;
        print('cropping to ',cropped_rows_and_cols)

        # we export the first stack without changes (apart cropping)
        num = export_files(filenames_fixed,start_slice,end_slice,angle,0,0,num,export_dir,cropped_rows_and_cols)
        #num = end_slice

        # now we check all stacks
        for i in range(1,len(dirnames)):
            print(dirnames[i-1],'->',dirnames[i])
        
            if reverse_files:
                filenames_moving=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[i])))[::-1]
            else:
                filenames_moving=np.sort(getListOfFiles(os.path.join(input_dir,dirnames[i])))
            
            ##
            #detect overlap and rotation
            I1 = cropped(np.array(io.imread(filenames_fixed[end_slice]),dtype='float'),cropped_rows_and_cols)
            res = detect_corr(I1,filenames_moving,overlap_range,cropped_rows_and_cols)
            start_slice = res[0]+1
            angle = res[1]+angle
            print('resulting angle =',angle)
            if (i==len(dirnames)-1):
                end_slice = len(filenames_moving)-1
            num = export_files(filenames_moving,start_slice,end_slice,angle,0,0,num,export_dir,cropped_rows_and_cols)
        
            #prepare for next round
            filenames_fixed = filenames_moving
        print("File copied successfully.")
    
#except shutil.SameFileError:
#    print("Source and destination represents the same file.")
except PermissionError:
    print("Permission denied.")
# except:
#     print("Error occurred while copying file.")
except:
    if compute_minmax == True:
        print("Breaking the code.")
        exit(1)



