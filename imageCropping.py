import tifffile as tiff
import os
import warnings
warnings.filterwarnings('ignore')


def imageCropping(indir, outdir, imdata):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    BX, BY, Width, Height = imdata
    imlist = os.listdir(indir)
    fin_im = imlist[-1]
    print('Cropping dataset: ', indir)
    for i in imlist: #[2798:]:
        im = tiff.imread(os.path.join(indir, i))
        im = im[BY : (BY+Height), BX : (BX+Width)]
        tiff.imsave(os.path.join(outdir, i), im)
        #print('saving cropped ', i)
        if i == fin_im:
            print('Spavada finita per ', indir)
            break

# indir = 'Z:\\giulia.saccomano\\stitched-2\\B023.0-012823-6E\\'
# outdir = 'E:\\LungTissue_nov23\\B023.0-012823-6E-crop\\'
# imdata = [0, 1632, 3872, 1458]
# imageCropping(indir, outdir, imdata)

# indir = 'Z:\\giulia.saccomano\\stitched-2\\B023.0-013033-A\\'
# outdir = 'E:\\LungTissue_nov23\\B023.0-013033-A-crop\\'
# imdata = [0, 1398, 3887, 1248]
# imageCropping(indir, outdir, imdata)

indir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12634_5G\12634_5G-right-stitched'
outdir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12634_5G\12634_5G-right-cropped'
imdata = [0, 1398, 3887, 1158]
imageCropping(indir, outdir, imdata)

indir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12634_5E\12634_5E\12634-5E-stitched'
outdir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12634_5E\12634_5E\12634-5E-cropped'
imdata = [0, 1134, 3889, 954]
imageCropping(indir, outdir, imdata)

indir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12941_1B\12941_1B-stitched'
outdir = r'Z:\giulia.saccomano\20240604_InHouse_GS\12941_1B\12941_1B-cropped'
imdata = [0, 1548, 3882, 1098]
imageCropping(indir, outdir, imdata)
