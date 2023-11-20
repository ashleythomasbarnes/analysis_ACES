import os 

inputfile = '/Users/abarnes/Dropbox/work/Smallprojects/aces/data/alma/12m7mtp_lowres/HNCO_7m12mTP_CubeMosaic_downsample9.fits'

importfits(inputfile, inputfile.replace('.fits', '.img'), overwrite=True)
imrebin(inputfile.replace('.fits', '.img'), inputfile.replace('.fits', '_downsamplespectrally.img'), factor=[1,1,10], overwrite=True)
exportfits(inputfile.replace('.fits', '_downsamplespectrally.fits'), velocity=True, overwrite=True)