##Stich Fields for 1-1,2-2,3-3,4-4,5-5 for IRS2 & W51 Main##

#FITS to casa.image

#importfits(fitsimage='W51_NH3_11_freq.fits', imagename='W51_NH3_11_.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='W51_NH3_22_freq.fits', imagename='W51_NH3_22.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='W51_NH3_33_freq.fits', imagename='W51_NH3_33.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='W51_NH3_44_freq.fits', imagename='W51_NH3_44.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='W51_NH3_55_freq.fits', imagename='W51_NH3_55.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)

#importfits(fitsimage='IRS2_NH3_11_freq.fits', imagename='IRS2_NH3_11.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='IRS2_NH3_22_freq.fits', imagename='IRS2_NH3_22.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='IRS2_NH3_33_freq.fits', imagename='IRS2_NH3_33.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='IRS2_NH3_44_freq.fits', imagename='IRS2_NH3_44.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)
#importfits(fitsimage='IRS2_NH3_55_freq.fits', imagename='IRS2_NH3_55.im', defaultaxes=True, defaultaxesvalues=['','','','I'], overwrite=True)

#Mosaic Fields

lm.defineoutputimage(nx=300, outputimage='NH3_11_mosaic.linmos', imagecenter='19h23m42.276 14d30m24.998' )
lm.makemosaic(images=['IRS2_NH3_11.im', 'W51_NH3_11.im'], weightimages=['/lustre/aoc/observers/nm-10487/data/run_pipeline/IRS2_11_cleaned_03.pb', '/lustre/aoc/observers/nm-10487/data/run_pipeline/W51_NH3_11_02.pb'])

lm.defineoutputimage(nx=300, outputimage='NH3_22_mosaic.linmos', imagecenter='19h23m42.276 14d30m24.998' )
lm.makemosaic(images=['IRS2_NH3_22.im', 'W51_NH3_22.im'], weightimages=['/lustre/aoc/observers/nm-10487/data/run_pipeline/IRS2_22_cleaned_03.pb', '/lustre/aoc/observers/nm-10487/data/run_pipeline/W51_NH3_22.pb'])

lm.defineoutputimage(nx=300, outputimage='NH3_33_mosaic.linmos', imagecenter='19h23m42.276 14d30m24.998' )
lm.makemosaic(images=['IRS2_NH3_33.im', 'W51_NH3_33.im'], weightimages=['/lustre/aoc/observers/nm-10487/data/run_pipeline/IRS2_33_cleaned_03.pb', '/lustre/aoc/observers/nm-10487/data/run_pipeline/W51_NH3_33.pb'])

lm.defineoutputimage(nx=300, outputimage='NH3_44_mosaic.linmos', imagecenter='19h23m42.276 14d30m24.998' )
lm.makemosaic(images=['IRS2_NH3_44.im', 'W51_NH3_44.im'], weightimages=['/lustre/aoc/observers/nm-10487/data/run_pipeline/IRS2_44_custom_mask_03.pb', '/lustre/aoc/observers/nm-10487/data/run_pipeline/NH3_44_02.pb'])

lm.defineoutputimage(nx=300, outputimage='NH3_55_mosaic.linmos', imagecenter='19h23m42.276 14d30m24.998' )
lm.makemosaic(images=['IRS2_NH3_55.im', 'W51_NH3_55.im'], weightimages=['/lustre/aoc/observers/nm-10487/data/run_pipeline/IRS2_55_custom_mask_02.pb', '/lustre/aoc/observers/nm-10487/data/run_pipeline/NH3_55_01.pb'])


exportfits(imagename='NH3_11_mosaic.linmos', fitsimage='NH3_11_mosaic.fits', dropstokes=True, overwrite=True)
exportfits(imagename='NH3_22_mosaic.linmos', fitsimage='NH3_22_mosaic.fits', dropstokes=True, overwrite=True)
exportfits(imagename='NH3_33_mosaic.linmos', fitsimage='NH3_33_mosaic.fits', dropstokes=True, overwrite=True)
exportfits(imagename='NH3_44_mosaic.linmos', fitsimage='NH3_44_mosaic.fits', dropstokes=True, overwrite=True)
exportfits(imagename='NH3_55_mosaic.linmos', fitsimage='NH3_55_mosaic.fits', dropstokes=True, overwrite=True)
