tclean(vis='IRS2_11.ms.contsub',
    imagename='IRS2_11',
    spw='0',
    field='0',
    specmode='cube',
    start='0.0km/s',
    width='0.1km/s',
    restfreq='23.69447GHz',
    deconvolver='multiscale',
    scales=[0,3,9],
    pbmask=0.2,
    threshold='5mJy', #Determined from lowest 3-sigma RMS value
    imsize=[200,200],
    restoringbeam='common',
    niter=int(1e6),
    #AUTOMASK PARAMETERS BELOW
    usemask='auto-multithresh',
    noisethreshold=4.0,
    sidelobethreshold=1.25, 
    minbeamfrac=0.3, 
    lownoisethreshold=2.0, 
    negativethreshold=3.0) 

tclean(vis='IRS2_22.ms.contsub',
    imagename='IRS2_22',
    spw='0',
    field='0',
    specmode='cube',
    start='0.0km/s',
    width='0.2km/s',
    restfreq='23.72260GHz',
    deconvolver='multiscale',
    scales=[0,3,9],
    pbmask=0.2, 
    threshold='5mJy', #Determined from lowest 3-sigma RMS value
    imsize=[200,200],
    restoringbeam='common',
    niter=int(1e6),
    #AUTOMASK PARAMETERS BELOW
    usemask='auto-multithresh',
    noisethreshold=3.0,
    sidelobethreshold=1.25, 
    minbeamfrac=0.3,
    lownoisethreshold=2.0, 
    negativethreshold=3.0)

tclean(vis='IRS2_33.ms.contsub',
    imagename='IRS2_33',
    spw='0', 
    field='0',
    specmode='cube',
    start='0.0km/s',
    width='0.1km/s',
    restfreq='23.87008GHz',
    deconvolver='multiscale',
    scales=[0,3,9],
    pbmask=0.2,
    threshold='6mJy', #Determined from lowest 3-sigma RMS value
    imsize=[200,200],
    restoringbeam='common',
    niter=int(1e6),
    #AUTOMASK PARAMETERS BELOW
    usemask='auto-multithresh',
    noisethreshold=4.0,
    sidelobethreshold=1.25, 
    minbeamfrac=0.3, 
    lownoisethreshold=2.0, 
    negativethreshold=3.0) 

tclean(vis='IRS2_44.ms.contsub',
    imagename='IRS2_44',
    spw='0',
    field='0',
    specmode='cube',
    start='0.0km/s',
    width='0.1km/s',
    restfreq='24.13935GHz',
    deconvolver='multiscale',
    scales=[0,3,9],
    pbmask=0.2,
    threshold='3.9mJy', #Determined from lowest 3-sigma RMS value
    imsize=[200,200],
    restoringbeam='common',
    niter=int(1e6),
    #AUTOMASK PARAMETERS BELOW
    usemask='auto-multithresh',
    noisethreshold=4.5,
    sidelobethreshold=1.25, 
    minbeamfrac=0.3, 
    lownoisethreshold=2.0, 
    negativethreshold=5.0) 

tclean(vis='IRS2_55.ms.contsub',
    imagename='IRS2_55_cleaned_02',
    spw='0',
    field='0',
    specmode='cube',
    start='0.0km/s',
    width='0.1km/s',
    restfreq='24.53292GHz',
    deconvolver='multiscale',
    scales=[0,3,9],
    pbmask=0.2,
    threshold='3mJy', #Determined from lowest 3-sigma RMS value
    imsize=[200,200],
    restoringbeam='common',
    niter=int(1e6),
    #AUTOMASK PARAMETERS BELOW
    usemask='auto-multithresh',
    noisethreshold=4.5,
    sidelobethreshold=1.25,
    minbeamfrac=0.3, 
    lownoisethreshold=2.0, 
    negativethreshold=5.0)

exit