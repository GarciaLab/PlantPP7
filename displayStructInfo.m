function displayStructInfo

disp('MeanFluoAll = the mean spot fluorescence across all cells as a function of time. Undetected cells are asigned a 0.')
disp('MeanOn = the mean spot fluorescence across active cells as a function of time. Undetected cells are asigned a NaN.')

disp('InstFractionON = the number of cells with spots at time t  divided by the total number of cells at time t.')

disp('AllParticles = a 2-dimensional array of particles x frames where the entries are the spot fluorescence values.')

disp('IntegralSoFarWithOff = a 2-dimensional array of particles x frames where the entries are the integrated spot fluorescence up to that frame. Inactive transcription spots are asigned a 0.')

disp('MeanAccumulatedFluoOn = mean across particles in the integrated spot fluorescence as a function of time. Inactive cells are not counted.')

disp('MeanAccumulatedFluoAll = mean across particles in the integrated spot fluorescence as a function of time. Inactive cells are counted by asigning them a value of 0.')

disp('CellsPerFrame = the median number of nuclei per frame')
