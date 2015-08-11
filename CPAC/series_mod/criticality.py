#in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')

# THIS SCRIPT USES Nvar * Ntimepoints LIKE MATRIX STRUCTURES

# Point process analysis for a signal. Values equal to 1 when the original value 
# is higher than the threshold (1*SD)
def point_process(signal):

    import numpy as np

    pp_signal = np.zeros(signal.shape[0])
    th = np.std(signal)
    
    # Extreme events over SD = 1
    pp_signal[signal > th] = 1
         
    # Just one marker each time the extreme event happens. Even if it last more timepoints,
    # we take 1.   For example: 00101110000100111011 =>   00101000000100100010  
    prev_value = 0     
    for index, value in enumerate(pp_signal):
        if prev_value == 1 and value == 1:
            pp_signal[index] = 0
        prev_value = value
              
    return pp_signal
    
    
# Conditional Rate Map. Given an fMRI, extract timeseries, calculate Point Process
# and then the Rate Map for each voxel given a seed   
def cond_rm(in_file, seed_location):

    import numpy as np
    import os
    import nibabel as nb
    # Treat fMRI image
    img = nb.load(in_file)
    #print img.shape
    data = img.get_data()
    
    (n_x, n_y, n_z, n_t) = data.shape
    
    K = np.zeros((n_x, n_y, n_z))
    # Extract each voxel
    seed_data = data[seed_location[0], seed_location[1], seed_location[2],:]
    # Extraction of PP signal
    pp_seed_data = point_process(seed_data) 
    # Count how many extreme events happen. This is needed for later calculation of the CRM ratio
    r = np.count_nonzero(pp_seed_data) 
        
    # As we have to compare the extreme events in the seed up to 2 time steps later,
    # we roll the series 2 times, ensuring the 1st value = 0. It could happen that 
    # comparing with the target, 2 extreme events counted as 1 if seed[idx]=extreme 
    # event and seed[idx+2]=extreme event, but it is very unlikely to happen.
    pp_seed_data_1 = np.roll(pp_seed_data, 1) 
    pp_seed_data_1[0] = 0
    pp_seed_data_1 = np.logical_or(pp_seed_data,pp_seed_data_1)
    pp_seed_data_2 = np.roll(pp_seed_data_1, 1) 
    pp_seed_data_2[0] = 0
    pp_seed_data_2 = np.logical_or(pp_seed_data_1,pp_seed_data_2)
    # example: 0100010010001000101001 => 0111011111101110111111
    
    # Calculate each PP signal
    for i_ in range(n_x):
        for j_ in range(n_y):
            for k_ in range(n_z):
                
                target_data = data[i_,j_,k_,:] 
                pp_target_data = point_process(target_data)
                
                # LOGIC AND (target/seed) and count(signal == 1), that will give you the X/r parameter [0,1]
                K[i_,j_,k_] = np.count_nonzero(np.logical_and(pp_seed_data_2,pp_target_data))/float(r)
    
    #create img with K values
    img_new = nb.Nifti1Image(K, header=img.get_header(), affine=img.get_affine())


    # Reconstruct the 3D volume
    cond_rm_file = os.path.join(os.getcwd(), 'cond_rm.nii.gz')
    img_new.to_filename(cond_rm_file)

    return cond_rm_file
  