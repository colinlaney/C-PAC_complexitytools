# THIS SCRIPT USES Nvar * Ntimepoints LIKE MATRIX STRUCTURES


def point_process(signal):
    """
    Point process analysis for a signal. Values equal to 1 when the original value 
    is higher than the threshold (1*SD)  
    
    Parameters
    ----------

    signal :  timeseries
    
    Returns
    -------

    pp_signal  :  Point Processed signal. No repetitions if the signal goes behind
    the threshold (activates), as explained in http://journal.frontiersin.org/article/10.3389/fphys.2012.00015/abstract
    """

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
    
    

def cond_rm(in_file, seed_location):
    """
    Conditional Rate Map. Given an fMRI, extract timeseries, calculate Point Process
    and then the Rate Map for each voxel given a seed   
    
    Parameters
    ----------

    in_file : 4D Nifti file
    seed_location  :  voxel to analyze as seed

    Returns
    -------

    cond_rm_img  :  3D Volume
    """


    import numpy as np
    import os
    import nibabel as nb
    from CPAC.criticallity import point_process  
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
    cond_rm_img = os.path.join(os.getcwd(), 'cond_rm.nii.gz')
    img_new.to_filename(cond_rm_img)

    return cond_rm_img
 
 
# 
def cluster_detection(in_file): 
    """
    Detects clusters after Point Processing a Brain 
    as described in http://journal.frontiersin.org/article/10.3389/fphys.2012.00015/abstract
    
    Parameters
    ----------

    in_file : 4D Nifti file

    Returns
    -------

    cluster_graph_img  :  4D Nifti file with an id for each cluster in each timestep
    """

    import numpy as np
    import os
    import nibabel as nb
    from CPAC.criticallity import point_process  

    # Treat fMRI image
    img = nb.load(in_file)
    data = img.get_data()
    
    (n_x, n_y, n_z, n_t) = data.shape
    
    # Get the PP data
    pp_data = np.zeros((n_x, n_y, n_z, n_t))
    for i_ in range(n_x):
        for j_ in range(n_y):
            for k_ in range(n_z):
                voxel_data = data[i_,j_,k_,:] 
                pp_data[i_,j_,k_,:] = point_process(voxel_data)
    
    cluster_graph_data_total = np.zeros((n_x, n_y, n_z, n_t))            
    for t_ in range(n_t):
        time_slice = pp_data[:,:,:,t_]
        cluster_graph_data = np.zeros((n_x, n_y, n_z))  
        cluster_number = 1
        
        for i_ in range(n_x):
            for j_ in range(n_y):
                for k_ in range(n_z):
                    
                    if time_slice[i_,j_,k_] == 1: # is active, check if it has active neighboours
                        if time_slice[i_-1,j_,k_] or time_slice[i_+1,j_,k_] \
                        or time_slice[i_,j_-1,k_] or time_slice[i_,j_+1,k_] \
                        or time_slice[i_,j_,k_-1] or time_slice[i_,j_,k_+1]:
                            
                            if cluster_graph_data[i_,j_,k_] == 0: # if is not in any previous cluster
                                this_cluster = (cluster_graph_data[i_-1,j_,k_] or cluster_graph_data[i_+1,j_,k_] \
                                or cluster_graph_data[i_,j_-1,k_] or cluster_graph_data[i_,j_+1,k_] \
                                or cluster_graph_data[i_,j_,k_-1] or cluster_graph_data[i_,j_,k_+1])
                                
                                if this_cluster == 0: #no neighbours in any previous cluster neither
                                    this_cluster = cluster_number
                                    cluster_graph_data[i_,j_,k_] = this_cluster
                                    cluster_number = cluster_number + 1
                                else: 
                                    #check cluster union
                                    merge_clusters = np.unique([cluster_graph_data[i_-1,j_,k_], cluster_graph_data[i_+1,j_,k_] \
                                , cluster_graph_data[i_,j_-1,k_], cluster_graph_data[i_,j_+1,k_] \
                                , cluster_graph_data[i_,j_,k_-1], cluster_graph_data[i_,j_,k_+1]])
                                    merge_clusters = merge_clusters[1:] #quit first value = 0
                                    
                                    this_cluster = merge_clusters[0]
                                    cluster_graph_data[i_,j_,k_] = this_cluster
                                    for cluster_to_merge in merge_clusters[1:]:
                                        cluster_graph_data[cluster_graph_data == cluster_to_merge] = this_cluster
                                    
                                    
                            else:
                                this_cluster = cluster_graph_data[i_,j_,k_]
                                
                            #find neighbours and give cluster_number
                            if time_slice[i_-1,j_,k_] == 1:
                                cluster_graph_data[i_-1,j_,k_] = this_cluster
                            elif time_slice[i_+1,j_,k_] == 1:
                                cluster_graph_data[i_+1,j_,k_] = this_cluster
                            elif time_slice[i_,j_-1,k_] == 1:
                                cluster_graph_data[i_,j_-1,k_] = this_cluster
                            elif time_slice[i_,j_+1,k_] == 1:
                                cluster_graph_data[i_,j_+1,k_] = this_cluster
                            elif time_slice[i_,j_,k_-1] == 1:
                                cluster_graph_data[i_,j_,k_-1] = this_cluster
                            elif time_slice[i_,j_,k_+1] == 1:
                                cluster_graph_data[i_,j_,k_+1] = this_cluster                                
                          
                                #find neighbours and give this_cluster
                                
                    # if not == 1¡, keep the search 
                        # if not neighbours, keep the search
                       
        cluster_graph_data_total[:,:,:,t_] = cluster_graph_data 
        
        img_new = nb.Nifti1Image(cluster_graph_data_total, header=img.get_header(), affine=img.get_affine())
        # Reconstruct the 4D volume
        cluster_graph_img = os.path.join(os.getcwd(), 'cluster_1N.nii.gz')
        img_new.to_filename(cluster_graph_img)
                
    return cluster_graph_img    
    
# Detects clusters after Point Processing a Brain 
# having 2 neighbouring distance in cross-sense
def cluster_detection_mod2(in_file):  

    import numpy as np
    import os
    import nibabel as nb
    from CPAC.criticallity import point_process  

    # Treat fMRI image
    img = nb.load(in_file)
    data = img.get_data()
    
    (n_x, n_y, n_z, n_t) = data.shape
    
    # Get the PP data
    pp_data = np.zeros((n_x, n_y, n_z, n_t))
    for i_ in range(n_x):
        for j_ in range(n_y):
            for k_ in range(n_z):
                voxel_data = data[i_,j_,k_,:] 
                pp_data[i_,j_,k_,:] = point_process(voxel_data)
    
    cluster_graph_data_total = np.zeros((n_x, n_y, n_z, n_t))            
    for t_ in range(n_t):
        time_slice = pp_data[:,:,:,t_]
        cluster_graph_data = np.zeros((n_x, n_y, n_z))  
        cluster_number = 1
        
        for i_ in range(n_x):
            for j_ in range(n_y):
                for k_ in range(n_z):
                    
                    if time_slice[i_,j_,k_] == 1: # is active, check if it has active neighboours
                        if time_slice[i_-1,j_,k_] or time_slice[i_+1,j_,k_] \
                        or time_slice[i_,j_-1,k_] or time_slice[i_,j_+1,k_] \
                        or time_slice[i_,j_,k_-1] or time_slice[i_,j_,k_+1] \
                        or time_slice[i_-2,j_,k_] or time_slice[i_+2,j_,k_] \
                        or time_slice[i_,j_-2,k_] or time_slice[i_,j_+2,k_] \
                        or time_slice[i_,j_,k_-2] or time_slice[i_,j_,k_+2]:
                            
                            if cluster_graph_data[i_,j_,k_] == 0: # if is not in any previous cluster
                                this_cluster = (cluster_graph_data[i_-1,j_,k_] or cluster_graph_data[i_+1,j_,k_] \
                                or cluster_graph_data[i_,j_-1,k_] or cluster_graph_data[i_,j_+1,k_] \
                                or cluster_graph_data[i_,j_,k_-1] or cluster_graph_data[i_,j_,k_+1] \
                                or cluster_graph_data[i_-2,j_,k_] or cluster_graph_data[i_+2,j_,k_] \
                                or cluster_graph_data[i_,j_-2,k_] or cluster_graph_data[i_,j_+2,k_] \
                                or cluster_graph_data[i_,j_,k_-2] or cluster_graph_data[i_,j_,k_+2])
                                
                                if this_cluster == 0: #no neighbours in any previous cluster neither
                                    this_cluster = cluster_number
                                    cluster_graph_data[i_,j_,k_] = this_cluster
                                    cluster_number = cluster_number + 1
                                else: 
                                    #check cluster union
                                    merge_clusters = np.unique([cluster_graph_data[i_-1,j_,k_], cluster_graph_data[i_+1,j_,k_] \
                                , cluster_graph_data[i_,j_-1,k_], cluster_graph_data[i_,j_+1,k_] \
                                , cluster_graph_data[i_,j_,k_-1], cluster_graph_data[i_,j_,k_+1] \
                                , cluster_graph_data[i_-2,j_,k_] , cluster_graph_data[i_+2,j_,k_] \
                                , cluster_graph_data[i_,j_-2,k_] , cluster_graph_data[i_,j_+2,k_] \
                                , cluster_graph_data[i_,j_,k_-2] , cluster_graph_data[i_,j_,k_+2]])
                                    merge_clusters = merge_clusters[1:] #quit first value = 0
                                    
                                    this_cluster = merge_clusters[0]
                                    cluster_graph_data[i_,j_,k_] = this_cluster
                                    for cluster_to_merge in merge_clusters[1:]:
                                        cluster_graph_data[cluster_graph_data == cluster_to_merge] = this_cluster
                                    
                                    
                            else:
                                this_cluster = cluster_graph_data[i_,j_,k_]
                                
                            #find neighbours and give cluster_number
                            if time_slice[i_-1,j_,k_] == 1:
                                cluster_graph_data[i_-1,j_,k_] = this_cluster
                            elif time_slice[i_+1,j_,k_] == 1:
                                cluster_graph_data[i_+1,j_,k_] = this_cluster
                            elif time_slice[i_,j_-1,k_] == 1:
                                cluster_graph_data[i_,j_-1,k_] = this_cluster
                            elif time_slice[i_,j_+1,k_] == 1:
                                cluster_graph_data[i_,j_+1,k_] = this_cluster
                            elif time_slice[i_,j_,k_-1] == 1:
                                cluster_graph_data[i_,j_,k_-1] = this_cluster
                            elif time_slice[i_,j_,k_+1] == 1:
                                cluster_graph_data[i_,j_,k_+1] = this_cluster
                            elif time_slice[i_-2,j_,k_] == 1:
                                cluster_graph_data[i_-1,j_,k_] = this_cluster
                            elif time_slice[i_+2,j_,k_] == 1:
                                cluster_graph_data[i_+1,j_,k_] = this_cluster
                            elif time_slice[i_,j_-2,k_] == 1:
                                cluster_graph_data[i_,j_-1,k_] = this_cluster
                            elif time_slice[i_,j_+2,k_] == 1:
                                cluster_graph_data[i_,j_+1,k_] = this_cluster
                            elif time_slice[i_,j_,k_-2] == 1:
                                cluster_graph_data[i_,j_,k_-1] = this_cluster
                            elif time_slice[i_,j_,k_+2] == 1:
                                cluster_graph_data[i_,j_,k_+1] = this_cluster     
                          
                                #find neighbours and give this_cluster
                                
                    # if not == 1¡, keep the search 
                        # if not neighbours, keep the search
                       
        cluster_graph_data_total[:,:,:,t_] = cluster_graph_data 
        
        img_new = nb.Nifti1Image(cluster_graph_data_total, header=img.get_header(), affine=img.get_affine())
        # Reconstruct the 4D volume
        cluster_graph_img = os.path.join(os.getcwd(), 'cluster_2N.nii.gz')
        img_new.to_filename(cluster_graph_img)   
        
    return cluster_graph_img
    
    
#def cluster_detection_mod3(in_file):  
#
#    import numpy as np
#    import os
#    import nibabel as nb
#    from CPAC.criticallity import point_process  
#
#    # Treat fMRI image
#    img = nb.load(in_file)
#    data = img.get_data()
#    
#    (n_x, n_y, n_z, n_t) = data.shape
#    
#    # Get the PP data
#    pp_data = np.zeros((n_x, n_y, n_z, n_t))
#    for i_ in range(n_x):
#        for j_ in range(n_y):
#            for k_ in range(n_z):
#                voxel_data = data[i_,j_,k_,:] 
#                pp_data[i_,j_,k_,:] = point_process(voxel_data)
#    
#    cluster_graph_data_total = np.zeros((n_x, n_y, n_z, n_t))            
#    for t_ in range(n_t):
#        time_slice = pp_data[:,:,:,t_]
#        cluster_graph_data = np.zeros((n_x, n_y, n_z))  
#        cluster_number = 1
#        
#        for i_ in range(n_x):
#            for j_ in range(n_y):
#                for k_ in range(n_z):
#                    
#                    if time_slice[i_,j_,k_] == 1: # is active, check if it has active neighboours
#                        if time_slice[i_-1,j_,k_] or time_slice[i_+1,j_,k_] \
#                        or time_slice[i_,j_-1,k_] or time_slice[i_,j_+1,k_] \
#                        or time_slice[i_,j_,k_-1] or time_slice[i_,j_,k_+1] \
#                        or time_slice[i_-2,j_,k_] or time_slice[i_+2,j_,k_] \
#                        or time_slice[i_,j_-2,k_] or time_slice[i_,j_+2,k_] \
#                        or time_slice[i_,j_,k_-2] or time_slice[i_,j_,k_+2] \
#                        or time_slice[i_-3,j_,k_] or time_slice[i_+3,j_,k_] \
#                        or time_slice[i_,j_-3,k_] or time_slice[i_,j_+3,k_] \
#                        or time_slice[i_,j_,k_-3] or time_slice[i_,j_,k_+3]:
#                            
#                            if cluster_graph_data[i_,j_,k_] == 0: # if is not in any previous cluster
#                                this_cluster = (cluster_graph_data[i_-1,j_,k_] or cluster_graph_data[i_+1,j_,k_] \
#                                or cluster_graph_data[i_,j_-1,k_] or cluster_graph_data[i_,j_+1,k_] \
#                                or cluster_graph_data[i_,j_,k_-1] or cluster_graph_data[i_,j_,k_+1] \
#                                or cluster_graph_data[i_-2,j_,k_] or cluster_graph_data[i_+2,j_,k_] \
#                                or cluster_graph_data[i_,j_-2,k_] or cluster_graph_data[i_,j_+2,k_] \
#                                or cluster_graph_data[i_,j_,k_-2] or cluster_graph_data[i_,j_,k_+2] \
#                                or cluster_graph_data[i_-3,j_,k_] or cluster_graph_data[i_+3,j_,k_] \
#                                or cluster_graph_data[i_,j_-3,k_] or cluster_graph_data[i_,j_+3,k_] \
#                                or cluster_graph_data[i_,j_,k_-3] or cluster_graph_data[i_,j_,k_+3])
#                                
#                                if this_cluster == 0: #no neighbours in any previous cluster neither
#                                    this_cluster = cluster_number
#                                    cluster_graph_data[i_,j_,k_] = this_cluster
#                                    cluster_number = cluster_number + 1
#                                else: 
#                                    #check cluster union
#                                    merge_clusters = np.unique([cluster_graph_data[i_-1,j_,k_], cluster_graph_data[i_+1,j_,k_] \
#                                , cluster_graph_data[i_,j_-1,k_], cluster_graph_data[i_,j_+1,k_] \
#                                , cluster_graph_data[i_,j_,k_-1], cluster_graph_data[i_,j_,k_+1] \
#                                , cluster_graph_data[i_-2,j_,k_] , cluster_graph_data[i_+2,j_,k_] \
#                                , cluster_graph_data[i_,j_-2,k_] , cluster_graph_data[i_,j_+2,k_] \
#                                , cluster_graph_data[i_,j_,k_-2] , cluster_graph_data[i_,j_,k_+2] \
#                                , cluster_graph_data[i_-3,j_,k_] , cluster_graph_data[i_+3,j_,k_] \
#                                , cluster_graph_data[i_,j_-3,k_] , cluster_graph_data[i_,j_+3,k_] \
#                                , cluster_graph_data[i_,j_,k_-3] , cluster_graph_data[i_,j_,k_+3]])
#                                    merge_clusters = merge_clusters[1:] #quit first value = 0
#                                    
#                                    this_cluster = merge_clusters[0]
#                                    cluster_graph_data[i_,j_,k_] = this_cluster
#                                    for cluster_to_merge in merge_clusters[1:]:
#                                        cluster_graph_data[cluster_graph_data == cluster_to_merge] = this_cluster
#                                    
#                                    
#                            else:
#                                this_cluster = cluster_graph_data[i_,j_,k_]
#                                
#                            #find neighbours and give cluster_number
#                            if time_slice[i_-1,j_,k_] == 1:
#                                cluster_graph_data[i_-1,j_,k_] = this_cluster
#                            elif time_slice[i_+1,j_,k_] == 1:
#                                cluster_graph_data[i_+1,j_,k_] = this_cluster
#                            elif time_slice[i_,j_-1,k_] == 1:
#                                cluster_graph_data[i_,j_-1,k_] = this_cluster
#                            elif time_slice[i_,j_+1,k_] == 1:
#                                cluster_graph_data[i_,j_+1,k_] = this_cluster
#                            elif time_slice[i_,j_,k_-1] == 1:
#                                cluster_graph_data[i_,j_,k_-1] = this_cluster
#                            elif time_slice[i_,j_,k_+1] == 1:
#                                cluster_graph_data[i_,j_,k_+1] = this_cluster
#                            elif time_slice[i_-2,j_,k_] == 1:
#                                cluster_graph_data[i_-1,j_,k_] = this_cluster
#                            elif time_slice[i_+2,j_,k_] == 1:
#                                cluster_graph_data[i_+1,j_,k_] = this_cluster
#                            elif time_slice[i_,j_-2,k_] == 1:
#                                cluster_graph_data[i_,j_-1,k_] = this_cluster
#                            elif time_slice[i_,j_+2,k_] == 1:
#                                cluster_graph_data[i_,j_+1,k_] = this_cluster
#                            elif time_slice[i_,j_,k_-2] == 1:
#                                cluster_graph_data[i_,j_,k_-1] = this_cluster
#                            elif time_slice[i_,j_,k_+2] == 1:
#                                cluster_graph_data[i_,j_,k_+1] = this_cluster     
#                            elif time_slice[i_-3,j_,k_] == 1:
#                                cluster_graph_data[i_-1,j_,k_] = this_cluster
#                            elif time_slice[i_+3,j_,k_] == 1:
#                                cluster_graph_data[i_+1,j_,k_] = this_cluster
#                            elif time_slice[i_,j_-3,k_] == 1:
#                                cluster_graph_data[i_,j_-1,k_] = this_cluster
#                            elif time_slice[i_,j_+3,k_] == 1:
#                                cluster_graph_data[i_,j_+1,k_] = this_cluster
#                            elif time_slice[i_,j_,k_-3] == 1:
#                                cluster_graph_data[i_,j_,k_-1] = this_cluster
#                            elif time_slice[i_,j_,k_+3] == 1:
#                                cluster_graph_data[i_,j_,k_+1] = this_cluster 
#                                #find neighbours and give this_cluster
#                                
#                    # if not == 1¡, keep the search 
#                        # if not neighbours, keep the search
#                       
#        cluster_graph_data_total[:,:,:,t_] = cluster_graph_data 
#       
#        img_new = nb.Nifti1Image(cluster_graph_data_total, header=img.get_header(), affine=img.get_affine())
#        # Reconstruct the 4D volume
#        cluster_graph = os.path.join(os.getcwd(), 'cluster_3N.nii.gz')
#        img_new.to_filename(cluster_graph)
#    
#    return cluster_graph_data_total  
    
  
def avalanche_detec(cluster_file):
    """
    Detects avalanches if a cluster file is given to it
    as described in http://journal.frontiersin.org/article/10.3389/fphys.2012.00015/abstract  
    
    Parameters
    ----------

    cluster_file : Previously calculated cluster_graph file

    Returns
    -------

    avalanche_img  :  Nifti file: 4D with an id for each avalanche
    """
    
    import numpy as np
    import nibabel as nb
    import os

    # Treat fMRI image
    img = nb.load(cluster_file)
    cluster_data = img.get_data()
  
    (n_x, n_y, n_z, n_t) = cluster_data.shape
    
    avalanche_id_total = np.zeros((n_x, n_y, n_z, n_t))
    avalanche_id_num = 1   
            
    for t_ in range(n_t):
     
        if t_ == 0: #if first timestep, all are candidates
            time_slice = cluster_data[:,:,:,t_]          
            time_slice_fut = cluster_data[:,:,:,t_+1]
            
            avalanche_id_now = avalanche_id_total[:,:,:,t_]   
            avalanche_id_fut = avalanche_id_total[:,:,:,t_+1]   
            
            for cluster in np.unique(cluster_data[:,:,:,t_])[1:]: #iterate over clusters
                # NEW AVALANCHE CASE
                if np.count_nonzero(time_slice_fut[(time_slice==cluster)]) >= 1 :                
                    avalanche_id_now[(time_slice==cluster)] = avalanche_id_num # assign cluster a aval_id
                    
                    for value in np.unique(time_slice_fut[(time_slice==cluster)])[1:]: #assing used clusters of t+1
                        avalanche_id_fut[(time_slice_fut==value)] = avalanche_id_num
                    
                    avalanche_id_num = avalanche_id_num +1 # Ivan: Would be interesting to define the length of sustaining avalanches
                    
                    avalanche_id_total[:,:,:,t_] = avalanche_id_now
                    avalanche_id_total[:,:,:,t_+1] = avalanche_id_fut
                    
                    
        elif t_ < (n_t-1):  #if not first timestep, check previous
            time_slice = cluster_data[:,:,:,t_]          
            time_slice_fut = cluster_data[:,:,:,t_+1]
            
            avalanche_id_now = avalanche_id_total[:,:,:,t_]   
            avalanche_id_fut = avalanche_id_total[:,:,:,t_+1] 
            
            for cluster in np.unique(cluster_data[:,:,:,t_])[1:]:
                # PREVIOUS AVALANCHE CASE
                if np.count_nonzero(avalanche_id_now[(time_slice==cluster)]) != 0: 
                    if np.count_nonzero(time_slice_fut[(time_slice==cluster)]) >= 1 :
                        
                        this_avalanche = avalanche_id_now[(time_slice==cluster)][0]
                        
                        for value in np.unique(time_slice_fut[(time_slice==cluster)])[1:]: #assing used clusters of the t+1
                            avalanche_id_fut[(time_slice_fut==value)] = this_avalanche
                        
                        avalanche_id_total[:,:,:,t_+1] = avalanche_id_fut 
                
                # Maybe this is wrong, since if it has a positive intersect with the past, already is part of the avalanches.
                # NEW AVALANCHE CASE ## HERE MAYBE PROBLEMS
                elif np.count_nonzero(avalanche_id_now[(time_slice==cluster)]) == 0: #and np.count_nonzero(time_slice_past[(time_slice==cluster)]) == 0:
                    if np.count_nonzero(time_slice_fut[(time_slice==cluster)]) >= 1 :
                        
                        avalanche_id_now[(time_slice==cluster)] = avalanche_id_num # assign cluster a aval_id
                        
                        for value in np.unique(time_slice_fut[(time_slice==cluster)])[1:]: #assing the used clusters of the t+1
                            avalanche_id_fut[(time_slice_fut==value)] = avalanche_id_num
                            
                        avalanche_id_num = avalanche_id_num + 1
                   
                        avalanche_id_total[:,:,:,t_] = avalanche_id_now
                        avalanche_id_total[:,:,:,t_+1] = avalanche_id_fut    
    
    
    img_new = nb.Nifti1Image(avalanche_id_total, header=img.get_header(), affine=img.get_affine())
    # Reconstruct the 4D volume
    avalanche_img = os.path.join(os.getcwd(), 'avalanche.nii.gz')
    img_new.to_filename(avalanche_img)  
    
    return avalanche_img    
    
    
    
    
    
    
    
    
    
    
    
    
    