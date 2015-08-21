# -*- coding: utf-8 -*-

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def create_nltsa(wf_name = 'nltsa_wf'):
  
    """  
    # Give which method we're doing (0 - PRE, 1 - IT, 2 - SFA)
    nltsa.inputs.inputspec.method_option = \
    methodOption
    # List of booleans, a measure in each one
    nltsa.inputs.inputspec.measures = \
    measures
    """    
    from CPAC.series_mod import calc_nltsa
    
     # Instantiate workflow with input name
    nltsa = pe.Workflow(wf_name)
    
    # Instantiate inputspec node
    inputspec = pe.Node(util.IdentityInterface(fields=['timeseries_one_d',
                                                       'method_option',
                                                       'measures']),
                        name='inputspec')
    
    # Instantiate calculate_centrality main function node
    calculate_nltsa = pe.Node(util.Function(input_names = ['timeseries_one_d',
                                                       'method_option',
                                                       'measures'],
                                                 output_names = ['out_list'],
                                                 function = calc_nltsa),
                                   name = 'calc_nltsa')
    
    # Connect inputspec node to main function node
    nltsa.connect(inputspec, 'timeseries_one_d', 
               calc_nltsa, 'timeseries_one_d')
    nltsa.connect(inputspec, 'method_option',
               calc_nltsa, 'method_option')
    nltsa.connect(inputspec, 'measures', 
               calc_nltsa, 'measures')
        
    # Instantiate outputspec node
    outputspec = pe.Node(util.IdentityInterface(fields=['nltsa_outputs']),
                         name = 'outputspec')
    
    # Connect function node output list to outputspec node
    nltsa.connect(calculate_nltsa, 'out_list',
               outputspec, 'nltsa_outputs')

    return nltsa               
               
               
def calc_nltsa(timeseries_one_d,
                    method_option,
                    measures):   
                        
    from CPAC.series_mod.utils import compute_corr
    from CPAC.series_mod.utils import compute_pcorr  
    from CPAC.series_mod.utils import compute_MI  
    from CPAC.series_mod.utils import compute_TE
    from CPAC.series_mod.gc import compute_pwcgc     
    #from CPAC.series_mod.criticality import compute_avalanche                  
               
        
    output_list = []    
        
    # Give which method we're doing (0 - PRE, 1 - IT, 2 - SFA)
    if method_option == 0:
        if measures[0] == True: # Give which measure to calc (0 - CORR, 1 - PCORR, 2 - PSI, 3 - PLV)  
            corr = compute_corr(timeseries_one_d)
            output_list.append(corr)                     
        if measures[1] == True:    #PCORR
            pcorr = compute_pcorr(timeseries_one_d)
            output_list.append(pcorr)            
            
            
            
        #if measures[2] == True:    #PSI
           #
        #if measures[3] == True:    #PLV
           #


                
    elif method_option == 1:  # Give which measure to calc (0 - ENT, 1 - CONDENT, 2 - MI, 3 - TE, 4 - ECC)
#        if measures[0] == True:    #ENT
#            
#        if measures[1] == True:    #CONDENT    
            
        if measures[2] == True:     #MI
            MI = compute_MI(timeseries_one_d)
            output_list.append(MI) 
        if measures[3] == True:    #TE 
            TE = compute_TE(timeseries_one_d)
            output_list.append(TE)            
            
        
#        if measures[4] == True:    #ECC
        
        if measures[5] == True:    #PWGC
            PWCGC = compute_pwcgc(timeseries_one_d)
            output_list.append(PWCGC) 
            
#    elif method_option == 2:  # Give which measure to calc (0 - AVALANCHE DETECTION, 1 - ?¿)
#        if measures[0] == True:    #AVALANCHE
#            AVALANCHE = compute_avalanche(fMRI) # This needs an fMRI file
#            output_list.append(AVALANCHE)     
#            
#        if measures[1] == True:    #FRACTALITY ?¿    
           
            
    #else:


    return output_list
    
    
def create_avalanche(wf_name = 'avalanche_wf'):
  
    """  
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_avalanche()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/rest_res_filt.nii.gz'
    >>> wf.run()
    """    
    from CPAC.series_mod import compute_avalanche      
    # Instantiate workflow with input name
    avalanche = pe.Workflow(wf_name)
    
    # Instantiate inputspec node
    inputspec = pe.Node(util.IdentityInterface(fields=['in_file']),
                        name='inputspec')
    
    # Instantiate calculate_centrality main function node
    calculate_avalanche = pe.Node(util.Function(input_names = ['timeseries_one_d'],
                                                 output_names = ['out_list'],
                                                 function = compute_avalanche),
                                   name = 'calculate_avalanche')
    
    # Connect inputspec node to main function node
    avalanche.connect(inputspec, 'timeseries_one_d', 
               calculate_avalanche, 'timeseries_one_d')

        
    # Instantiate outputspec node
    outputspec = pe.Node(util.IdentityInterface(fields=['avalanche_outputs']),
                         name = 'outputspec')
    
    # Connect function node output list to outputspec node
    avalanche.connect(calculate_avalanche, 'out_list',
               outputspec, 'avalanche_outputs')

    return avalanche               
      