
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
                                   name = 'calculate_centrality')
    
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
               outputspec, 'centrality_outputs')

    return nltsa               
               
               
def calc_nltsa(timeseries_one_d,
                    method_option,
                    measures):   
                        
    from CPAC.series_mod.utils import compute_corr
    from CPAC.series_mod.utils import compute_pcorr  
    from CPAC.series_mod.utils import compute_MI  
    from CPAC.series_mod.utils import compute_TE                    
        
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
            
    # Otherwise, CRITICALLITY [MAYBE ANOTHER WF]
    #else:


    return output_list