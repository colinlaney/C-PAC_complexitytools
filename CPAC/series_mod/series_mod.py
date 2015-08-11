# -*- coding: utf-8 -*-

# coding: utf-8
#import os
#import sys
#import re
#import commands
import nipype.pipeline.engine as pe
#import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
#from CPAC.series_mod.utils import compute_ApEn


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

    
#def create_ApEn():
#
#    """
#    ApEn calculation for fMRI file
#
#    Parameters
#    ----------
#
#    None
#
#    Returns
#    -------
#    corr : workflow
#        Correlation Workflow
#
#    Notes
#    -----
#
#    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_
#
#    Workflow Inputs: ::
#
#        inputspec.in_file : string (existing nifti file)
#            Input EPI 4D Volume
#
#	inputspec.m_param : parameter m, time series size
#
#        inputspec.r_param : parameter r - factor?¿
#
#    Workflow Outputs: ::
#
#        outputspec.result_vector : ApEn value
#        
#    References
#    ---------- 
#
#    Examples
#    --------
#    >>> from CPAC import series_mod
#    >>> wf = series_mod.create_ApEn()
#    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
#    >>> wf.inputs.inputspec.m_param = 30
#    >>> wf.inputs.inputspec.m_param = 4
#    >>> wf.run()
#
#    """
#
#
#
#    ApEn = pe.Workflow(name='ApEn')
#    inputNode = pe.Node(util.IdentityInterface(fields=[
#                                                'in_file',
#                                                'm_param',
#						'r_param'
#                                                ]),
#                        name='inputspec')
#
#
#    outputNode = pe.Node(util.IdentityInterface(fields=[
#                                                    'result_vector']),
#                        name='outputspec')
#
#
#
#    ApEn_calc_Node = pe.Node(util.Function(input_names=['in_file', 'm_param','r_param'],
#                                   output_names=['result_vector'],
#                     function=compute_ApEn),
#                     name='ApEn_calc')
#
#
#    ApEn.connect(inputNode, 'in_file',
#                    ApEn_calc_Node, 'in_file')
#    ApEn.connect(inputNode, 'm_param',
#                    ApEn_calc_Node, 'm_param')  
#    ApEn.connect(inputNode, 'r_param',
#                    ApEn_calc_Node, 'r_param')  
#                    
#    ApEn.connect(ApEn_calc_Node, 'result_vector',
#                 outputNode, 'result_vector')
#
#
#
#    return ApEn    