# -*- coding: utf-8 -*-

# coding: utf-8
#import os
#import sys
#import re
#import commands
import nipype.pipeline.engine as pe
#import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from CPAC.series_mod.utils import compute_ROI_corr
from CPAC.series_mod.utils import compute_ROI_pcorr  
from CPAC.series_mod.utils import compute_MI  
from CPAC.series_mod.utils import compute_TE
#from CPAC.series_mod.utils import compute_ApEn


def create_nltsa(wf_name = 'nltsa_wf'):
  
    """  
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_ROI_corr()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()
    
    
    
    
    workflow.connect(resample_functional_to_template, 'out_file',
                                 nltsa, 'inputspec.subject')
                # Subject mask/parcellation image
                workflow.connect(template_dataflow, 'outputspec.out_file',
                                 nltsa, 'inputspec.template')
                # Give which method we're doing (0 - PRE, 1 - IT, 2 - SFA)
                nltsa.inputs.inputspec.method_option = \
                methodOption
                # List of booleans, a measure in each one
                nltsa.inputs.inputspec.measures = \
                measures




    """
    if method_option == 0:
        calc_degree = True
    elif method_option == 1:
        calc_eigen = True
    elif method_option == 2:
        calc_lfcd = True
    # Weighting
    out_binarize = weight_options[0]
    out_weighted = weight_options[1]
    
    
    
    nltsa = pe.Workflow(wf_name)
    
    
    # Give which method we're doing (0 - PRE, 1 - IT, 2 - SFA)
    if method_option == 0:
        if measures[0] == True:
            
            ROI_corr = pe.Workflow(wf_name)
            inputNode = pe.Node(util.IdentityInterface(fields=[
                                                        'in_file',
                                                        'mask_file',
                                                        'voxelwise'
                                                        ]),
                                name='inputspec')
        
        
            outputNode = pe.Node(util.IdentityInterface(fields=[
                                                            'corr_mat']),
                                name='outputspec')
        
            
        
            corr_matNode = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                           output_names=['corr_mat'],
                             function=compute_ROI_corr),
                             name='corr_calc')
        
        
            nltsa.connect(inputNode, 'in_file',
                            corr_matNode, 'in_file')
            nltsa.connect(inputNode, 'mask_file',
                            corr_matNode, 'mask_file')  
                            
            nltsa.connect(corr_matNode, 'corr_mat',
                         outputNode, 'corr_mat')
     elif method_option == 1:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=True)
    # Otherwise, compute blocksize with regards to available memory
    else:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=False)



    '''
    Method to calculate centrality and map them to a nifti file
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file)
        path to mask/parcellation unit
    method_option : integer
        0 - degree centrality calculation, 1 - eigenvector centrality calculation, 2 - lFCD calculation
    threshold_option : an integer
        0 for probability p_value, 1 for sparsity threshold, 
        2 for actual threshold value, and 3 for no threshold and fast approach
    threshold : a float
        pvalue/sparsity_threshold/threshold value
    weight_options : list (boolean)
        list of booleans, where, weight_options[0] corresponds to binary counting 
        and weight_options[1] corresponds to weighted counting (e.g. [True,False]) 
    allocated_memory : string
        amount of memory allocated to degree centrality
    
    Returns
    -------
    out_list : list
        list containing out mapped centrality images
    '''
    
    # Import packages
    from CPAC.network_centrality import load,\
                                        get_centrality_by_rvalue,\
                                        get_centrality_by_sparsity,\
                                        get_centrality_fast,\
                                        map_centrality_matrix,\
                                        calc_blocksize,\
                                        convert_pvalue_to_r
    from CPAC.cwas.subdist import norm_cols
    
    # Check for input errors
    if weight_options.count(True) == 0:
        raise Exception("Invalid values in weight options" \
                        "At least one True value is required")
    # If it's sparsity thresholding, check for (0,1]
    if threshold_option == 1:
        if threshold <= 0 or threshold > 1:
            raise Exception('Threshold value must be a positive number'\
                            'greater than 0 and less than or equal to 1.'\
                            '\nCurrently it is set at %d' % threshold)
    if method_option == 2 and threshold_option != 2:
        raise Exception('lFCD must use correlation-type thresholding.'\
                         'Check the pipline configuration has this setting')
    import time
    start = time.clock()
    
    # Init variables
    out_list = []
    ts, aff, mask, t_type, scans = load(datafile, template)
    
    # If we're doing eigenvectory centrality, need entire correlation matrix
    if method_option == 0 and threshold_option == 1:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    sparsity_thresh=threshold)
    elif method_option == 1:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=True)
    # Otherwise, compute blocksize with regards to available memory
    else:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=False)
    # Normalize the timeseries for easy dot-product correlation calc.
    ts_normd = norm_cols(ts.T)
    
    # P-value threshold centrality
    if threshold_option == 0:
        r_value = convert_pvalue_to_r(scans, threshold)
        centrality_matrix = get_centrality_by_rvalue(ts_normd, 
                                                     mask, 
                                                     method_option, 
                                                     weight_options, 
                                                     r_value, 
                                                     block_size)
    # Sparsity threshold
    elif threshold_option == 1:
        centrality_matrix = get_centrality_by_sparsity(ts_normd, 
                                                       method_option, 
                                                       weight_options, 
                                                       threshold, 
                                                       block_size)
    # R-value threshold centrality
    elif threshold_option == 2:
        centrality_matrix = get_centrality_by_rvalue(ts_normd, 
                                                     mask, 
                                                     method_option, 
                                                     weight_options, 
                                                     threshold, 
                                                     block_size)
    # For fast approach (no thresholding)
    elif threshold_option == 3:
        centrality_matrix = get_centrality_fast(ts, method_option)
    # Otherwise, incorrect input for threshold_option
    else:
        raise Exception('Option must be between 0-3 and not %s, check your '\
                        'pipeline config file' % str(threshold_option))
    
    # Print timing info
    print 'Timing:', time.clock() - start
 
    # Map the arrays back to images
    for mat in centrality_matrix:
        centrality_image = map_centrality_matrix(mat, aff, mask, t_type)
        out_list.append(centrality_image)
    
    # Finally return
    return out_list


    return ROI_corr


def create_ROI_corr():

    """
    Simple Network Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    corr : workflow
        Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.mask_file : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.corr_matrix : ROI_number * ROI_number array


    Corr Workflow Procedure:

    1. Generate Correlation Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_ROI_corr()
    >>> wf.inputs.inputspec.voxelwise = 0
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_ROI_corr()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    ROI_corr = pe.Workflow(name='ROI_corr')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file',
                                                'voxelwise'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'corr_mat']),
                        name='outputspec')



    corr_matNode = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['corr_mat'],
                     function=compute_ROI_corr),
                     name='corr_calc')


    ROI_corr.connect(inputNode, 'in_file',
                    corr_matNode, 'in_file')
    ROI_corr.connect(inputNode, 'mask_file',
                    corr_matNode, 'mask_file')  
                    
    ROI_corr.connect(corr_matNode, 'corr_mat',
                 outputNode, 'corr_mat')



    return ROI_corr


def create_ROI_pcorr():

    """
    Simple Network Partial Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    pcorr : workflow
        Partial Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.mask_file : string (existing nifti file)
            Input Whole Brain Mask of EPI 4D Volume

    Workflow Outputs: ::

        outputspec.pcorr_matrix : ROI_number * ROI_number array


    Corr Workflow Procedure:

    1. Generate Partial Correlation Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_ROI_corr()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_ROI_corr()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    ROI_pcorr = pe.Workflow(name='ROI_pcorr')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file',
                                                'voxelwise'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'pcorr_mat']),
                        name='outputspec')


    pcorr_matNode = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['pcorr_mat'],
                     function=compute_ROI_pcorr),
                     name='pcorr_calc')


    ROI_pcorr.connect(inputNode, 'in_file',
                    pcorr_matNode, 'in_file')
    ROI_pcorr.connect(inputNode, 'mask_file',
                    pcorr_matNode, 'mask_file')  
                    
    ROI_pcorr.connect(pcorr_matNode, 'pcorr_mat',
                 outputNode, 'pcorr_mat')



    return ROI_pcorr    
    
def create_MI():

    """
    Simple Network Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    corr : workflow
        Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.bins : integer
            How many states do you want to be the data transformed to

    Workflow Outputs: ::

        outputspec.MI_matrix : [ROI_number * ROI_number] array

    Corr Workflow Procedure:

    1. Generate MI Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_MI()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    from CPAC import series_mod
    wf = series_mod.create_MI()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    MI = pe.Workflow(name='MI_comp')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file',
                                                'voxelwise'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'MI_mat']),
                        name='outputspec')


    MI_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['MI_mat'],
                     function=compute_MI),
                     name='MI_mat')


    MI.connect(inputNode, 'in_file',
                    MI_mat, 'in_file')
    MI.connect(inputNode, 'mask_file',
                    MI_mat, 'mask_file') 
                    
    MI.connect(MI_mat, 'MI_mat',
                 outputNode, 'MI_mat')

    return MI
    
    
def create_TE():

    """
    Simple Network Correlation Matrix calculation for ROIs in the mask file

    Parameters
    ----------

    None

    Returns
    -------
    corr : workflow
        Correlation Workflow

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/series_mod/series_mod.py>`_

    Workflow Inputs: ::

        inputspec.in_file : string (existing nifti file)
            Input EPI 4D Volume

        inputspec.bins : integer
            How many states do you want to be the data transformed to

    Workflow Outputs: ::

        outputspec.TE_matrix : [ROI_number * ROI_number] array

    Corr Workflow Procedure:

    1. Generate TE Matrix from the input EPI 4D volume and EPI mask.


    Workflow Graph:


    Detailed Workflow Graph:

        
    References
    ---------- 

    Examples
    --------
    >>> from CPAC import series_mod
    >>> wf = series_mod.create_TE()
    >>> wf.inputs.inputspec.in_file = '/home/data/Project/subject/func/in_file.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/Project/subject/func/mask.nii.gz'
    >>> wf.run()
    
    in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    wf = series_mod.create_TE()
    wf.inputs.inputspec.in_file = in_file
    wf.inputs.inputspec.mask_file = mask_file
    wf.base_dir = '/home/asier/git/C-PAC/CPAC/series_mod'
    wf.run()

    """



    TE = pe.Workflow(name='TE_comp')
    inputNode = pe.Node(util.IdentityInterface(fields=[
                                                'in_file',
                                                'mask_file',
                                                'voxelwise'
                                                ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                    'TE_mat']),
                        name='outputspec')


    TE_mat = pe.Node(util.Function(input_names=['in_file', 'mask_file'],
                                   output_names=['MI_mat'],
                     function=compute_TE),
                     name='MI_mat')


    TE.connect(inputNode, 'in_file',
                    TE_mat, 'in_file')
    TE.connect(inputNode, 'mask_file',
                    TE_mat, 'mask_file') 
                    
    TE.connect(TE_mat, 'TE_mat',
                 outputNode, 'TE_mat')

    return TE
    
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