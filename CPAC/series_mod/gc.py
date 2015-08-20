def compute_pwcgc(in_file, lag = 1):
    """
    Computes Pairwise Conditional Granger Causality from a 1D datafile and returns a np.array.
    
    Parameters
    ----------

    in_file : 1D file

    Returns
    -------

    data_array =  voxel(x,y,z) * voxel(x,y,z)

    """

    from CPAC.series_mod import pwcgc  
    import numpy as np

    # treatment 1D matrix
    data = np.genfromtxt(in_file, skip_header = 1)[:,2:].T
    # TEST IF THE STRUCTURE OF THE DATA IS ADECUATE. I DON'T KNOW HOW 1D FILES FROM CPAC ARE
    
    pwcgc_mat = pwcgc(data, lag)
    
    np.save(in_file[:-7]+'_pwcgc.npy', pwcgc_mat)

    return pwcgc_mat










def pwcgc(tsdata, p):
    
    import numpy as np
    from series_mod.gc import tsdata_to_var


    n = tsdata.shape[0]
    F = np.zeros([n,n])
    
    # full regression   
    [AF,SIG,E] = tsdata_to_var(tsdata, p)    
    LSIG = np.log(abs(np.diag(SIG)))
    
    for j_ in range(n):
    
        # reduced regression
        jo = np.arange(n) # omit j
        jo = np.delete(jo,j_)

        [AF,SIGj,E] = tsdata_to_var(tsdata[jo], p)  
        LSIGj = np.log(abs(np.diag(SIGj)))
    
        F[jo,j_] = LSIGj-LSIG[jo]
    
#        for ii_ in range(n-1):
#            i_ = jo[ii_]
#            F[i_,j_] = LSIGj[ii_]-LSIG[i_]
        
    return F
    

# autocov_to_mvgc
#
# Calculate conditional time-domain MVGC (multivariate Granger causality)
#
# <matlab:open('autocov_to_mvgc.m') code>
#
# Syntax
#
#     F = autocov_to_mvgc(G,x,y)
#
# Arguments
#
# See also <mvgchelp.html#4 Common variable names and data structures>.
#
# _input_
#
#     G          autocovariance sequence
#     x          vector of indices of target (causee) multi-variable
#     y          vector of indices of source (causal) multi-variable
#
# _output_
#
#     F          Granger causality
#
# Description
#
# Returns the time-domain MVGC
#
# <<eq_mvgc.png>>
#
# from the variable |Y| (specified by the vector of indices |y|) to the
# variable |X| (specified by the vector of indices |x|), conditional on all
# other variables |Z| represented in |G|, for a stationary VAR process with
# autocovariance sequence |G|. See ref. [1] for details.
#
# The caller should take note of any warnings issued by this function and test
# results with a call <isbad.html |isbad|>|(F,false)|.
#
# References
#
# [1] L. Barnett and A. K. Seth,
# <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
#     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
# Inference>, _J. Neurosci. Methods_ 223, 2014
# [ <matlab:open('mvgc_preprint.pdf') preprint> ].
#
# See also
#
# <autocov_to_var.html |autocov_to_var|> |
# <isbad.html |isbad|>
#
# (C) Lionel Barnett and Anil K. Seth, 2012.

def mvgc(tsdata, x, y, p):
    
    import numpy as np
    from series_mod import tsdata_to_var

    # Local Variables: G, F, xzy, n, xz, SIGR, SIG, y, x, z
    # Function calls: autocov_to_var, log, det, NaN, length, autocov_to_mvgc, size
    n = tsdata.shape[0]
    
    
     #WARNING PART 1!!
#    x = x.flatten(0).conj() #be carefull for multi variate case
#    # vectorise
#    y = y.flatten(0).conj()
#    # vectorise
    
    
    z = np.arange(n)
    z = np.delete(z,[np.array(np.hstack((x, y)))])
    # indices of other variables (to condition out)
    xz = np.array(np.hstack((x, z)))
    xzy = np.array(np.hstack((xz, y)))
    F = 0
    # full regression
    ixgrid1 = np.ix_(xzy,xzy)
    [AF,SIG,E] = tsdata_to_var(tsdata[ixgrid1], p) #G[ixgrid1,:])
    # reduced regression
    ixgrid2 = np.ix_(xz,xz)
    [AF,SIGR,E] = tsdata_to_var(tsdata[ixgrid2], p) #G[ixgrid2,:])
    # reduced regression
    
    #WARNING PART 2!!
    #x = np.arange(np.size(x,axis=1)+1)
    
    ###########
    ixgrid3 = np.ix_(x,x)
    F = np.log(abs(np.linalg.det(SIGR[ixgrid3])))-np.log(abs(np.linalg.det(SIG[ixgrid3]))) 
    #####  not tested
    return F
    
    
  

def tsdata_to_var(X, p):

    import numpy as np
    from matplotlib import pylab

    X = X[:,:,np.newaxis]
    [n, m, N] = X.shape
    N=1
    # assert(p < m,'too many lags');
    p1 = p+1
    
    A = 0 #nan
    SIG = 0 #nan
    E = 0 #nan
    
    X = pylab.demean(X, axis=1)
    # no constant term

    M = np.dot(N, m-p)
    n_p = np.dot(n, p)
    # stack lags
    X0 = np.reshape(X[:,p1-1:m], (n, M))
    # concatenate trials for unlagged observations
    XL = np.zeros((n, p, M))
    for k in range(1, (p)+1): #if lag = 1, only 1 iteration
        XL[:,k-1,:] = np.reshape(X[:,p1-k-1:m-k,:], (n, M))
        # concatenate trials for k-lagged observations
        
    XL = np.reshape(XL, (n_p, M))
    # stack lags
    
    A = np.linalg.lstsq(XL.T,X0.T)[0].T
    # so A(:,:,k) is the k-lag coefficients matrix
    
    # OLS using QR decomposition
    
    E = X0-np.dot(A, XL)
    # residuals
    
    SIG = np.dot(E, E.T)/(M-1)
    # residuals covariance matrix
    
    
    return [A, SIG, E]  
    