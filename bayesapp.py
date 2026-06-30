#!/usr/bin/python3

import sys
import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import time
import argparse
from sys import argv
from bayesapp_helpfunctions import *

if __name__=='__main__':

    ### start timing
    start_time = time.time()

    ## current version of BayesApp
    version = 2.1

    ### welcome message
    printt('#######################################################################################')
    printt('RUNNING bayesapp.py version %s \n - for instructions type: python bayesapp.py -h' % version)
    command = "python bayesapp.py"
    for aa in argv[1:]:
        if os.path.sep in aa:
            display_aa = os.path.basename(aa) 
        else:
            display_aa = aa
        if ' ' in display_aa:
            command += " \"%s\"" % display_aa
        else:
            command += " %s" % display_aa
    printt('command used: %s' % command)
    printt('#######################################################################################')

    ### input values 
    parser = argparse.ArgumentParser(description='BayesApp - calculates the pair distance distribution from SAXS or SANS data.')
      
    # mandatory inputs
    parser.add_argument('-f', '--data_file', required=True, help='path to data')

    # optional inputs
    parser.add_argument('-qmin', '--qmin', default='0', help='Minimum q-value for the scattering curve.')
    parser.add_argument('-qmax', '--qmax', default='100', help='Maximum q-value for the scattering curve.')
    parser.add_argument('-units', '--units', default='auto', help='units of data: A or nm')

    parser.add_argument('-G', '--Guinier', action='store_true', default=False, help='Guinier plot')
    parser.add_argument('-qmaxRg', '--Guinier_qmaxRg', default='0.7', help='qmax*Rg in Guinier analysis')
    parser.add_argument('-G_skip', '--Guinier_skip', default='0', help='skip first points in Guinier analysis')

    parser.add_argument('-K', '--Kratky', action='store_true', default=False, help='Kratky plot')
    parser.add_argument('-K_dim', '--Kratky_dim', action='store_true', default=False, help='Dimensionless Kratky plot')
    parser.add_argument('-Mw', '--Kratky_Mw', action='store_true', default=False, help='Mw from Kratky plot')

    parser.add_argument('-P', '--Porod', action='store_true', default=False, help='Porod plot')
    parser.add_argument('-P_lim', '--Porod_limit', default='0.0', help='Porod limit start')

    parser.add_argument('-d', '--dmax', default='', help='dmax (add f in front to fix, e.g. f100)')
    parser.add_argument('-T', '--transformation', default='A', help='transformation')
    parser.add_argument('-b', '--nrebin', default='500', help='rebin to approximately this number of points')
    parser.add_argument('-a', '--alpha', default='', help='alpha (add f in front to fix, e.g. f100)')
    parser.add_argument('-s', '--smear', default='', help='smear')
    parser.add_argument('-n', '--prpoints', default='', help='prpoints')
    parser.add_argument('-c', '--noextracalc', default='100', help='number of extra calculations')
    parser.add_argument('-rescale_mode', '--rescale_mode', default='C', help='error rescale mode, if wrong')
    parser.add_argument('-nbin', '--nbin', default=' ', help='nbin')
    parser.add_argument('-Bg', '--Bg', default='0', help='Bg')

    parser.add_argument('-fitbackground', '--fitbackground', default='Y', help='fitbackground: Y or N, default is Y')
    parser.add_argument('-logx', '--logx', action='store_true', default=False, help='logarithmic log axis')
    parser.add_argument('-make_pr_bin', '--make_pr_bin', action='store_true', default=False, help='make_pr_bin')
    parser.add_argument('-pr_binsize', '--pr_binsize', default='1', help='pr_binsize')
    parser.add_argument('-skip', '--skip_first', default='', help='skip first points from analysis')
    parser.add_argument('-outlier_ite', '--outlier_ite', action='store_true', default=False, help='outlier_ite')
    parser.add_argument('-fast', '--fast_run', action='store_true', default=False, help='fast_run')
    
    parser.add_argument('-show', '--show', action='store_true', default=False, help='show plots') 
    parser.add_argument('-z', '--zip_compress', action='store_true', default=False, help='compress output files in zip file')

    args = parser.parse_args()

    #################################################################################
    ## adjust input
    #################################################################################

    # read data name
    data = os.path.basename(args.data_file)

    ## remove spaces and brackets from name
    data = data.replace(" ","_").replace("(","").replace(")","").replace("[","").replace("]","")
    
    ## naming bug fix: fortran77 cannot use long file names
    if len(data)>48:
        data = 'data_name_too_long_for_fortran77.dat'

    ## import skip first, qmin and qmax
    try:
        skip_first = int(args.skip_first) # skip first points
    except:
        skip_first = ''
    try:
        qmin = float(args.qmin)
    except:
        qmin = 0.0
    try:
        qmax = float(args.qmax)
    except:
        qmax = 100.0

    ## read number of header and footer lines in data file
    header,footer = get_header_footer(args.data_file)

    ## read q from data, and update qmin and qmax if needed
    try:
        if skip_first:
            q_check = np.genfromtxt(args.data_file,skip_header=header+skip_first,skip_footer=footer,usecols=[0],unpack=True)
        else:
            q_check = np.genfromtxt(args.data_file,skip_header=header,skip_footer=footer,usecols=[0],unpack=True)
    except:
        if skip_first:
            q_check = np.genfromtxt(args.data_file,skip_header=header+skip_first,skip_footer=footer,usecols=[0],unpack=True,encoding='cp855')
        else:
            q_check = np.genfromtxt(args.data_file,skip_header=header,skip_footer=footer,usecols=[0],unpack=True,encoding='cp855')
    if q_check[0] > qmin:
        qmin = q_check[0]
    if q_check[-1] < qmax:
        qmax = q_check[-1]

    ## check validity of q range
    q_diff = qmax - qmin
    if q_diff < 0.0:
        printt('\n\n!!!ERROR!!!\nqmin should be smaller than qmax.\n')
        sys.exit()

    ## automatically detect units
    if args.units == 'auto':
        if np.max(q_check) > 2.0:
            units = 'nm'
            if qmax > 7:
                qmax = 7
        else:
            units = 'A'
            if qmax > 0.7:
                qmax = 0.7

    ## rename and convert input
    dmax = args.dmax
    transformation = args.transformation
    pr_binsize = float(args.pr_binsize)
    prpoints = args.prpoints
    qmaxRg_in = float(args.Guinier_qmaxRg)
    Guinier_skip_in = int(args.Guinier_skip)
    Porod_limit = float(args.Porod_limit)

    ## remove path from filename (for output)
    filename = os.path.basename(args.data_file)

    ## ensure data is at the current location, and naming is correct
    if not (os.path.exists(data)):
        shutil.copy2(args.data_file, data) # copy data file to current location and change name
    
    ## ensure BIFT executable (bift) and Fortran source code (bift.f) are at the current location
    path = os.path.dirname(os.path.realpath(__file__))
    exe = 'bift.exe' if os.name == 'nt' else 'bift'
    if not (os.path.exists(exe) and os.path.samefile(os.path.join(path, exe), exe)):
        shutil.copy2(os.path.join(path, exe), '.') # copy bift executable to current location
    if not (os.path.exists('bift.f') and os.path.samefile(os.path.join(path, 'bift.f'), 'bift.f')):
        shutil.copy2(os.path.join(path, 'bift.f'), '.') # copy Fortran source code to current location

    ## print start bayesapp message
    printt("=================================================================================")
    printt('    Reading data:                   %s' % filename)
    printt('        header lines in datafile:   %d' % header)
    printt('        footer lines in datafile:   %d' % footer)
    printt('        q min in data:              %f' % q_check[0])
    printt('        q min used                  %f' % qmin)
    printt('        q max in data:              %f' % q_check[-1])
    printt('        q max used                  %f' % qmax)
    printt('        number of points in data:   %d' % len(q_check))
    printt('        data rebinned to around:    %s points' % args.nrebin)
    printt('        units (of q) assumed to be: 1/%s' % units)
    printt("=================================================================================")

    ##################################
    # beginning of outlier while loop 
    ##################################


    # the idea of the next section is first to run a fast run with BIFT to reasonalbe input params, then run a more accurate run

    # if all parameters are set, skip the fast run to avoid an error
    if not dmax == '' and not transformation == 'A' and not skip_first == '' and not prpoints == '' and args.fast_run:
        # making fast run false, because else, in this special case:
        # fast run will first be skipped as all the above are provided (dmax, transformation etc), 
        # and then the normal run will also be skipped, so Bayesapp is not run, which gives errors
        printt('\n    WARNING: changing --fast_run to False, since both dmax, tranformation, prpoints, and skip_first are provided, and these are the numbers that should be determined by the initial fast run\n')
        args.fast_run = False
        if not args.fast_run:
            printt('fast run is now False')

    # the following outlier loop can be used to remove outliers from the data iteratively (max 20 iterations) - if opted for
    CONTINUE_OUTLIER = True
    count_ite,max_ite,Noutlier_prev,outliers_removed = 0,20,1e3,0
    while CONTINUE_OUTLIER:
        count_auto = 0
        #############################################################
        # beginning of auto dmax/transformation/skip_first while loop
        #############################################################
        # run up to 3 times (count_auto) to automatically determine good input params
        while (dmax == '' or transformation == 'A' or skip_first == '' or prpoints == '') and count_auto < 3:
            # write an inputfile to be given to bift executable (bift)
            f = open("inputfile.dat",'w')
            f.write('%s\n' % data)
            f.write('%f\n' % qmin)
            f.write('%f\n' % qmax)
            f.write('%s\n' % args.Bg)
            f.write('200\n') # nrebin data to 200 points in the fast run
            f.write('%s\n' % dmax)
            f.write('\n')
            if args.alpha:
                f.write('%s\n' % args.alpha)
            else:
                f.write('f5\n') # fix log(alpha) to 5 in the fast run if no user input is given
            f.write('%s\n' % args.smear)
            f.write('\n')
            f.write('\n')
            if prpoints == '':
                f.write('70\n') # set pr points (number of points in r of p(r)) in the fast run - if no user input is given
            else:
                f.write('%s\n' % prpoints)
            f.write('\n') # 0 extra error calculations in the fast run - no need for good error estimates
            if transformation == 'A':
                f.write('D\n') # use Debye transformations (positivity constraint), if nothing is opted for by the user
            else:
                f.write('%s\n' % transformation)
            f.write('%s\n' % args.fitbackground)
            f.write('%s\n' % args.rescale_mode) # rescale method. N: non-constant, C: constant
            if args.rescale_mode == 'N':
                f.write('%s\n' % args.nbin)
            else:
                f.write('\n')
            f.write('\n')
            f.close()

            ## run bayesfit (fast run)
            printt("=================================================================================")
            printt("    Fast run for input parameter estimation")
            printt("=================================================================================")
            with open('inputfile.dat', 'r') as input_file:
                result = subprocess.run([os.path.join(os.getcwd(), exe)], stdin=input_file)

            ## estimate dmax from fast run
            if dmax == '':
                ## retrive dmax from parameter file (parameters.dat) using the read_params() function from bayesapp_helpfunctions.py
                dmax_value = read_params()[1]
                dmax = '%s' % dmax_value

            ## set transformation. If thre were negative points in the p(r) from the fast run, set to negative (N)
            if transformation == 'A':
                # A for auto (not selected by the user)
                r,pr,d_pr = np.genfromtxt('pr.dat',skip_header=0,usecols=[0,1,2],unpack=True)
                threshold_pr = np.max(pr)*1E-5
                contain_negative_entries = np.any(pr<-threshold_pr)
                if contain_negative_entries:
                    transformation = 'N' # N for negative (allow negative values)
                    if dmax[0] == 'f':
                        # if dmax if forced to a value by the user by prefix f, do not change dmax, else change it to use input (or default value)
                        pass
                    else:
                        dmax = args.dmax # reset dmax
                        if dmax == '':
                            printt("\n\n\nChanged transformation - new fast run with this tranformation\n\n\n")
                else:
                    transformation = 'D'

            ## set prpoints from dmax in fast run
            # if dmax is large, use more points in pr (slower but needed for large particles)
            if (prpoints == '' and dmax): 
                if dmax[0] == 'f':
                    dmax_aa = float(dmax[1:])
                else:
                    dmax_aa = float(dmax)

                if units == 'nm':
                    dmax_aa *= 10 # convert dmax from nm to angstrom

                threshold_dmax_aa = 180
                if dmax_aa > threshold_dmax_aa:
                    tmp = int(np.amin([dmax_aa/3,150])) # set prpoints to dmax/3, but max 150
                    prpoints = '%d' % tmp
                else:
                    prpoints = '60' # min value of prpoints

            ## determine skip_first
            if (skip_first == '' and dmax):
                # read Rg value from the fast run
                Rg_value = read_params()[2]
                # read data
                try:
                    q,I,dI = np.genfromtxt(data,skip_header=header,skip_footer=footer,usecols=[0,1,2],unpack=True)
                except:
                    q,I,dI = np.genfromtxt(data,skip_header=header,skip_footer=footer,usecols=[0,1,2],unpack=True,encoding='cp855')
                ## remove entries with zero for the error
                idx_nonzero = np.where(dI!=0)
                q,I,dI = q[idx_nonzero],I[idx_nonzero],dI[idx_nonzero]

                ## use Guinier analysis to determine how many first points to skip (if any)
                xx = q**2
                yy = np.log(abs(I))
                dyy = dI/abs(I)
                M = len(q)
                qmax_Guinier = 1.25/Rg_value
                if qmax_Guinier > qmin:
                    last = M-np.where(q<qmax_Guinier)[0][-1]
                    skip_first_max = 25 # never skip more than 25 points
                    n = np.min([len(q[:-last]),skip_first_max+2])
                    skip_first_array = range(0,n)
                    i = 0
                    chi2r_prev = 1e99
                    YELLOW_CARDS = 0
                    while i<n and YELLOW_CARDS<2:
                        s = skip_first_array[i]
                        x = xx[s:-last]
                        y = yy[s:-last]
                        dy = dyy[s:-last]
                        N = len(x)
                        degree = 2
                        c = np.polyfit(x,y,degree)
                        poly2_function = np.poly1d(c)
                        yfit = poly2_function(x)
                        R = (yfit-y)/dy
                        chi2r = np.sum(R**2)/(N-3)
                        chi2r_dif = (chi2r_prev - chi2r)/chi2r
                        # condition for convergence - relative chi2r difference should be above 30%
                        if abs(R[0]) < 1.5 and chi2r_dif<0.3:
                            # 2 yellow cards (red card) is convergence
                            YELLOW_CARDS += 1
                        i +=1
                        chi2r_prev = chi2r
                    try:
                        skip_first = np.max([s-2,0])
                    except:
                        skip_first = 0
                    #check that slope is negative
                    a,b = np.polyfit(xx[skip_first:-last],yy[skip_first:-last],1)
                    if a>=0:
                        skip_first = 0
                        Rg_calc = Rg_value
                    else:
                        Rg_calc = np.sqrt(-3*a)
                    #check qmaxRg is below 1.4  
                    qmaxRg = qmax_Guinier*Rg_value
                    qmaxRg_calc = qmax_Guinier*Rg_calc
                    if qmaxRg > 1.4:
                        skip_first = 0
                    elif qmaxRg_calc > 1.4:
                        skip_first = 0
                else:
                    skip_first = 0
                if skip_first>5:
                    # convert skip_first to new qmin
                    try:
                        q_check = np.genfromtxt(data,skip_header=header+skip_first,skip_footer=footer,usecols=[0],unpack=True)
                    except:
                        q_check = np.genfromtxt(data,skip_header=header+skip_first,skip_footer=footer,usecols=[0],unpack=True,encoding='cp855')
                    if q_check[0] > qmin:
                        qmin = q_check[0]
                    # downscale prpoints (likely smaller dmax)
                    if prpoints: 
                        if float(prpoints) > 70:
                            prpoints = '70'
                    # find new dmax
                    if dmax:
                        if dmax[0] == 'f':
                            pass
                        else:
                            dmax = '' # reset dmax

        ## summarize results from fast run
        printt("\n\n\n")
        printt("    Estimated input parameters from fast run:")
        printt("        dmax:                           %s" % dmax)
        printt("        transformation:                 %s" % transformation)
        printt("        number of first points skipped: %d" % skip_first)
        printt("        qmin:                           %e" % qmin)
        printt("        number of points in p(r):       %s" % prpoints)
        printt("\n\n\n")

        ###############################################################
        ## make actual run with determined (or manually set) values
        ###############################################################
        # if it fails, change the transformation and try again (SECOND TRY), if if fails again, then exit

        CONTINUE_Trans = True
        SECOND_TRY = False
        while CONTINUE_Trans and not args.fast_run:
                
            ## make input file for running bift (accurate run)
            f = open("inputfile.dat",'w')
            f.write('%s\n' % data)
            f.write('%f\n' % qmin)
            f.write('%f\n' % qmax)
            f.write('%s\n' % args.Bg)
            f.write('%s\n' % args.nrebin)
            f.write('%s\n' % dmax)
            f.write('\n')
            f.write('%s\n' % args.alpha)
            f.write('%s\n' % args.smear)
            f.write('\n')
            f.write('\n')
            f.write('%s\n' % prpoints)
            f.write('%s\n' % args.noextracalc)
            f.write('%s\n' % transformation)
            f.write('%s\n' % args.fitbackground)
            f.write('%s\n' % args.rescale_mode) # rescale method. N: non-constant, C: constant, I: intensity-dependent
            if args.rescale_mode == 'N':
                f.write('%s\n' % args.nbin)
            else:
                f.write('\n')
            f.write('\n')
            f.close()

            ## run BIFT
            try:
                os.remove('parameters.dat') # remove parameters file from the initial fast run
            except FileNotFoundError:
                pass # file does not exist → nothing to do
            printt("=================================================================================")
            printt("    Running BayesApp with estimated input parameters")
            printt("=================================================================================") 
            with open('inputfile.dat', 'r') as input_file:
                result = subprocess.run([os.path.join(os.getcwd(), exe)], stdin=input_file)

            ## import params data to check that BIFT was running ok (if not, algorithm will change transformation and try again)
            try:
                dmax_value = read_params()[1] # if there is no parameters.dat file, this will give error
                int(dmax_value) # if dmax_valule is nan, this will give error
                if dmax[0] == 'f':
                    dmax = 'f%f' % dmax_value # set dmax
                else:
                    dmax = '%f' % dmax_value # set dmax
                CONTINUE_Trans = False
            except:
                if SECOND_TRY:
                    # second try also failed - give error, suggestions, and quit
                    CONTINUE_Trans = False
                    printt("=================================================================================")
                    printt("    Could not find a solution. Try changing Maximum distance, Transformation, alpha, Number of points in pr, or Skip first points")
                    printt("=================================================================================")
                    exit()
                elif transformation in ['N','M']:
                    # change transformation from N or M to D (positivity constraint)
                    transformation = 'D'
                    printt("=================================================================================")
                    printt("    No solution - running again with new transformation: %s" % transformation)
                    printt("=================================================================================")
                    SECOND_TRY = True
                elif transformation == 'D':
                    # change transformation from D (positivity constraint) to N (allow negative values)
                    transformation = 'N'
                    printt("=================================================================================")
                    printt("    No solution - running again with new transformation: %s" % transformation)
                    printt("=================================================================================")
                    SECOND_TRY = True
                
        ## import data and fit
        qdat,Idat,sigma = np.genfromtxt('data.dat',skip_header=0,usecols=[0,1,2],unpack=True)
        sigma_rs = np.genfromtxt('rescale.dat',skip_header=3,usecols=[2],unpack=True)
        qfit,Ifit = np.genfromtxt('fit.dat',skip_header=1,usecols=[0,1],unpack=True)

        ## interpolate fit on q-values from data
        Ifit_interp = np.interp(qdat,qfit,Ifit)
        with open('fit_q.dat','w') as f:
            for x,y in zip(qdat,Ifit_interp):
                f.write('%10.10f %10.10f\n' % (x,y))

        ## calculate residuals
        R = (Idat-Ifit_interp)/sigma
        maxR = np.ceil(np.amax(abs(R)))
        R_rs = (Idat-Ifit_interp)/sigma_rs
        maxR_rs = np.ceil(np.amax(abs(R_rs)))

        ## outlier analysis
        # calculate p-value
        x = np.linspace(-10,10,1000)
        pdx = np.exp(-x**2/2)
        norm = np.sum(pdx)
        p = np.zeros(len(R))  
        for i in range(len(R)):
            idx_i = np.where(x>=abs(R[i]))
            p[i] = np.sum(pdx[idx_i])
        p /= norm # normalize value
        p *= len(R) # correction for multiple testing
        idx = np.where(p<0.01) # find statistical outliers
        Noutlier = len(idx[0])
        idx_max = np.argmax(abs(R)) # find the worst outlier
        filename_outlier = 'outlier_filtered.dat'
        if Noutlier:
            with open(filename_outlier,'w') as f:
                f.write('# data, with worst outlier filtered out\n')
                for i in range(len(R))          :
                    if i!=idx_max:
                        f.write('%e %e %e\n' % (qdat[i],Idat[i],sigma[i]))
    
        ## retrive output from parameter file
        I0,dmax_out,Rg,chi2r,background,alpha_out,Ng,Ns,evidence,Prob,Prob_str,assessment,beta,Run_max,Run_max_expect,dRun_max_expect,p_Run_max_str,NR,NR_expect,dNR_expect,p_NR,prpoints_float,d_dmax_out,d_Rg = read_params()

        ## if there are many outliers, then try to gradually (in steps of 50) increase number of points in p(r) and rerun (until prpoints is above 180)
        if Noutlier > 1 and Noutlier < int(Noutlier_prev*0.8) and not args.fast_run:
            Noutlier_prev = Noutlier
            if dmax[0] == 'f':
                dmax = 'f%f' % dmax_out
            else:
                dmax = '%f' % dmax_out # update dmax value
            if prpoints_float < 190:
                printt("=================================================================================")
                printt("    number of outliers: %d" % Noutlier)
                printt("        trying to improve fit by increasing number of points in p(r)")
                printt("        old number of points in p(r): %.0f" % prpoints_float)
                prpoints_float += 50
                printt("        new number of points in p(r): %.0f" % prpoints_float)
                prpoints = '%.0f' % prpoints_float
            else:
                CONTINUE_OUTLIER = False
            printt("=================================================================================")
        # remove worst outlier and run again
        # this is then done iterative untill all outliers are removed
        elif args.outlier_ite and Noutlier:
            data = filename_outlier
            printt("=================================================================================")
            printt("number of outliers: %d" % Noutlier)
            printt("removing worst oulier and rerunning")
            printt("=================================================================================")
            CONTINUE_OUTLIER = Noutlier
            outliers_removed += 1 # count number of removed outliers
        else:
            CONTINUE_OUTLIER = False
        count_ite += 1
        if count_ite >= max_ite:
            CONTINUE_OUTLIER = False
            printt('max iterations in outlier removal reached (=%d). There could be something wrong with error estimates in data' % max_ite)

    ###########################
    # end of oulier while loop 
    ###########################
    
    ## import p(r)
    r,pr,d_pr = np.genfromtxt('pr.dat',skip_header=0,usecols=[0,1,2],unpack=True)

    if args.make_pr_bin:
        ## if opted for, intepolate pr on grid with binsize of pr_binsize
        if units == 'nm':
            pr_binsize /= 10
        r_bin = np.arange(0,r[-1],pr_binsize)
        pr_bin = np.interp(r_bin,r,pr)
        n = len(r)/len(r_bin)
        pr_bin_max = np.interp(r_bin,r,pr+d_pr)
        pr_bin_min = np.interp(r_bin,r,pr-d_pr)
        d_pr_bin = ((pr_bin_max-pr_bin_min)/2)/np.sqrt(n)
        with open('pr_bin.dat','w') as f:
            for x,y,z in zip(r_bin,pr_bin,d_pr_bin):
                f.write('%10.10f %10.10e %10.10e\n' % (x,y,z))

    ###########################
    # plot results!
    ###########################

    ## general plotting settings
    markersize = 4
    linewidth = 1

    ## plot p(r)
    plt.plot(r,np.zeros(len(r)),linestyle='--',color='grey',zorder=0)
    # plt.errorbar(r,pr,yerr=d_pr,marker='.',markersize=markersize,linestyle='None',color='black')
    plt.plot(r,pr,marker='.',markersize=markersize,linewidth=linewidth,color='black')
    if args.make_pr_bin:
        # plt.errorbar(r_bin,pr_bin,d_pr_bin,marker='.',markersize=markersize,linewidth=linewidth,color='green',label='p(r), fixed binsize')
        plt.plot(r_bin,pr_bin,marker='.',markersize=markersize,linewidth=linewidth,color='green',label='p(r), fixed binsize')
        plt.legend(frameon=False)
    # plt.plot(r_dense,pr_dense,linewidth=linewidth,color='black')
    # plt.plot(r_dense,spline2(r_dense),linewidth=linewidth,color='red',label='p(r) 0.1')
    # plt.plot(r_dense,spline3(r_dense),linewidth=linewidth,color='green',label='p(r) 0.2')
    # plt.plot(r_dense,spline4(r_dense),linewidth=linewidth,color='blue',label='p(r) 0.5')
    plt.xlabel(r'$r$ [%s]' % units)
    plt.ylabel(r'$p(r)$')
    plt.title('p(r)')
    # plt.legend()
    plt.tight_layout()
    plt.savefig('pr.png',dpi=200)

    ## plot data, fit and residuals, not rescaled 
    TRUNCATE_ANALYSIS = 0
    if TRUNCATE_ANALYSIS:
        f,(p0,p1,p2) = plt.subplots(3,1,gridspec_kw={'height_ratios': [4,1,5]},sharex=True)
    else:
        f,(p0,p1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},sharex=True)
    p0.errorbar(qdat,Idat,yerr=sigma,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0,label='data')
    if args.logx:
        p0.set_xscale('log')
        p0.plot(qdat,Ifit_interp,color='black',linewidth=linewidth,label='p(r) fit')
    else:
        p0.plot(qfit,Ifit,color='black',linewidth=linewidth,zorder=1,label='p(r) fit') 
    p0.set_ylabel(r'$I(q)$')
    p0.set_yscale('log')
    p0.set_title('p(r) fit to data')

    p1.plot(qdat,R,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0)
    if args.logx:
        p1.set_xscale('log')
        p1.plot(qdat,Idat*0,linewidth=linewidth,color='black',zorder=1)
        if Noutlier:
            p1.plot(qdat,-3*np.ones(len(Idat)),linewidth=linewidth,linestyle='--',color='grey',zorder=2,label=r'$\pm 3\sigma$')
            p1.plot(qdat,3*np.ones(len(Idat)),linewidth=linewidth,linestyle='--',color='grey',zorder=3)
    else:
        p1.plot(qfit,Ifit*0,linewidth=linewidth,color='black',zorder=1)
        if Noutlier:
            p1.plot(qfit,-3*np.ones(len(Ifit)),linewidth=linewidth,linestyle='--',color='grey',zorder=2,label=r'$\pm 3\sigma$')
            p1.plot(qfit,3*np.ones(len(Ifit)),linewidth=linewidth,linestyle='--',color='grey',zorder=3)
    if TRUNCATE_ANALYSIS:
        p2.plot(qdat,Idat/sigma,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0)
        p2.plot(qdat,Idat*0,linewidth=linewidth,color='black',zorder=1)
        p2.plot(qdat,Idat/Idat*2,linewidth=linewidth,linestyle='--',color='grey',zorder=1)
        p2.set_ylim(-5,10)
        #p2.set_yscale('log')
        p2.set_xlabel(r'$q$ [%s$^{-1}$]' % units)
    else:
        p1.set_xlabel(r'$q$ [%s$^{-1}$]' % units)
        
    ## plot outliers
    if Noutlier>1:
        p0.plot(qdat[idx],Idat[idx],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='grey',zorder=4,label='potential outliers')
        p1.plot(qdat[idx],R[idx],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='grey',zorder=4)
        p0.plot(qdat[idx_max],Idat[idx_max],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='black',zorder=4,label='worst outlier')
        p1.plot(qdat[idx_max],R[idx_max],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='black',zorder=4)
    elif Noutlier == 1:
        p0.plot(qdat[idx_max],Idat[idx_max],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='black',zorder=4,label='potential outlier')
        p1.plot(qdat[idx_max],R[idx_max],linestyle='none',marker='o',markerfacecolor='none',markeredgecolor='black',zorder=4)
    p1.set_ylabel(r'$\Delta I(q)/\sigma$')
    try:
        p1.set_ylim(-maxR,maxR)
        if Noutlier:
            p1.set_yticks([-maxR,-3,0,3,maxR])
        else:
            p1.set_yticks([-maxR,0,maxR])
    except:
        printt("    WARNING: Some residuals are either NaN or inf - bad fit?")
        printt("         probably just a numerical instability")
        printt("         try changing the number of points in p(r)")
    p0.legend(frameon=False)
    plt.tight_layout()
    plt.savefig('Iq.png',dpi=200)
    
    ## Guinier analysis
    if args.Guinier:
        try:
            qmaxRg = float(qmaxRg_in)
        except:
            qmaxRg = 1.25
        Rg_Guinier = Rg
        try:
            Guinier_skip = int(Guinier_skip_in) 
        except:
            idx = np.where(qdat*Rg_Guinier<=qmaxRg)
            q2 = qdat[idx]**2
            lnI = np.log(Idat[idx])
            dlnI = sigma[idx]/Idat[idx]
            
            Guinier_skip = 0
            CONTINUE_GUINIER = True
            while CONTINUE_GUINIER:
                try:
                    a,b = np.polyfit(q2[Guinier_skip:],lnI[Guinier_skip:],1,w=1/dlnI[Guinier_skip:])
                    fit = b+a*q2[Guinier_skip:]
                    R = (lnI[Guinier_skip:]-fit)
                    Rmean = np.mean(abs(R))
                    if abs(R[0]) > Rmean*2:
                        # if upturn, skip more points
                        Guinier_skip += 1
                    else:
                        CONTINUE_GUINIER = False
                except:
                    CONTINUE_GUINIER = False
        
        if qdat[Guinier_skip]*Rg<=qmaxRg:
            # determine Rg from Guinier iteratively
            for i in range(7):
                idx = np.where(qdat*Rg_Guinier<=qmaxRg)
                q2 = qdat[idx]**2
                lnI = np.log(Idat[idx])
                dlnI = sigma[idx]/Idat[idx]
                n = len(idx[0])-Guinier_skip
                # forgot the motivation for these lines?:
                #while (Guinier_skip > 0) and (n<10):
                #    Guinier_skip = Guinier_skip-1
                #    n = n+1
                try:
                    if len(q2[Guinier_skip:]) > 4:
                        a,b = np.polyfit(q2[Guinier_skip:],lnI[Guinier_skip:],1,w=1/dlnI[Guinier_skip:])
                        fit = b+a*q2[Guinier_skip:]
                        Rg_Guinier = (Rg_Guinier + np.sqrt(-3*a))/2
                        Error_Guinier = False   
                    else: 
                        Error_Guinier = True
                except: 
                    Error_Guinier = True                 
            if Error_Guinier:
                error_message = '\nERROR in Guinier fit\n - do you have a defined Guinier region?\n - maybe try to skip some of the first points?\n - interparticle interactions may lead to a negative slope at low q\n - SANS contrast match may lead to a negative slope at low q'
                printt(error_message)
                f,(p0,p1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},sharex=True)
                p0.text(0.1,0.7,error_message,transform=p0.transAxes)
                plt.tight_layout()
                plt.savefig('Guinier.png',dpi=200)
            else:
                qmaxRg = np.sqrt(q2[-1])*Rg_Guinier
                R = (lnI[Guinier_skip:]-fit)/dlnI[Guinier_skip:]
                Rmax = np.ceil(np.amax(abs(R)))
                chi2r_Guinier = np.sum(R**2)/(n-2)

                f,(p0,p1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},sharex=True)
                p0.errorbar(q2,lnI,yerr=dlnI,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0)
                p0.plot(q2[Guinier_skip:],fit,color='black',linewidth=linewidth,zorder=1,label=r'$R_g$=%1.2f, $q_{max}R_g$=%1.2f, $\chi^2_r$=%1.1f, skipped_points=%d' % (Rg_Guinier,qmaxRg,chi2r_Guinier,Guinier_skip))
                p1.plot(q2[Guinier_skip:],R,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0)
                p1.plot(q2,q2-q2,color='black',linewidth=linewidth,zorder=1)
                p0.set_ylabel(r'$ln(I)$')
                p1.set_xlabel(r'$q^2$ [%s$^{-2}$]' % units)
                # p1.set_ylabel(r'$\Delta lnI/\sigma_{lnI}$')
                p1.set_ylabel(r'$\Delta lnI\cdot I/\sigma$') # inserted the expression for the uncertainty of lnI
                p1.set_ylim([-Rmax,Rmax])
                p1.set_yticks([-Rmax,0,Rmax])
                p0.set_title('Guinier plot')
                p0.legend(frameon=False)
                plt.tight_layout()
                plt.savefig('Guinier.png',dpi=200)
        else:
            Rg_Guinier = float('nan')
            error_message = '\nERROR in Guinier fit\n - do you have a defined Guinier region?\n - maybe you skipped too many points?\n - maybe your sample is large (>hundreds of nm)?'
            printt(error_message)
            f,(p0,p1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},sharex=True)
            p0.text(0.1,0.7,error_message,transform=p0.transAxes)
            plt.tight_layout()
            plt.savefig('Guinier.png',dpi=200)


    ## Kratky and Mw calculation (assuming pure protein)

    # subtract constant from data (background estimated by BIFT)
    I_sub,I0_sub = Idat-background,I0-background
    
    # calculate Mw from integrating the Kratky plot to get Porod volume and convert to Mw
    qRg = qdat*Rg
    if units == 'nm':
        qdat_aa = qdat*0.1
        Rg_aa = Rg*10
    else:
        qdat_aa = qdat
        Rg_aa = Rg
    qm = np.amin([8.0/Rg_aa,np.amax(qdat_aa)])
    relative_uncertainty = np.max([Rg_aa/300,0.1]) # Ficher et al 2010 J. Appl. Cryst. (2010). 43, 101-109
    idx = np.where(qRg <= 8.0)
    yy = qdat_aa**2*I_sub
    dq_aa = (np.amax(qdat_aa[idx])-np.amin(qdat_aa[idx]))/len(idx[0])
    Qt = np.sum(yy[idx])*dq_aa # scattering invariant
    Vt = 2*np.pi**2*I0_sub/Qt # A^3
    Vt_nm = Vt/1000
    MwP = 0.625 * Vt_nm # Petoukhov et al 2012, 0.625 kDa/nm3
    
    # Piiadov et al 2018 Protein Science  https://doi.org/10.1002/pro.3528
    qm2,qm3,qm4 = qm**2,qm**3,qm**4
    A = -2.114e6*qm**4 + 2.920e6*qm3 - 1.472e6*qm2 + 3.349e5*qm - 3.577e4
    B =                  12.09*qm3   - 9.39*qm2    + 3.03*qm    + 0.29
    Vm = A+B*Vt # A^3
    Vm_nm = Vm/1000 # nm^3
    MwF = 0.83 * Vm_nm # Squire and Himmel 1979, 0.83 kDa/nm3 --> 0.83/1000 kDa/A3
    dMwF = MwF * relative_uncertainty

    ## Kratky
    if args.Kratky or args.Kratky_dim:
        # Mw in label?
        if args.Kratky_Mw:
            label = 'Mw = %1.1f +/- %1.1fkDa' % (MwF,dMwF)
        else:
            label = ''

        # Kratky plot
        plt.figure()
        if args.Kratky_dim:
            # dimensionless Kratky plot
            x, y, dy = qRg, qRg**2*I_sub/I0_sub, qRg*qRg*sigma/I0_sub
            plt.ylabel(r'$I/I(0) (q R_G)^2$')
            plt.xlabel(r'$q R_g$')
        else:
            # standard Kratky plot
            x, y, dy = qdat, qdat**2*I_sub, qdat*qdat*sigma
            plt.ylabel(r'$I q^2$')
            plt.xlabel(r'$q$ [%s$^{-1}$]' % units)
        plt.errorbar(x,y,yerr=dy,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0,label=label)
        plt.plot(x,np.zeros_like(x),linestyle='--',color='grey',zorder=1) 
        if label: 
            plt.legend(frameon=False)
        plt.title('Kratky plot')
        plt.tight_layout()
        plt.savefig('Kratky.png',dpi=200)
    
    ## Porod
    if args.Porod:
        # transform data
        y = qdat**4 * (Idat - background)
        dy = qdat**4 * sigma

        # constant fit at Porod limit
        if Porod_limit:
            qm_Porod = Porod_limit
        else:
            useful_qmax = np.pi*Ng/dmax_out
            qm_Porod = useful_qmax*0.95
        if np.amax(qdat) <= qm_Porod:
            qm_Porod = 0.9*np.amax(qdat)
        idx = np.where(qdat>qm_Porod)
        a = np.polyfit(qdat[idx],y[idx],0,w=1/dy[idx])

        # Porod plot
        f,p0 = plt.subplots(1,1)
        p0.errorbar(qdat,y,yerr=dy,linestyle='none',marker='.',markersize=markersize,color='red',zorder=0)
        p0.plot([qm_Porod,qm_Porod],[np.amin(y),np.amax(y)],color='grey',linestyle='--')
        p0.plot(qdat[idx],qdat[idx]/qdat[idx]*a,color='black',label='fit with constant')
        p0.set_title('Porod plot')
        p0.set_ylabel(r'$I q^4$')
        p0.set_xlabel(r'$q$ [%s$^{-1}$]' % units)
        p0.legend(frameon=False)
        plt.tight_layout()
        plt.savefig('Porod.png',dpi=200)

    ## import and plot data with rescaled errors
    if Prob < 0.003:
        f,(p0,p1) = plt.subplots(2,1,gridspec_kw={'height_ratios': [4,1]},sharex=True)

        # rescaled data
        p0.errorbar(qdat,Idat,yerr=sigma_rs,linestyle='none',marker='.',markersize=markersize,color='blue',zorder=0,label='data with rescaled errors')
        if args.logx:
            p0.set_xscale('log')
            p0.plot(qdat,Ifit_interp,color='black',linewidth=linewidth,zorder=1,label='p(r) fit')
        else:
            p0.plot(qfit,Ifit,color='black',linewidth=linewidth,zorder=1,label='p(r) fit')
        p0.set_ylabel(r'$I(q)$')
        p0.set_yscale('log')
        p0.set_title('p(r) fit to data with rescaled errors')
        p0.legend(frameon=False)

        # residuals
        p1.plot(qdat,R_rs,linestyle='none',marker='.',markersize=markersize,color='blue',zorder=0)
        if args.logx:
            p1.set_xscale('log')
            p1.plot(qdat,np.zeros_like(qdat),linewidth=linewidth,color='black',zorder=1)
        else:
            p1.plot(qfit,np.zeros_like(qfit),linewidth=linewidth,color='black',zorder=1)
        p1.set_xlabel(r'$q$ [%s$^{-1}$]' % units)
        p1.set_ylabel(r'$\Delta I(q)/\sigma_\mathrm{rescale}$')
        try:
            p1.set_ylim(-maxR_rs,maxR_rs)
            p1.set_yticks([-maxR_rs,0,maxR_rs])
        except:
            printt("WARNING: Some residuals are either NaN or inf - bad fit?")
            printt("         probably just a numerical instability")
            printt("         try changing the number of points in p(r)")
        plt.tight_layout()
        plt.savefig('Iq_rs.png',dpi=200)
        
    ## output values
    printt("\n\n\n")
    printt("    Estimated parameters after runnning BayesApp:")
    # printt("        dmax:                       %s" % dmax)
    printt("        dmax:                       %.4f +- %.4f %s" % (dmax_out,d_dmax_out,units))
    printt("        Rg from p(r):               %.4f +- %.4f %s" % (Rg,d_Rg,units))
    if args.Guinier:
        printt("        Rg from Guinier analysis:   %.4f" % Rg_Guinier)
    # if args.Kratky_Mw:
    printt("        Mw (assuming protein):      %1.1f +/- %1.1f kDa" % (MwF,dMwF))
    printt("        transformation:             %s" % transformation)
    printt("        skip first points:          %d" % skip_first)
    printt("        number of points in p(r):   %s" % prpoints)
    printt("        number of outliers:         %d" % Noutlier)
    printt("        number of outliers removed: %d" % outliers_removed)

    ## export SASBDB-readable p(r) file
    with open('pr.dat','w') as f:
        f.write('# p(r) from BayesApp version %s\n' % version)
        f.write('# data = %s\n' % data)
        f.write('# Radius of gyration from p(r), Rg = %.4f +- %.4f %s\n' % (Rg,d_Rg,units))
        if args.Guinier:
            f.write('# Radius of gyration from Guinier, Rg = %.4f %s\n' % (Rg_Guinier,units))
        else:
            f.write('# Radius of gyration from Guinier, Rg = %.4f %s\n' % (float('nan'),units))
        f.write('# Maximum distance, dmax = %.4f +- %.4f %s\n' % (dmax_out,d_dmax_out,units))
        if units == 'A':
            f.write('# apparent Porod Volume, Vp = %.2f %s^3\n' % (Vt,units))
            f.write('# refined particle Volume, Vm = %.2f %s^3\n' % (Vm,units))
        elif units == 'nm':
            f.write('# truncatated Porod Volume, Vp = %.2f %s^3\n' % (Vt_nm,units))
            f.write('# estimated particle Volume, Vm = %.2f %s^3\n' % (Vm_nm,units))
        else:
            printt('ERROR: no units defined')
            exit()
        f.write('# Molecular weight, Mw = %.2f +- %.2f kDa\n' % (MwF,dMwF))
        f.write('#\n')
        f.write('# Notes and references:\n')
        f.write('# Vp is a optimized for Mw estimation, not a best estimate of particle volume (Petoukhov et al., 2007, DOI 10.1107/S0021889812007662)\n')
        f.write('# Vm is a refined estimate of the particle volume (Piiadov et al., 2018, DOI 10.1002/pro.3528)\n')
        f.write('# Mw was calculated from Vm assuming protein density 83 kDa/nm3 (Squire and Himmel, 1979, DOI: 10.1016/0003-9861(79)90563-0, PMID: 507801)\n')
        f.write('# Vp(in nm^3)*0.625 kDa/nm3 gives a Mw of %.2f kDa +- 20%s (Petoukhov et al., 2007, DOI 10.1107/S0021889812007662)\n' % (MwP,'%'))
        f.write('# If BayesApp was useful for you work, please cite the most recent publication: Larsen and Pedersen, 2021, DOI 10.1107/S1600576721006877\n')
        f.write('#\n')
        f.write('# r p(r) sigma_p(r)\n')
        for r_i,pr_i,d_pr_i in zip(r,pr,d_pr):
            f.write('%f %f %f\n' % (r_i,pr_i,d_pr_i))
    
    ### end timing
    end_time = time.time()-start_time
    printt('\n    total time:                   %0.1f seconds' % end_time)
    printt("\n\n\n")

    ## compress output files to zip file
    if args.zip_compress:
        import zipfile
        import glob 
        # zip_filename = f'results_{prefix}.zip'
        zip_filename = f'results_{data}.zip'
        printt('\n    compressing output to zip file: %s' % zip_filename)
        # files_to_zip = ['filename', 'bift.f', exe, 'pr.dat', 'pr_bin.dat', 'pr_smooth.dat','data.dat', 'fit.dat', 'fit_q.dat', 'parameters.dat', 'rescale.dat','outlier_filtered.dat', 'scale_factor.dat', 'bayesapp.log', 'inputfile.dat']
        files_to_zip = ['filename', 'bift.f', exe, 'pr.dat', 'pr_bin.dat','data.dat', 'fit.dat', 'fit_q.dat', 'parameters.dat', 'rescale.dat','outlier_filtered.dat', 'scale_factor.dat', 'bayesapp.log', 'inputfile.dat']
        files_to_zip.extend(glob.glob('*.png')) # add png images
        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file in files_to_zip:
                if os.path.exists(file):
                    zipf.write(file)
    
    ## show plots
    if args.show:
        plt.show()
