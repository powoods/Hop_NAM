#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:50:46 2024

@author: patrick.woods
"""
#v6#

def est_nam_H2_means(pheno, rep=False, get_emmeans = False, get_binary_means = False, qn_transform = False):
    """
    Uses an Input Phenotype Data File To Estimate the Broad Sense Heritability for Colony Counts or Estimated Marginal Means.


    Parameters
    ----------
    pheno : .csv file containing Bine 1 and Bine 2 colony counts
    rep : boolean. Default is False. Assigning to True will include Rep in the model when calculating heritability.
    get_emmeans : boolean. Default is False. Assigning to True will cause function to return a pandas dataframe with EMMEANS.
    get_binary_means : boolean. Default is False. Assigning to True will cause funcion to return a pandas dataframe with EMMEANS converted to the binary 'resistant' or 'susceptible' phenotype scale.
    qn_transform : boolean. Default is False. Assigning to True will cause function to set mean colony counts less than 1 to missing and quantile normalize all mean colony counts greater than or equal to 1.
    Returns
    -------
    Broad sense heritability as a percentage of total phenotypic variance.
    
    OR
    
    Estimated Marginal Means (Emmeans).
    
    OR
    
    Quantile Estimated Marginal Means (Emmeans_qn).
    
    OR
    
    Binary Estimated Marginal Means (0, 1).

    """
    import pandas as pd  # package for tabular dataframes
    import pymer4 as pymer
    import numpy as np  # package for working with complex data arrays
    import matplotlib.pyplot as plt  # package for plotting
    from pymer4.models import Lmer  # package for fitting linear mixed models
    import rpy2.robjects as robjects #needed to access R functionality
    import rpy2.robjects.packages as rpackages #needed to access R functionality
    from rpy2.robjects.vectors import StrVector #needed to access R functionality
    from rpy2.robjects import pandas2ri #needed to convert R objects to Pandas objects
    pandas2ri.activate() #activates the conversion functionality of the previously imported code
    if not rep and get_emmeans == False:
        print('Calculating heritability without rep included in the model...')
    
    elif rep == True and get_emmeans == False:
        print('Calculating heritability with rep included in the model...')
        
    elif rep == False and get_emmeans == True and qn_transform == True:
        print('Calculating quantile normalized emmeans without rep included in the model...')
        
    elif rep == True and get_emmeans == True and qn_transform == True:
        print('Calculating quantile normalized emmeans with rep included in the model...')
    
    elif rep == True and get_emmeans == True and get_binary_means == True:
         print('Calculating Binary Emmeans with rep included in the model...')   
        
    elif rep == False and get_emmeans==True and get_binary_means == True:
        print('Calculating Binary Emmeans without rep included in the model...')
    
    elif rep == False and get_emmeans == True:
        print('Calculating Emmeans without rep included in the model...')
    
    elif rep and get_emmeans == True:
        print('Calculating Emmeans with rep included in the model...')
    
    else:
        print('Function not provided enough parameters.')
    

    # reading in the raw phenotype data for the NAM population
    nam = pd.read_csv(pheno)

    ### Extract and Analyze Bench A ###
    # Creating a conditional filter to match all rows where the Bench column = A
    bencha_cond = nam['Bench'] == 'A'
    # Using the conditional filter to extract only Bench A rows
    bencha = nam[bencha_cond]

    # Aggregating all samples by their Plant_ID factor
    ga = bencha.groupby(['Plant_ID'])
    # Calculating the simple means for Bine 1 across all Plant_ID
    bine1a_mean = ga[['Bine 1 Colony counts']].mean()
    bine1a_mean.rename(columns={'Bine 1 Colony counts': 'mean_counts'}, inplace=True)  # Renaming the Bine 1 Colony counts column to mean_counts
    
    bine1a_mean = bine1a_mean.reset_index()
    
    # Calculating the simple means for Bine 2 across all Plant_ID
    bine2a_mean = ga[['Bine 2 Colony counts']].mean()
    bine2a_mean.rename(columns={'Bine 2 Colony counts': 'mean_counts'}, inplace=True)  # Renaming the Bine 2 Colony counts column to mean_counts (needs to match the previously renamed column)
    
    bine2a_mean = bine2a_mean.reset_index()

    # rbinding Bine 1 and Bine 2 means together
    bine1a_2a_means = pd.concat([bine1a_mean, bine2a_mean])

    # Aggregating by common Plant_ID
    g1a_2a = bine1a_2a_means.groupby(['Plant_ID'])
    # Calculating the simple means across Bines 1 and 2
    bencha_overall_mean = g1a_2a[['mean_counts']].mean()

    bencha_overall_mean['Bench'] = 'A'
    bencha_overall_mean['Bench'] = bencha_overall_mean['Bench'].astype('category')
    
    bencha_overall_mean = bencha_overall_mean.reset_index()
    
    ### Extract and Analyze Bench B ###
    # Creating a conditional filter to match all rows where the Bench column = A
    benchb_cond = nam['Bench'] == 'B'
    # Using the conditional filter to extract only Bench A rows
    benchb = nam[benchb_cond]

    # Aggregating all samples by their Plant_ID factor
    gb = benchb.groupby(['Plant_ID'])
    # Calculating the simple means for Bine 1 across all Plant_ID
    bine1b_mean = gb[['Bine 1 Colony counts']].mean()
    bine1b_mean.rename(columns={'Bine 1 Colony counts': 'mean_counts'}, inplace=True)  # Renaming the Bine 1 Colony counts column to mean_counts

    bine1b_mean = bine1b_mean.reset_index()

    # Calculating the simple means for Bine 2 across all Plant_ID
    bine2b_mean = gb[['Bine 2 Colony counts']].mean()
    bine2b_mean.rename(columns={'Bine 2 Colony counts': 'mean_counts'}, inplace=True)  # Renaming the Bine 2 Colony counts column to mean_counts (needs to match the previously renamed column)

    bine2b_mean = bine2b_mean.reset_index()

    # rbinding Bine 1 and Bine 2 means together
    bine1b_2b_means = pd.concat([bine1b_mean, bine2b_mean])

    # Aggregating by common Plant_ID
    g1b_2b = bine1b_2b_means.groupby(['Plant_ID'])
    # Calculating the simple means across Bines 1 and 2
    benchb_overall_mean = g1b_2b[['mean_counts']].mean()

    benchb_overall_mean['Bench'] = 'B'
    benchb_overall_mean['Bench'] = benchb_overall_mean['Bench'].astype('category')

    benchb_overall_mean = benchb_overall_mean.reset_index()

    ### Combine Bench A and Bench B into a Single DataFrame for Model Building ###
    means_for_heritability = pd.concat([bencha_overall_mean, benchb_overall_mean])
    means_for_heritability[['Plant_ID', 'Bench']] = means_for_heritability[['Plant_ID', 'Bench']].astype('category')
    means_for_heritability = means_for_heritability.reset_index()
    ### Build Model Based on User Defined Condition ###
    emmeans = rpackages.importr('emmeans')
    lme4 = rpackages.importr('lme4')
    nlme = rpackages.importr('nlme')
    base = rpackages.importr('base')
    stats = rpackages.importr('stats')
    bn = rpackages.importr('bestNormalize')
    set_column = getattr(base, '$<-')


    if rep == False:
        mod_r = lme4.lmer('mean_counts ~ (1|Plant_ID)',data=means_for_heritability)
        # extrating variance components and converting them to R dataframe
        var_comp = base.as_data_frame(nlme.VarCorr(mod_r))
        pd_df = pandas2ri.rpy2py_dataframe(var_comp)
        G = pd_df.iloc[0, 3]
        error = pd_df.iloc[1, 3]
        total_var = G + error
        heritability = G/total_var * 100

    else:
        mod_r = lme4.lmer('mean_counts ~ (1|Plant_ID) + (1|Bench)',
                        data=means_for_heritability)
        # extrating variance components and converting them to R dataframe
        var_comp = base.as_data_frame(nlme.VarCorr(mod_r))
        pd_df = pandas2ri.rpy2py_dataframe(var_comp)
        G = pd_df.iloc[0, 3]
        Rep = pd_df.iloc[1, 3]
        error = pd_df.iloc[2, 3]
        total_var = G + Rep + error
        heritability = G/total_var * 100

    if get_emmeans == False:
        return print('The broad sense heritability is: ' + str(round(heritability, 2)) + '%')
    
    else:
        if get_binary_means == False:
            if rep == False:
                mod_m_no_rep = stats.lm('mean_counts ~ Plant_ID', data = means_for_heritability)
                em_out_no_rep = emmeans.emmeans(mod_m_no_rep, 'Plant_ID')
                em_df_no_rep = base.as_data_frame(em_out_no_rep)
                em_df_no_rep[0] = base.as_character(em_df_no_rep[0])
                em_df_no_rep_pd = pandas2ri.rpy2py_dataframe(em_df_no_rep)
                em_df_no_rep_pd.loc[em_df_no_rep_pd['emmean'] < 0, 'emmean'] = 0
                em_df_no_rep_pd['emmean'] = round(em_df_no_rep_pd['emmean'], 2) #rounding the emmeans to two decimal places
                if qn_transform == False:
                
                    return em_df_no_rep_pd
                else:
                    em_df_no_rep_pd.loc[em_df_no_rep_pd['emmean'] < 1, 'emmean'] = np.NaN
                    em_df_no_rep_r = robjects.conversion.py2rpy(em_df_no_rep_pd)
                    qn = bn.orderNorm(em_df_no_rep_r[1])
                    transformed_df_no_rep = base.as_data_frame(qn[0])
                    transformed_df_no_rep = set_column(transformed_df_no_rep, 'Plant_ID', em_df_no_rep_r[0])
                    transformed_df_no_rep_pd = pandas2ri.rpy2py_dataframe(transformed_df_no_rep)
                    transformed_df_no_rep_pd.rename(columns={transformed_df_no_rep_pd.columns[0]: 'mean_counts_qn'}, inplace=True)
                    
                    return transformed_df_no_rep_pd
            else:
                mod_m_rep = lme4.lmer('mean_counts ~ Plant_ID + (1|Bench)', data = means_for_heritability)
                em_out_rep = emmeans.emmeans(mod_m_rep, 'Plant_ID')
                em_df_rep = base.as_data_frame(em_out_rep)
                em_df_rep[0] = base.as_character(em_df_rep[0])
                em_df_rep_pd = pandas2ri.rpy2py_dataframe(em_df_rep)
                em_df_rep_pd.loc[em_df_rep_pd['emmean'] < 0, 'emmean'] = 0
                em_df_rep_pd['emmean'] = round(em_df_rep_pd['emmean'],2) #rounding the emmeans to two decimal places
                if qn_transform == False:
                
                    return em_df_rep_pd
                else:
                    em_df_rep_pd.loc[em_df_rep_pd['emmean'] < 1, 'emmean'] = np.NaN
                    em_df_rep_r = robjects.conversion.py2rpy(em_df_rep_pd)
                    qn = bn.orderNorm(em_df_rep_r[1])
                    transformed_df_rep = base.as_data_frame(qn[0])
                    transformed_df_rep = set_column(transformed_df_rep, 'Plant_ID', em_df_rep_r[0])
                    transformed_df_rep_pd = pandas2ri.rpy2py_dataframe(transformed_df_rep)
                    transformed_df_rep_pd.rename(columns={transformed_df_rep_pd.columns[0]: 'mean_counts_qn'}, inplace=True)
                
                    return transformed_df_rep_pd
        else:
            if rep == False:
                mod_m_no_rep = stats.lm('mean_counts ~ Plant_ID', data = means_for_heritability)
                em_out_no_rep = emmeans.emmeans(mod_m_no_rep, 'Plant_ID')
                em_df_no_rep = base.as_data_frame(em_out_no_rep)
                em_df_no_rep[0] = base.as_character(em_df_no_rep[0])
                em_df_no_rep_pd_bin = pandas2ri.rpy2py_dataframe(em_df_no_rep)
                em_df_no_rep_pd_bin['emmean'] = round(em_df_no_rep_pd_bin['emmean'], 2)
                
                em_df_no_rep_pd_bin.loc[em_df_no_rep_pd_bin['emmean'] >= 1, 'emmean'] = 1 
                em_df_no_rep_pd_bin.loc[em_df_no_rep_pd_bin['emmean'] < 1, 'emmean'] = 0 
                
                return em_df_no_rep_pd_bin
            else:
                mod_m_rep = lme4.lmer('mean_counts ~ Plant_ID + (1|Bench)', data = means_for_heritability)
                em_out_rep = emmeans.emmeans(mod_m_rep, 'Plant_ID')
                em_df_rep = base.as_data_frame(em_out_rep)
                em_df_rep[0] = base.as_character(em_df_rep[0])
                em_df_rep_pd_bin = pandas2ri.rpy2py_dataframe(em_df_rep)
                em_df_rep_pd_bin['emmean'] = round(em_df_rep_pd_bin['emmean'],2)
                
                em_df_rep_pd_bin.loc[em_df_rep_pd_bin['emmean'] >= 1, 'emmean'] = 1
                em_df_rep_pd_bin.loc[em_df_rep_pd_bin['emmean'] < 1, 'emmean'] = 0
                
                return em_df_rep_pd_bin