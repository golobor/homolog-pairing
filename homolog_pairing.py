#!/usr/bin/env python
import warnings                                                                  
import numpy as np                                                               
import scipy.stats as st                                                         
import pandas as pd
    
import cooler                                                                    
                                                                                 
import cooltools                                                                 
import cooltools.num                                                             
from cooltools.num import numutils                                               

import click

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('COOLER_PATH', nargs=1, type=click.Path(exists=True))
@click.option(
    '-o',
    '--output',
    type=click.File('w'),
    default='-',
    show_default=True,
    help='path to the output bedgraph file (by default, printed into stdout)')
@click.option(
    '--window-bp',
    type=int, 
    default=20000, 
    show_default=True,
    help='window size to aggregate the contact frequency between homologous loci')
@click.option(
    '--homolog-suffixes', 
    nargs=2, 
    default = ('_057', '_439'),    
    show_default=True,
    help='suffixes used to distinguish the pair of homolog chromosomes')

@click.option(
    '--ignore-diags',
    type=int, 
    default=0, 
    show_default=True,
    help='the number of diagonals of the Hi-C matrix to ignore due to '
         'contamination by short-distance Hi-C ')
@click.option(
    '--poisson-perc',
    type=float, 
    default=None, 
    show_default=True,
    help='if provided, report the N-th percentile of the sampling distribution '
        'of the pairing score at each location, assuming the trans-homolog '
        'interaction counts are Poisson-distributed')
@click.option(
    '--balance',
    type=bool, 
    default=True, 
    show_default=True,
    help='use interaction frequencies normalized by iterative correction')
@click.option(
    '--normalize-by-cis',
    type=bool, 
    default=True, 
    show_default=True,
    help='divide the pairing score by the cis interaction frequency in the corresponding window')
@click.option(
    '--normalize-by-median',
    type=bool, 
    default=True, 
    show_default=True,
    help='divide the pairing score by its genome-wide median')
@click.option(
    '--transform',
    type=click.Choice(['linear', 'log2', 'log10']),
    default='log2', 
    show_default=True,
    help='report either log-transformed or non-transformed (linear) pairing scores')
@click.option(
    '--report-per-homolog',
    type=bool,
    default=True,
    show_default=True,
    help='if True, duplicate the pairing scores for each of the homologs '
         'otherwise, report the pairing score once for each homolog pair')

def make_pairing_bedgraph(
    cooler_path,
    output,
    window_bp,
    homolog_suffixes,
    ignore_diags,
    poisson_perc,
    balance,
    normalize_by_cis,
    normalize_by_median,
    transform,
    report_per_homolog, 
):
    """Generate a bedgraph file with a homolog pairing score.
    """
    clr = cooler.Cooler(cooler_path)
    
    df = get_homolog_pairing_score(
        clr,
        window_bp,
        homolog_suffixes,    
        ignore_diags,
        poisson_perc,
        balance,
        normalize_by_cis,
        report_per_homolog,
        )

    if transform == 'log2':                                                      
        df['pairing'] = np.log2(df['pairing'])
        df['pairing'][~np.isfinite(df['pairing'])] = np.nan
    if transform == 'log10':                                                     
        df['pairing'] = np.log10(df['pairing'])
        df['pairing'][~np.isfinite(df['pairing'])] = np.nan

    gw_median_signal = np.nanmedian(df['pairing'])
    if normalize_by_median:                                                          
        if transform.startswith('log'): 
            df['pairing'] -= gw_median_signal
        else:
            df['pairing'] /= gw_median_signal

    df['pairing'] = np.nan_to_num(df['pairing'])
    #df.rename({'pairing':'dataValue', 'start':'chromStart', 'end':'chromEnd'})
        
    bedgraph_header = (
        """track type=bedGraph
        name="{transform}pairing{percentile}"
        description="{transform}pairing{percentile}, {window_bp}bp window at {res}bp resolution{norm_by_cis}{normalize_by_median}{ignore_diags}. Summary: median: {median}"
        visibility=full color=0,0,0 altColor=50,50,50 graphType=points yLineMark=0 yLineOnOff=on
        group="{transform}pairing"
        """.replace('\n', ' ').format(
            transform = (transform +'_') if transform else '',
            percentile = ('_'+str(poisson_perc)+'th_perc') if poisson_perc else '',
            window_bp=window_bp,
            res = clr.info['bin-size'],
            norm_by_cis = ', normalized by cis diagonal' if normalize_by_cis else '',
            normalize_by_median = ', normalized by the GW-median' if normalize_by_median else '',
            ignore_diags = ', ignored {} diagonals'.format(ignore_diags) if ignore_diags else '',
            median = gw_median_signal
            ) + '\n')
    
    output.write(bedgraph_header)           

    df.to_csv(
        output,
        sep='\t',
        header=None,
        index=False)
        
                                                                                 
def get_homolog_pairing_score(
    clr,
    window_bp=10000,
    homolog_suffixes = ('_057', '_439'),    
    ignore_diags = 0,
    poisson_perc = None,
    balance = True,
    normalize_by_cis = True,
    report_per_homolog = True):

    pairing_dfs = []                              

    hom_chroms = clr.chromnames                                                        
    base_chroms = []   
    for chrom in hom_chroms:                                                         
        base_chrom = ''
        
        if chrom.endswith(homolog_suffixes[0]):                                          
            base_chrom = chrom[:-len(homolog_suffixes[0])]
        if chrom.endswith(homolog_suffixes[1]):                                          
            base_chrom = chrom[:-len(homolog_suffixes[1])]
        if (                                                                         
            ((base_chrom+homolog_suffixes[0]) in hom_chroms)                                  
            and ((base_chrom+homolog_suffixes[1]) in hom_chroms)                              
            and base_chrom not in base_chroms):                                      

            base_chroms.append(base_chrom)                                           

    window = int(np.floor(window_bp / 2 / clr.info['bin-size']))
    bins = clr.bins()[:]
    for chrom in base_chroms:                                                        

        if (not balance) or poisson_perc is not None:
            trans_mat_raw = clr.matrix(balance=False).fetch(
                chrom + homolog_suffixes[0], chrom+homolog_suffixes[1]).astype(np.float)
            trans_diag_raw = _take_big_diagonal_pixel(
                trans_mat_raw, window, ignore_diags, agg_func=np.nansum)
        if balance:
            trans_mat_balanced = clr.matrix(balance=True).fetch(
                chrom + homolog_suffixes[0], chrom+homolog_suffixes[1]).astype(np.float)
            trans_diag_balanced = _take_big_diagonal_pixel(
                trans_mat_balanced, window, ignore_diags, agg_func=np.nansum)
           
        trans_diag = trans_diag_balanced if balance else trans_diag_raw
        
        if poisson_perc is not None:
            trans_diag *= st.poisson.isf(poisson_perc/100, trans_diag_raw) / trans_diag_raw

        if normalize_by_cis:
            cis_hom1_mat = clr.matrix(balance=balance).fetch(
                chrom + homolog_suffixes[0], chrom+homolog_suffixes[0]).astype(np.float)
            cis_hom2_mat = clr.matrix(balance=balance).fetch(
                chrom + homolog_suffixes[1], chrom+homolog_suffixes[1]).astype(np.float)
            cis_hom1_diag = _take_big_diagonal_pixel(
                cis_hom1_mat, window, ignore_diags, agg_func=np.nansum)
            cis_hom2_diag = _take_big_diagonal_pixel(
                cis_hom2_mat, window, ignore_diags, agg_func=np.nansum)

            with warnings.catch_warnings():                                              
                warnings.simplefilter("ignore")                                          
                pairing = trans_diag / ((cis_hom1_diag + cis_hom2_diag)/2)                
        else:                                                                        
            pairing = trans_diag                                                     

        chrom_bins = bins[bins.chrom == chrom+homolog_suffixes[0]]                       
        pairing_df = pd.DataFrame(dict(binid=chrom_bins.index.values))
        pairing_df['chrom'] = chrom
        pairing_df['start'] = chrom_bins.start.values                             
        pairing_df['end'] = chrom_bins.end.values                                 
        pairing_df['pairing'] = pairing

        if report_per_homolog:
            pairing_df_1 = pairing_df.copy()
            pairing_df_2 = pairing_df.copy()
            pairing_df_1['chrom'] = chrom + homolog_suffixes[0]
            pairing_df_2['chrom'] = chrom + homolog_suffixes[1]
            pairing_df_1['binid'] = bins[
                bins['chrom'] == chrom + homolog_suffixes[0]].index.values
            pairing_df_2['binid'] = bins[
                bins['chrom'] == chrom + homolog_suffixes[1]].index.values

            if (np.searchsorted(bins['chrom'], chrom + homolog_suffixes[0], 'left') <
                np.searchsorted(bins['chrom'], chrom + homolog_suffixes[1], 'left')):

                pairing_dfs.append(pairing_df_1)                                            
                pairing_dfs.append(pairing_df_2)                                            
            else:
                pairing_dfs.append(pairing_df_2)                                            
                pairing_dfs.append(pairing_df_1)                                            
        else:
            pairing_dfs.append(pairing_df)


    pairing_dfs = pd.concat(pairing_dfs, ignore_index=True)
    pairing_dfs.sort_values('binid', inplace=True)
    pairing_dfs.drop('binid', axis=1, inplace=True)

    return pairing_dfs                                                                 

def _take_big_diagonal_pixel(
    mat, pad_bins=10, ignore_diags=3, agg_func =np.nanmean):
    """                                                                          
    """                                                                          
    if (ignore_diags):                                                           
        mat = mat.copy()                                                         
        for i in range(-ignore_diags, ignore_diags+1):                           
            numutils.set_diag(mat, np.nan, i)                                    
                                                                                 
    with warnings.catch_warnings():                                              
        warnings.simplefilter("ignore")                                          
                                                                                 
        N = mat.shape[0]                                                         
        score = np.nan * np.ones(N)                                              
        for i in range(0, N):                                                    
            lo = max(0, i-pad_bins)                            
            hi = min(i+pad_bins+1, N)                                              
            # nanmean of interactions to reduce the effect of bad bins           
            score[i] = agg_func(mat[lo:hi, lo:hi])                               
                                                                                 
    return score                                                                 

if __name__ == '__main__':
    make_pairing_bedgraph()

