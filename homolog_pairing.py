import warnings                                                                  
import numpy as np                                                               
import scipy.stats as st                                                         
import pandas as pd
    
import cooler                                                                    
                                                                                 
import cooltools                                                                 
import cooltools.num                                                             
from cooltools.num import numutils                                               
                                                                                                                                                                                                                                  
def _take_big_diagonal_pixel(mat, pad_bins=10, ignore_diags=3, agg_func =np.nanmean):
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
                                                                                 
                                                                                 
def get_homolog_pairing_score(
    clr,
    window_bp=10000,
    hom_suffixes = ('_057', '_439'),    
    ignore_diags = 0,
    balance = True,
    poisson_perc = None,
    normalize_by_cis = True,
):
    
    pairing_dfs = []                                                                   

    hom_chroms = clr.chromnames                                                        
    base_chroms = []   
    for chrom in hom_chroms:                                                         
        base_chrom = ''
        
        if chrom.endswith(hom_suffixes[0]):                                          
            base_chrom = chrom[:-len(hom_suffixes[0])]
        if chrom.endswith(hom_suffixes[1]):                                          
            base_chrom = chrom[:-len(hom_suffixes[1])]
        print(chrom, base_chrom)
        if (                                                                         
            ((base_chrom+hom_suffixes[0]) in hom_chroms)                                  
            and ((base_chrom+hom_suffixes[1]) in hom_chroms)                              
            and base_chrom not in base_chroms):                                      

            base_chroms.append(base_chrom)                                           

    window = int(np.floor(window_bp / 2 / clr.info['bin-size']))
    bins = clr.bins()[:]
    for chrom in base_chroms:                                                        
        chrom_bins = bins[bins.chrom == chrom+hom_suffixes[0]]                       

        if (not balance) or poisson_perc is not None:
            trans_mat_raw = clr.matrix(balance=False).fetch(chrom + hom_suffixes[0], chrom+hom_suffixes[1]).astype(np.float)
            trans_diag_raw = _take_big_diagonal_pixel(trans_mat_raw, window, ignore_diags, agg_func=np.nansum)
        if balance:
            trans_mat_balanced = clr.matrix(balance=True).fetch(chrom + hom_suffixes[0], chrom+hom_suffixes[1]).astype(np.float)
            trans_diag_balanced = _take_big_diagonal_pixel(trans_mat_balanced, window, ignore_diags, agg_func=np.nansum)
           
        trans_diag = trans_diag_balanced if balance else trans_diag_raw
        
        if poisson_perc is not None:
            trans_diag *= st.poisson.isf(poisson_perc/100, trans_diag_raw) / trans_diag_raw

        if normalize_by_cis:
            cis_hom1_mat = clr.matrix(balance=balance).fetch(chrom + hom_suffixes[0], chrom+hom_suffixes[0]).astype(np.float)
            cis_hom2_mat = clr.matrix(balance=balance).fetch(chrom + hom_suffixes[1], chrom+hom_suffixes[1]).astype(np.float)
            cis_hom1_diag = _take_big_diagonal_pixel(cis_hom1_mat, window, ignore_diags, agg_func=np.nansum)
            cis_hom2_diag = _take_big_diagonal_pixel(cis_hom2_mat, window, ignore_diags, agg_func=np.nansum)

            pairing = trans_diag / ((cis_hom1_diag + cis_hom2_diag)/2)                
        else:                                                                        
            pairing = trans_diag                                                     

        pairing_df = pd.DataFrame(dict(binid=chrom_bins.index.values))
        pairing_df['chrom'] = chrom
        pairing_df['start'] = chrom_bins.start.values                             
        pairing_df['end'] = chrom_bins.end.values                                 
        pairing_df['pairing_score'] = pairing

        pairing_dfs.append(pairing_df.copy())                                            

    pairing_dfs = pd.concat(pairing_dfs, ignore_index=True)
    pairing_dfs.sort_values('binid', inplace=True)
    pairing_dfs.drop('binid', axis=1, inplace=True)

    return pairing_dfs                                                                 


def make_pairing_bedgraph(
    clr_path,
    out_bedgraph_path,
    window_bp=10000,
    hom_suffixes = ('_057', '_439'),    
    ignore_diags = 0,
    balance = True,
    poisson_perc = None,
    normalize_by_cis = True,
    transform = 'log2',
    subtract_median = False,
):
    clr = cooler.Cooler(clr_path)
    
    df = get_homolog_pairing_score(
        clr,
        window_bp,
        hom_suffixes,    
        ignore_diags,
        balance,
        poisson_perc,
        normalize_by_cis,
        transform,
        subtract_median
        )

    if transform == 'log2':                                                      
        df['pairing'] = np.log2(df['pairing'])
        df['pairing'][~np.isfinite(df['pairing'])] = np.nan
    if transform == 'log':                                                       
        df['pairing'] = np.log(df['pairing'])
        df['pairing'][~np.isfinite(df['pairing'])] = np.nan
    if transform == 'log10':                                                     
        df['pairing'] = np.log10(df['pairing'])
        df['pairing'][~np.isfinite(df['pairing'])] = np.nan

    if subtract_median:                                                          
        df['pairing'] -= np.nanmedian(df['pairing'])                                         

    df['pairing'] = np.nan_to_num(df['pairing'])
    df.rename({'pairing':'dataValue'})
        
    bedgraph_header = (
        """track type=bedGraph
        name="{transform}pairing{percentile}"
        description="{transform}pairing{percentile}, {window_bp}bp window at {res}bp resolution{norm_by_cis}{subtract_median}{ignore_diags}"
        visibility=full color=0,0,0 altColor=50,50,50 graphType=points yLineMark=0 yLineOnOff=on
        group="{transform}pairing"
        """.replace('\n', ' ').format(
            transform = (transform +'_') if transform else '',
            percentile = ('_'+str(poisson_perc)+'th_perc') if poisson_perc else '',
            window_bp=window_bp,
            res = clr.info['bin-size'],
            norm_by_cis = ', normalize by cis diagonal' if normalize_by_cis else '',
            subtract_median = ', subtract median' if subtract_median else '',
            ignore_diags = ', ignore {} diagonals'.format(ignore_diags) if ignore_diags else '') 
            + '\n')
    return bedgraph_header
    
    with open(out_bedgraph_path, 'w') as out_f:
        
        out_f.write(bedgraph_header)           

        df.to_csv(
            out_f,
            sep='\t',
            header=None,
            index=False)
        
