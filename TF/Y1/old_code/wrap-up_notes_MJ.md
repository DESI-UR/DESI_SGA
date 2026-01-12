Hi! This is a quick list of wrap-up notes, things that are currently not working, and explanations for the end of 2023. As it stands, `Iron_Abell2151-Virgo_ITFR_MJ.ipynb` is the most up-to-date notebook, where I was trying to fit the TFR for Abell-2151 and Virgo clusters. Below is a breakdown of what's going on in that notebook, as well as what the other notebooks in `Y1` are.<p>
    
- Helpful Things to Look At 
    - `DESI_SGA/TF/Y1/jointlinefit.ipynb` is the original write-up of how to do `MultiLinFit`. The cells in `Iron_Abell2151-Virgo_ITFR_MJ.ipynb` were pretty much copied from here. 
    - The notebooks in `TF/SV/` were all used for the Coma cluster, and the most-recently-edited ones are good reading IMO. 
    - The two notebooks in `Y1` that say "explicit" are from doing the cluster membership without using the membership function that's written in `Iron_Abell2151-Virgo_ITFR_MJ.ipynb`. They're also (we think) incorrect, and whatever is wrong with the cluster membership affects those, as well, but they're done just one at a time. 
    

- `Iron_Abell2151-Virgo_ITFR_MJ.ipynb` (Main Notebook)<p>
    - The most up-to-date notebook. 
    - The section <b>Cluster Membership</b> follows the same methods as used for the Coma cluster, but is giving strange outputs for galaxy locations. 
    - Both Abell-2151 and Virgo need their visual inspection (VI) done. There are no generated cutouts for Abell-2151, as the number of galaxies fell below 20 so I didn't make them. Virgo has VI cutouts in `\cache`. Other than VI, the <b>Filtering</b> sections are pretty straightforward if you read cell-by-cell. 
    - Starting with <b>Fitting the TFR</b> is where things get wonky. 
        - You'll see a graph of "Real Data" and another of "incorrect cov2". The Real Data graph should be deleted once `cov2` is resolved. 
        - `cov2` is not working correctly, so most of this section has a provision for "when you get `cov2` working correctly" and the workaround cells can be deleted eventually. (Those should have notes that say so.) Performing the fit and plotting results does not work without `scipy.optimize.minimize` working with the correct `cov2` calculation. 
    - The multi-line fitting (<b>Encapsulate Multi-Dataset Fits in a Single Class</b>) also currently doesn't work due to the `cov2` problem. Everything up to "Run an MCMC" works, but there's an unexpected <i>NaN</i> in `tau` that's keeping it from going on (due to the `cov2` situation). 

That's about everything that's currently wrong with that notebook. Bottom line: cluster membership is still wonky, and there's something up with the covariances that's maybe based on that, maybe not. 
    
- `Iron_Cluster_Membership_MJ.ipynb`
    - This notebook is <i>all</i> of the cluster membership possibilities, with the NEST number (Cluster ID number) and the number of galaxies in that cluster. This is... definitely wrong right now. Probably for similar reasons as to the cluster membership wonkiness from the other notebook. 
    
If there's anything in those notebooks that doesn't make sense please feel free to reach out! Find me at <i>mkell19@u.rochester.edu<i>. <p> 
    
Cheers, and good luck! —MJ