import numpy as np
from tqdm import tqdm

def select_galaxy(table, sga_table, column_classification):
    """ Remove all galaxies that do not have 
    1) center point and two unique points
    or
    2) a symmetric pair and one unique point

    Parameters
    ----------
    table:
        the database you want to sort through

    sga_table:
        SGA fits

    column_classification (str):
        All good galaxies are marked as a 1, this creates a new column to mark

    Returns
    -------
    New table with only galaxies marked as 1
    """
    sga_dict = {}
    for i in range(len(sga_table)):
        sga_dict[sga_table['SGA_ID'][i]] = i
    
    #target distance that can be considered in center
    #units R26
    center_dist_lim = 0.001

    #minimum distance between two targets to be considered unique
    #units R26
    unique_dist_lim = 0.01

    #empty column to classify galaxy
    table[column_classification] = 0
    
    #for each unique SGA ID
    for i in tqdm(np.unique(table['SGA_ID'])):
    
        #identify all the galaxies and targets
        obs_id = np.logical_and(table['SGA_ID'] == i, table['TARGETID'] > 0)
        #logical_and returns if both statements are true
        
        #makes a table of the targets corresponding to this galaxy
        obs = table[obs_id]
    
        sga_id = sga_dict[i] 
    
    #if the galaxy has more than three observations
        if len(obs) >= 3:
    
    #-----
    # a
    #-----
            #check to see if there is a center observation
            if np.any(obs['DIST_R26'] < center_dist_lim):
    
                #since there is a center observation, check if there are two unique points
                not_center = obs[obs['DIST_R26'] > unique_dist_lim]
    
                if len(not_center) >= 2 and (np.max(not_center['DIST_R26']) - np.min(not_center['DIST_R26'])) >= unique_dist_lim:
    
                    #since there 2 unique points, classify it as viable galaxy
                    table['{}'.format(column_classification)][obs_id] = 1
    #----
    # b
    #----
            #since there isn't a center observation
            else:
                
                #split targets into above and below SGA declination
                if (sga_table['PA'][sga_id] < 45) or (sga_table['PA'][sga_id] > 135):
    
                    left_index = obs['TARGET_DEC'] - sga_table['DEC'][sga_id] > 0
    
                else:
                    left_index = obs['TARGET_RA'] - sga_table['RA'][sga_id] > 0
                    
                left = obs[left_index]
                right = obs[~left_index]
    
                if len(left) > 0 and len(right) > 0:
    
                    for j in range(len(left)):
                        
                        # check that there are 2 symmetric observations
                        if np.any(np.abs(right['DIST_R26'] - left['DIST_R26'][j]) < unique_dist_lim):
    
                        #check if there is a third point
                            if (np.any(np.abs(right['DIST_R26'] - left['DIST_R26'][j]) > unique_dist_lim)
                            or np.any(np.abs(left['DIST_R26'] - left['DIST_R26'][j]) > unique_dist_lim)):
    
                                #viable galaxy
                                table[column_classification][obs_id] = 1
    
    return (table[table[column_classification]==1])

def velocity_selection(loa_new, column_classification):
    """ Remove all galaxies that have < 3 observations that are all less than 1000 km/s

    Parameters
    ----------
    loa_new:
        the database you want to sort through

    column_classification (str):
        All good galaxies are marked as a 1, this creates a new column to mark

    Returns
    -------
    New table with only galaxies marked as 0
    """
    
    # create a new table that removes all galaxies with < 3 observations that are all less than 1000 km/s
    loa_new[column_classification] = 0 
    for sga_id in tqdm(np.unique(loa_new['SGA_ID'])):
        
        #identify all targets in the galaxy
        targ_list = loa_new[loa_new['SGA_ID']==sga_id]
    
        #get velocities for each point
        velocity = np.array(targ_list['Velocity'])
    
        #------------------------------------------------------------------
        # axe any galaxies that have < 3 points with valid velocities
        #------------------------------------------------------------------
        counter = 0
    
        for i in range(len(velocity)):
            if abs(velocity[i]) < 1000:
                counter += 1
                
        if counter < 3:
            #print (sga_idx, ':not enough points')
            continue
    
        elif counter >= 3:
            loa_new[column_classification][loa_new['SGA_ID']==sga_id] = 1
        
    return (loa_new[loa_new[column_classification]==1])