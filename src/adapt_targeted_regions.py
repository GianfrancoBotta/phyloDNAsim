import pandas as pd

def adapt_insertion(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['pos'])
    # Shift all lines after the insertion and only the end of the index where the insertion is (if the insertion is in a targeted interval)
    mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
    regions_df.loc[mask, ['start', 'end']] += len(info['insertion'])
    if(both):
        regions_df.loc[mask, ['start", "end']] += len(info['insertion'])
    else:
        regions_df.loc[mask, ['end']] += len(info['insertion'])
    
    return(regions_df)

def adapt_deletion(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['start'])
    # Shift all lines after the deletion and only the end of the index where the deletion is (if the deletion is in a targeted interval)
    mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
    regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
    if(both):
        regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
    else:
        regions_df.loc[mask, ['end']] -= info['end'] - info['start']
    
def adapt_CNV(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['start'])
    # Shift all lines after the CNV and only the end of the index where the CNV is (if the CNV is in a targeted interval)
    mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
    regions_df.loc[mask, ['start', 'end']] += (info['end'] - info['start']) * info['rep_num']
    if(both):
        regions_df.loc[mask, ['start', 'end']] += (info['end'] - info['start']) * info['rep_num']
    else:
        regions_df.loc[mask, ['end']] += (info['end'] - info['start']) * info['rep_num']
        
def adapt_aneuploidy(regions_df, info):
    # Rows to duplicate
    block = regions_df[regions_df["chr"] == info['chrom_name']]
    block_dups = []
    for i in range(1, info['rep_num']):
        b = block.copy()
        b[["start", "end"]] += i * info['chrom_length']
        block_dups.append(b)
    
    # Insert duplicated rows into the dataframe
    mask = regions_df["chr"] == info['chrom_name']
    last_idx = regions_df.index[mask].max()
    regions_df = pd.concat(
        [
            regions_df.loc[:last_idx],
            *block_dups,
            regions_df.loc[last_idx + 1:],
        ],
        ignore_index=True,
    )
    return(regions_df)

def adapt_translocation(regions_df, info):
    idx1, both1 = find_idx(regions_df, info['chrom_name1'], info['bkpt1'])
    idx2, both2 = find_idx(regions_df, info['chrom_name2'], info['bkpt2'])
    
    # Split a region if the translocation happens in its middle
    if(both1):
        idx1 = idx1 - 1
    else:
        row = regions_df.iloc[idx1-1].copy()
        line1 = row.copy()
        line1['end'] = info['bkpt1']
        line2 = row.copy()
        line2['start'] = info['bkpt1']

        # Drop the old row and insert the new ones
        regions_df = pd.concat([
            regions_df.iloc[:idx1-1],
            pd.DataFrame([line1, line2]),
            regions_df.iloc[idx1:]
        ]).reset_index(drop=True)
    
    if(both2):
        idx2 = idx2 - 1
    else:
        row = regions_df.iloc[idx2-1].copy()
        line1 = row.copy()
        line1['end'] = info['bkpt2']
        line2 = row.copy()
        line2['start'] = info['bkpt2']

        # Drop the old row and insert the new ones
        regions_df = pd.concat([
            regions_df.iloc[:idx2-1],
            pd.DataFrame([line1, line2]),
            regions_df.iloc[idx2:]
        ]).reset_index(drop=True)
    idx2 = idx2 - 1
    
    mask1 = regions_df["chr"] == info['chrom_name1']
    last_idx1 = regions_df.index[mask1].max()
    mask2 = regions_df["chr"] == info['chrom_name2']
    last_idx2 = regions_df.index[mask2].max()
    
    # Translocate the regions
    if(info['normal']):
        slice1 = regions_df.loc[idx1:last_idx1].copy()
        slice2 = regions_df.loc[idx2:last_idx2].copy()
        
        # Change slice 1
        slice1['chr'] = info['chrom_name2']
        slice1[['start', 'end']] += info['bkpt2'] - info['bkpt1']
        
        # Change slice 2
        slice2['chr'] = info['chrom_name1']
        slice2[['start', 'end']] +=  info['bkpt1'] - info['bkpt2']
    
        # Compose the new dataframe
        regions_df = pd.concat(
            [
                regions_df.loc[:idx1],
                slice2,
                regions_df.loc[idx1+slice1.shape[0]+1:idx2],
                slice1,
                regions_df.loc[idx2+slice2.shape[0]+1:],
            ],
            ignore_index=True
        )
        
    return(regions_df)

def adapt_chromothripsis(regions_df, info):
    df_list = []
    mask = regions_df["chr"] == info['chrom_name']
    previous_idx = regions_df.index[mask].min()
    for bkpt in info['bkpts']:
        current_idx, both = find_idx(regions_df, info['chrom_name'], bkpt)
        
        # Split a region if the translocation happens in its middle
        if(both):
            current_idx = current_idx - 1
        else:
            row = regions_df.iloc[current_idx-1].copy()
            line1 = row.copy()
            line1['end'] = info['bkpt']
            line2 = row.copy()
            line2['start'] = info['bkpt']

            # Drop the old row and insert the new ones
            regions_df = pd.concat([
                regions_df.iloc[:current_idx-1],
                pd.DataFrame([line1, line2]),
                regions_df.iloc[current_idx:]
            ]).reset_index(drop=True)
        
        df_list.append(regions_df.loc[previous_idx:current_idx])
        previous_idx = current_idx + 1
    
    # Shuffle the order of the df according to the mutational event
    df_list_shuffled = [None] * len(df_list)
    for i, idx in enumerate(info['order']):
        df = df_list[idx]
        if i == 0 and idx != 0:
            offset = -info['bkpts'][idx-1]
        elif i != 0:
            offset = info['bkpts'][info['order'][i-1]] - info['bkpts'][info['order'][i-1] - 1]
        else:
            offset = 0

        df[['start', 'end']] += offset
        
        df_list_shuffled[i] = df
        
    # Compose the new dataframe
    regions_df = pd.concat(
        [
            regions_df.loc[:info['start']],
            *df_list_shuffled,
            regions_df.loc[info['end']+1:]
        ],
        ignore_index=True
        )
    
    return(regions_df)
    
def adapt_BFB(regions_df, info):
    for i in range(len(info['bkpts'])):
        idx, both = find_idx(regions_df, info['chrom_name'], info['bkpt'])
        # Shift all lines after the breakpoint and only the end of the index where the breakpoint is (if the breakpoint is in a targeted interval)
        mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
        if i == 0:
            regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i]
            if(both):
                regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i]
            else:
                regions_df.loc[mask, ['end']] += info['bkpts'][i]
        else:
            regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i] - info['bkpts'][i-1]
            if(both):
                regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i] - info['bkpts'][i-1]
            else:
                regions_df.loc[mask, ['end']] += info['bkpts'][i] - info['bkpts'][i-1]
    
    return(regions_df)

def find_idx(regions_df, chr, loc):
    mask = (
        (regions_df['chr'] == chr) &
        (regions_df['start'] >= loc)
    )

    idx = regions_df.index[mask][0]
    both = True
    
    if(regions_df.iloc[idx-1, :]['end'] >= loc):
        idx = idx - 1
        both = False
    
    return(idx, both)