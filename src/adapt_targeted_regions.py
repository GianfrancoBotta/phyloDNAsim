import pandas as pd

def adapt_insertion(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['pos'])
    if(idx is not None):
        # Shift all lines after the insertion and only the end of the index where the insertion is (if the insertion is in a targeted interval)
        mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
        regions_df.loc[mask, ['start', 'end']] += len(info['insertion'])
        if(both):
            regions_df.loc[mask, ['start', 'end']] += len(info['insertion'])
        else:
            regions_df.loc[mask, ['end']] += len(info['insertion'])
    
    return(regions_df)

def adapt_deletion(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['start'])
    if(idx is not None):
        # Shift all lines after the deletion and only the end of the index where the deletion is (if the deletion is in a targeted interval)
        mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
        regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
        if(both):
            regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
        else:
            regions_df.loc[mask, ['end']] -= info['end'] - info['start']
    
    return(regions_df)
    
def adapt_CNV(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['start'])
    if(idx is not None):
        # Shift all lines after the CNV and only the end of the index where the CNV is (if the CNV is in a targeted interval)
        mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
        regions_df.loc[mask, ['start', 'end']] += (info['end'] - info['start']) * info['rep_num']
        if(both):
            regions_df.loc[mask, ['start', 'end']] += (info['end'] - info['start']) * info['rep_num']
        else:
            regions_df.loc[mask, ['end']] += (info['end'] - info['start']) * info['rep_num']
    
    return(regions_df)
        
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
    # Adapt dataframes at breakpoints
    idx1, both1 = find_idx(regions_df, info['chrom_name1'], info['bkpt1'])
    if(idx1 is not None):
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
                regions_df.iloc[idx1+1:]
            ]).reset_index(drop=True)
        mask1 = regions_df["chr"] == info['chrom_name1']
        last_idx1 = regions_df.index[mask1].max()
        # Take the df slice
        slice1 = regions_df.loc[idx1:last_idx1].copy()
        slice1['chr'] = info['chrom_name2']
        slice1[['start', 'end']] += info['bkpt2'] - info['bkpt1']
    else:
        mask1 = regions_df["chr"] == info['chrom_name1']
        idx1 = regions_df.index[mask1].max() + 1
        slice1 = pd.DataFrame(columns=regions_df.columns)

    idx2, both2 = find_idx(regions_df, info['chrom_name2'], info['bkpt2'])
    if(idx2 is not None):
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
        mask2 = regions_df["chr"] == info['chrom_name2']
        last_idx2 = regions_df.index[mask2].max()
        # Take the df slice
        slice2 = regions_df.loc[idx2:last_idx2].copy()
        slice2['chr'] = info['chrom_name1']
        slice2[['start', 'end']] +=  info['bkpt1'] - info['bkpt2']
    else:
        mask2 = regions_df["chr"] == info['chrom_name2']
        idx2 = regions_df.index[mask2].max() + 1
        slice2 = pd.DataFrame(columns=regions_df.columns)

    # Normal translocation
    if(info['normal']):
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
    
    # Anormal translocation
    else:
        # Compose the new dataframe
        regions_df = pd.concat(
            [
                regions_df.loc[:idx1],
                regions_df.loc[idx1+slice1.shape[0]+1:idx2],
                slice1,
                slice2,
                regions_df.loc[idx2+slice2.shape[0]+1:],
            ],
            ignore_index=True
        )
        
    return(regions_df)

def adapt_chromothripsis(regions_df, info):
    df_list = []
    previous_idx, _ = find_idx(regions_df, info['chrom_name'], info['start'])
    if(previous_idx is not None):
        for bkpt in info['bkpts']:
            current_idx, both = find_idx(regions_df, info['chrom_name'], info['start'] + bkpt)
            if(current_idx is not None):
                # Split a region if the translocation happens in its middle
                if(both):
                    current_idx = current_idx - 1
                else:
                    row = regions_df.iloc[current_idx-1].copy()
                    line1 = row.copy()
                    line1['end'] = bkpt
                    line2 = row.copy()
                    line2['start'] = bkpt

                    # Drop the old row and insert the new ones
                    regions_df = pd.concat([
                        regions_df.iloc[:current_idx-1],
                        pd.DataFrame([line1, line2]),
                        regions_df.iloc[current_idx:]
                    ]).reset_index(drop=True)
            
                df_list.append(regions_df.loc[previous_idx:current_idx])
                previous_idx = current_idx + 1
                
        # Append last interval
        df_list.append(regions_df.loc[previous_idx:previous_idx+1])
        
        # Shuffle the order of the df according to the mutational event
        order = [x for x in info['order'] if x < len(df_list)]
        df_list_shuffled = [None] * len(df_list)
        for i, idx in enumerate(order):
            df = df_list[idx].copy()
            if i == 0 and idx != 0:
                offset = -info['bkpts'][idx-1]
            elif i != 0 and order[i-1] != 0:
                offset = info['bkpts'][order[i-1]] - info['bkpts'][order[i-1] - 1]
            elif i != 0 and order[i-1] == 0:
                offset = info['bkpts'][order[i-1]]
            else:
                offset = 0

            df.loc[:, ['start', 'end']] += offset
            
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
        idx, both = find_idx(regions_df, info['chrom_name'], info['bkpts'][i])
        if(idx is not None):
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
    
    if(sum(mask)):
        idx = regions_df.index[mask][0]
        both = True
        
        if(regions_df.iloc[idx-1, :]['end'] >= loc):
            idx = idx - 1
            both = False
    else:
        idx = None
        both = False
    
    return(idx, both)