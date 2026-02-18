import pandas as pd

def adapt_insertion(regions_df, info):
    idx, both = find_idx(regions_df, info['chrom_name'], info['pos'])
    existence = info['chrom_name'] in regions_df['chr'].values
    if(idx is not None and existence):
        # Shift all lines after the insertion and only the end of the index where the insertion is (if the insertion is in a targeted interval)
        mask = (regions_df.index > idx) & (regions_df['chr'] == info['chrom_name'])
        regions_df.loc[mask, ['start', 'end']] += len(info['insertion'])
        mask_idx = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
        if(both):
            regions_df.loc[mask_idx, ['start', 'end']] += len(info['insertion'])
        else:
            regions_df.loc[mask_idx, ['end']] += len(info['insertion'])
    
    return(regions_df)

def adapt_deletion(regions_df, info):
    idx_start, both_start = find_idx(regions_df, info['chrom_name'], info['start'])
    idx_end, both_end = find_idx(regions_df, info['chrom_name'], info['end'])
    existence = info['chrom_name'] in regions_df['chr'].values
    
    if(idx_start is not None and idx_end is not None and existence):
        
        # Case 1: deletion starts and ends between two targeted regions, covering zero, one, or more regions
        if(both_start and both_end):
            mask = (regions_df.index > idx_end) & (regions_df['chr'] == info['chrom_name'])
            regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
            
            if idx_end > idx_start:
                regions_df = regions_df.drop(
                    regions_df.index[idx_start:idx_end+1] # .index methods is end index exclusive
                ).reset_index(drop=True)
                
                
        # Case 2: deletion starts or ends within a targeted region, covering zero, one, or more regions
        if(not both_start):
            if idx_start == idx_end:
                mask = (regions_df.index > idx_end) & (regions_df['chr'] == info['chrom_name'])
                regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
                mask_idx = (regions_df.index == idx_end) & (regions_df['chr'] == info['chrom_name'])
                regions_df.loc[mask_idx, ['end']] -= info['end'] - info['start']
            
            if idx_end > idx_start:
                if(both_end):
                    mask = (regions_df.index >= idx_end) & (regions_df['chr'] == info['chrom_name'])
                    regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
                    
                    regions_df.loc[idx_start, ['end']] = info['start']
                    regions_df = regions_df.drop(
                        regions_df.index[idx_start+1:idx_end+1] # .index methods is end index exclusive
                    ).reset_index(drop=True)
                    
                    
                else:
                    mask = (regions_df.index > idx_end) & (regions_df['chr'] == info['chrom_name'])
                    regions_df.loc[mask, ['start', 'end']] -= info['end'] - info['start']
                    
                    regions_df.loc[idx_start, ['end']] = info['start']
                    regions_df.loc[idx_end, ['end']] = info['end']
                    regions_df = regions_df.drop(
                        regions_df.index[idx_start+1:idx_end] # .index methods is end index exclusive
                    ).reset_index(drop=True)
        
    if(idx_start is not None and idx_end is None and existence):
        # Set idx_end to last idx in the chromosome
        mask = (regions_df['chr'] == info['chrom_name'])
        idx_end = regions_df.index[mask].max()
        
        if(not both_start):
            regions_df.loc[idx_start, ['end']] = info['start']
            idx_start = idx_start + 1
            
        regions_df = regions_df.drop(
                    regions_df.index[idx_start:idx_end+1] # .index methods is end index exclusive
                ).reset_index(drop=True)
    
    return(regions_df)
    
def adapt_CNV(regions_df, info):

    existence = info['chrom_name'] in regions_df['chr'].values
    
    # Case 1: CNV starts and ends between two targeted regions, covering zero, one, or more regions
    # Nothing to do
    
    
    # Case 2: CNV starts or ends within a targeted regions, covering one or more regions
    idx_start, both_start = find_idx(regions_df, info['chrom_name'], info['start'])
    if(not both_start and idx_start is not None and existence):
        # Split the regions in the middle of the CNV
        row = regions_df.loc[idx_start].copy()
        line1 = row.copy()
        line1['end'] = info['start']
        line2 = row.copy()
        line2['start'] = info['start']

        # Drop the old row and insert the new ones
        regions_df = pd.concat([
            regions_df.loc[:idx_start-1],
            pd.DataFrame([line1, line2]),
            regions_df.loc[idx_start+1:]
        ]).reset_index(drop=True)
        
    idx_start, both_start = find_idx(regions_df, info['chrom_name'], info['start'])
    
    idx_end, both_end = find_idx(regions_df, info['chrom_name'], info['end'])
    if(not both_end and idx_end is not None and existence):
        # Split the regions in the middle of the CNV
        row = regions_df.loc[idx_end].copy()
        line1 = row.copy()
        line1['end'] = info['end']
        line2 = row.copy()
        line2['start'] = info['end']

        # Drop the old row and insert the new ones
        regions_df = pd.concat([
            regions_df.loc[:idx_end-1],
            pd.DataFrame([line1, line2]),
            regions_df.loc[idx_end+1:]
        ]).reset_index(drop=True)
        
    idx_end, both_end = find_idx(regions_df, info['chrom_name'], info['end'])
        
    if(idx_start is not None and existence):
        # Shift all lines after the CNV
        if(idx_end is not None):
            mask = (regions_df.index > idx_end) & (regions_df['chr'] == info['chrom_name'])
            regions_df.loc[mask, ['start', 'end']] += (info['end'] - info['start']) * (info['rep_num'] - 1)
            mask_idx = (regions_df.index == idx_end) & (regions_df['chr'] == info['chrom_name'])
            if(both_end):
                regions_df.loc[mask_idx, ['start', 'end']] += (info['end'] - info['start']) * (info['rep_num'] - 1)
            else:
                regions_df.loc[mask_idx, ['end']] += (info['end'] - info['start']) * (info['rep_num'] - 1)
        else:
            # Set idx_end to last idx in the chromosome
            mask = (regions_df['chr'] == info['chrom_name'])
            idx_end = regions_df.index[mask].max()

        if(idx_start < idx_end):
            # Rows to duplicate
            block = regions_df.loc[idx_start:idx_end]
            block_dups = []
            for i in range(1, info['rep_num']):
                b = block.copy()
                b[["start", "end"]] += i * (info['end'] - info['start'])
                block_dups.append(b)
            
            # Insert duplicated rows into the dataframe
            mask = regions_df["chr"] == info['chrom_name']
            regions_df = pd.concat(
                [
                    regions_df.loc[:idx_end],
                    *block_dups,
                    regions_df.loc[idx_end + 1:],
                ],
                ignore_index=True,
            )
        
    return(regions_df)
        
def adapt_aneuploidy(regions_df, info):
    existence = info['chrom_name'] in regions_df['chr'].values
    if existence:
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
    existence1 = info['chrom_name1'] in regions_df['chr'].values
    if existence1:
        # Adapt dataframes at breakpoints
        idx1, both1 = find_idx(regions_df, info['chrom_name1'], info['bkpt1'])
        if(idx1 is not None):
            # Split a region if the translocation happens in its middle
            if(both1):
                idx1 = idx1 - 1
            else:
                row = regions_df.loc[idx1].copy()
                line1 = row.copy()
                line1['end'] = info['bkpt1']
                line2 = row.copy()
                line2['start'] = info['bkpt1']
    
                # Drop the old row and insert the new ones
                regions_df = pd.concat([
                    regions_df.loc[:idx1-1],
                    pd.DataFrame([line1, line2]),
                    regions_df.loc[idx1+1:]
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

    else: # if chrom_name1 is not in the panel
        idx1 = find_previous_idx(info['chrom_name1'], regions_df)
        slice1 = pd.DataFrame(columns=regions_df.columns)
        
    
    existence2 = info['chrom_name2'] in regions_df['chr'].values
    if existence2:
        # Adapt dataframes at breakpoints
        idx2, both2 = find_idx(regions_df, info['chrom_name2'], info['bkpt2'])
        if(idx2 is not None):
            if(both2):
                idx2 = idx2 - 1
            else:
                row = regions_df.loc[idx2].copy()
                line1 = row.copy()
                line1['end'] = info['bkpt2']
                line2 = row.copy()
                line2['start'] = info['bkpt2']
    
                # Drop the old row and insert the new ones
                regions_df = pd.concat([
                    regions_df.loc[:idx2-1],
                    pd.DataFrame([line1, line2]),
                    regions_df.loc[idx2+1:]
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
            
    else: # if chrom_name2 is not in the panel
        idx2 = find_previous_idx(info['chrom_name2'], regions_df)
        slice2 = pd.DataFrame(columns=regions_df.columns)
        
    # Normal translocation
    if(info['normal']):
        # Compose the new dataframe
        min_idx = min(idx1, idx2)
        max_idx = max(idx1, idx2)
        regions_df = pd.concat(
            [
                regions_df.loc[:min_idx],
                slice2,
                regions_df.loc[min_idx+slice1.shape[0]+1:max_idx],
                slice1,
                regions_df.loc[max_idx+slice2.shape[0]+1:],
            ],
            ignore_index=True
        )
    
    # Anormal translocation
    else:
        # Compose the new dataframe
        min_idx = min(idx1, idx2)
        max_idx = max(idx1, idx2)
        regions_df = pd.concat(
            [
                regions_df.loc[:min_idx],
                regions_df.loc[min_idx+slice1.shape[0]+1:max_idx],
                slice1,
                slice2, # You should change also the offset for the start and end
                regions_df.loc[max_idx+slice2.shape[0]+1:],
            ],
            ignore_index=True
        )
        
    return(regions_df)

def adapt_chromothripsis(regions_df, info):
    existence = info['chrom_name'] in regions_df['chr'].values
    df_list = []
    previous_idx, _ = find_idx(regions_df, info['chrom_name'], info['start'])
    if(previous_idx is not None and existence):
        for bkpt in info['bkpts']:
            current_idx, both = find_idx(regions_df, info['chrom_name'], info['start'] + bkpt)
            if(current_idx is not None):
                # Split a region if the translocation happens in its middle
                if(both and current_idx != 0):
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
                if(current_idx != 0):
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
        start_idx, _ = find_idx(regions_df, info['chrom_name'], info['start'])
        end_idx, _ = find_idx(regions_df, info['chrom_name'], info['end'])
        regions_df = pd.concat(
            [
                regions_df.loc[:start_idx],
                *df_list_shuffled,
                regions_df.loc[end_idx:]
            ],
            ignore_index=True
            )
    
    return(regions_df)
    
def adapt_BFB(regions_df, info):
    existence = info['chrom_name'] in regions_df['chr'].values
    if existence:
        for i in range(len(info['bkpts'])):
            idx, both = find_idx(regions_df, info['chrom_name'], info['bkpts'][i])
            if(idx is not None):
                # Shift all lines after the breakpoint and only the end of the index where the breakpoint is (if the breakpoint is in a targeted interval)
                mask = (regions_df.index > idx) & (regions_df['chr'] == info['chrom_name'])
                if i == 0:
                    regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i]
                else:
                    regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i] - info['bkpts'][i-1]
                    
                mask = (regions_df.index == idx) & (regions_df['chr'] == info['chrom_name'])
                if i == 0:
                    if(both):
                        regions_df.loc[mask, ['start', 'end']] += info['bkpts'][i]
                    else:
                        regions_df.loc[mask, ['end']] += info['bkpts'][i]
                else:
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
        
        if(idx != 0 and regions_df.iloc[idx-1, :]['end'] >= loc and regions_df.iloc[idx-1, :]['chr'] == chr):
            idx = idx - 1
            both = False
    else:
        idx = None
        both = False
    
    return(idx, both)

def find_previous_idx(chr_name, 
                      regions_df, 
                      chr_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                                   'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                                   'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']):
    
    present_chr = [c for c in chr_order if c in regions_df['chr'].values]
    
    if chr_name not in chr_order:
        raise ValueError(f"Chromosome {chr_name} not found in chr_order")
    
    chr_idx = chr_order.index(chr_name) - 1

    if chr_idx < 0:
        idx = 0
    
    elif chr_order[chr_idx] not in present_chr:
        idx = find_previous_idx(chr_order[chr_idx], regions_df)
        
    elif(chr_order[chr_idx] in present_chr):
        mask = regions_df['chr'] == chr_order[chr_idx]
        idx = regions_df.index[mask].max() + 1
         
    return idx