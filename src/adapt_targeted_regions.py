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
        if info['rep_num'] > 0:
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
        else:
            # Rows to delete
            regions_df = regions_df[regions_df["chr"] != info['chrom_name']].reset_index(drop=True)
    
    return(regions_df)


def adapt_translocation(regions_df, info):
    # Extract chromosome dataframes
    chr1_df = regions_df[regions_df['chr'] == info['chrom_name1']].copy()
    chr2_df = regions_df[regions_df['chr'] == info['chrom_name2']].copy()

    # Offsets
    offset1 = info['bkpt2'] - info['bkpt1']
    offset2 = info['bkpt1'] - info['bkpt2']

    # Split and shift chromosomes
    slice1 = split_chromosome(chr1_df, info['bkpt1'], info['chrom_name2'], offset1)
    slice2 = split_chromosome(chr2_df, info['bkpt2'], info['chrom_name1'], offset2)

    # Remove original chromosomes from regions_df
    regions_df = regions_df[~regions_df['chr'].isin([info['chrom_name1'], info['chrom_name2']])].copy()

    # Concatenate back
    if info['normal']:
        # normal translocation: slice2 before slice1
        regions_df = pd.concat([regions_df, slice2, slice1], ignore_index=True)
    else:
        # abnormal translocation: slice1 then slice2
        regions_df = pd.concat([regions_df, slice1, slice2], ignore_index=True)

    # Sort by chromosome and start to maintain order
    regions_df = regions_df.sort_values(['chr','start', 'end']).reset_index(drop=True)
    
    return regions_df

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
                    line1['end'] = info['start'] + bkpt
                    line2 = row.copy()
                    line2['start'] = info['start'] + bkpt

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


def find_idx(regions_df, chr_name, loc):
    chr_df = regions_df[regions_df['chr'] == chr_name]

    if chr_df.empty:
        return None, False

    # Case 1 — inside a region
    inside = chr_df[
        (chr_df['start'] <= loc) &
        (chr_df['end'] > loc)
    ]

    if not inside.empty:
        return inside.index[0], False  # False = inside region

    # Case 2 — between regions
    after = chr_df[chr_df['start'] > loc]

    if not after.empty:
        return after.index[0], True  # True = between regions

    # Case 3 — after last region
    return None, False

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

# Split and shift function for a single chromosome
def split_chromosome(chrom_regions, bkpt, new_chr, offset):
    if chrom_regions.empty:
        return pd.DataFrame(columns=chrom_regions.columns)

    # Case 1: breakpoint inside a region
    inside_mask = (chrom_regions['start'] <= bkpt) & (chrom_regions['end'] > bkpt)
    if inside_mask.any():
        split_idx = chrom_regions.index[inside_mask][0]

        row = chrom_regions.loc[split_idx].copy()

        # Split the row
        top = row.copy()
        top['end'] = bkpt
        bottom = row.copy()
        bottom['start'] = bkpt

        # Regions staying in original chromosome
        before = chrom_regions.loc[:split_idx-1]
        top_df = pd.DataFrame([top])

        # Regions to move
        after = pd.concat([pd.DataFrame([bottom]), chrom_regions.loc[split_idx+1:]], ignore_index=True)
    
    else:
        # Case 2 & 3: breakpoint before first region or between regions
        after_mask = chrom_regions['start'] > bkpt
        if after_mask.any():
            split_idx = chrom_regions.index[after_mask][0]
            before = chrom_regions.loc[:split_idx-1] if split_idx > 0 else pd.DataFrame(columns=chrom_regions.columns)
            after = chrom_regions.loc[split_idx:].copy()
            top_df = pd.DataFrame(columns=chrom_regions.columns)  # nothing stays in original
        else:
            # Breakpoint after last region → nothing moves
            return chrom_regions.copy()

    # Apply offset and new chromosome only to the moved part
    after[['start','end']] += offset
    after['chr'] = new_chr

    # Combine regions that stay + moved regions
    new_df = pd.concat([before, top_df, after], ignore_index=True)
    
    return new_df