import pandas as pd

def find_idx(regions_df, chr_name, loc):
    """
    Find the index of the region that contains or follows `loc` on `chr_name`.

    Returns
    -------
    (index, between) where:
      - index   : integer label in regions_df, or None if loc is past all regions
      - between : True  → loc falls *between* two regions (or before the first)
                  False → loc falls *inside* a region
    """
    chr_df = regions_df[regions_df['chr'] == chr_name]

    if chr_df.empty:
        return None, False

    # Case 1 — inside a region [start, end)
    inside = chr_df[(chr_df['start'] <= loc) & (chr_df['end'] > loc)]
    if not inside.empty:
        return inside.index[0], False

    # Case 2 — between regions (loc < start of next region)
    after = chr_df[chr_df['start'] > loc]
    if not after.empty:
        return after.index[0], True

    # Case 3 — past the last region on this chromosome
    return None, False


def _safe_concat(parts):
    """pd.concat that ignores empty DataFrames and always resets the index."""
    non_empty = [p for p in parts if p is not None and not p.empty]
    if not non_empty:
        return pd.DataFrame()
    return pd.concat(non_empty, ignore_index=True)


def _split_region_at(regions_df, chr_name, pos):
    """
    If `pos` falls inside a region on `chr_name`, split that region into two
    rows ([start, pos) and [pos, end)) and return the updated dataframe together
    with the index of the *second* half.  If `pos` is already a boundary the
    dataframe is returned unchanged and the index of the next region is returned.
    """
    idx, between = find_idx(regions_df, chr_name, pos)

    if idx is None or between:
        # pos is past all regions or already on a boundary — nothing to split
        return regions_df, idx

    row = regions_df.loc[idx]
    if row['start'] == pos:
        # pos is exactly the start of this region — already a boundary
        return regions_df, idx

    # Split
    top = row.copy()
    top['end'] = pos
    bottom = row.copy()
    bottom['start'] = pos

    # Safe slicing: avoid negative indices when idx == 0
    before = regions_df.loc[:idx - 1] if idx > 0 else pd.DataFrame(columns=regions_df.columns)
    after  = regions_df.loc[idx + 1:]

    regions_df = _safe_concat([before, pd.DataFrame([top]), pd.DataFrame([bottom]), after])

    # The bottom half is now at the position right after `top`
    new_idx, _ = find_idx(regions_df, chr_name, pos)
    return regions_df, new_idx


def _split_and_shift(chrom_regions, bkpt, new_chr, offset):
    """
    Split a single-chromosome DataFrame at `bkpt`.
    Regions before the breakpoint keep their chromosome name.
    Regions from the breakpoint onward get `new_chr` and are shifted by `offset`.
    """
    if chrom_regions.empty:
        return chrom_regions.copy()

    cols = chrom_regions.columns

    # Breakpoint inside a region → split it
    inside_mask = (chrom_regions['start'] <= bkpt) & (chrom_regions['end'] > bkpt)
    if inside_mask.any():
        split_idx = chrom_regions.index[inside_mask][0]
        row    = chrom_regions.loc[split_idx]
        top    = row.copy(); top['end']   = bkpt
        bottom = row.copy(); bottom['start'] = bkpt

        before = chrom_regions.loc[:split_idx - 1] if split_idx > chrom_regions.index[0] else pd.DataFrame(columns=cols)
        after  = _safe_concat([pd.DataFrame([bottom]), chrom_regions.loc[split_idx + 1:]])
    else:
        after_mask = chrom_regions['start'] > bkpt
        if after_mask.any():
            split_idx = chrom_regions.index[after_mask][0]
            before = chrom_regions.loc[:split_idx - 1] if split_idx > chrom_regions.index[0] else pd.DataFrame(columns=cols)
            after  = chrom_regions.loc[split_idx:].copy()
            top    = None
        else:
            # Breakpoint is after all regions — nothing moves
            return chrom_regions.copy()

    after = after.copy()
    after[['start', 'end']] += offset
    after['chr'] = new_chr

    parts = [before]
    if inside_mask.any():
        parts.append(pd.DataFrame([top]))
    parts.append(after)

    return _safe_concat(parts)


def _sanitize(regions_df):
    # Remove regions that have been fully consumed
    regions_df = regions_df[regions_df['end'] > regions_df['start']].reset_index(drop=True)
    return regions_df


def find_previous_idx(
    chr_name,
    regions_df,
    chr_order=('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
               'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
               'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'),
):
    chr_order = list(chr_order)
    present_chr = [c for c in chr_order if c in regions_df['chr'].values]

    if chr_name not in chr_order:
        raise ValueError(f"Chromosome {chr_name} not found in chr_order")

    chr_idx = chr_order.index(chr_name) - 1

    if chr_idx < 0:
        return 0

    if chr_order[chr_idx] not in present_chr:
        return find_previous_idx(chr_order[chr_idx], regions_df, chr_order)

    mask = regions_df['chr'] == chr_order[chr_idx]
    return int(regions_df.index[mask].max()) + 1


# ---------------------------------------------------------------------------
# Insertion
# ---------------------------------------------------------------------------

def adapt_insertion(regions_df, info):
    """
    Shift coordinates after a point insertion.

    info keys: chrom_name, pos, insertion (sequence string)
    """
    if info['chrom_name'] not in regions_df['chr'].values:
        return regions_df

    idx, between = find_idx(regions_df, info['chrom_name'], info['pos'])

    if idx is None:
        return regions_df

    ins_len = len(info['insertion'])
    chrom = info['chrom_name']

    # All regions strictly after the insertion point shift entirely
    mask_after = (regions_df.index > idx) & (regions_df['chr'] == chrom)
    regions_df.loc[mask_after, ['start', 'end']] += ins_len

    # The region that *contains* the insertion: only the end shifts
    # (if between==True the insertion is before this region → shift both)
    mask_at = (regions_df.index == idx) & (regions_df['chr'] == chrom)
    if between:
        regions_df.loc[mask_at, ['start', 'end']] += ins_len
    else:
        regions_df.loc[mask_at, 'end'] += ins_len

    return regions_df


# ---------------------------------------------------------------------------
# Deletion
# ---------------------------------------------------------------------------

def adapt_deletion(regions_df, info):
    """
    Remove a deleted segment from the coordinate space.

    info keys: chrom_name, start, end
    """
    chrom   = info['chrom_name']
    del_start = info['start']
    del_end   = info['end']
    del_len   = del_end - del_start

    if chrom not in regions_df['chr'].values:
        return regions_df

    idx_start, both_start = find_idx(regions_df, chrom, del_start)
    idx_end,   both_end   = find_idx(regions_df, chrom, del_end)

    # ---- deletion start and end are both *between* targeted regions ----------
    if both_start and (both_end or idx_end is None):
        # Shift everything after the deletion
        if idx_end is not None:
            mask = (regions_df.index >= idx_end) & (regions_df['chr'] == chrom)
        else:
            # deletion extends past all regions on this chromosome
            mask = pd.Series(False, index=regions_df.index)

        regions_df.loc[mask, ['start', 'end']] -= del_len

        # Drop any whole regions that fall entirely inside the deletion
        if idx_end is not None and idx_end > idx_start:
            drop_slice = regions_df.index[idx_start : idx_end]
        elif idx_end is None:
            chr_mask  = regions_df['chr'] == chrom
            all_idxs  = regions_df.index[chr_mask]
            drop_slice = all_idxs[all_idxs >= idx_start]
        else:
            drop_slice = []

        regions_df = regions_df.drop(index=drop_slice).reset_index(drop=True)
        return regions_df

    # ---- deletion start falls *inside* a region -----------------------------
    if not both_start and idx_start is not None:

        if idx_end is None:
            # Deletion extends past all remaining regions on this chromosome
            chr_mask = regions_df['chr'] == chrom
            last_idx = regions_df.index[chr_mask].max()

            # Truncate the region that contains del_start
            regions_df.loc[idx_start, 'end'] = del_start

            # Drop everything after it on this chromosome
            drop_slice = regions_df.index[idx_start + 1 : last_idx + 1]
            regions_df = regions_df.drop(index=drop_slice).reset_index(drop=True)
            return regions_df

        if idx_start == idx_end:
            # Deletion contained within a single region
            mask = (regions_df.index > idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask, ['start', 'end']] -= del_len
            mask_at = (regions_df.index == idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask_at, 'end'] -= del_len
            return regions_df

        # idx_end > idx_start
        if both_end:
            # Deletion ends between regions: shift everything from idx_end onward
            mask = (regions_df.index >= idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask, ['start', 'end']] -= del_len

            # Truncate the region containing del_start
            regions_df.loc[idx_start, 'end'] = del_start

            # Drop fully-covered regions between idx_start and idx_end
            drop_slice = regions_df.index[idx_start + 1 : idx_end]
            regions_df = regions_df.drop(index=drop_slice).reset_index(drop=True)
        else:
            # Deletion ends inside another region
            mask = (regions_df.index > idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask, ['start', 'end']] -= del_len

            regions_df.loc[idx_start, 'end'] = del_start
            regions_df.loc[idx_end,   'end'] -= del_len

            # Drop fully-covered regions between the two partial ones
            drop_slice = regions_df.index[idx_start + 1 : idx_end]
            regions_df = regions_df.drop(index=drop_slice).reset_index(drop=True)

    return  _sanitize(regions_df)


# ---------------------------------------------------------------------------
# Copy-Number Variation (CNV / tandem duplication)
# ---------------------------------------------------------------------------

def adapt_CNV(regions_df, info):
    """
    Duplicate a segment `rep_num` times in tandem.

    info keys: chrom_name, start, end, rep_num
    """
    chrom   = info['chrom_name']
    seg_len = info['end'] - info['start']
    extra   = seg_len * (info['rep_num'] - 1)   # net coordinate increase

    if chrom not in regions_df['chr'].values:
        return regions_df

    # 1. Ensure clean boundaries at both CNV edges
    regions_df, idx_start = _split_region_at(regions_df, chrom, info['start'])
    regions_df, idx_end   = _split_region_at(regions_df, chrom, info['end'])

    if idx_start is None:
        return regions_df

    # 2. Determine the last index involved in the duplicated block
    if idx_end is not None:
        # Shift everything after the CNV block
        _, end_between = find_idx(regions_df, chrom, info['end'])
        mask_after = (regions_df.index >= idx_end) & (regions_df['chr'] == chrom)
        if end_between:
            regions_df.loc[mask_after, ['start', 'end']] += extra
        else:
            mask_after_strict = (regions_df.index > idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask_after_strict, ['start', 'end']] += extra
            mask_at = (regions_df.index == idx_end) & (regions_df['chr'] == chrom)
            regions_df.loc[mask_at, 'end'] += extra

        block_end = idx_end - 1 if end_between else idx_end
    else:
        # CNV extends to the end of the chromosome
        chr_mask  = regions_df['chr'] == chrom
        block_end = int(regions_df.index[chr_mask].max())

    # 3. Duplicate the block (rep_num - 1) times
    if idx_start <= block_end:
        block = regions_df.loc[idx_start : block_end].copy()
        duplicates = []
        for i in range(1, info['rep_num']):
            dup = block.copy()
            dup[['start', 'end']] += i * seg_len
            duplicates.append(dup)

        after = regions_df.loc[block_end + 1 :]
        regions_df = _safe_concat([regions_df.loc[:block_end], *duplicates, after])

    return regions_df


# ---------------------------------------------------------------------------
# Aneuploidy (whole-chromosome gain / loss)
# ---------------------------------------------------------------------------

def adapt_aneuploidy(regions_df, info):
    """
    Duplicate or delete a whole chromosome.

    info keys: chrom_name, rep_num (>0 = copies to add, 0 = deletion), chrom_length
    """
    chrom = info['chrom_name']

    if chrom not in regions_df['chr'].values:
        return regions_df

    if info['rep_num'] > 0:
        block    = regions_df[regions_df['chr'] == chrom].copy()
        last_idx = int(regions_df.index[regions_df['chr'] == chrom].max())
        duplicates = []
        for i in range(1, info['rep_num']):
            dup = block.copy()
            dup[['start', 'end']] += i * info['chrom_length']
            duplicates.append(dup)

        regions_df = _safe_concat([
            regions_df.loc[:last_idx],
            *duplicates,
            regions_df.loc[last_idx + 1:],
        ])
    else:
        regions_df = regions_df[regions_df['chr'] != chrom].reset_index(drop=True)

    return regions_df


# ---------------------------------------------------------------------------
# Translocation
# ---------------------------------------------------------------------------

def adapt_translocation(regions_df, info):
    """
    Reciprocal translocation between two chromosomes at their respective breakpoints.

    info keys: chrom_name1, chrom_name2, bkpt1, bkpt2, normal (bool)
    """
    chr1_df = regions_df[regions_df['chr'] == info['chrom_name1']].copy()
    chr2_df = regions_df[regions_df['chr'] == info['chrom_name2']].copy()

    offset1 = info['bkpt2'] - info['bkpt1']
    offset2 = info['bkpt1'] - info['bkpt2']

    slice1 = _split_and_shift(chr1_df, info['bkpt1'], info['chrom_name2'], offset1)
    slice2 = _split_and_shift(chr2_df, info['bkpt2'], info['chrom_name1'], offset2)

    # Remove the two original chromosomes
    rest = regions_df[
        ~regions_df['chr'].isin([info['chrom_name1'], info['chrom_name2']])
    ].copy()

    regions_df = _safe_concat([rest, slice2, slice1])

    # Sort using a proper chromosome order to avoid lexicographic issues
    chr_order = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
                 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
                 'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
    chr_rank = {c: i for i, c in enumerate(chr_order)}
    regions_df['_rank'] = regions_df['chr'].map(lambda c: chr_rank.get(c, 999))
    regions_df = (regions_df
                  .sort_values(['_rank', 'start', 'end'])
                  .drop(columns='_rank')
                  .reset_index(drop=True))

    return regions_df


# ---------------------------------------------------------------------------
# Chromothripsis
# ---------------------------------------------------------------------------

def adapt_chromothripsis(regions_df, info):
    """
    Rearrange segments produced by chromothripsis.

    info keys: chrom_name, start, end, bkpts (list of offsets from start),
               order (permutation list)
    """
    chrom = info['chrom_name']

    if chrom not in regions_df['chr'].values:
        return regions_df

    # 1. Split at every breakpoint so segment boundaries are clean
    absolute_bkpts = sorted(set(
        [info['start']] +
        [info['start'] + b for b in info['bkpts']] +
        [info['end']]
    ))

    for pos in absolute_bkpts:
        regions_df, _ = _split_region_at(regions_df, chrom, pos)

    # 2. Collect the resulting segments between consecutive breakpoints
    df_list = []
    for i in range(len(absolute_bkpts) - 1):
        seg_start = absolute_bkpts[i]
        seg_end   = absolute_bkpts[i + 1]
        seg = regions_df[
            (regions_df['chr'] == chrom) &
            (regions_df['start'] >= seg_start) &
            (regions_df['end']   <= seg_end)
        ].copy()
        df_list.append(seg)

    if not df_list:
        return regions_df

    # 3. Reorder segments according to `info['order']`
    order = [x for x in info['order'] if x < len(df_list)]
    if len(order) != len(df_list):
        order = order + [i for i in range(len(df_list)) if i not in order]

    new_pos = absolute_bkpts[0]
    shuffled = []
    for idx in order:
        seg = df_list[idx].copy()
        seg_len = absolute_bkpts[idx + 1] - absolute_bkpts[idx]
        offset  = new_pos - absolute_bkpts[idx]
        seg[['start', 'end']] += offset
        shuffled.append(seg)
        new_pos += seg_len

    # 4. Rebuild the dataframe around the rearranged window
    before = regions_df[
        (regions_df['chr'] == chrom) & (regions_df['end'] <= info['start'])
    ]
    after = regions_df[
        (regions_df['chr'] == chrom) & (regions_df['start'] >= info['end'])
    ]
    other_chrs = regions_df[regions_df['chr'] != chrom]

    regions_df = _safe_concat([other_chrs, before, *shuffled, after])

    chr_order = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
                 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
                 'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
    chr_rank = {c: i for i, c in enumerate(chr_order)}
    regions_df['_rank'] = regions_df['chr'].map(lambda c: chr_rank.get(c, 999))
    regions_df = (regions_df
                  .sort_values(['_rank', 'start', 'end'])
                  .drop(columns='_rank')
                  .reset_index(drop=True))

    return regions_df


# ---------------------------------------------------------------------------
# Breakage-Fusion-Bridge (BFB)
# ---------------------------------------------------------------------------

def adapt_BFB(regions_df, info):
    """
    Apply BFB cycle coordinate shifts.

    info keys: chrom_name, bkpts (list of breakpoint positions, in original coords)

    Each breakpoint triggers a fold-back: regions after bkpt[i] are shifted by
    the length of the segment [bkpt[i-1], bkpt[i]] (or bkpt[0] for the first).
    The breakpoints are in *original* coordinates, so we must translate them
    before each lookup.
    """
    chrom = info['chrom_name']

    if chrom not in regions_df['chr'].values:
        return regions_df

    cumulative_shift = 0  # track how much we've already shifted the coordinate space

    for i, bkpt_orig in enumerate(info['bkpts']):
        bkpt_current = bkpt_orig + cumulative_shift  # position in the *current* df

        idx, between = find_idx(regions_df, chrom, bkpt_current)
        if idx is None:
            continue

        seg_len = info['bkpts'][i] if i == 0 else info['bkpts'][i] - info['bkpts'][i - 1]

        # Shift all rows strictly after the breakpoint
        mask_after = (regions_df.index > idx) & (regions_df['chr'] == chrom)
        regions_df.loc[mask_after, ['start', 'end']] += seg_len

        # Shift the boundary row
        mask_at = (regions_df.index == idx) & (regions_df['chr'] == chrom)
        if between:
            regions_df.loc[mask_at, ['start', 'end']] += seg_len
        else:
            regions_df.loc[mask_at, 'end'] += seg_len

        cumulative_shift += seg_len

    return regions_df