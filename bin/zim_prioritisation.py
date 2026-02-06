#!/usr/bin/env python3
import pandas as pd
import argparse

def load_zim_db(path):
    """
    Zimbabwe cohort frequency table (TSV).

    Required columns (case-insensitive):
      - Chr
      - Start
      - Ref
      - Alt

    And either:
      - zim_AF
    or:
      - zim_AC and zim_AN (AF will be computed as AC/AN)
    """
    df = pd.read_csv(path, sep='\t')
    cols = {c.lower(): c for c in df.columns}

    required = ['chr', 'start', 'ref', 'alt']
    for r in required:
        if r not in cols:
            raise SystemExit(
                f"Zimbabwe DB file must contain a '{r}' column "
                f"(case-insensitive). Got: {df.columns.tolist()}"
            )

    # Ensure zim_AF exists, or compute from AC/AN
    if 'zim_af' not in cols:
        if 'zim_ac' in cols and 'zim_an' in cols:
            c_ac = cols['zim_ac']
            c_an = cols['zim_an']
            # make sure numeric
            df[c_ac] = pd.to_numeric(df[c_ac], errors='coerce')
            df[c_an] = pd.to_numeric(df[c_an], errors='coerce')
            df['zim_AF'] = df[c_ac] / df[c_an]
        else:
            raise SystemExit(
                "Zimbabwe DB must have zim_AF or zim_AC + zim_AN "
                "(case-insensitive)."
            )
    else:
        c_af = cols['zim_af']
        if c_af != 'zim_AF':
            df.rename(columns={c_af: 'zim_AF'}, inplace=True)

    # Standardise coordinate columns to match ANNOVAR
    df.rename(columns={
        cols['chr']: 'Chr',
        cols['start']: 'Start',
        cols['ref']: 'Ref',
        cols['alt']: 'Alt'
    }, inplace=True)

    return df[['Chr', 'Start', 'Ref', 'Alt', 'zim_AF']]


def _make_gnomad_max_col(df, candidates, new_name):
    """
    From a list of candidate gnomAD AF columns, create a single column
    that's the row-wise max across them (ignoring NaN).
    """
    if not candidates:
        return None

    # ensure numeric
    for c in candidates:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    df[new_name] = df[candidates].max(axis=1, skipna=True)
    return new_name


def prioritise(multianno_path, zim_db_path, out_prefix,
               max_gnomad_af=0.01,
               max_gnomad_afr_af=0.01,
               max_zim_af=0.05):

    print(f"Reading ANNOVAR multianno file: {multianno_path}")
    ann = pd.read_csv(multianno_path, sep='\t', low_memory=False)
    print(f"  n = {len(ann)} variants before filtering")

    # Core columns
    for col in ['Chr', 'Start', 'End', 'Ref', 'Alt']:
        if col not in ann.columns:
            raise SystemExit(f"Missing required column '{col}' in multianno file.")

    # Your multianno has gnomAD exome + genome and duplicated AF columns.
    # We'll detect AF columns and create two meta-columns:
    #   - gnomad_AF_max     (max over global AF columns)
    #   - gnomad_AFR_AF_max (max over AFR AF columns)
    gnomad_all_candidates = [
        c for c in ann.columns
        if c.startswith('AF') and ('afr' not in c.lower())
    ]
    gnomad_afr_candidates = [
        c for c in ann.columns
        if 'afr' in c.lower() and c.startswith('AF')
    ]

    gnomad_all = _make_gnomad_max_col(ann, gnomad_all_candidates, 'gnomad_AF_max')
    gnomad_afr = _make_gnomad_max_col(ann, gnomad_afr_candidates, 'gnomad_AFR_AF_max')

    print("Using gnomAD columns:")
    print(f"  global AF (max over {gnomad_all_candidates}) -> {gnomad_all}")
    print(f"  AFR AF   (max over {gnomad_afr_candidates}) -> {gnomad_afr}")

    # Load Zimbabwe cohort DB
    print(f"Reading Zimbabwe cohort DB: {zim_db_path}")
    zim = load_zim_db(zim_db_path)
    print(f"  n = {len(zim)} distinct cohort variants")

    # Merge on genomic position + alleles
    merged = ann.merge(
        zim,
        on=['Chr', 'Start', 'Ref', 'Alt'],
        how='left'
    )

    # Variants not in Zim DB → zim_AF = 0
    merged['zim_AF'] = merged['zim_AF'].fillna(0.0)

    # Impact categories
    func_col = 'Func.refGene'
    exonic_func_col = 'ExonicFunc.refGene'

    if func_col not in merged.columns or exonic_func_col not in merged.columns:
        print("WARNING: Func.refGene / ExonicFunc.refGene not found; "
              "skipping consequence-based filters.")
        merged['impact_tier'] = 'unknown'
    else:
        lof_terms = {
            'frameshift_deletion', 'frameshift_insertion',
            'frameshift_variant', 'stopgain', 'stoploss',
            'splicing', 'splice_site'
        }
        missense_terms = {'nonsynonymous SNV', 'missense_variant'}

        cond_lof = merged[exonic_func_col].isin(lof_terms) | \
                   merged[func_col].str.contains('splicing', na=False)
        cond_missense = merged[exonic_func_col].isin(missense_terms)
        cond_synonymous = merged[exonic_func_col].isin(
            {'synonymous SNV', 'synonymous_variant'}
        )

        merged['impact_tier'] = 'other'
        merged.loc[cond_synonymous, 'impact_tier'] = 'synonymous'
        merged.loc[cond_missense, 'impact_tier'] = 'missense'
        merged.loc[cond_lof, 'impact_tier'] = 'LoF'

    # Frequency filters: gnomAD global + AFR
    freq_mask = True
    if gnomad_all is not None:
        freq_mask = freq_mask & (
            merged[gnomad_all].isna() | (merged[gnomad_all] <= max_gnomad_af)
        )
    if gnomad_afr is not None:
        freq_mask = freq_mask & (
            merged[gnomad_afr].isna() | (merged[gnomad_afr] <= max_gnomad_afr_af)
        )

    # Zimbabwe cohort filter (this is the Zi-mendelian twist)
    zim_mask = merged['zim_AF'] <= max_zim_af

    # Impact filter: keep LoF + missense by default
    impact_mask = merged['impact_tier'].isin(['LoF', 'missense'])

    final_mask = freq_mask & zim_mask & impact_mask

    filtered = merged[final_mask].copy()
    print(f"  n = {len(filtered)} variants after global+AFR+Zimbabwe AF and impact filters")

    # Flags for debugging / reporting
    filtered['flag_high_zim_af'] = merged['zim_AF'] > 0.01
    if gnomad_all is not None:
        filtered['flag_common_gnomad'] = merged[gnomad_all] > 0.01
    if gnomad_afr is not None:
        filtered['flag_common_gnomad_afr'] = merged[gnomad_afr] > 0.01

    # Outputs
    out_all = f"{out_prefix}.zim_annotated.tsv"
    out_filtered = f"{out_prefix}.zim_prioritised.tsv"

    merged.to_csv(out_all, sep='\t', index=False)
    filtered.to_csv(out_filtered, sep='\t', index=False)

    print(f"Wrote annotated file (with zim_AF)  : {out_all}")
    print(f"Wrote prioritised variant list      : {out_filtered}")


def main():
    ap = argparse.ArgumentParser(
        description="Zimbabwe-specific prioritisation on ANNOVAR multianno output."
    )
    ap.add_argument("--multianno", required=True,
                    help="ANNOVAR *hg38_multianno.txt file")
    ap.add_argument("--zim_db", required=True,
                    help="Zimbabwe Mendelian cohort frequency table (TSV)")
    ap.add_argument("--out_prefix", required=True,
                    help="Prefix for output files")

    ap.add_argument("--max_gnomad_af", type=float, default=0.01,
                    help="Max global gnomAD AF (default 0.01)")
    ap.add_argument("--max_gnomad_afr_af", type=float, default=0.01,
                    help="Max AFR gnomAD AF (default 0.01)")
    ap.add_argument("--max_zim_af", type=float, default=0.05,
                    help="Max AF in Zimbabwe cohort (default 0.05)")

    args = ap.parse_args()

    prioritise(
        args.multianno,
        args.zim_db,
        args.out_prefix,
        max_gnomad_af=args.max_gnomad_af,
        max_gnomad_afr_af=args.max_gnomad_afr_af,
        max_zim_af=args.max_zim_af
    )


if __name__ == "__main__":
    main()
