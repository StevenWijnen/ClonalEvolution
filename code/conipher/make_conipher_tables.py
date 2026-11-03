from pathlib import Path
import pandas as pd
import numpy as np 
import os 
import pysam 
import pyranges as pr

BASE_PATH = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/annotated_mafs")
# BASE_FACETS_LN_MET = Path("/DATA/s.wijnen/Projects/ClonalEvolution/output/ln_met/facets/results")
BASE_FACETS_LN_MET = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux/ln_met/results")
#BASE_FACETS_PRIMARY = Path("/DATA/s.wijnen/Projects/ClonalEvolution/output/primary/facets/results")
BASE_FACETS_PRIMARY = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux/primary/results")

OUTPUT_PATH = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/conipher_tables_full_oncoanalyser")
BAM_FILES = Path("/data_pool_1/ClonalEvolution/oncoanalyser/output_full_runs")

def get_copy_number_from_segments(mut_df, segments_df):    
    maf_subset = mut_df[['Chromosome', 'vcf_pos', "tcn", "lcn"]].copy()
    maf_subset['Start'] = maf_subset['vcf_pos'] - 1
    maf_subset['End'] = maf_subset['vcf_pos']
    
    segments_df = segments_df.rename(columns={
        'start': 'Start', 
        'end': 'End', 
        'chrom': 'Chromosome'
    }).copy()
    
    segments_df['Start'] = segments_df['Start'] - 1
    
    mut_pr = pr.PyRanges(maf_subset)
    seg_pr = pr.PyRanges(segments_df)
    
    joined = mut_pr.join(seg_pr, how="left").as_df()

    if 'tcn_b' in joined.columns:
        joined['tcn'] = joined['tcn'].fillna(joined['tcn.em'])
    joined['tcn'] = joined['tcn'].replace(-1, np.nan).fillna(2)
    
    if 'lcn_b' in joined.columns:
        joined['lcn'] = joined['lcn'].fillna(joined['lcn.em'])
    
    joined.loc[joined['lcn'].replace(-1,np.nan).isna(), 'lcn'] = 1
    
    if ((joined['lcn'] == -1 ).any()):
        print(joined[joined['lcn'] == -1])
        input()
    
        
    return joined[['Chromosome', 'vcf_pos', 'tcn', 'lcn']]


def get_bam_counts(bam_path, chrom, pos, ref_allele, alt_allele):
    """
    Returns (ref_count, alt_count) at a 1-based genomic position.
    """
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        ref_count = 0
        alt_count = 0
        
        # VCF is 1-based; Pysam is 0-based. 
        target_pos = pos - 1 

        for pileupcolumn in samfile.pileup("chr" + str(chrom), target_pos, target_pos + 1, truncate=True):
            if pileupcolumn.pos == target_pos:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        aln = pileupread.alignment
                        if aln.is_duplicate or aln.is_secondary or aln.is_supplementary:
                            continue
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == ref_allele:
                            ref_count += 1
                        elif base == alt_allele:
                            alt_count += 1

        return ref_count, alt_count


def make_conipher_table(case_id, primary_id, ln_met_id, normal_id): 
    primary_path = BASE_PATH / f"{primary_id}_facets_annotated_filtered.maf"
    ln_met_path = BASE_PATH / f"{ln_met_id}_facets_annotated_filtered.maf"

    if not (os.path.exists(primary_path) and os.path.exists(ln_met_path) ):
        print("No bueno 1")
        print(primary_path)
        print(ln_met_path)
        return
    
    primary_maf = pd.read_csv(primary_path, delimiter="\t")
    ln_met_maf = pd.read_csv(ln_met_path, delimiter="\t")
    
    # Force Chromosome to string (and strip whitespace) and vcf_pos to integer
    #Somehow this was neccessary as for example on PRADO_164 it was not droping duplicates for chromsome in LN met maf was an object?

    for df in [primary_maf, ln_met_maf]:
        df['Chromosome'] = df['Chromosome'].astype(str).str.strip()
        df['vcf_pos'] = df['vcf_pos'].astype(int)
    
    facets_file_primary = BASE_FACETS_PRIMARY / f"{primary_id}/{primary_id}_{normal_id}.txt"
    facets_file_ln = BASE_FACETS_LN_MET / f"{ln_met_id}/{ln_met_id}_{normal_id}.txt"


    if not (os.path.exists(facets_file_primary) and os.path.exists(facets_file_ln) ):
        print("No bueno 2")
        return
    
    data_primary = pd.read_csv(facets_file_primary, delimiter="\t")
    data_ln_met =  pd.read_csv(facets_file_ln, delimiter="\t")

    purity_primary = data_primary.iloc[1, 2]
    purity_ln_met = data_ln_met.iloc[1, 2]
    
    
    ploidy_primary = data_primary.iloc[1, 3]
    ploidy_ln_met = data_ln_met.iloc[1, 3]
    
    
    segments_primary = pd.read_csv(BASE_FACETS_PRIMARY /f"{primary_id}/cncf.csv")
    segments_ln_met = pd.read_csv(BASE_FACETS_LN_MET /  f"{ln_met_id}/cncf.csv")
    
    cols_to_keep = ["Chromosome", "vcf_pos", "Reference_Allele", "Tumor_Seq_Allele2"]
    
    master_muts = pd.concat([
            primary_maf[cols_to_keep], 
            ln_met_maf[cols_to_keep]
        ]).drop_duplicates(subset=["Chromosome", "vcf_pos"])    

    # 2. Use the streamlined function
    def prepare_conipher_data(master_muts, maf, sample_id, bam_path, purity, segments_df):
       
        df = master_muts.merge(
                maf.drop(columns=["Reference_Allele", "Tumor_Seq_Allele2"]),
                on=["Chromosome", "vcf_pos"],
                how="left"
            )
        cn_df = get_copy_number_from_segments(df, segments_df)
        df = df.merge(cn_df, on=["Chromosome", "vcf_pos"], how="left")
        df['tcn'] = df['tcn_y']
        df['lcn'] = df['lcn_y']

        df = df.drop(columns=['tcn_x', 'tcn_y', 'lcn_x', 'lcn_y'], errors='ignore')


        missing_mask = df["t_alt_count"].isna()
        if missing_mask.any():
            print(f"Querying BAM for {missing_mask.sum()} missing mutations in {sample_id}...")
            results = df[missing_mask].apply(
                lambda row: get_bam_counts(bam_path, row['Chromosome'], row['vcf_pos'], 
                                        row['Reference_Allele'], row['Tumor_Seq_Allele2']),
                axis=1
            )
            
            # Unpack the (ref, alt) tuples into the dataframe
            df.loc[missing_mask, "t_ref_count"] = [r[0] for r in results]
            df.loc[missing_mask, "t_alt_count"] = [r[1] for r in results]
            
        # Standard cleanup for the rest
        df["t_depth"] = df["t_ref_count"] + df["t_alt_count"]

        df["Tumor_Sample_Barcode"] = sample_id
        df["purity"] = df["purity"].fillna(purity)
        #Old version I just set these as a default. Now I look them up in the facets output
        # df["tcn"] = df["tcn"].fillna(2)
        # df["lcn"] = df["lcn"].fillna(1)
        return df        
        
    def make_table(maf_file, ploidy):
        rc = maf_file.shape[0]
        conipher_df = pd.DataFrame({
                "CASE_ID": [case_id] * rc,
                "SAMPLE": case_id + "_" + maf_file['Tumor_Sample_Barcode'], 
                "CHR": maf_file["Chromosome"],
                "POS": maf_file["vcf_pos"],
                "REF": maf_file["Reference_Allele"],
                "ALT": maf_file["Tumor_Seq_Allele2"],
                "REF_COUNT": maf_file["t_ref_count"],
                "VAR_COUNT": maf_file["t_alt_count"],
                "DEPTH": maf_file["t_depth"],
                # Ensure column names used in the calculation exist in maf_file
                "COPY_NUMBER_A": maf_file["tcn"] - maf_file["lcn"], 
                "COPY_NUMBER_B": maf_file["lcn"],
                "ACF": np.where(maf_file["purity"].isna(), 0.25, maf_file["purity"]),
                "PLOIDY": [ploidy] * rc 
            })

        conipher_df = conipher_df.dropna(subset=['COPY_NUMBER_A', 'COPY_NUMBER_B'])
            

            
        return conipher_df
    
    prim_bam_path = BAM_FILES / "primary" / case_id / "alignments/dna" /  f"{primary_id}.redux.bam"
    ln_met_bam_path = BAM_FILES / "ln_met" / case_id / "alignments/dna" /  f"{ln_met_id}.redux.bam"

    primary_filled = prepare_conipher_data(master_muts, primary_maf, primary_id, prim_bam_path , purity_primary, segments_primary)
    ln_met_filled = prepare_conipher_data(master_muts, ln_met_maf, ln_met_id, ln_met_bam_path, purity_ln_met, segments_ln_met)
    
    # 3. Use your existing make_table logic on these "filled" dataframes
    table_primary = make_table(primary_filled, ploidy_primary)
    table_ln_met = make_table(ln_met_filled, ploidy_ln_met)
    
    
    combined_df = pd.concat([table_primary, table_ln_met], axis=0)
    
    
    combined_df.to_csv(f"{OUTPUT_PATH}/{case_id}.tsv", sep="\t")
    

sample_mapper = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")

for index, row in sample_mapper.iterrows():
    case_id = row["Sample"]
    primary_id = row["DNA_primary_cegat"]
    ln_met_id = row["DNA_ln_cegat"]
    normal_id = row["DNA_normal_cegat"]
    if(case_id == "PRADO_117"):
        make_conipher_table(case_id, primary_id, ln_met_id, normal_id)