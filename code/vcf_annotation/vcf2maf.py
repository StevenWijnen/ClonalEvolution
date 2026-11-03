import pandas as pd 
import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

samples = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/cleaned_full_df.csv")

def path_builder(source, group_id, tumor_id):
    path = f"/data_pool_1/ClonalEvolution/oncoanalyser/output_full_runs/{source}/{group_id}/purple/{tumor_id}.purple.somatic.vcf.gz"
    return path    

def process_sample(site, row):
    group_id = str(row['Sample']).strip()
    if site == "primary":
        tumor_id = str(row["DNA_primary_cegat"]).strip()
    else:
        tumor_id = str(row["DNA_ln_cegat"]).strip()
    normal_id = str(row["DNA_normal_cegat"]).strip()
    tumor_vcf = path_builder(source=site, group_id=group_id, tumor_id=tumor_id)        
    
    if not os.path.isfile(tumor_vcf):
        print(f"VCF not found, skipping: {tumor_vcf}")
        return

    output_dir = f"/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/{site}"
    os.makedirs(output_dir, exist_ok=True)
    
    uncompressed_vcf = os.path.join(output_dir, f"{tumor_id}.purple.somatic.vcf")
    if not os.path.exists(uncompressed_vcf):
        print(f"Unzipping {tumor_vcf} to {uncompressed_vcf}")
        try:
            subprocess.run(
                f"gunzip -c {tumor_vcf} > {uncompressed_vcf}",
                shell=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Unzipping failed for {tumor_id}: {e}")
            return

    output_maf = os.path.join(output_dir, f"{tumor_id}.purple.maf")

    command = [
        "perl", "/DATA/s.wijnen/Projects/SAREK_Additional_Modules/packages/vcf2maf-1.6.22/vcf2maf.pl",
        "--input-vcf", uncompressed_vcf,
        "--output-maf", output_maf,
        "--tumor-id", tumor_id,
        "--normal-id", normal_id,
        "--vep-path", "/home/s.wijnen/miniforge3/envs/vcf2maf_clean/bin",
        "--ref-fasta", "/DATA/s.wijnen/reference/human/oncoanalyser/GRCh38_masked_exclusions_alts_hlas.fasta", 
        "--vep-data", "/DATA/s.wijnen/reference/vep_cache",
        "--retain-info", "ROQ",
        "--ncbi-build=GRCh38"
    ]

    print(f"Running vcf2maf for {tumor_id}")
    try: 
        subprocess.run(command, check=True)
        if os.path.exists(uncompressed_vcf):
            print(f"Success! Removing uncompressed VCF: {uncompressed_vcf}")
            os.remove(uncompressed_vcf)
    except subprocess.CalledProcessError as e:
        print(f"vcf2maf failed for {tumor_id}: {e}")

# Build list of all jobs
jobs = [
    (site, row)
    for site in ["primary", "ln_met"]
    for _, row in samples.iterrows()
]

MAX_WORKERS = 34

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(process_sample, site, row): (site, row) for site, row in jobs}
    for future in as_completed(futures):
        site, row = futures[future]
        try:
            future.result()
        except Exception as e:
            print(f"Unexpected error for {site} / {row['Sample']}: {e}")