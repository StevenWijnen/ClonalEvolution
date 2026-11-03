import subprocess
import os
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

BASE_PATH = Path("/data_pool_1/ClonalEvolution/oncoanalyser/output_full_runs/primary")

def annotate_vep(case_id : str, tumor_id:str, normal_id:str):
    ss = f"{tumor_id}_vs_{normal_id}"
    vcf_file = BASE_PATH / case_id / "purple" / f"{tumor_id}.purple.somatic.vcf.gz"
    output_file = BASE_PATH / case_id / "maf" /f"{case_id}.purple.pvacseq.vep.vcf"

    command = [
        "/home/s.wijnen/miniforge3/envs/pvactools/bin/vep",
        "--input_file", str(vcf_file),
        "--output_file", str(output_file),
        "--format", "vcf", 
        "--vcf",
        "--tsl",
        "--hgvs",
        "--plugin", "Wildtype",
        "--plugin", "Frameshift",
        "--cache", 
        "--offline",
        "--assembly", "GRCh38",
        "--fasta", "/DATA/s.wijnen/reference/human/oncoanalyser/GRCh38_masked_exclusions_alts_hlas.fasta",
        "--dir_cache", "/DATA/s.wijnen/reference/vep_cache",
        "--dir_plugins", "/home/s.wijnen/.vep/Plugins", # Use absolute path
        "--no_stats",
        "--fork", "2"
    ]

    subprocess.run(command, check=True)
    
    

def build_HLA_string(case_id, tumor_id):
    path = BASE_PATH / f"{case_id}/lilac/{tumor_id}.lilac.tsv"
    lilac_output = pd.read_csv(path,  sep="\t")
    
    unique_alleles = lilac_output["Allele"].unique().astype(str)
    allele_string = ",".join(f"HLA-{a}" for a in unique_alleles)
    
    return allele_string

def run_pvacseq(case_id:str, tumor_id:str, normal_id:str, hla_string:str, tumor_purity:str):
    input_vcf = BASE_PATH / case_id / "maf" /f"{case_id}.purple.pvacseq.vep.vcf"
    output_path = Path(f"/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/pvacseq/primary/{case_id}")
    output_path.mkdir(parents=True, exist_ok=True)
    
    command = [
        "pvacseq", "run",
        str(input_vcf),           
        tumor_id,  
        hla_string,
        "NetMHCpan", "NetMHCpanEL",
        str(output_path),         
        "--normal-sample-name",normal_id,
        "--tumor-purity", tumor_purity,
        "--minimum-fold-change", "1",
        "--binding-threshold", "500",
        "--tdna-vaf", "0",
        "--pass-only",
        "-e1", "8,9,10,11", 
        "-c", "0", #Foldchange difference
        "--iedb-install-directory", "/DATA/s.wijnen/reference/human/MHC/"
    ]
    print(command)
    subprocess.run(command, check=True)

def get_purity(case_id, tumor_id):
    path = BASE_PATH / f"{case_id}/purple/{tumor_id}.purple.qc"
    purple_qc = pd.read_csv(path,  sep="\t")
    return purple_qc.iloc[3,1]

def process_patient(row):
    case_id = row["case_id"]
    ln_tumor_id = row["primary_tumor_id"]
    normal_id = row["normal_id"]
    annotate_vep(case_id, ln_tumor_id, normal_id)
    hla_string = build_HLA_string(case_id, ln_tumor_id)
    purity = get_purity(case_id, ln_tumor_id)
    run_pvacseq(case_id, ln_tumor_id, normal_id, hla_string, purity)
    
    return f"Finished {case_id}"

def parallel_main():
    df = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/df_good_quality_patient_samples.csv")
    patient_data = df.to_dict(orient="records")
    with ProcessPoolExecutor(max_workers=16) as executor:
        results = list(executor.map(process_patient, patient_data))
    print("All patients processed.")

def main():
    df = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/df_good_quality_patient_samples.csv")
    patient_data = df.to_dict(orient="records")
    for x in patient_data:
        process_patient(x)
if __name__ == "__main__":
    parallel_main()        

