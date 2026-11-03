import os
import pandas as pd
from pathlib import Path
import subprocess

BASE_BAM_PATH = Path("/data_pool_1/ClonalEvolution/oncoanalyser/output_full_runs")
OUTPUT_DIR = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux")

def file_builder(site:str, case_id:str, sample_id:str) -> Path:
    return BASE_BAM_PATH / site / case_id / "alignments/dna" / f"{sample_id}.redux.bam"

def runfacets(site:str, case_id:str, tumor_id:str, normal_id:str)-> None:
    print("lets goooo")
    normal_bam_file = file_builder(site, case_id, normal_id)
    tumor_bam_file = file_builder(site, case_id, tumor_id)
    output_dir = OUTPUT_DIR / site
    
    output_prefix = os.path.join(output_dir, "pileup", f"{tumor_id}")    
    results_output = os.path.join(output_dir, "results", f"{tumor_id}")

    facets_outdir = os.path.join(results_output)
    os.makedirs(results_output, exist_ok=True)
    

    command1=[]
    command2=[]
    snp_pileup_exe = "/DATA/s.wijnen/Projects/SAREK_Additional_Modules/packages/FACETS/snp-pileup"
    
    VCF = "/DATA/s.wijnen/reference/human/clonalEvolution/Homo_sapiens_assembly38_doublefiltered.dbsnp138.vcf"
    print(normal_bam_file)
    print(tumor_bam_file)
    if os.path.isfile(normal_bam_file) and os.path.isfile(tumor_bam_file):    
        # First command: SNP-pileup
        # if not os.path.isfile(f"{output_prefix}.snp_pileup.gz"):
        command1 = [
            "Rscript", "/DATA/s.wijnen/Projects/SAREK_Additional_Modules/packages/FACETS/snp-pileup-wrapper.R", 
            "--snp-pileup-path", snp_pileup_exe, 
            "--vcf-file", VCF, 
            "--pseudo-snps", "100", #default 50 
            "--normal-bam", normal_bam_file,
            "--tumor-bam", tumor_bam_file,
            "--output-prefix", output_prefix, 
            "--min-read-counts", "30,10" #default 0,0
        ]
        # Second command: Run FACETS
        # if not os.path.isfile(f"{facets_outdir}/{tumor_id}_{normal_id}_hisens.rds"):
        #https://github.com/taylor-lab/facets-suite/wiki/1.-Running-FACETS
        #https://github.com/dariober/cnv_facets?tab=readme-ov-file#citation--getting-help
        command2 = [
                "Rscript", "/DATA/s.wijnen/Projects/facets/code/run-facets-wrapper.R", 
                "--counts-file", f"{output_prefix}.snp_pileup.gz",
                "--sample-id", f"{tumor_id}_{normal_id}",
                "--purity-cval", "1000", 
                "--cval", "500", 
                "--min-nhet", "25",
                "--purity-min-nhet", "25", 
                "--directory", facets_outdir,
                "--facets-lib-path", "/home/s.wijnen/R/x86_64-pc-linux-gnu-library/4.4",
                "--normal-depth", "50",
                "--snp-window-size", "250",
                "--genome", "hg38"
            ]

        # Run both commands
        for command in [command1, command2]:
            try:
                if len(command) >0:
                    subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                print(f"❌ Error processing {tumor_id}: {e}")
    else:
        print(f"⚠️ Missing BAM file for {tumor_id} or {normal_id}")



def main(): 
    samples = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/cleaned_full_df.csv")
    for _, row in samples.iterrows():
        case_id = str(row["Sample"]).strip()
        prim_tumor_id = str(row["DNA_primary_cegat"]).strip()
        ln_met_tumor_id = str(row["DNA_ln_cegat"]).strip()
        normal_id = str(row["DNA_normal_cegat"]).strip()
        
        runfacets("primary", case_id, prim_tumor_id, normal_id)
        runfacets("ln_met", case_id, ln_met_tumor_id, normal_id)
        

main()