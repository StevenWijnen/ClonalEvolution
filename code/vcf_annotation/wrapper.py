import os
import pandas as pd
from pathlib import Path
import subprocess


PATH_TO_FACETS_PRIMARY = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux/primary/results")
PATH_TO_FACETS_LN = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux/ln_met/results")
MAF_PATH = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/")


OUTPUT_PATH = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/annotated_mafs")


def read_normal_maf(path):
    if not os.path.isfile(path):
        print(f"Missing file: {path}")
        raise Exception("Maf output not found")
    df = pd.read_csv(path, delimiter="\t", comment="#")

    return df


def filter_pass(maf_data):
    return maf_data[(maf_data["FILTER"] == "PASS")]

def filter_t_alt(maf_data):
    return maf_data[(maf_data["t_alt_count"] >= 5)]

def filter_ROQ(maf_data):
    return maf_data[(maf_data["ROQ"] >= 70)]

def filter_microsec(maf_data, tumor_id, normal_id,green_list):
    # Load microsec file
    microsec_file = Path(
        f"/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/microsec_cleaned/"
        f"{tumor_id}_vs_{normal_id}.tsv"
    )
    if os.path.isfile(microsec_file):
        print("lets gooo")
        microsec = pd.read_csv(microsec_file, sep="\t")
    else:
        return maf_data
    
    # Keep only rows where msec_filter_all_called == True
    microsec["msec_filter_all_called"]  =  microsec["msec_filter_all"] == "Artifact suspicious"
    microsec_filtered = microsec[microsec["msec_filter_all_called"]]

    # Create a key column for easy matching
    microsec_filtered["key"] = (
        microsec_filtered["Chr"].astype(str) + "_" + microsec_filtered["Pos"].astype(str)
    )

    maf_data["key"] = (
        maf_data["Chromosome"].astype(str) + "_" + maf_data["vcf_pos"].astype(str)
    )

    # Build a set for fast lookup
    to_filter_set = set(microsec_filtered["key"])

    # Flag rows that should be filtered
    maf_data["microsec_flag"] = maf_data["key"].isin(to_filter_set)

    maf_data = maf_data[~maf_data["microsec_flag"] | maf_data["snv_id"].isin(green_list[0])]


    return maf_data

def facets_annotate_maf(maf_file, tumor_id, normal_id,  ffpe=False ):
    
    maf_file = OUTPUT_PATH /  f"{tumor_id}_filtered.maf"
    if ffpe:          
       facets_output = PATH_TO_FACETS_PRIMARY / tumor_id / f"{tumor_id}_{normal_id}_hisens.rds"
    else:
        facets_output = PATH_TO_FACETS_LN /tumor_id / f"{tumor_id}_{normal_id}_hisens.rds"

    if not os.path.isfile(facets_output):
        print(f"Missing file: {facets_output}")
        raise Exception("Facets output not found")
    output_dir = OUTPUT_PATH / f"{tumor_id}_facets_annotated_filtered.maf"
    command = [
            "Rscript", "/DATA/s.wijnen/Projects/facets/code/annotate-maf-wrapper.R", 
            "--maf-file", maf_file, 
            "--facets-output", facets_output, 
            "--output", output_dir,
        ]
        
    subprocess.run(command, check=True)
    subprocess.run(["rm", maf_file])

def write_maf(maf, tumor_id):
    file = OUTPUT_PATH / f"{tumor_id}_filtered.maf"
    maf.to_csv(file, sep="\t", index=False)
    return file
    

def annotate_ffpe_mafs(maf, tumor_id, normal_id ):    
    maf = filter_pass(maf)
    maf = filter_t_alt(maf)
    #TODO Check for duplicates?  Keep the highest impact variant
    # maf = filter_ROQ(maf, green_list)
    # maf = filter_microsec(maf, tumor_id, normal_id, green_list)
    if maf.shape[0] == 0:
        print(tumor_id)
        return
    file = write_maf(maf, tumor_id)
    
    facets_annotate_maf(file, tumor_id, normal_id, ffpe=True)
        
def annotate_ff_mafs(maf,tumor_id, normal_id):

    maf = filter_pass(maf)
    maf = filter_t_alt(maf)
    
    #TODO Check for duplicates? Keep the highest impact variant
    if maf.shape[0] == 0:
        print(tumor_id)
        return
    file = write_maf(maf, tumor_id)    
    facets_annotate_maf(file, tumor_id, normal_id, ffpe=False) # fix this one
    


df = pd.read_csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")

for _, row in df.iterrows():
    case_id = row["Sample"]
    primary_id = row["DNA_primary_cegat"]
    ln_met_id = row["DNA_ln_cegat"]
    normal_id = row["DNA_normal_cegat"]
    
    primary_path = MAF_PATH / f"primary/mafs/{primary_id}.purple.maf"
    ln_met_path = MAF_PATH / f"ln_met/mafs/{ln_met_id}.purple.maf"
    
    # green_list_path = Path("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/greenlists")
    
    
    if not os.path.isfile(primary_path):
        print(f"Missing file: {primary_path}")
        continue
    
    if not os.path.isfile(ln_met_path):
        print(f"Missing file: {ln_met_path}")
        continue
    
    
    primary_maf = pd.read_csv(primary_path, delimiter="\t", comment="#")
    # primary_maf["snv_id"] = (
    #     primary_maf.loc[:, "Chromosome"].astype(str) + "_" + primary_maf.loc[:, "vcf_pos"].astype(str)
    # )
    
    ln_maf = pd.read_csv(ln_met_path, delimiter="\t", comment="#")
    
    # ln_maf["snv_id"] = (
    #     ln_maf.loc[:,"Chromosome"].astype(str) + "_" + ln_maf.loc[:,"vcf_pos"].astype(str)
    # )
    # try:
    #     green_list = pd.read_csv(green_list_path / f"{case_id}_greenlist.txt", header=None)
    # except pd.errors.EmptyDataError:
    #     green_list = pd.DataFrame(columns=[0])
         
    annotate_ffpe_mafs(primary_maf,primary_id, normal_id)
    annotate_ff_mafs(ln_maf,ln_met_id, normal_id)
    