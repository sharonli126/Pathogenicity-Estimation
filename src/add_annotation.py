# this file have all the function for adding the REVEL, MGI, gnomad_4_1 gene constraint, am_pathgencity, and gene_ontonology

### IMPORTANT please look at download.ymal to download all required database before adding anntation with functions in this file
### IMPORTANT please create a duckdb file to run the am_pathogencity annotation, since it is too huge

import pandas as pd
import duckdb as db

### REVEL
# obtaining the revel dictionary
def create_revel_dictionary(revel_file_path): # a text file

   revel_dictionary = {}
   try:
       with open(revel_file_path,'r') as rct:
           for line in rct:
               tmp = line.strip().split(',')
               chr="chr"+tmp[0]
               pos=":"+tmp[2]
               ref= tmp[3]
               var=">"+tmp[4]
               revel= tmp[7]
               dict_key= chr+pos+ref+var
               dict_score= revel
               revel_dictionary[dict_key]=dict_score

       return (revel_dictionary )     
    
   except:
       return {} # a dictionary with {genomic_nomencalture: revel score}

# please obtain ids_dictionary from create_revel_dictionary 
# adding revel value
def add_revel (ids_dictionary, dataset): # the revel dictionary and patient data (a data fream variable)

    dataset['revel'] = pd.NA

    for i in range(len(dataset)):
        value = dataset.at[i, 'genomic_nomencalture']
        if value in ids_dictionary:
            dataset.at[i, "revel"] = ids_dictionary[value]
        else:
            dataset.at[i, "revel"] = pd.NA

    return dataset # a dataset with revel scores


### gene ontonlogy 
def read_gmt(file_path): # a text file
    # Initialize lists to store data
    gmt_dict = {}

    # Read the GMT file line by line
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            description = parts[0]
            url = parts[1]
            genes = parts[2:]
            for i in genes:
                if i in gmt_dict:
                    gmt_dict[i].append(description)
                else:
                    gmt_dict[i] = [description]

    return gmt_dict # a dictionary with {gene: description}


def filter_gene_set(gmt, pattern): # pattern are keyword would like to obtain
    filtered_dict = {}
    
    for gene_set, data in gmt.items():
        if pattern.match(gene_set) or pattern.match(data['Description']):
            filtered_dict[gene_set] = data

    return filtered_dict

# obtain the gene_dict by read_gmt first
# adding gene_ontonology (term) into the dataset
def add_gmt (data, gene_dict):
    data["term"] = pd.NA
    for i in range(len(data)):
        symbol = data.at[i, 'SYMBOL']
        if symbol in gene_dict:
            data.at[i, 'term'] = gene_dict[symbol]
        else:
            pass
    return data # data with added gene ontonology term

        
### mgi
def mgi_data(file_path_mouse_marker, file_path_human_mouse_symbol): # two csv files, one with mouse marker + notes, another one with human to mouse gene symbol
    mouse_marker = pd.read_csv(file_path_mouse_marker, sep=',')
    mgi_df = mouse_marker[['symbol', 'name', 'note']]
    column_names = ["human_marker_symbol", "human_gene_id", "mouse_marker_symbol", "mgi_marker_id", "mammalian_phenotype_id"]
    human_mouse_symbol = pd.read_csv(file_path_human_mouse_symbol, sep='\t', header=None, index_col=False, names = column_names)

    return mgi_df, human_mouse_symbol # return filtered mouse marker and notes, and unchanged human to mouse gene symbol  csv

# please obtain mgi_df and human_mouse_symbol by mgi_data first
def add_mgi (data, mgi_df, human_mouse_symbol): # patient data, filtered mouse marker and ntoes, and human to mouse symbol csv
    result = human_mouse_symbol.merge(mgi_df[['symbol', 'note']], 
                                  left_on='mouse_marker_symbol', right_on='symbol', 
                                  how='left')

    result = result.dropna(subset=['symbol'])
    result.drop(columns=['symbol'], inplace=True)

    data["mouse_symbol"] = pd.NA
    data["mgi_notes"] = pd.NA

    for index, value in enumerate(data["SYMBOL"]):
        # Find matching rows in mouse_marker_filtered (mgi_df)
        matching_rows = result[result["human_marker_symbol"] == value]
        
        if not matching_rows.empty:
            # If a match is found, update the corresponding row in data
            data.at[index, 'mouse_symbol'] = matching_rows.iloc[0]["mouse_marker_symbol"]
            data.at[index, 'mgi_notes'] = []
            for i in range(len(matching_rows)):
                data.at[index, 'mgi_notes'].append(matching_rows.iloc[i]["note"])
    return data

### am_pathogencity
def am_df(con, file_path):
    # since alpha missense dataset is huge, using duck db is more efficient
    con.execute(f"""
        CREATE TABLE alpha_missense AS 
        SELECT * FROM read_csv_auto('{file_path}', sep='\t');
    """)

    # creating genomic_nomencalture in am dataset
    con.execute("""
        ALTER TABLE alpha_missense 
        ADD COLUMN genomic_nomencalture STRING;

        UPDATE alpha_missense
        SET genomic_nomencalture = "#CHROM" || ':' || "POS" || REF || '>' || ALT; 
    """)    

# please run am_df first
# add am_pathogenecity
def add_am_path(con,tablename):
    con.execute(f"""
    UPDATE {tablename}
    SET am_pathogenicity = am.am_pathogenicity
    FROM alpha_missense AS am
    WHERE {tablename}.genomic_nomencalture = am.genomic_nomencalture;
    """)  # am score will be added when the patient data is created in duckdb


# Gnomad 4.1
# please download the gnomad version 4.1 file
def filted_gnomad_4_1 (file_path, file_path_to_store = None):
    # this function is to obtain valid gomad_4_1 records with observed patterns
    gnomad_4_1 = pd.read_csv(file_path, sep = '\t')

    # Filter rows where canonical and mane_select are both TRUE
    gene_with_TT = gnomad_4_1[
        ((gnomad_4_1['canonical'] == True) & (gnomad_4_1['mane_select'] == True)) 
    ]
    genes_with_multiple_rows = gene_with_TT.groupby('gene').filter(lambda group: len(group) > 1)
    genes_with_multiple_rows_ENSG = genes_with_multiple_rows[genes_with_multiple_rows['gene_id'].str.contains('ENSG')]
    genes_with_single_row = gene_with_TT.groupby('gene').filter(lambda group: len(group) == 1)
    gene_with_missing = gene_with_TT[(gene_with_TT['gene']).isna()] 
    filtered_tt = pd.concat([genes_with_multiple_rows_ENSG, genes_with_single_row, gene_with_missing])


    # Filter rows where canonical and mane_select are either TRUE
    gene_in_filtered_tt = filtered_tt['gene_id'].unique()
    remaining_genes_df = gnomad_4_1[~gnomad_4_1['gene_id'].isin(gene_in_filtered_tt)]
    gene_with_TorT = remaining_genes_df[
        ((remaining_genes_df['canonical'] == True) | (remaining_genes_df['mane_select'] == True)) &
        (remaining_genes_df['gene_id'].str.contains('ENSG'))
    ]

    filtered_gnomad_4_1 = pd.concat([filtered_tt, gene_with_TorT])
    # Drop rows where 'gene_id' is either 29970 or 7730
    filtered_gnomad_4_1 = filtered_gnomad_4_1[~filtered_gnomad_4_1['gene_id'].isin([29970, 7730])]

    # Store the filtered file
    if file_path_to_store:
        filtered_gnomad_4_1.to_csv(file_path_to_store, sep='\t', index=False)

    return filtered_gnomad_4_1 

def add_gnomad_4_1 (data, filtered_gnomad_4_1): # patient data and the filtered gnomad dataset
    selected_gnomad_cols = ['gene_id', 'lof.oe_ci.upper', 'lof.pLI', 'lof.z_score', 'mis.oe', 'mis.z_score', 'syn.oe', 'syn.z_score']

    # Rename the columns to match your target columns in data
    gnomad_file_renamed = filtered_gnomad_4_1[selected_gnomad_cols].rename(columns={
        'lof.oe_ci.upper': 'LOEUF_upper',
        'lof.pLI': 'pLI',
        'lof.z_score': 'lof_z_4_1',
        'mis.oe': 'mis_oe_4_1',
        'mis.z_score': 'mis_z_4_1',
        'syn.oe': 'syn_oe_4_1',
        'syn.z_score': 'syn_z_4_1'
    })

    data = pd.merge(data, gnomad_file_renamed, how='left', left_on='Gene', right_on='gene_id')
    data = data.drop(columns='gene_id')

    return data # patient data with added gene constraint

def adding_all_an(dataset, mgi_df, human_mouse_symbol, gene_dict, ids_dictionary, filtered_gnomad_4_1, dataset_name, file_path):
    dataset = add_revel(ids_dictionary = ids_dictionary, dataset= dataset)
    dataset = add_gmt (data = dataset, gene_dict = gene_dict)
    dataset = add_mgi (data = dataset, mgi_df = mgi_df, human_mouse_symbol = human_mouse_symbol)
    dataset = add_gnomad_4_1 (data = dataset, filtered_gnomad_4_1 = filtered_gnomad_4_1 )
    return dataset # patient dat with all annotation except alpha missense
    
def extract_data(dataset):
    # to extract all the reqruied attributes
    combined = dataset[["genomic_nomencalture", "Consequence", "IMPACT", 
                        "Amino_acids", "Codons", "SYMBOL", "AF", "gnomAD3_AF", 
                        "gnomAD_AF", "CADD_PHRED", "SpliceAI_pred_DS_AG", 
                        "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                        "SpliceAI_pred_DS_DL", "revel", "term", 
                        'mouse_symbol', 'mgi_notes', 'LOEUF_upper', 
                        'lof_z_4_1', 'mis_oe_4_1', 'mis_z_4_1', 'syn_oe_4_1', 
                        'syn_z_4_1', 'PolyPhen', 'ClinVar_CLNSIG']]
    
    return combined