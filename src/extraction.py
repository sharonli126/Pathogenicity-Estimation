# This is for when you want to look at record of specific genetic variants

# IMPORTANT run with ductdb

import duckdb as db

def pathogenic_variant_table(con, table_and_variant_dict, table_name):
    # Create the table with an extra "from" column of type VARCHAR to store dataset names
    sample_table = next(iter(table_and_variant_dict))
    con.execute(f"""
        CREATE OR REPLACE TABLE {table_name} AS 
        SELECT *, CAST(NULL AS VARCHAR) AS "from" FROM {sample_table} WHERE 1=0;
    """)
    
    # Iterate over each dataset and its variants
    for table in table_and_variant_dict:
        for variant in table_and_variant_dict[table]:
            # Insert matching variants into the result table with the dataset name
            con.execute(f"""
                INSERT INTO {table_name}
                SELECT *, '{table}' AS "from" 
                FROM {table}
                WHERE genomic_nomencalture ILIKE '%{variant}%'
            """)
    
    # Fetch the final results as a DataFrame
    return con.execute(f"SELECT * FROM {table_name}").fetchdf() # obtain a table with all the records (annotations and score) of tagted genetic variants