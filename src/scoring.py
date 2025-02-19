# This file contains the algorithms and the calculation of the total score. The final dataset will be obtained the final_score function.

# IMPORTANT please install duckdb and run all the function with duckdb

def add_col(con, tablename):
    # Add all the scoring columns for later scoring
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN am_pathogenicity REAL;")
    con.execute(f"""
        CREATE OR REPLACE TABLE {tablename}  AS
        SELECT * 
        FROM {tablename}  c
        LEFT JOIN data_freq d
        ON c.genomic_nomencalture = d.uniqueValue;
    """)
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN consequence_score INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN splice_score INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN maf_score INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN caddscore INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN am_score INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN mgiscore INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN gnomadscore INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN revelscore INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN genoterm_presence INTEGER;")
    con.execute(f"ALTER TABLE {tablename} ADD COLUMN total REAL;")

def scoring(con, tablename):
# minor allele frequency of Ensembl
    con.execute(f"""
        UPDATE {tablename}
        SET maf_score = 
            CASE 
                WHEN AF IS NULL THEN 3 
                WHEN AF < 0.0001 THEN 3
                WHEN AF < 0.005 THEN 2
                WHEN AF < 0.01 THEN 1
                ELSE 0
            END;
    """) 


# gnomAD maf
    data_limit = con.execute(f"SELECT gnomAD3_AF FROM {tablename} LIMIT 10;").fetchdf()
    data_type = data_limit["gnomAD3_AF"].dtype

    if data_type == 'O':  # 'O' indicates object type, usually strings
        # If 'gnomAD3_AF' contains string values, run this query
        con.execute(f"""
            UPDATE {tablename}
            SET maf_score = maf_score +
                CASE 
                    WHEN gnomAD3_AF = '.' 
                         AND (gnomAD_AF < 0.0001 OR gnomAD_AF IS NULL) 
                    THEN 3 -- when gnomAD3_AF is invalid, use gnomAD_AF, if it is missing or < 0.0001, score 3
                    WHEN gnomAD3_AF IS NULL THEN 3 -- if gnomAD3_AF is missing score 3
                    WHEN gnomAD3_AF != '.' AND CAST(gnomAD3_AF AS DOUBLE) < 0.0001 THEN 3 -- if gnomAD3_AF <0.0001 score 3
                    WHEN gnomAD3_AF != '.' AND CAST(gnomAD3_AF AS DOUBLE) < 0.005 THEN 2 -- if gnomAD3_AF <0.005 score 2
                    WHEN gnomAD3_AF != '.' AND CAST(gnomAD3_AF AS DOUBLE) < 0.01 THEN 1 -- if gnomAD3_AF < 0.01 score 2
                    ELSE 0
                END;
        """)
    else:  # If 'gnomAD3_AF' is numeric
        con.execute(f"""
            UPDATE {tablename}
            SET maf_score = maf_score +
                CASE 
                    WHEN gnomAD3_AF IS NULL THEN 3 -- if gnomAD3_AF is missing score 3
                    WHEN gnomAD3_AF < 0.0001 THEN 3 -- if gnomAD3_AF <0.0001 score 3
                    WHEN gnomAD3_AF < 0.005 THEN 2 -- if gnomAD3_AF <0.005 score 2
                    WHEN gnomAD3_AF  < 0.01 THEN 1 -- if gnomAD3_AF < 0.01 score 2
                    ELSE 0
                END;
        """)

# AF False positive 
    con.execute(f"""
            UPDATE {tablename}
            SET maf_score = maf_score +
                CASE
                    WHEN distinct_patientDatabaseId_count > 139 THEN -3 -- if allele frequency calculated from APF database > 139 person score -3
                    ELSE 0
                END;
        """)   


# caddscore
    con.execute(f"""
            UPDATE {tablename}
            SET caddscore =
                CASE
                    WHEN CADD_PHRED > 20 THEN 3 -- if CADD_PHRED >20 score 3
                    WHEN CADD_PHRED > 10 THEN 2 -- if CADD_PHRED >10 score 2
                    WHEN CADD_PHRED > 5 THEN 1 -- if CADD_PHRED >5 score 1
                    ELSE 0 -- else no score
                END;
        """)   

# mgi
    con.execute(f"""
    UPDATE {tablename}
    SET mgiscore =
        CASE
            WHEN mgi_notes IS NULL THEN 0
                WHEN mgi_notes ~~* '%immune%' 
                OR mgi_notes ~~* '%apopto%'
                OR mgi_notes ~~* '% lethal%' 
                OR mgi_notes ~~* '%fatal%' 
                OR mgi_notes ~~* '%severe %' 
                OR mgi_notes ~~* '%die%' 
                OR mgi_notes ~~* '%death%' 
                OR mgi_notes ~~* '%cytotoxic%' 
                OR mgi_notes ~~* '%null mutation%'
            THEN 1 -- for mgi_note with any of these keywords, score 1
            ELSE -1 -- else score -1
        END;
    """)


    
# revel
    con.execute(f"""
        UPDATE {tablename}
        SET revelscore =
            CASE
                WHEN revel > 0.75 THEN 2 -- if revel > 0.75 score 2
                WHEN revel > 0.5 THEN 1 -- if revel > 0.5 score 1
                ELSE 0 -- else no score
            END;
    """)

# alpha missense
    con.execute(f"""
        UPDATE {tablename}
        SET am_score =
            CASE
                WHEN am_pathogenicity > 0.75 THEN 2 -- if am_pathogenicity > 0.75 score 2
                WHEN am_pathogenicity > 0.5 THEN 1 -- if am_pathogenicity > 0.5 score 1
                ELSE 0 -- else no score
            END;
    """)

# geno_term
    con.execute(f"""
    UPDATE {tablename}
    SET genoterm_presence = 0;
    """)

    con.execute(f"""
    UPDATE {tablename}
    SET genoterm_presence = 
        CASE 
            WHEN term ~~* '%INTERFERON%'::text 
                OR term ~~* '%CYTOKINE%'::text 
                OR term ~~* '%CHEMOKINE%'::text 
                OR term ~~* '%TOLL%'::text 
                OR term ~~* '%LYMPHO%'::text 
                OR term ~~* '%MACROPHAGE%'::text 
                OR term ~~* '%LEUCOCY%'::text 
                OR term ~~* '%IMMUN%'::text 
                OR term ~~* '%NF-KAP%'::text 
                OR term ~~* '%B CELL RECEPTOR%'::text 
                OR term ~~* '%APOPTO%'::text 
            THEN 1 -- if any of the keyword present show 1
            ELSE 0 -- else show 0
        END;
""")

# Consequence

    con.execute(f"""

    UPDATE {tablename}
    SET consequence_score = 
        CASE 
            WHEN IMPACT ~~* '%HIGH%' THEN 2 -- if impact is high score 2
            WHEN IMPACT ~~* '%MODERATE%' THEN 1 -- if impact is moderate score 1
            WHEN IMPACT ~~* '%LOW%' THEN 0 -- if impact is low score 0
            ELSE 0
        END;
    """)
    
# splice
    con.execute(f"""
    UPDATE {tablename}
    SET splice_score = 
        CASE 
            WHEN (SpliceAI_pred_DS_AG > 0.8 OR SpliceAI_pred_DS_AL > 0.8 OR SpliceAI_pred_DS_DG >0.8 OR SpliceAI_pred_DS_DL > 0.8) THEN 2  -- if any SpliceAI_pred > 0.8 score 2
            WHEN (SpliceAI_pred_DS_AG > 0.6 OR SpliceAI_pred_DS_AL > 0.6 OR SpliceAI_pred_DS_DG >0.6 OR SpliceAI_pred_DS_DL > 0.6) THEN 1-- if any SpliceAI_pred > 0.6 score 1
            ELSE 0 -- else score 0
        END;

    """)

# for gnomad_4_1
    con.execute(f"""
    UPDATE {tablename}
    SET gnomadscore = 
    CASE
        WHEN Consequence ~~* '%missense_variant%' THEN CASE
            WHEN mis_z_4_1 > 3 THEN 2 -- if mis_z_4_1 > 3 score 2
            WHEN mis_z_4_1 > 1.7 THEN 1 -- if mis_z_4_1 > 2 score 1
            WHEN mis_z_4_1 < 0 THEN -1 -- if mis_z_4_1 <0 score -1
            ELSE 0 -- else score 0
            END
        WHEN Consequence ~~* '%synonymous_variant%' THEN CASE
            WHEN syn_z_4_1 > 3 THEN 2 -- if syn_z_4_1 > 3 score 2
            WHEN syn_z_4_1 > 1.7 THEN 1 -- if syn_z_4_1 > 2 score 1
            WHEN syn_z_4_1 < 0 THEN -1 -- if syn_z_4_1 <0 score -1
            ELSE 0 -- else score 0
            END
        ELSE CASE 
            WHEN LOEUF_upper < 0.5 THEN 2  -- if LOEUF_upper < 0.5 score 2
            WHEN LOEUF_upper < 0.67 THEN 1  -- if LOEUF_upper < 0.67 score 1
            WHEN LOEUF_upper > 1.2 THEN - 1  -- if LOEUF_upper > 1.2 score -1
            ELSE 0 -- else score 0
            END
        END;
    
                """)



def final_score(con, tablename):
    # to scale the total score of all different type of variant with different annotation to the same range
    con.execute(f"""
        UPDATE {tablename}
        SET total = CASE
            WHEN Consequence ~~* '%missense_variant%' THEN -- for missense variants, the total score is 18, with specific annotation
                (maf_score + caddscore + mgiscore + revelscore + am_score + consequence_score + gnomadscore) / 18
            WHEN Consequence ~~* '%splice%' THEN -- for splice variant, the total score is 16, with specific annotation
                (maf_score + caddscore + mgiscore + splice_score + consequence_score + gnomadscore) / 16
            ELSE -- for other variants, the total score is 14, with all general annoations
                (maf_score + caddscore + mgiscore  + consequence_score + gnomadscore) / 14
            END;
    """)

    data = con.execute(f"SELECT * FROM {tablename} ORDER BY total DESC").fetchdf()

    return data

    

