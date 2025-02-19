import pandas as pd

def validate_table(table_and_variant_dict):
    df = []
    
    for table_name, (table, variants) in table_and_variant_dict.items():
        
        for variant in variants:
            # Filtering the table rows for the specific variant
            filtered_rows = table[table['genomic_nomencalture'] == variant]
            
            if not filtered_rows.empty:
                conseq = filtered_rows['Consequence'].iloc[0]
                specific_variant_score = filtered_rows['total'].iloc[0]

                # Calculate metrics
                top_num = (table['total'] >= specific_variant_score).sum()
                total = len(table)
                percentage = round(top_num / total * 100, 2)
                above_count = (table['total'] > specific_variant_score).sum()
                below_count = (table['total'] < specific_variant_score).sum()

                # Append the result for this variant to the result list
                df.append({
                    'Genetic Variant': variant,
                    'Consequence': conseq,
                    'from': table_name,
                    'Score': specific_variant_score,
                    'Top': top_num,
                    'Percentage': percentage,
                    '# of record Above': above_count,
                    '# of record below': below_count,
                    'Total records': total
                })
            else:
                # Handling case where variant is not found
                df.append({
                    'Genetic Variant': variant,
                    'Consequence': 'Variant not found',
                    'from': table_name,
                    'Score': None,
                    'Top': None,
                    'Percentage': None,
                    '# of record Above': None,
                    '# of record below': None,
                    'Total records': len(table)
                })
    
    # Convert the result list into a DataFrame and return
    result_df = pd.DataFrame(df)
    return result_df
