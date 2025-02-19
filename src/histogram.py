# To show the tuning result in the patient dataset with an identified pathogenic variant
# You are able to plot multiple plots, but each plot can only show one targeted genetic variant position
# Plots can be from different dataset
''' 
datasets should be in a format of tuples in a list with [(dataset variable name, the genogenomic_nomencalture of targed variant, the patient id)]
datasets = [
    (test_data_result, 'chr4:54295165A>T', '1'),  
    (test_data_result, 'chr19:1231143T>C', '2'), 
    (test_data_result, 'chr17:50375415A>C', '3'),
    (test_data_result, 'chr1:89186388T>C', '4')
    ]

FYI test data_result can be diferent dataset.

'''

import matplotlib.pyplot as plt
import math

def histogram (datasets):

    # Create a grid of subplots
    num_plots = len(datasets)
    
    if num_plots > 2:
        # Create a grid of subplots
        cols = math.ceil(num_plots/2)  
        rows = 2  # only 2 rows is always allowed
        fig, axs = plt.subplots(rows, cols, figsize=(9 * cols, 6 * rows))

        for i, (dataset, variant, cohort) in enumerate(datasets):
            row = i // cols
            col = i % cols
            
            # Plot the Distribution of 'total' Scores
            axs[row, col].hist(dataset['total'], bins=30, color='blue', alpha=0.7)
            axs[row, col].set_xlabel('Total Score', fontsize = 13)
            axs[row, col].set_ylabel('Variant Count', fontsize = 13)
            axs[row, col].set_title(f'Distribution for {variant} in Patient {cohort}', fontsize = 15)

            
            # Get the score for the identified pathogenic variant
            specific_variant_score = dataset.loc[dataset['genomic_nomencalture'] == variant, 'total'].values[0]
            
            # Plot a vertical line at the score of the specific variant
            axs[row, col].axvline(specific_variant_score, color='red', linestyle='dashed', linewidth=2, label=f'{variant} (Score = {specific_variant_score})')
            axs[row, col].legend()
            
            # Count Scores Above and Below the Specific Variant's Score
            above_count = (dataset['total'] > specific_variant_score).sum()
            precentage = round((dataset['total'] >= specific_variant_score).sum() / len(dataset) * 100, 2)
            below_count = (dataset['total'] < specific_variant_score).sum()
            top_num = (dataset['total'] >= specific_variant_score).sum()
            
            # Display the counts on the plot
            axs[row, col].text(specific_variant_score + 0.05, axs[row, col].get_ylim()[1]*0.75, f'Top: {top_num} ({precentage}%) \nAbove: {above_count}\nBelow: {below_count}', color='red', ha='left', fontsize = 15)

    else: # when only 1 or 2 variants are required
        if num_plots == 1:
            rows = 1
            cols = 1
        else:
            rows = 1
            cols = 2
        fig, axs = plt.subplots(rows, cols, figsize=(10, 5 * rows), squeeze=False)

        # Flatten the axs array for easier iteration
        axs = axs.flatten()
        # Plot each dataset in its respective subplot
        for i, (dataset, variant, patient) in enumerate(datasets):
            axs[i].hist(dataset['total'], bins=30, color='blue', alpha=0.7)
            axs[i].set_title(f'Patient {patient}: {variant}', fontsize=14)
            axs[i].set_xlabel('Total Score', fontsize=13)
            axs[i].set_ylabel('Variant Count', fontsize=13)

        # Get the score for the identified pathogenic variant
            specific_variant_score = dataset.loc[dataset['genomic_nomencalture'] == variant, 'total'].values[0]
            
            # Plot a vertical line at the score of the specific variant
            axs[i].axvline(specific_variant_score, color='red', linestyle='dashed', linewidth=2, label=f'{variant} (Score = {specific_variant_score})')
            axs[i].legend()
            
            # Count Scores Above and Below the Specific Variant's Score
            above_count = (dataset['total'] > specific_variant_score).sum()
            precentage = round((dataset['total'] >= specific_variant_score).sum() / len(dataset) * 100, 2)
            below_count = (dataset['total'] < specific_variant_score).sum()
            top_num = (dataset['total'] >= specific_variant_score).sum()

            axs[i].text(specific_variant_score + 0.05, axs[i].get_ylim()[1]*0.75, f'Top: {top_num} ({precentage}%) \nAbove: {above_count}\nBelow: {below_count}', color='red', ha='left', fontsize = 15)
    
    plt.tight_layout()

    return plt.show()