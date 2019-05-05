import pandas as pd
import numpy as np
import sys


# Sample argv[1]: filenames_birna2.txt
df_files = pd.read_csv(sys.argv[1], delimiter='\t', index_col=0, names=['Q1', 'Q2'])
print(df_files.head())

# Sample argv[2]: raw_human_data_aug_info_interactions.csv
df_interact = pd.read_csv(sys.argv[2], delimiter=',', index_col=0, names=['len1', 'len2', 'start1', 'end1', 'start2', 'end2'])
print(df_interact.head())


df = pd.concat([df_files, df_interact], axis=1, join='inner')
print(df.head())
print(len(df))

df['outof'] = '.'
df['rank'] = '.'
df['rank_count'] = '.'

df['half'] = '.'
df['i_half'] = '.'

df['half_count'] = '.'
df['i_half_count'] = '.'

for idx, row in df.iterrows():
    name = str(idx)

    curr_name = name + "_bi_pair_0.5_1"
    curr_df = pd.read_csv(curr_name, delimiter='\t')
    sorted_curr_df = curr_df.sort_values(by=['unpaired_score'], ascending=False)
    sorted_curr_df = sorted_curr_df.reset_index(drop=True)
    print(curr_name)
    print(sorted_curr_df.head())

    flag_rank = False
    flag_half = False
    flag_i_half = False

    rank_count = 0
    half_count = 0
    i_half_count = 0

    for idx_2, row_2 in sorted_curr_df.iterrows():

        _from_1 = row_2['start1']
        _till_1 = min(row_2['start1'] + row_2['win_size1'] - 1, row['len1'] - 1)
        
        _from_2 = row_2['start2']
        _till_2 = min(row_2['start2'] + row_2['win_size2'] - 1, row['len2'] - 1)


        i_from_1 = row['start1']
        i_till_1 = min(row['end1'], row['len1']-1)
        
        i_from_2 = row['start2']
        i_till_2 = min(row['end2'], row['len2']-1)


        common_section_1 = min(_till_1, i_till_1) - max(_from_1, i_from_1) + 1
        common_section_2 = min(_till_2, i_till_2) - max(_from_2, i_from_2) + 1

        if common_section_1 <= 0 or common_section_2 <= 0:
            continue
        

        rank_count += 1
        if not flag_rank:
            df.ix[idx, 'rank'] = idx_2 + 1
            df.ix[idx, 'outof'] = len(sorted_curr_df)
            flag_rank = True


        if common_section_1 >= (_till_1 - _from_1 + 1)/2 and common_section_2 >= (_till_2 - _from_2 + 1)/2:
            half_count += 1
            if not flag_half:
                df.ix[idx, 'half'] = idx_2 + 1
                flag_half = True


        if common_section_1 >= (i_till_1 - i_from_1 + 1)/2 and common_section_2 >= (i_till_2 - i_from_2 + 1)/2:
            i_half_count += 1
            if not flag_i_half:
                df.ix[idx, 'i_half'] = idx_2 + 1
                flag_i_half = True

    df.ix[idx, 'rank_count'] = rank_count
    df.ix[idx, 'half_count'] = half_count
    df.ix[idx, 'i_half_count'] = i_half_count


outfile = str(sys.argv[1]).split('.')[0] + '_pair_ranks.csv'
df.to_csv(outfile)


