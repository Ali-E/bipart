import pandas as pd
import numpy as np
import sys



df_files = pd.read_csv(sys.argv[1], delimiter='\t', index_col=0, names=['Q1', 'Q2'])
# df_files = pd.read_csv(sys.argv[1], delimiter='\t', names=['id_name', 'Q1', 'Q2'])
print(df_files.head())

df_interact = pd.read_csv(sys.argv[2], delimiter=',', index_col=0, names=['len1', 'len2', 'start1', 'end1', 'start2', 'end2'])
# df_interact = pd.read_csv(sys.argv[2], delimiter=',', names=['id_name', 'len1', 'len2', 'start1', 'end1', 'start2', 'end2'])
print(df_interact.head())


df = pd.concat([df_files, df_interact], axis=1, join='inner')
print(df.head())
print(len(df))

df['outof1'] = '.'
df['outof2'] = '.'

df['rank1'] = '.'
df['rank2'] = '.'

df['rank_count1'] = '.'
df['rank_count2'] = '.'

df['half1'] = '.'
df['half2'] = '.'
df['i_half1'] = '.'
df['i_half2'] = '.'

df['half_count1'] = '.'
df['half_count2'] = '.'
df['i_half_count1'] = '.'
df['i_half_count2'] = '.'

for idx, row in df.iterrows():
    # name = row['id_name']
    name = str(idx)
    for i in range(2):
        curr_name = name + "_bi_single_" + str(i) + ".0.5_1"
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

            _from = row_2['start']
            _till = min(row_2['start'] + row_2['win_size'] - 1, row['len' + str(i+1)] - 1)
            
            i_from = row['start' + str(i+1)]
            i_till = min(row['end' + str(i+1)], row['len' + str(i+1)]-1)

            common_section = min(_till, i_till) - max(_from, i_from) + 1
            if common_section <= 0:
                continue
            

            rank_count += 1
            if not flag_rank:
                df.ix[idx, 'rank'+str(i+1)] = idx_2 + 1
                df.ix[idx, 'outof'+str(i+1)] = len(sorted_curr_df)
                flag_rank = True


            if common_section >= (_till-_from+1)/2:
                half_count += 1
                if not flag_half:
                    df.ix[idx, 'half'+str(i+1)] = idx_2 + 1
                    flag_half = True


            if common_section >= (i_till-i_from+1)/2:
                i_half_count += 1
                if not flag_i_half:
                    df.ix[idx, 'i_half'+str(i+1)] = idx_2 + 1
                    flag_i_half = True

        df.ix[idx, 'rank_count'+str(i+1)] = rank_count
        df.ix[idx, 'half_count'+str(i+1)] = half_count
        df.ix[idx, 'i_half_count'+str(i+1)] = i_half_count


outfile = str(sys.argv[1]).split('.')[0] + '_ranks.csv'
df.to_csv(outfile)


