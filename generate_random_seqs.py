import sys
import random
import os
import pandas


alph_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
alph_list = [0, 1, 2, 3]

max_num = int(sys.argv[2])
len_list = range(5, max_num, 1)


def search_file(file_name, l1, l2):
    with open(file_name) as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split('=')
            keyword = 'FTable(0,' + str(l1-1) + ',0,' + str(l2-1) + ')'
            print("keyword: ", keyword)
            if parts[0] == keyword:
                return float(parts[1][:-1])
    return 'Does not exist!'


alpha_values = []

with open('random_seqs.fa', 'w') as f:
    for dummy_i in range(int(sys.argv[1])):
        f.write('>test_'+str(dummy_i)+'\n')
            
        l1 = random.choice(len_list)
        l2 = random.choice(len_list)
           
        alpha_file = 'random_seq_number_' + str(dummy_i) + '.fa'
        with open(alpha_file, 'w') as fn:
            for i in range(l1):
                r = random.choice(alph_list)
                f.write(str(alph_dict[r]))
                fn.write(str(r)+'\n')
            
            f.write('&')

            for i in range(l2):
                r = random.choice(alph_list)
                f.write(str(alph_dict[r]))
                fn.write(str(r)+'\n')
            
            fn.write('\n')

        os.system('./bpmax.check ' + str(l1) + ' ' + str(l2) + ' < ' + alpha_file + ' > ' + 'temp_out.txt')
        
        final_alpha_value = search_file('temp_out.txt', l1, l2)

        print(final_alpha_value)
        alpha_values.append(final_alpha_value)
 
        f.write('\n')


# os.system('rm temp_out.txt')
os.system('./bpmax ' + 'random_seqs.fa')

df = pandas.read_csv('random_seqs.fa.scores.1_1', delimiter='\t')

print(df.head())

orig_values = list(df['full'])

for i in range(len(orig_values)):
    if orig_values[i] != alpha_values[i]:
        print('Test ' + str(i) + ': Failure!')
    else:
        print('Test ' + str(i) + ': Passed!')













