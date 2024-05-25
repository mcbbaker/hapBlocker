#!/usr/bin/env python

# purpose: 
# usage: python3 hapBlocker.py vcf.gz index.gz.[csi,tbi] samples min_snps percent-threshold comparators
# notes: 'samples' should be the path to the file with the samples (named as they are in the vcf file) in the order the user 
#        wants the samples to appear in the plots. also assuming that 'Chr' is somewhere in the vcf ID's for the chromosomes
#        'comparator' refers to the first X rows in the sample files that you want to compare/color the samples with. for example,
#        the NAM populations would have the first 2 rows in the sample file as comparators. if using a population (e.g. MAGIC or depletion)
#        where you want samples to be compared with every other sample above it, use the number of samples minus 1


import sys
from pysam import VariantFile

index = sys.argv[2]
vcf = VariantFile(sys.argv[1],mode='r',index_filename=index)
sample_file = open(sys.argv[3],'r') 
min_snps = int(sys.argv[4])
threshold = int(sys.argv[5])
comparator = int(sys.argv[6])

chroms = []
lengths = {}
block_info = {}
sample_list = []
parents = []
results = {}

## FUNCTIONS

# takes in a list of tuples and determines how many tuples contain the same genotype (% identity) 
#   tuple_list: list of tuples. each tuple contains two genotypes, e.g. ((0/0),(1/1))
#   returns: the % identity
def get_percent_identical(tuple_list):
    identical = 0
    total = len(tuple_list)
    for i in tuple_list:
        if i[0] == i[1]:
            identical += 1
    return (identical/total)*100

# clears the data from the dictionary for a specific sample and resets it so this sample
# is placed back in the 'init' state
#   block_info: dictionary holding state/block information
#   i: sample to clear
#   start: start of the new block
#   parents: list of parents/comparators
def clear_dict_for_sample(block_info, i, start, parents):
    block_info[i]['State'] = 'Init'
    block_info[i]['Start'] = start
    block_info[i]['Block Num'] += 1
    block_info[i]['Hap Sample'] = ''
    for j in parents:
        if j in block_info[i].keys():
            block_info[i][j]['Identical'] = 0
            block_info[i][j]['Total'] = 0
            block_info[i][j]['Recent'].clear()

def merge_and_print(sample, results):
    prev = None
    for num in results[sample].keys():
        # first block for sample
        if prev == None:
            prev = results[sample][num]
        # remainder of blocks for sample
        else:
            # if same chromosome and coloring -- update end
            if results[sample][num][0] == prev[0] and results[sample][num][5] == prev[5]:
                prev[2] = results[sample][num][2]
            else:
                print(*prev, sep='\t')
                results[sample][num][3] = prev[3]+1
                prev = results[sample][num]
    print(*prev, sep='\t')



## MAIN

# PARSE THROUGH SAMPLE FILE
cur = 0
for line in sample_file:
    line = line.strip('\n')
    sample_list.append(line)
    results[line] = {}
    if cur < comparator:
        parents.append(line)
    cur += 1

# PREPROCESS HEADER FOR CHROM INFO
for record in vcf.header.records:
    # only parse through contig records
    if (record.type == 'CONTIG'):
        name = record.get('ID')
        # if chromosome, store name and length
        if ('Chr' in name):
            chroms.append(name)
            lengths[name] = int(record.get('length'))
        # ignore unitigs
        else:
            pass

# CREATE DICT TO STORE BLOCK/STATE INFO
# add first sample separately (no comparisons)
block_info[sample_list[0]] = {'State': 'Top', 'GT':'', 'Block Num': 1}
# loop through remainder of samples
cur = 1
for i in sample_list[1::]:
    block_info[i] = {'State': 'Init', 'Start': 1, 'GT':'', 'Block Num': 1, 'Hap Sample': ''}
    # loop through all samples above 
    for j in sample_list[0:cur:1]:
        # add identical/total info if it should be compared to/is a parent
        if j in parents:
            block_info[i][j] = {'Identical':0, 'Total':0, 'Recent': []}
    cur += 1

# PROCESS CALLS ONE CHROM AT A TIME
for chrom in chroms:
    # add entire chrom block for top sample to results
    results[sample_list[0]][block_info[sample_list[0]]['Block Num']] = list((chrom, 1, lengths[chrom], block_info[sample_list[0]]['Block Num'], sample_list[0], sample_list[0]))
    block_info[sample_list[0]]['Block Num'] += 1
    
    # loop through calls for each chrom
    interval = chrom + ":1-" + str(lengths[chrom])
    for call in vcf.fetch(region=interval):
        cur = 0
        for i in sample_list:
            # update genotype
            block_info[i]['GT'] = call.samples.get(i).get('GT')
            
            # Top state that only applies to the first sample -- already stored block for entire chrom, do nothing
            if block_info[i]['State'] == 'Top':
                pass

            # Init state for trying to initalize a block, keep adding snps until reached min snp count
            elif block_info[i]['State'] == 'Init':
                num_snp_pass = True
                # loop through all samples above
                for j in sample_list[0:cur:1]:
                    if j in parents:
                        # only compare if both genotypes are not missing
                        if block_info[i]['GT'][0] != None and block_info[j]['GT'][0] != None:
                            block_info[i][j]['Total'] += 1
                            block_info[i][j]['Recent'].append((block_info[i]['GT'],block_info[j]['GT']))
                            if block_info[i]['GT'] == block_info[j]['GT']:
                                block_info[i][j]['Identical'] += 1
                        if block_info[i][j]['Total'] < min_snps:
                            num_snp_pass = False
                # enough snps to see if there is a sample that passes threshold
                if num_snp_pass == True:
                    passing_sample = None
                    for j in sample_list[0:cur:1]:
                        if j in parents:
                            if passing_sample == None and (block_info[i][j]['Identical']/block_info[i][j]['Total']) * 100 >= threshold:
                                passing_sample = j
                    # there is a passing sample, so extend
                    if passing_sample != None:
                        block_info[i]['State'] = 'Extend'
                        block_info[i]['Hap Sample'] = passing_sample
                    # there is no passing sample, so fail and store block in results
                    else:
                        results[i][block_info[i]['Block Num']] = list((call.chrom, block_info[i]['Start'], call.pos, block_info[i]['Block Num'], i, i))
                        clear_dict_for_sample(block_info, i, call.pos+1, parents)


            # a block has been initialized - extend until % identical in most recent snps goes below threshold
            elif block_info[i]['State'] == 'Extend':
                hap_sample = block_info[i]['Hap Sample']
                if block_info[i]['GT'][0] != None and block_info[hap_sample]['GT'][0] != None:
                    # update total
                    block_info[i][hap_sample]['Total'] += 1
                    # check if identical
                    if block_info[i]['GT'] == block_info[hap_sample]['GT']:
                        block_info[i][hap_sample]['Identical'] += 1
                    # pop and append
                    block_info[i][hap_sample]['Recent'].pop(0)
                    block_info[i][hap_sample]['Recent'].append((block_info[i]['GT'],block_info[hap_sample]['GT']))
                    if get_percent_identical(block_info[i][hap_sample]['Recent']) < threshold:
                        results[i][block_info[i]['Block Num']] = list((call.chrom, block_info[i]['Start'], call.pos, block_info[i]['Block Num'], i, hap_sample))
                        clear_dict_for_sample(block_info, i, call.pos+1, parents)
            
            cur += 1

    # AFTER ALL CALLS FOR CHROM ARE PROCESSED
    # top sample already fully stored, store remainder of blocks to reach the end of the chrom and clear data for next chrom
    for i in sample_list[1::]:
        if block_info[i]['State'] == 'Extend':
            results[i][block_info[i]['Block Num']] = list((chrom, block_info[i]['Start'], lengths[chrom], block_info[i]['Block Num'], i, block_info[i]['Hap Sample']))
        else:
            # would only occur if last block ended on the last coordinate of the chromosome
            if block_info[i]['Start'] > lengths[chrom]:
                pass
            else:
                results[i][block_info[i]['Block Num']] = list((chrom, block_info[i]['Start'], lengths[chrom], block_info[i]['Block Num'], i, i))
        clear_dict_for_sample(block_info, i, 1, parents)

for sample in results.keys():
    merge_and_print(sample, results)

vcf.close()
sample_file.close()
