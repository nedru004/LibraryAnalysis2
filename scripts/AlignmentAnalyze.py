#!/usr/bin/env python
import copy
import csv
import os
import numpy
import pandas as pd
import re
import itertools

codon = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I',
         'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
         'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
         'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
         'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D',
         'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
         'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
         'GGG': 'G'}


def merge_pairs(in1, in2, out, bbmerge_location=r'java -cp ../bbmap/current/ jgi.BBMerge'):
    #bbmerge_location = r'java -cp ' + os.path.abspath(os.path.join(__file__, "..",'bbmap/current')) + ' jgi.BBMerge'
    bbmerge_command = bbmerge_location
    merge_out = os.path.splitext(in1)[0]+'_merge'+os.path.splitext(in1)[1]
    unmerge_out1 = os.path.splitext(in1)[0]+'_unmerge'+os.path.splitext(in1)[1]
    unmerge_out2 = os.path.splitext(in2)[0]+'_unmerge'+os.path.splitext(in2)[1]
    for file in [out, merge_out, unmerge_out1, unmerge_out2]:
        if os.path.exists(file):
            os.remove(file)
    bbmerge_command += f' in1={in1} in2={in2} out={merge_out} outu1={unmerge_out1} outu2={unmerge_out2}'
    print(f'{bbmerge_command}')
    ret = os.system(bbmerge_command)
    assert ret == 0
    filenames = [merge_out, unmerge_out1, unmerge_out2]
    with open(out, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    for tempfile in [merge_out, unmerge_out1, unmerge_out2]:
        os.remove(tempfile)


def correct_pairs(in1, in2, bbmerge_location=r'java -cp ../bbmap/current/ jgi.BBMerge'):
    bbmerge_command = bbmerge_location
    if 'fastq.gz' in in1:
        out1 = in1.split('.fastq.gz')[0]+'_corrected.fastq.gz'
        out2 = in2.split('.fastq.gz')[0] + '_corrected.fastq.gz'
    else:
        out1 = os.path.splitext(in1)[0]+'_corrected'+os.path.splitext(in1)[1]
        out2 = os.path.splitext(in2)[0]+'_corrected'+os.path.splitext(in2)[1]
    bbmerge_command += f' in1={in1} in2={in2} out1={out1} out2={out2} ecco mix'
    print(bbmerge_command)
    ret = os.system(bbmerge_command)
    assert ret == 0
    return out1, out2


def david_call_variants(sam_file, wt, outfile, app, root):
    quality = 25
    quality_nt = 30
    # Analyze Codons
    loop_read = open(sam_file, 'r')
    mutation_dict = {}
    variants = {'WT': 0}
    wt_count = [0]*int(len(wt)/3)
    file = open(outfile + '.csv', 'w')
    wt_file = open(outfile +'_wt.csv', 'w')
    variant_check = app.variant_check.get()
    read_count = 0
    percentage_reads = list(range(0, app.reads, int(app.reads/99)+(app.reads % 99 > 0)))
    while True:  # loop through each read from ngs
        try:
            tmp_r1 = next(loop_read).split('\t')
        except:
            break
        if '@' in tmp_r1[0]:
            continue
        # if read is paired to the next file
        if tmp_r1[6] == '=':
            tmp_r2 = next(loop_read).split('\t')
            if tmp_r1[0] == tmp_r2[0]:  # paired reads have the same name
                # identify forward strand
                if int(f'{int(tmp_r1[1]):012b}'[::-1][4]):
                    r1 = tmp_r2
                    r2 = tmp_r1
                else:
                    r1 = tmp_r1
                    r2 = tmp_r2
                read_list = [r1, r2]
            else:
                raise Exception('Paired read not found')
        else:  # only if there is one read
            r1 = tmp_r1
            read_list = [r1]
        tmpcounts = {r1[0]}  # initialize name
        indel = 0  # move this inside the loop if you want to ignore indels from one of the reads
        tmp_wt_count = [0]*int(len(wt)/3+1)
        read_length = [100000, 0]
        for r3 in read_list:
            if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and int(r3[4]) > quality:  # if forward read is aligned and higher than quality threshold
                # set index
                ref = int(r3[3])-1  # subtract 1 for python numbering
                read_length[0] = min(read_length[0], ref)
                read = 0
                # cycle through each codon until end of read or end of reference sequence
                fragment_cigar = list(filter(None, re.split(r'(\d+)', r3[5])))
                for cigar in [fragment_cigar[ii:ii+2] for ii in range(0, len(fragment_cigar), 2)]:
                    if 'X' == cigar[1]:
                        # this is a substitution
                        for i_sub in range(int(cigar[0])):  # not sure how this range works
                            if ord(r3[10][read+i_sub].upper()) - 33 > quality_nt:
                                codon_diff = (ref+i_sub) % 3  # find the correct frame
                                # record codon if you can see the entire codon
                                i_codon = r3[9][read+i_sub-codon_diff:read+i_sub-codon_diff+3].upper()
                                if read+i_sub-codon_diff+3 < len(r3[9]) and read+i_sub-codon_diff >= 0 and 'N' not in i_codon:
                                    tmpcounts.add(str(int((ref+i_sub-codon_diff)/3+1))+'_' + i_codon)
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    elif '=' == cigar[1]:
                        # this is a perfect match
                        codon_start = ref + ref % 3
                        codon_end = ref + int(cigar[0]) - (int(cigar[0]) % 3)
                        if codon_end - codon_start >= 3:  # if the perfect match spans an entire codon then count it
                            # make sure this doesnt exceed the length of the gene using min
                            wt_range = range(int(codon_start/3+1), min(int((codon_end-3)/3+1)+1, int(len(wt)/3)), 1)  # need to add 1 for range
                            for wt_i in wt_range:
                                tmp_wt_count[wt_i] = 1
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    elif 'D' == cigar[1]:
                        # this is a deletion
                        indel = 1
                        # only record if the deletion is a multiple of 3 (no frameshift)
                        if int(cigar[0]) % 3 == 0:
                            tmpcounts.add(str(int(ref/3+1))+'_'+str(int(ref/3+1+int(cigar[0])/3))+'del')
                        ref += int(cigar[0])
                    elif 'I' == cigar[1]:
                        indel = 1
                        # this is an insertion
                        # only record if the insertion is a multiple of 3 (no frameshift)
                        if int(cigar[0]) % 3 == 0:
                            insertion = r3[9][read:read+int(cigar[0])].upper()
                            tmpcounts.add(str(int(ref/3+1)) + '_' + str(int(ref/3+1+int(cigar[0])/3)) + 'ins'+insertion)
                        read += int(cigar[0])
                    elif 'S' == cigar[1]:
                        # soft clipping
                        read += int(cigar[0])
                    elif 'M' == cigar[1]:
                        # bad read (N) - ignore
                        ref += int(cigar[0])
                        read += int(cigar[0])
                    else:
                        raise TypeError('Cigar format not found')
                read_length[1] = max(read_length[1], ref)
        if len(tmpcounts) > 1: # and not indel:  # only record if mutations found # only recording if no indels were found
            # append data
            sorted_counts = sorted(tmpcounts, reverse=True)
            #totalcounts.append(sorted_counts)  # this is really slow. Instead write to file
            file.write(','.join(sorted_counts) + '\n')
            mutation = sorted_counts[1:]
            for mut in mutation:
                try:
                    mutation_dict[mut] += 1  # simple dictionary of mutations
                except KeyError:
                    mutation_dict[mut] = 1
            if variant_check and read_length[0] == 0 and read_length[1] >= len(wt):
                try:
                    variants[','.join(mutation)] += 1
                except KeyError:
                    variants[','.join(mutation)] = 1
        elif len(tmpcounts) == 1 and variant_check and read_length[0] == 0 and read_length[1] >= len(wt):
            variants['WT'] += 1
        if sum(tmp_wt_count) > 1:  # and not indel:
            wt_positions = [xi for xi, x in enumerate(tmp_wt_count) if x]  # all positions wt codon was found
            for wt_idx in wt_positions:
                wt_count[wt_idx] += 1  # only add to total count after both reads have been processed
            wt_file.write(','.join(str(x) for x in wt_positions) + '\n')
        read_count += 1
        if read_count in percentage_reads:
            root.update_idletasks()
            app.progress['value'] += 1
            root.update()
    # Add wt to mutation_dict
    for i, wt_c in enumerate(wt_count[1:]):
        mutation_dict[str(i+1)+'_'+str(wt.seq[i*3:i*3+3])] = wt_c
    loop_read.close()
    file.close()
    wt_file.close()
    with open(outfile + '_mutlist.csv', 'w') as f:
        for mut in mutation_dict.keys():
            f.write(','.join(mut.split('_'))+',')
            f.write(str(mutation_dict[mut])+'\n')
    mutations = pd.DataFrame({'name': list(mutation_dict.keys()), 'count': list(mutation_dict.values())})
    mutations[['position', 'AA']] = mutations['name'].str.split('_', expand=True)
    mutations = mutations.replace({'AA': codon})
    mutations['position'] = pd.to_numeric(mutations['position'])
    matrix = mutations.pivot_table(index='AA', columns='position', values='count', aggfunc=sum)
    matrix.to_csv(outfile + '_matrix.csv')
    with open(outfile + '_variants.csv', 'w') as f:
        for mut in variants.keys():
            f.write(str(variants[mut]) + ',')
            f.write(mut + '\n')

    #os.remove(sam_file)
    #return mutation_table, correlation_matrix, correlation_matrix_wt


def david_paired_analysis(mut_files, wt_files, app, root):
    percentage_reads = list(range(0, app.reads, int(app.reads/99)))
    root.update_idletasks()
    app.progress['value'] = 0
    root.update()
    outfile = os.path.splitext(mut_files)[0]
    mutation_file = open(mut_files, 'r')
    ## Mutation Pairs ##
    mutation_pair_dict = {}
    for mutation_list in mutation_file:
        mutations = mutation_list.strip('\n').split(',')[1:]
        if len(mutations) > 1:
            for mut_1, mut_2 in itertools.combinations(mutations, 2):
                try:
                    mutation_pair_dict[mut_1+'+'+mut_2] += 1  # simple dictionary of mutations
                except KeyError:
                    mutation_pair_dict[mut_1+'+'+mut_2] = 1
    with open(outfile+'_mut_pairs.csv', 'w') as f:
        for mut in mutation_pair_dict.keys():
            f.write(mut+',')
            f.write(str(mutation_pair_dict[mut])+'\n')
    ## WT pairs ##
    wt_positions = open(wt_files, 'r')
    correlations_wt = {}
    read_count = 0
    for positions in wt_positions:
        for wt_pos in itertools.combinations(positions.strip('\n').split(',')[1:], 2):
            try:
                correlations_wt[','.join(wt_pos)] += 1
            except KeyError:
                correlations_wt[','.join(wt_pos)] = 1
        read_count += 1
        if read_count in percentage_reads:
            root.update_idletasks()
            app.progress['value'] += 1
            root.update()
    wt_matrix = pd.DataFrame({'pair': list(correlations_wt.keys()), 'count': list(correlations_wt.values())})
    wt_matrix[['pos1', 'pos2']] = wt_matrix['pair'].str.split(',', expand=True)
    wt_matrix['pos1'] = pd.to_numeric(wt_matrix['pos1'])
    wt_matrix['pos2'] = pd.to_numeric(wt_matrix['pos2'])
    matrix = wt_matrix.pivot_table(index='pos1', columns='pos2', values='count', aggfunc=sum)
    # matrix = pd.DataFrame(correlation_matrix_wt)
    matrix.to_csv(outfile + '_correlation_wt_table.csv')

def align_all_bbmap(sequencing_file, reference, sam_file, bbmap_script=r'java -cp ../bbmap/current/ align2.BBMap',
                    max_gap=500, paired_sequencing_file=None):
    if paired_sequencing_file:
        command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                  f"in2={paired_sequencing_file} maxindel={max_gap} local outm={sam_file}"
    else:
        command = f"{bbmap_script} ref={reference} in={sequencing_file} " \
                  f"maxindel={max_gap} local outm={sam_file}"
    print(command)
    ret = os.system(command)
    assert ret == 0