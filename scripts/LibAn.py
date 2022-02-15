# python
import AlignmentAnalyze
import os
import time
import Bio
import Bio.Entrez
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

def run(wt_seq, seq_file1, seq_file2):
    # Sanity checks
    #app.run.config(bg='green', activebackground='green', relief=tk.SUNKEN, state='disabled')
    #assert os.path.isfile(f'{args.bowtie}bowtie2'), f"Can't find bowtie2 at {args.bowtie}"
    #assert os.path.isfile(f'{args.bowtie}bowtie2-build'), f"Can't find bowtie2-build at {args.bowtie}"
    assert os.path.isfile(seq_file1), f"given sequencing file, '{seq_file1}', does not exist"
    assert Bio.SeqIO.read(wt_seq, 'fasta').seq.translate(), f"given refrence/wildtype sequence file " \
                                                                f"'{wt_seq}' is not a valid FASTA file " \
                                                                f"containing one unambiguous DNA sequence!"

    # variables from arguments
    wt_seq = Bio.SeqIO.read(wt_seq, 'fasta')
    sequencing_file = seq_file1
    paired_sequencing_file = seq_file2
    quiet = 0

    # output files
    rootname = os.path.splitext(seq_file1)[0]
    log_file = f'{rootname}_logfile.log'

    # initialize some variables
    programstart = time.time()

    ### Process Sequencing Records ###

    # merge paired reads
    if paired_sequencing_file and not os.path.exists(rootname.split('.fastq')[0]+'_corrected.fastq.gz'):
        message = f'Merging paired reads from {sequencing_file} and {paired_sequencing_file} using bbmerge.'
        if not quiet:
            print(message)
        print('Merging paired reads.\n')
        sequencing_file, paired_sequencing_file = AlignmentAnalyze.correct_pairs(sequencing_file, paired_sequencing_file)
        # mergedfile = f"{rootname}_merged.{seqreadtype}"
        # AlignmentAnalyze.merge_pairs(sequencing_file, paired_sequencing_file, mergedfile)
        # sequencing_file = mergedfile

    # align reads (using bbmap)
    if not os.path.exists(f'{rootname}.sam'):  # and app.aamuts_file:
        message = f'Aligning all sequences from {sequencing_file} to {wt_seq} using bbmap.'
        print(message)
        print('Aligning sequencing reads to reference.\n')
        AlignmentAnalyze.align_all_bbmap(sequencing_file, wt_seq, f'{rootname}.sam', max_gap=len(wt_seq), paired_sequencing_file=paired_sequencing_file)

    # determine the number of reads
    lines = sum(1 for i in open(f'{rootname}.sam', 'rb'))
    reads = int((lines-3) / 2)
    print('Total number of reads to analyze: '+str(reads))
    print(f'Total number of reads to analyze: {str(reads)}\n')

    # call variants / find mutations
    print(f'Calling mutations/variants\n')
    # David's work
    AlignmentAnalyze.david_call_variants(f"{rootname}.sam", wt_seq, rootname, app, root)
    print(f'Finding paired mutations\n')
    if run_correlation:
        AlignmentAnalyze.david_paired_analysis(rootname + '.csv', rootname + '_wt.csv', app, root)
    # os.system(f'java -cp ../bbmap/current/ var2.CallVariants2 in={rootname}.sam ref={app.wt_file} ploidy=1 out={rootname}.vcf 32bit')

    seq_analyze_time = time.time() - programstart
    time_per_seq = round(seq_analyze_time / reads * 1000, 3)
    print(f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    print(f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    print('Finished Analysis')
