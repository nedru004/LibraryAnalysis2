# python
import AlignmentAnalyze
import os
import time
import Bio
import Bio.Entrez
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import Progressbar

def run():
    # Sanity checks
    #app.run.config(bg='green', activebackground='green', relief=tk.SUNKEN, state='disabled')
    #assert os.path.isfile(f'{args.bowtie}bowtie2'), f"Can't find bowtie2 at {args.bowtie}"
    #assert os.path.isfile(f'{args.bowtie}bowtie2-build'), f"Can't find bowtie2-build at {args.bowtie}"
    assert os.path.isfile(app.seq_file), f"given sequencing file, '{app.seq_file}', does not exist"
    assert os.path.isfile(app.wt_file), f"given refrence/wildtype file name '{app.wtseq}' does not exist!"
    assert Bio.SeqIO.read(app.wt_file, "fasta").seq.translate(), f"given refrence/wildtype sequence file " \
                                                                f"'{app.wt_file}' is not a valid FASTA file " \
                                                                f"containing one unambiguous DNA sequence!"
    assert app.domains_file is None or os.path.isfile(app.domains_file), f"given domains file, '{app.domains_file}', does not exist"
    assert app.muts_file is None or os.path.isfile(app.muts_file), f"given domains file, '{app.muts_file}', does not exist"
    assert app.aamuts_file is None or os.path.isfile(app.aamuts_file), f"given domains file, '{app.aamuts_file}', does not exist"
    assert (app.muts_file or app.domains_file or app.aamuts_file or not app.deep_check.get()), \
        "Targeted analysis cannot be done if neither mutations file nor domains files are provided"
    #assert not app.correlations_check.get() or app.aamuts_file or app.muts_file or app.domains_file, \
    #    "Cannot analyze correlations between designed mutations if designed mutations are not provided in either " \
    #    "mutations file, amino acid mutation file, and/or domains file"
    assert app.control_correlations_file is None or os.path.isfile(app.control_correlations_file), \
        f"given refrence correlations file, '{app.control_correlations_file}', does not exist"

    # variables from arguments
    wt_seq = Bio.SeqIO.read(app.wt_file, "fasta")
    sequencing_file = app.seq_file
    paired_sequencing_file = app.paired_file
    quiet = 0

    # output files
    rootname = os.path.splitext(app.seq_file)[0]
    log_file = f'{rootname}_logfile.log'

    # initialize some variables
    programstart = time.time()

    message = 'Reference sequence read\n'
    app.output.insert('end', message)
    root.update()

    ### Process Sequencing Records ###

    # merge paired reads
    if paired_sequencing_file and not os.path.exists(rootname.split('.fastq')[0]+'_corrected.fastq.gz'):
        message = f'Merging paired reads from {sequencing_file} and {paired_sequencing_file} using bbmerge.'
        if not quiet:
            print(message)
        app.output.insert('end', 'Merging paired reads.\n')
        root.update()
        sequencing_file, paired_sequencing_file = AlignmentAnalyze.correct_pairs(sequencing_file, paired_sequencing_file)
        # mergedfile = f"{rootname}_merged.{seqreadtype}"
        # AlignmentAnalyze.merge_pairs(sequencing_file, paired_sequencing_file, mergedfile)
        # sequencing_file = mergedfile

    # align reads (using bbmap)
    if not os.path.exists(f'{rootname}.sam'):  # and app.aamuts_file:
        message = f'Aligning all sequences from {sequencing_file} to {wt_seq} using bbmap.'
        print(message)
        app.output.insert('end', 'Aligning sequencing reads to reference.\n')
        root.update()
        AlignmentAnalyze.align_all_bbmap(sequencing_file, app.wt_file, f'{rootname}.sam', max_gap=len(wt_seq), paired_sequencing_file=paired_sequencing_file)

    # determine the number of reads
    lines = sum(1 for i in open(f'{rootname}.sam', 'rb'))
    reads = int((lines-3) / 2)
    print('Total number of reads to analyze: '+str(reads))
    app.reads = reads
    app.output.insert('end', f'Total number of reads to analyze: {str(reads)}\n')
    root.update()

    # call variants / find mutations
    app.output.insert('end', f'Calling mutations/variants\n')
    root.update()
    # David's work
    AlignmentAnalyze.david_call_variants(f"{rootname}.sam", wt_seq, rootname, app, root)
    app.output.insert('end', f'Finding paired mutations\n')
    root.update()
    if app.correlations_check.get():
        AlignmentAnalyze.david_paired_analysis(rootname + '.csv', rootname + '_wt.csv', app, root)
    # os.system(f'java -cp ../bbmap/current/ var2.CallVariants2 in={rootname}.sam ref={app.wt_file} ploidy=1 out={rootname}.vcf 32bit')

    seq_analyze_time = time.time() - programstart
    time_per_seq = round(seq_analyze_time / reads * 1000, 3)
    print(f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    app.reads = reads
    app.output.insert('end', f'Analyzed {reads} in {round(seq_analyze_time, 1)} seconds. {time_per_seq} ms per read.')
    app.output.insert('end', 'Finished Analysis')
    root.update()

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.winfo_toplevel().title("Library Alignment and Analysis")
        self.quiet_check = tk.IntVar()
        self.verbose_check = tk.IntVar()
        self.deep_check = tk.IntVar()
        self.correlations_check = tk.IntVar()
        self.parallel_check = tk.IntVar()
        self.bryan = tk.IntVar()
        self.domain_check = tk.IntVar()
        self.AA_check = tk.IntVar()
        self.variant_check = tk.IntVar()

        self.seq = tk.Button(self, text='Input sequencing reads (fastq)', command=self.browse_seq)
        self.seq.pack()
        self.seq_file = None

        self.paired = tk.Button(self, text='Input paired sequencing reads (fastq)', command=self.browse_paired)
        self.paired.pack()
        self.paired_file = None

        self.wtseq = tk.Button(self, text='Input wildtype sequence (fasta)', command=self.browse_wt)
        self.wtseq.pack()
        self.wt_file = None

        self.run = tk.Button(self, text='Run analysis', command=run)
        self.run.pack()

        self.mut = tk.Button(self, text='designed mutation NT list (csv)', command=self.browse_mut)
        self.mut.pack()
        self.muts_file = None

        self.mutaa = tk.Button(self, text='designed mutation AA list (csv)', command=self.browse_mutaa)
        self.mutaa.pack()
        self.aamuts_file = None

        self.control_correlations = tk.Button(self, text='Correlation file - baseline (csv)', command=self.browse_correlation)
        self.control_correlations.pack()
        self.control_correlations_file = None

        self.domains = tk.Button(self, text='Expected mutant domains (csv)', command=self.browse_domains)
        self.domains.pack()
        self.domains_file = None

        self.variant = tk.Checkbutton(self, text='variant analysis', variable=self.variant_check)
        self.variant.pack()

        self.correlation = tk.Checkbutton(self, text='Correlation analysis', variable=self.correlations_check)
        self.correlation.pack()

        self.output = tk.Text(self, height=10, width=60)
        self.output.pack()

        self.progress = Progressbar(self, orient='horizontal', length=500, mode='determinate')
        self.progress.pack()

    def browse_seq(self):
        self.seq_file = filedialog.askopenfilename(title="Select a File")
        if self.seq_file != '':
            self.seq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.seq_file = None
    def browse_paired(self):
        self.paired_file = filedialog.askopenfilename(title="Select a File")
        if self.paired_file != '':
            self.paired.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.paired_file = None
    def browse_wt(self):
        self.wt_file = filedialog.askopenfilename(title="Select a File")
        if self.wt_file != '':
            self.wtseq.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.wt_file = None
    def browse_mut(self):
        self.muts_file = filedialog.askopenfilename(title="Select a File")
        if self.muts_file != '':
            self.mut.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.muts_file = None
    def browse_mutaa(self):
        self.aamuts_file = filedialog.askopenfilename(title="Select a File")
        if self.aamuts_file != '':
            self.mutaa.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.aamuts_file = None
    def browse_domains(self):
        self.domains_file = filedialog.askopenfilename(title="Select a File")
        if self.domains_file != '':
            self.domains.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.domains_file = None
    def browse_correlation(self):
        self.control_correlations_file = filedialog.askopenfilename(title="Select a File")
        if self.control_correlations_file != '':
            self.control_correlations.config(bg='green', activebackground='green', relief=tk.SUNKEN)
        else:
            self.control_correlations_file = None

if __name__ == "__main__":
    root = tk.Tk()
    root.geometry('700x600')
    app = Application(master=root)
    app.mainloop()
