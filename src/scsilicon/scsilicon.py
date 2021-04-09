import pandas as pd
import numpy as np
import sys
import os
import tasklogger
import wget
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
from . import utils

def download_ref_data(params):
    ref_path = os.path.join(params.data_dir, params.ref)
    if not os.path.exists(ref_path):
        os.makedirs(ref_path)
    tasklogger.log_info('downloading reference files...')
    #download ref fasta
    url = 'http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.gz'.format(params.ref)
    wget.download(url, out=ref_path)

    # download each chrom fasta file
    chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    for chrom in chroms:
        url = 'http://hgdownload.cse.ucsc.edu/goldenpath/{0}/chromosomes/{1}.fa.gz'.format(params.ref, chrom)
        wget.download(url, out=ref_path)

    #download chrome size file
    url = 'http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.chrom.sizes'.format(params.ref)
    wget.download(url, out=ref_path)

    #download dbsnp file
    command = 'curl -s http://hgdownload.cse.ucsc.edu/goldenPath/{0}/database/snp151Common.txt.gz'.format(params.ref) + ''' |zcat|grep -v "^#" | grep "single" | grep "exact" | awk '{print $5 "\t" $2 "\t" $3+1 "\t" $10 "\t" $7 "\t" $8}' > ''' + '{0}/dbsnp.{1}.bed'.format(ref_path, params.ref)
    code = os.system(command)
    if code != 0:
        tasklogger.log_error('Errors happend when downloading dbsnp file...')

class SCSiliconParams:

    def __init__(self, out_dir='./', data_dir='./data', ref='hg19', chrom='chr22', layout='SE', coverage=5, isize=260, threads=1, verbose=1):
        self.out_dir = out_dir
        self.data_dir = data_dir
        self.ref = ref
        self.chrom = chrom
        self.layout = layout
        self.coverage = coverage
        self.isize = isize
        self.threads = threads
        if self.chrom != 'all':
            self.ref_file = os.path.join(self.data_dir, self.ref, self.chrom + '.fa.gz')
        else:
            self.ref_file = os.path.join(self.data_dir, self.ref, self.ref + '.fa.gz')
        self.chrom_size_file = os.path.join(self.data_dir, self.ref, self.ref + '.chrom.sizes')
        self.dbsnp_file = os.path.join(self.data_dir, self.ref, 'dbsnp.'+ self.ref +'.bed')
        self.profile_file = os.path.join(utils.root_path(), 'data/normal.profile')
        self._check_params()
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_in(['hg19', 'hg38'], ref=self.ref)
        utils.check_in(['SE', 'PE'], layout=self.layout)
        utils.check_positive(coverage=self.coverage)
        utils.check_int(isize=self.isize)
        utils.check_positive(isze=self.isize)
        utils.check_int(threads=self.threads)
        utils.check_positive(threads=self.threads)


    def set_params(self, **params):
        """Set the parameters of SCSilicon.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        out_dir : string, optional, default: './'
            The output directory path
        
        ref : string, optional, default: hg19
            The reference genome version: hg19 or hg38

        chrom: string, optional, default: chr22
            The chromosome number for reads generation:  all or a specific chromosome

        layout: string, optional, default: 'PE'
            The reads laryout: PE or SE (PD for paired-end and SE for single-end)

        coverage: int, optional, default: 5
            The sequencing coverage
        
        isize: int, optional, default: 260
            The mean insert size for paired-end sequencing

        threads: int, optional, default: 1
            The number of threads to use for reads generation
       
        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'out_dir' in params and params['out_dir'] != self.out_dir:
            self.out_dir = params['out_dir']
            del params['out_dir']
        if 'ref' in params and params['ref'] != self.ref:
            self.ref = params['ref']
            del params['ref']
        if 'chrom' in params and params['chrom'] != self.chrom:
            self.chrom = params['chrom']
            del params['chrom']
        if 'layout' in params and params['layout'] != self.layout:
            self.layout = params['layout']
            del params['layout']
        if 'coverage' in params and params['coverage'] != self.coverage:
            self.coverage = params['coverage']
            del params['coverage']
        if 'isize' in params and params['isize'] != self.isize:
            self.isize = params['isize']
            del params['isize']
        if 'threads' in params and params['threads'] != self.threads:
            self.threads = params['threads']
            del params['threads']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))


class SNPSimulator:
    def __init__(self, cell_no=1, snp_no=1000, verbose=1):
        self.cell_no = cell_no
        self.snp_no = snp_no
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_int(cell_no=self.cell_no)
        utils.check_positive(cell_no=self.cell_no)
        utils.check_int(snp_no=self.snp_no)
        utils.check_positive(snp_no=self.snp_no)

    def set_params(self, **params):
        """Set the parameters of SCSilicon.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        cell_no : int, optional, default: 1
            The sample number for this simulation
        
        snp_no : int, optional, default: 1000
            The SNP number of each sample
        
        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.cell_no = params['cell_no']
            del params['cell_no']
        if 'snp_no' in params and params['snp_no'] != self.snp_no:
            self.snp_no = params['snp_no']
            del params['snp_no']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _generate_snp_file(self, params):
        tasklogger.log_info('start generating snp file...')
        all_snps = pd.read_csv(params.dbsnp_file, header=None, sep='\t')
        if params.chrom != 'all':
            all_snps = all_snps[all_snps[1]==params.chrom]

        # randomm select some ratio snps
        snps = all_snps.sample(n=int(self.snp_no * 1.2))

        #out snp file
        for sample in self.samples:
            snp_file =os.path.join(params.out_dir, sample + '-snps.txt')
            self.samples[sample]['snp_file'] = snp_file
            sample_snp = snps.sample(n=self.snp_no)
            sample_snp.to_csv(snp_file, index=False, header=None, sep='\t')
    
    def _generate_reads_for_snp(self, params):
        for sample in self.samples:
            tasklogger.log_info('start generating reads for '+sample+'...')
            fasta_file = os.path.join(params.out_dir, sample + '-ref.fasta')
            # self.samples[sample]['fasta_file'] = fasta_file
            command = 'scssim simuvars -r {0} -s {1} -o {2}'.format(params.ref_file, self.samples[sample]['snp_file'], fasta_file)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating varations with snp file...')
            
            #generate reads
            out_prefix = os.path.join(params.out_dir, sample)
            command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(fasta_file, params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating reads...')
            
            # clean sample fasta file
            os.remove(fasta_file)
            os.remove(fasta_file+'.fai')
    
    # def _generate_reads_for_snp(self, params):
    #     tasklogger.log_info('start genrating reads with fastq format...')
    #     for sample in self.samples:
    #         out_prefix = os.path.join(params.out_dir, sample)
    #         command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(self.samples[sample]['fasta_file'], params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
    #         code = os.system(command)
    #         if code != 0:
    #             tasklogger.log_error('Errors happend when simulating reads...')
            
    #         # clean sample fasta file
    #         os.remove(self.samples[sample]['fasta_file'])

    def sim_samples(self, params):
        self._generate_snp_file(params)
        # self._sim_fasta_for_snps(params)
        self._generate_reads_for_snp(params)

class CNVSimulator:
    def __init__(self, cell_no=1, bin_len=500000, seg_no=10, cluster_no=1, normal_frac=0.4, noise_frac=0.1, verbose=1):
        self.cell_no = cell_no
        self.bin_len = bin_len
        self.seg_no = seg_no
        self.cluster_no = cluster_no
        self.normal_frac = normal_frac
        self.noise_frac = noise_frac
        self._check_params()
        self.cnv = None
        self.segments = None
        self.clusters = None
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_int(cell_no=self.cell_no)
        utils.check_positive(cell_no=self.cell_no)
        utils.check_int(bin_len=self.bin_len)
        utils.check_positive(bin_len=self.bin_len)
        utils.check_int(seg_no=self.seg_no)
        utils.check_positive(seg_no=self.seg_no)
        utils.check_int(cluster_no=self.cluster_no)
        utils.check_positive(cluster_no=self.cluster_no)
        utils.check_between(0,1,normal_frac=self.normal_frac)
        utils.check_between(0,1,noise_frac=self.noise_frac)

    def set_params(self, **params):
        """Set the parameters of SCSilicon.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        cell_no : int, optional, default: 1
            The sample number for this simulation
        
        bin_len : int, optional, default: 500000
            The fixed bin length
        
        seg_no : int, optional, default: 10
            The segment number for each chromosome

        cluster_no : int, optional, default: 1
            The cell cluster number

        normal_frac : float, optional, default: 0.4
            The fraction of normal cells

        noise_frac : float, optional, default: 0.1
            The noise fraction for cnv matrix
       
        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.cell_no = params['cell_no']
            del params['cell_no']
        if 'bin_len' in params and params['bin_len'] != self.bin_len:
            self.bin_len = params['bin_len']
            del params['bin_len']
        if 'seg_no' in params and params['seg_no'] != self.seg_no:
            self.seg_no = params['seg_no']
            del params['seg_no']
        if 'cluster_no' in params and params['cluster_no'] != self.cluster_no:
            self.cluster_no = params['cluster_no']
            del params['cluster_no']
        if 'normal_frac' in params and params['normal_frac'] != self.normal_frac:
            self.normal_frac = params['normal_frac']
            del params['normal_frac']
        if 'noise_frac' in params and params['noise_frac'] != self.noise_frac:
            self.noise_frac = params['noise_frac']
            del params['noise_frac']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _split_chr_to_bins(self, params):
        """Split chromosomes to fixed-lenght bins

        Parameters
        ----------
        bin_len : int
            fixed-bin-length

        Returns
        -------
        ref: Dataframe of pandas
        """
        tasklogger.log_info('Splitting chromosome into bins with a fixed length of ' + str(self.bin_len) + '...')
        ref = pd.DataFrame(
            columns=['Chromosome', 'Start', 'End'])
        chrom_sizes = pd.read_csv(params.chrom_size_file, header=None, sep='\t')
        
        if params.chrom != 'all':
            chrom = params.chrom
            chrom_size = int(chrom_sizes[chrom_sizes[0]==chrom][1])
            tasklogger.log_info('Splitting ' + chrom + '...')
            start = 1
            end = self.bin_len
            count = 1
            while(start <= chrom_size):
                ref = ref.append({
                    'Chromosome': chrom,
                    'Start': start,
                    'End': min(end, chrom_size),
                }, ignore_index=True)
                count += 1
                start = end + 1
                end = self.bin_len * count     
        else:
            for index, row in chrom_sizes.iterrows():
                [chrom, chrom_size] = row
                tasklogger.log_info('Splitting ' + chrom + '...')
                start = 1
                end = self.bin_len
                count = 1
                while(start <= chrom_size):
                    ref = ref.append({
                        'Chromosome': chrom,
                        'Start': start,
                        'End': min(end, chrom_size),
                    }, ignore_index=True)
                    count += 1
                    start = end + 1
                    end = self.bin_len * count            
        return ref

    def _generate_cnv_matrix(self, params):
        ref = self._split_chr_to_bins(params)
        columns = ['cell' + str(i+1) for i in range(self.cell_no)]
        cnv = pd.DataFrame(columns=columns)
        segments = pd.DataFrame(columns=['Chromosome', 'Start_bin_no', 'End_bin_no'])
        clusters_file = pd.DataFrame(columns=['Cell', 'Cluster'])

        if params.chrom != 'all':
            all_chroms = [params.chrom]
        else:
            all_chroms = sorted(np.unique(ref['Chromosome']), key=lambda x: int(x[3:]))

        # generate cell clusters
        if self.cluster_no < 1 or self.cluster_no > self.cell_no*(1-self.normal_frac) - 1:
            self.cluster_no = 1
        if self.cluster_no == 1:
            clusters  = [0, self.cell_no]
        else:
            clusters = np.sort(np.random.choice(range(int(self.cell_no*self.normal_frac)+1, self.cell_no-1), self.cluster_no - 2, replace=False))
            clusters = [0] + [int(self.cell_no*self.normal_frac)] + list(clusters) + [self.cell_no]
        cluster_count = 1
        for i in range(len(clusters) - 1):
            start = clusters[i]
            end = clusters[i+1]
            for j in range(start, end):
                clusters_file = clusters_file.append([{
                    'Cell': 'cell' + str(j+1),
                    'Cluster': cluster_count,
                }], ignore_index=True)
            cluster_count += 1


        for chrom in all_chroms:
            chrom_ref = ref[ref['Chromosome'] == chrom]
            total_bins = chrom_ref.shape[0]
            sub_matrix = np.zeros((total_bins, self.cell_no),dtype=int)
            if self.seg_no < 1 or self.seg_no > total_bins:
                self.seg_no = 10
            if self.seg_no == 1:
                breaks = []
            else:
                breaks = np.sort(np.random.choice(range(1, total_bins-1), self.seg_no - 1, replace=False))
            
            breaks = [0] + list(breaks) + [total_bins]
            
            for i in range(len(breaks)-1):
                for j in range(len(clusters)-1):
                    sub_matrix[breaks[i]:breaks[i+1],clusters[j]:clusters[j+1]] = np.random.randint(1,10)
            index = chrom_ref['Chromosome'] + ':' + chrom_ref['Start'].astype(str) + '-' + chrom_ref['End'].astype(str)
            sub_matrix = pd.DataFrame(sub_matrix, columns=columns, index=index)
            cnv = cnv.append(sub_matrix)
            for i in range(len(breaks) - 1):
                start = breaks[i]
                end = breaks[i+1] - 1
                segments = segments.append([{
                    'Chromosome': chrom,
                    'Start_bin_no': start,
                    'End_bin_no': end
                }], ignore_index=True)
        cnv.iloc[:,0:max(0, round(self.cell_no*self.normal_frac))] = 2

        #add noise
        total_values = cnv.shape[0] * cnv.shape[1]
        noise = np.random.randint(0, total_values+1, size=round(total_values*self.noise_frac))
        for i in noise:
            cnv.iat[int((i-1)/cnv.shape[1]), (i-1) % cnv.shape[1]] = np.random.randint(1, 10)

        cnv = cnv.T
        cnv = cnv.astype('int32')
        self.cnv = cnv
        self.segments = segments
        self.clusters = clusters_file

        self.cnv.to_csv(os.path.join(params.out_dir, 'cnv.csv'))
        self.segments.to_csv(os.path.join(params.out_dir, 'segments.csv'))
        self.clusters.to_csv(os.path.join(params.out_dir, 'clusters.csv'))

    def visualize_cnv_matrix(self, out_prefix):
        cmap = sns.cubehelix_palette(start = 1.5, rot = 3, gamma=0.8, as_cmap = True)
        ax = sns.heatmap(self.cnv, cmap=cmap, xticklabels=False, yticklabels=False)
        plt.savefig(out_prefix + '.pdf', dpi=300)

    def _generate_cnv_bed_file(self, params):
        for sample, cnvs in self.cnv.iterrows():
            bed_file = os.path.join(params.out_dir, sample + '.bed')
            self.samples[sample]['bed_file'] = bed_file
            with open(bed_file, 'w') as out_bed:
                for index, value in enumerate(cnvs):
                    if value == 2:
                        continue
                    chrom = cnvs.index[index].split(':')[0]
                    start = cnvs.index[index].split(':')[1].split('-')[0]
                    end = cnvs.index[index].split(':')[1].split('-')[1]
                    out_bed.write('\t'.join(['c', chrom, start, end, str(value), str(np.random.randint(1, value+1))])+'\n')
    
    def _generate_reads_for_cnv(self, params):
        for sample in self.samples:
            tasklogger.log_info('start generating reads for '+sample+'...')
            fasta_file = os.path.join(params.out_dir, sample + '-ref.fasta')
            # self.samples[sample]['fasta_file'] = fasta_file
            command = 'scssim simuvars -r {0} -v {1} -o {2}'.format(params.ref_file, self.samples[sample]['bed_file'], fasta_file)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating varations with cnv file...')

            # generate reads
            out_prefix = os.path.join(params.out_dir, sample)
            command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(fasta_file, params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating reads...')

            # clean sample fasta and bed file
            os.remove(fasta_file)
            os.remove(fasta_file+'.fai')
            os.remove(self.samples[sample]['bed_file'])

    # def _generate_reads_for_cnv(self, params):
    #     tasklogger.log_info('start genrating reads with fastq format...')
    #     for sample in self.samples:
    #         out_prefix = os.path.join(params.out_dir, sample)
    #         command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(self.samples[sample]['fasta_file'], params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
    #         code = os.system(command)
    #         if code != 0:
    #             tasklogger.log_error('Errors happend when simulating reads...')

    def sim_samples(self, params):
        self._generate_cnv_matrix(params)
        self._generate_cnv_bed_file(params)
        # self._sim_fasta_for_cnv(params)
        self._generate_reads_for_cnv(params)

class IndelSimulator:
    def __init__(self, cell_no=1, in_no=100, del_no=100, verbose=1):
        self.cell_no = cell_no
        self.in_no = in_no
        self.del_no = del_no       
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_int(cell_no=self.cell_no)
        utils.check_positive(cell_no=self.cell_no)
        utils.check_int(in_no=self.in_no)
        utils.check_positive(in_no=self.in_no)
        utils.check_int(del_no=self.del_no)
        utils.check_positive(del_no=self.del_no)

    def set_params(self, **params):
        """Set the parameters of SCSilicon.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        cell_no : int, optional, default: 1
            The sample number for this simulation
        
        in_no : int, optional, default: 1000
            The insertion number of each sample
        
        del_no : int, optional, default: 1000
            The deletion number of each sample
        
        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.cell_no = params['cell_no']
            del params['cell_no']
        if 'in_no' in params and params['in_no'] != self.snp_no:
            self.snp_no = params['in_no']
            del params['in_no']
        if 'del_no' in params and params['del_no'] != self.snp_no:
            self.snp_no = params['del_no']
            del params['del_no']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _generate_random_seq(self):
        chrs = ['a', 't', 'c', 'g']
        seq_len = np.random.choice(range(4, 10))
        seq = ''
        for i in range(seq_len):
            seq += np.random.choice(chrs)
        return seq

    def _generate_indel_file(self, params):
        tasklogger.log_info('start generating indel profile file...')
        chrom_sizes = pd.read_csv(params.chrom_size_file, header=None, sep='\t')
        if params.chrom != 'all':
            chrom = params.chrom
            chrom_size = int(chrom_sizes[chrom_sizes[0]==chrom][1])
        else:
            chroms = {}
            for index, row in chrom_sizes.iterrows():
                [chrom, chrom_size] = row
                chroms[chrom] = chrom_size

        #out indel file
        for sample in self.samples:
            indel_file =os.path.join(params.out_dir, sample + '-indel.txt')
            self.samples[sample]['indel_file'] = indel_file
            
            if params.chrom != 'all':
                ins = np.sort(np.random.choice(range(1, chrom_size), self.in_no, replace=False))
                dels = np.sort(np.random.choice(range(1, chrom_size), self.del_no, replace=False))
                with open(indel_file, 'w') as output:
                    for insert in ins:
                        output.write('\t'.join(['i', chrom, str(insert), self._generate_random_seq(),np.random.choice(['homo', 'het'])])+'\n')
                    for delete in dels:
                        output.write('\t'.join(['d', chrom, str(delete), str(np.random.choice(range(4, 10))),np.random.choice(['homo', 'het'])])+'\n')
            else:
                with open(indel_file, 'w') as output:
                    for i in range(self.in_no):
                        chrom = np.random.choice(chroms.keys())
                        chrom_size = chroms[chrom]
                        insert = np.random.choice(range(1, chrom_size))
                        output.write('\t'.join(['i', chrom, str(insert), self._generate_random_seq(),np.random.choice(['homo', 'het'])])+'\n')
                    for j in range(self.del_no):
                        chrom = np.random.choice(chroms.keys())
                        chrom_size = chroms[chrom]
                        delete = np.random.choice(range(1, chrom_size))
                        output.write('\t'.join(['d', chrom, str(delete), str(np.random.choice(range(4, 10))),np.random.choice(['homo', 'het'])])+'\n')
    
    def _generate_reads_for_indel(self, params):
        for sample in self.samples:
            tasklogger.log_info('start generating reads for '+sample+'...')
            fasta_file = os.path.join(params.out_dir, sample + '-ref.fasta')
            # self.samples[sample]['fasta_file'] = fasta_file
            command = 'scssim simuvars -r {0} -v {1} -o {2}'.format(params.ref_file, self.samples[sample]['indel_file'], fasta_file)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating varations with indel file...')
            
            #generate reads
            out_prefix = os.path.join(params.out_dir, sample)
            command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(fasta_file, params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating reads...')
            
            # clean sample fasta file
            os.remove(fasta_file)
            os.remove(fasta_file+'.fai')
    
    # def _generate_reads_for_snp(self, params):
    #     tasklogger.log_info('start genrating reads with fastq format...')
    #     for sample in self.samples:
    #         out_prefix = os.path.join(params.out_dir, sample)
    #         command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(self.samples[sample]['fasta_file'], params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
    #         code = os.system(command)
    #         if code != 0:
    #             tasklogger.log_error('Errors happend when simulating reads...')
            
    #         # clean sample fasta file
    #         os.remove(self.samples[sample]['fasta_file'])

    def sim_samples(self, params):
        self._generate_indel_file(params)
        # self._sim_fasta_for_snps(params)
        self._generate_reads_for_indel(params)

class SNVSimulator:
    def __init__(self, cell_no=1, snv_no=100, verbose=1):
        self.cell_no = cell_no
        self.snv_no = snv_no
        self._check_params()
        self.samples = dict.fromkeys(['cell' + str(i+1) for i in range(self.cell_no)])
        for sample in self.samples:
            self.samples[sample] = {}
        tasklogger.set_level(verbose)


    def _check_params(self):
        """Check SCSilicon parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as threads='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_int(cell_no=self.cell_no)
        utils.check_positive(cell_no=self.cell_no)
        utils.check_int(snv_no=self.snv_no)
        utils.check_positive(snv_no=self.snv_no)

    def set_params(self, **params):
        """Set the parameters of SCSilicon.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        cell_no : int, optional, default: 1
            The sample number for this simulation
        
        snv_no : int, optional, default: 1000
            The SNV number of each sample
        
        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print log messages

        Returns
        -------
        self
        """

        # parameters
        if 'cell_no' in params and params['cell_no'] != self.cell_no:
            self.cell_no = params['cell_no']
            del params['cell_no']
        if 'snv_no' in params and params['snv_no'] != self.snv_no:
            self.snv_no = params['snv_no']
            del params['snv_no']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']
        self._check_params()
        self.get_params()
        return self

    def get_params(self):
        print(vars(self))

    def _generate_snv_file(self, params):
        tasklogger.log_info('start generating snv file...')
        chrom_sizes = pd.read_csv(params.chrom_size_file, header=None, sep='\t')
        if params.chrom != 'all':
            chrom = params.chrom
            chrom_size = int(chrom_sizes[chrom_sizes[0]==chrom][1])
        else:
            chroms = {}
            for index, row in chrom_sizes.iterrows():
                [chrom, chrom_size] = row
                chroms[chrom] = chrom_size

        #out indel file
        for sample in self.samples:
            snv_file =os.path.join(params.out_dir, sample + '-snv.txt')
            self.samples[sample]['snv_file'] = snv_file
            
            if params.chrom != 'all':
                snvs = np.sort(np.random.choice(range(1, chrom_size), self.snv_no, replace=False))
                with open(snv_file, 'w') as output:
                    for snv in snvs:
                        alleles = np.random.choice(['a','t','g','c'],2, replace=False)
                        output.write('\t'.join(['s', chrom, str(snv), alleles[0], alleles[1], np.random.choice(['homo', 'het'])])+'\n')
            else:
                with open(snv_file, 'w') as output:
                    for i in range(self.snv_no):
                        chrom = np.random.choice(chroms.keys())
                        chrom_size = chroms[chrom]
                        snv = np.random.choice(range(1, chrom_size))
                        alleles = np.random.choice(['a','t','g','c'],2, replace=False)
                        output.write('\t'.join(['s', chrom, str(snv), alleles[0], alleles[1], np.random.choice(['homo', 'het'])])+'\n')
    
    def _generate_reads_for_snv(self, params):
        for sample in self.samples:
            tasklogger.log_info('start generating reads for '+sample+'...')
            fasta_file = os.path.join(params.out_dir, sample + '-ref.fasta')
            # self.samples[sample]['fasta_file'] = fasta_file
            command = 'scssim simuvars -r {0} -v {1} -o {2}'.format(params.ref_file, self.samples[sample]['snv_file'], fasta_file)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating varations with snv file...')
            
            #generate reads
            out_prefix = os.path.join(params.out_dir, sample)
            command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(fasta_file, params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
            code = os.system(command)
            if code != 0:
                tasklogger.log_error('Errors happend when simulating reads...')
            
            # clean sample fasta file
            os.remove(fasta_file)
            os.remove(fasta_file+'.fai')
    
    # def _generate_reads_for_snp(self, params):
    #     tasklogger.log_info('start genrating reads with fastq format...')
    #     for sample in self.samples:
    #         out_prefix = os.path.join(params.out_dir, sample)
    #         command = 'scssim genreads -i {0} -m {1} -l {2} -c {3} -s {4} -t {5} -o {6}'.format(self.samples[sample]['fasta_file'], params.profile_file, params.layout, params.coverage, params.isize, params.threads, out_prefix)
    #         code = os.system(command)
    #         if code != 0:
    #             tasklogger.log_error('Errors happend when simulating reads...')
            
    #         # clean sample fasta file
    #         os.remove(self.samples[sample]['fasta_file'])

    def sim_samples(self, params):
        self._generate_snv_file(params)
        # self._sim_fasta_for_snps(params)
        self._generate_reads_for_snv(params)
