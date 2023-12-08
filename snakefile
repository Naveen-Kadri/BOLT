
localrules:prepare_phenotypes, prepare_pc,combine_LD_score,combine_LD_prune


#FILES
plink_binaries = expand('/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/TOSEQ/{{breed}}/CHR{{chr}}/imputed_{{sex}}.{ext}',ext=['bim', 'fam', 'bed'])
phenofile = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/GWAS/{pheno}/{breed}/{sex}/std5/gam1/phenotypes.txt'
OUT_DIR = '/cluster/work/pausch/naveen/BOLT'
pc_file='/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/GWAS/GRM/{breed}/{sex}/genome.eigenvec'


#WILDCARDS
chromosomes = range (1,30)
breeds = ['bv', 'fv']
sexes = ['male', 'female']
phenos = ['GRR']

resources = {
    'bv' : {
        'male' :   ['04:00:00', 4000, 10],
        'female' : ['04:00:00', 8000, 10]
    },
    'fv' : {
        'male' : ['04:00:00', 4000, 10],
        'female':['24:00:00', 8000, 20]
    }
}


#PROGRAMS
ldsc='/cluster/home/nkadri/PROGRAMS/ldsc/ldsc.py '

rule all:
    input:
        expand(OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/manhattan_{breed}_{sex}.tiff', pheno=phenos, breed=breeds,sex=sexes),
        expand(OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/bychr_{breed}_{sex}.pdf', pheno=phenos, breed=breeds,sex=sexes),
         
rule prepare_phenotypes:
    '''header is required'''
    input:
        pheno=phenofile
    output:
        outfile=OUT_DIR + '/{pheno}/{breed}/{sex}/phenotypes.txt'
    shell:
        '''
        echo -e "FID\tIID\t{wildcards.pheno}\tnobs" > {output.outfile}
        cat {input.pheno} >>{output.outfile}
        '''


rule prepare_pc:
    '''
    adding headers
    '''
    input:
        infile=pc_file
    output:
        outfile=OUT_DIR + '/PC/{breed}_{sex}.txt'
    params:
        npc=4
    script:
        'prepare_pc.py'

        
rule LD_score:
    input:
        plink_binaries=lambda wc : [myfile.format(breed=wc.breed, chr=wc.chr,sex='male') for myfile in plink_binaries],
    output:
        outfile=expand(OUT_DIR + '/LDSCORE/{{breed}}/chr{{chr}}.{ext}', ext=['log','l2.M','l2.M_5_50','l2.ldscore.gz'])
    conda:
        'ldsc'
    params:
        outfile=lambda wc, output:output.outfile[0][:-4],
        geno=lambda wc, input:input.plink_binaries[0][:-4],
        window_size=1000
    resources:
        mem_mb=16000,
        walltime='24:00:00'
    shell:
        '''
        python {ldsc}\
	--bfile {params.geno} \
	--l2\
	--ld-wind-kb {params.window_size}\
	--out {params.outfile} \
        '''

        

rule combine_LD_score:
    input:
        infiles=expand(OUT_DIR + '/LDSCORE/{{breed}}/chr{chr}.l2.ldscore.gz', chr=chromosomes)
    output:
        outfile=OUT_DIR + '/LDSCORE/{breed}/genome.l2.ldscore.gz'
    run:
        import gzip
        with gzip.open (output.outfile, 'wt') as out:
            for fnum,myfile in enumerate(input.infiles):
                with gzip.open (myfile, 'rt') as inf:
                    for lnum, line in enumerate (inf):
                        if lnum==0:
                            if fnum==0:
                                out.write (line)
                        else:
                            out.write (line)
            

rule LD_prune:
    input:
        plink_binaries=lambda wc : [myfile.format(breed=wc.breed, chr=wc.chr,sex='male') for myfile in plink_binaries]
    output:
        outfiles=expand(OUT_DIR + '/MODELSNPS/{{breed}}/chr{{chr}}.{ext}', ext=['prune.in','prune.out'])
    params:
        geno=lambda wc,input:input.plink_binaries[0][:-4],
        out=lambda wc,output:output.outfiles[0][:-9],
        window_size='1000kb',
        step_count=1,
        r2_thresh=0.8,
        maf_thresh = 0.005
    resources:
        mem_mb=16000,
        walltime='04:00:00'
    shell:
        '''
        module load gcc/8.2.0 plink/1.9-beta6.18
        plink --cow \
        --bfile {params.geno} \
        --maf {params.maf_thresh} \
        --indep-pairwise {params.window_size} {params.step_count} {params.r2_thresh} \
        --out {params.out}  
        '''

rule combine_LD_prune:
    input:
        infiles=expand(OUT_DIR + '/MODELSNPS/{{breed}}/chr{chr}.prune.in', chr=chromosomes)
    output:
        outfile=OUT_DIR + '/MODELSNPS/{breed}/genome.prune.in'
    shell:
        '''
        cat {input.infiles} >{output.outfile}
        '''
        
rule BOLT:
    input:
        binaries =lambda wc: expand (plink_binaries, chr=chromosomes, breed=wc.breed,sex=wc.sex),
        phenofile =rules.prepare_phenotypes.output.outfile,
        covariates = rules.prepare_pc.output.outfile,
        model_snps =ancient(rules.combine_LD_prune.output.outfile),
        LD_score =ancient(rules.combine_LD_score.output.outfile)
    output:
        outfile=OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/result.txt'
    params:
        bim = lambda wc, input : input.binaries[0][:-4].replace ("CHR1", f"CHR{{{chromosomes[0]}:{chromosomes[-1]}}}")+'.bim',
        bed = lambda wc, input : input.binaries[0][:-4].replace ("CHR1", f"CHR{{{chromosomes[0]}:{chromosomes[-1]}}}")+'.bed',
        fam=lambda wc, input : input.binaries[0][:-4]+'.fam',
        qcovar_col = "PC{1:4}",
        missing_snp = 0.1,
        missing_ind = 0.1
    resources:
        walltime=lambda wc : resources [wc.breed] [wc.sex] [0],
        mem_mb=lambda wc : resources [wc.breed] [wc.sex] [1]
    threads:
        lambda wc : resources [wc.breed] [wc.sex] [2]
    envmodules:
        'intel/19.1.0',
        'bolt-lmm/2.4.1'
    shell:
        '''
        bolt --lmm \
        --Nautosomes=29 \
        --bed={params.bed} \
        --bim={params.bim} \
        --fam={params.fam} \
        --phenoFile={input.phenofile} \
        --LDscoresFile={input.LD_score} \
        --LDscoresCol=L2 \
        --phenoCol={wildcards.pheno} \
        --covarFile={input.covariates} \
        --qCovarCol={params.qcovar_col} \
        --maxMissingPerSnp={params.missing_snp} \
        --maxMissingPerIndiv={params.missing_ind} \
        --statsFile={output.outfile} \
        --verboseStats \
        --numThreads={threads} \
        --modelSnps={input.model_snps} \
        --maxModelSnps=2000000 
        '''

rule format_result:
    input:
        infile=rules.BOLT.output.outfile
    output:
        outfile=OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/forplot.txt'
    shell:
        '''
        awk 'NR>1{{ print $2,$3,$7,$14 }}' {input.infile} > {output.outfile}
        '''

rule manhattan:
    input:
        infile=rules.format_result.output.outfile
    output:
        plot_file=OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/manhattan_{breed}_{sex}.tiff'
    params:
        maf_thresh = 0.005
    resources:
        mem_mb=8000,
        walltime='01:00:00'
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    script:
        'manhattan.R'    

rule plot_by_chr:
    input:
        infile=rules.format_result.output.outfile
    output:
        plot_file=OUT_DIR + '/BOLT/{pheno}/{breed}/{sex}/bychr_{breed}_{sex}.pdf'
    params:
        maf_thresh = 0.005
    resources:
        mem_mb=8000,
        walltime='01:00:00'
    envmodules:
        'gcc/8.2.0',
        'r/4.2.2'
    script:
        'plot_by_chr.R'    

        
