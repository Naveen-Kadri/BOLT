#infile='/cluster/work/pausch/naveen/RECOMBINATION/REFALT/GRM/male/genome.eigenvec'
#npc = 4
infile =snakemake.input.infile
outfile=snakemake.output.outfile
npc =snakemake.params.npc
out = open (outfile, 'w')

header ='\t'.join  (["FID", "IID"] + ["PC"+str(num) for num in range(1, npc+1)])
out.write (f'{header}\n')

with open (infile) as inf:
    for line in inf:
        spl=line.rstrip().split()
        tw = '\t'.join (spl [:npc+2])
        out.write (f'{tw}\n')


out.close()
