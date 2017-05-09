rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
           "variants.database",
           "sample_network.bmp",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/delly_out/deletions.bcf",
           "pindel_out/pindel_calls_INV",
           "cnvkit_out/germline_cnvkit.heatmap",
           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/svaba.log",
    output: 
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """