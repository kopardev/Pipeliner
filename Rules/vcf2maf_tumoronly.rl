rule vcf2maf_tumoronly_mutect2:
  input: vcf=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf"
  output: maf=config['project']['workpath']+"/mutect2_out/oncotator_out/{x}.maf",
          vcf=temp("mutect2_out/{x}.FINALmutect2.vep.vcf"),
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['units'][wildcards.x],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load VEP/92; module load vcf2maf/1.6.16; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 2 --vcf-tumor-id {params.tumorsample} --tumor-id {params.tumorsample} --ncbi-build {params.build} --species {params.species}"

rule vcf2maf_tumoronly_mutect:
  input: vcf=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf"
  output: maf=config['project']['workpath']+"/mutect_out/oncotator_out/{x}.maf",
          vcf=temp("mutect_out/{x}.FINAL.vep.vcf"),
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['units'][wildcards.x],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load VEP/92; module load vcf2maf/1.6.16; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 2 --vcf-tumor-id {params.tumorsample} --tumor-id {params.tumorsample} --ncbi-build {params.build} --species {params.species}"

rule vcf2maf_tumoronly_vardict:
  input: vcf=config['project']['workpath']+"/vardict_out/{x}.FINAL.vcf"
  output: maf=config['project']['workpath']+"/vardict_out/oncotator_out/{x}.maf",
          vcf=temp("vardict_out/{x}.FINAL.vep.vcf"),
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['units'][wildcards.x],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load VEP/92; module load vcf2maf/1.6.16; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 2 --vcf-tumor-id {params.tumorsample} --tumor-id {params.tumorsample} --ncbi-build {params.build} --species {params.species}"

rule vcf2maf_tumoronly_merged:
  input: vcf=config['project']['workpath']+"/merged_somatic_variants/{x}.merged.vcf",
  output: maf=config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{x}.maf",
          vcf=temp("merged_somatic_variants/{x}.merged.vep.vcf"),
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['units'][wildcards.x],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load VEP/92; module load vcf2maf/1.6.16; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --retain-info set --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 2 --vcf-tumor-id {params.tumorsample} --tumor-id {params.tumorsample} --ncbi-build {params.build} --species {params.species}"