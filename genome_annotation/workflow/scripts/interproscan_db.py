import gffutils

interpro_results = gffutils.create_db(snakemake.input["interpro_gff"], snakemake.output["interpro_results_db"])