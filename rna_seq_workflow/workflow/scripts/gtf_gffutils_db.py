import gffutils

gtf_db = gffutils.create_db(snakemake.input["gtf"], snakemake.output["gtf_database"])