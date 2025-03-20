process MAKE_JSONS { 

    label 'base'

    input : 
        val fasta1
        val fasta2
    
    output : 
        tuple path("*.json"), emit : jsons
    
    script : 
    """
    #!/bin/bash

    module load miniconda3/24.1.2-py310

    conda activate rnaseq_basic

    python3 ${params.bin}/make_json_from_fasta.py -fa1 ${fasta1} -fa2 ${fasta2}

    """
}

process AF3_ALIGNMENT {

    label 'af3_alignment'

    input :
        val json_file_path
    
    output : 
        tuple path("${json_basename}"), path("${json_basename}_data.json"), emit: af3_alignment

    script:
    """ 
    #!/bin/bash

    json_basename=\$(${json_file_path} .json)

    json_basename=\$(basename ${json_file_path} .json)
    json_dir=\$(dirname ${json_file_path})
    json_file=\$(basename ${json_file_path})

    proteinA=\$(echo "\${json_basename}" | cut -d'_' -f1)

    singularity exec \\
        --bind "\${json_dir}:/root/af_input" \\
        --bind "\$PWD:/root/af_output" \\
        --bind /fs/ess/PCON0160/ALPHAFOLD3/models:/root/models \\
        --bind /fs/project/pub_data/alphafold3/3.0.0:/root/public_databases \\
        docker://benpasto/alphafold3:latest \\
        python3 /app/alphafold/run_alphafold.py \\
        --norun_inference \\
        --json_path=/root/af_input/\${json_file} \\
        --model_dir=/root/models \\
        --db_dir=/root/public_databases \\
        --output_dir=/root/af_output

    A_lower=\$(echo "\$proteinA" | awk '{print tolower(\$0)}')

    dir=\$PWD
    mv \$A_lower ./\$json_basename
    cd ./\$json_basename
    mv \${A_lower}_data.json \${json_basename}_data.json
    cd \$dir
    mv ./\$json_basename/\${json_basename}_data.json .

    """
}

process AF3_INFERENCE {


    label 'af3_inference'

    input : 
        tuple val(json_dir), val(json_file)
    
    output : 
        path("*.pae_min.tsv"), emit : summary_conf

    script:
    """ 
    #!/bin/bash

    json_basename=\$(${json_file} .afm)

    singularity exec --nv\
        --bind ${json_dir}:/root/af_input \
        --bind \$PWD:/root/af_output \
        --bind /fs/ess/PCON0160/ALPHAFOLD3/models:/root/models \
        --bind /fs/project/pub_data/alphafold3/3.0.0:/root/public_databases \
        docker://benpasto/alphafold3:latest \
        python3 /app/alphafold/run_alphafold.py \
        --norun_data_pipeline \
        --json_path="/root/af_input/${json_file}" \
        --model_dir=/root/models \
        --db_dir=/root/public_databases \
        --output_dir=/root/af_output
    
    python3 ${params.bin}/get_sum_conf.py \$PWD/\${json_basename} \${json_basename}

    """
}

process COMBINED_CONFIDENCE_SUMMARY { 

    label 'base'

    publishDir "$params.results", mode : 'copy', pattern : 'merged_pae_scores.tsv'

    input : 
        val summary_confs
    
    output : 
        path "merged_pae_scores.tsv"
    
    script : 
    """
    #!/bin/bash

    cat ${summary_confs.join(' ')} > merged_pae_scores.tsv

    """

}