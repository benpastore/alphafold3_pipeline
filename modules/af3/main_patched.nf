process MAKE_JSONS { 

    label 'base'

    input : 
        val fasta1
        val fasta2
    
    output : 
        path("*.json"), emit : jsons
    
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

    // Publish all MSA outputs to the cache
    publishDir "$params.msa_cache", mode: 'copy', pattern: "*_data.json"

    input:
        val(json_file_path)

    output:
        path("*_data.json"), emit: af3_alignment

    script:
    """
    #!/bin/bash
    
    json_basename=\$(basename ${json_file_path} .json)
    json_dir=\$(dirname ${json_file_path})
    json_file=\$(basename ${json_file_path})

    proteinA=\$(echo "\${json_basename}" | cut -d'_' -f1)
    A_lower=\$(echo "\$proteinA" | awk '{print tolower(\$0)}')
    cached_file="${params.msa_cache}/\${json_basename}_data.json"

    if [[ -f "\$cached_file" ]]; then
        echo "Cached MSA found at \$cached_file — using cached result."
        cp \$cached_file ./\${json_basename}_data.json
    else
        echo "No cached MSA found — running alignment."

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

        # Rename output for consistent naming
        mv \$PWD/\${A_lower}/\${A_lower}_data.json ./\${json_basename}_data.json
    fi
    """
}


process AF3_INFERENCE {


    //errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' } 
    //maxRetries 3

    label 'af3_inference'

    input : 
        val(json_file_path)
    
    output : 
        path("*.pae_min.tsv")
        path("*_data")
        path("pae_min.tmp"), emit : summary_conf

    publishDir {
       def base = json_file_path.getFileName().toString().replaceFirst(/\.json$/, '')
       def prefix = base.replaceFirst(/_data$/, '')
       return "${params.inference_models}/${prefix}"
    }, mode: 'move', pattern: '!*.tmp'
    
    script:
    """
    #!/bin/bash

    MODEL_DB=/fs/ess/PCON0160/ALPHAFOLD3/models
    AF_DB=/fs/project/pub_data/alphafold3/3.0.0

    json_basename=\$(basename ${json_file_path} .json)
    json_dir=\$(dirname ${json_file_path})
    json_file=\$(basename ${json_file_path})
    proteinA=\$(echo "\${json_basename}" | cut -d'_' -f1)
    A_lower=\$(echo "\$proteinA" | awk '{print tolower(\$0)}')

    # move MODEL_DB, AF_DB, 


    singularity exec --nv\
        --bind \${json_dir}:/root/af_input \
        --bind \$PWD:/root/af_output \
        --bind /fs/ess/PCON0160/ALPHAFOLD3/models:/root/models \
        --bind /fs/project/pub_data/alphafold3/3.0.0:/root/public_databases \
        docker://benpasto/alphafold3:latest \
        python3 /app/alphafold/run_alphafold.py \
        --norun_data_pipeline \
        --json_path="/root/af_input/\${json_file}" \
        --model_dir=/root/models \
        --db_dir=/root/public_databases \
        --output_dir=/root/af_output
    
    python3 ${params.bin}/get_sum_conf.py \$PWD/\${A_lower} \${json_basename}
    echo \$PWD/\${json_basename}
    mv \$A_lower \$json_basename
    
    cp \$PWD/\${json_basename}.pae_min.tsv pae_min.tmp
    
    """
}

process COMBINED_CONFIDENCE_SUMMARY { 

    label 'base'

    publishDir "$params.results", mode : 'move', pattern : 'merged_pae_scores.tsv'

    input : 
        val summary_confs
    
    output : 
        path("*.tsv")
    
    script : 
    """
    #!/bin/bash

    echo -e "pae_min\tsample" > header 
    cat header ${summary_confs.join(' ')} | grep -v pae > tmp
    cat header tmp > merged_pae_scores.tsv

    """

}