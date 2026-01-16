process AF3_ALIGNMENT_PRECOMPUTE {


    // there are many problems with this process do not use until this message is removed.

    label 'af3_alignment'

    input:
        val json_file_path
        val msa_directory

    output:
        tuple path("${json_basename}"), path("${json_basename}.afm"), emit: af3_alignment

    script:
    """ 
    #!/bin/bash

    set -e  # Exit on error

    [ ! -d ${msa_directory} ] && mkdir -p ${msa_directory}

    json_basename=\$(basename ${json_file_path} .json)
    json_dir=\$(dirname ${json_file_path})
    json_file=\$(basename ${json_file_path})

    # Extract protein names from JSON filename (assuming format: ProteinA_ProteinB.json)
    proteinA=\$(echo "\${json_basename}" | cut -d'_' -f1)
    proteinB=\$(echo "\${json_basename}" | cut -d'_' -f2)

    # Define paths for precomputed individual MSAs
    MSA_CACHE_DIR="${msa_directory}"
    MSA_A="\${MSA_CACHE_DIR}/\${proteinA}.a3m"
    MSA_B="\${MSA_CACHE_DIR}/\${proteinB}.a3m"
    PAIRED_MSA="\${MSA_CACHE_DIR}/\${proteinA}_\${proteinB}.a3m"
    AFM_OUTPUT_DIR="\$PWD/afm_output"

    AFM_CACHE_DIR="${msa_directory}"
    AFM_FILE="\${AFM_CACHE_DIR}/\${proteinA}_\${proteinB}.afm"

    if [[ -f "\$AFM_FILE" ]]; then
        echo "AFM file already exists: \${AFM_FILE}. Skipping computation."
        cp "\$AFM_FILE" ./  # Copy AFM file to working directory
        exit 0
    fi

    if [[ -f "\$MSA_A" && -f "\$MSA_B" ]]; then
        echo "Using precomputed MSAs for \${proteinA} and \${proteinB} to compute paired MSA."
        
        # Generate the paired MSA with co-evolution constraints
        singularity exec \
            --bind "\$PWD:/root/af_output" \
            docker://benpasto/alphafold3:latest \
            python3 /app/alphafold/tools/pair_msa.py \
            --msa1 "\$MSA_A" \
            --msa2 "\$MSA_B" \
            --output_msa "\$PAIRED_MSA"
        
        # Run AlphaFold3 to compute AFM file from the paired MSA
        singularity exec \
            --bind "\$PWD:/root/af_output" \
            --bind "\${json_dir}:/root/af_input" \
            --bind "\$MSA_CACHE_DIR:/root/msas" \
            --bind "\$AFM_OUTPUT_DIR:/root/afm_output" \
            docker://benpasto/alphafold3:latest \
            python3 /app/alphafold/run_alphafold.py \
            --norun_inference \
            --json_path=/root/af_input/\${json_file} \
            --msa_path=/root/msas/\${proteinA}_\${proteinB}.a3m \
            --model_dir=/root/models \
            --db_dir=/root/public_databases \
            --output_dir=/root/afm_output
        
        mv /root/afm_output/features/\${proteinA}_\${proteinB}.afm "\$AFM_FILE"

    else
        echo "One or both MSAs are missing. Running full AF3 alignment (no inference) to compute all MSAs and AFM."

        singularity exec \
            --bind "\${json_dir}:/root/af_input" \
            --bind "\$PWD:/root/af_output" \
            --bind /fs/ess/PCON0160/ALPHAFOLD3/models:/root/models \
            --bind /fs/project/pub_data/alphafold3/3.0.0:/root/public_databases \
            docker://benpasto/alphafold3:latest \
            python3 /app/alphafold/run_alphafold.py \
            --norun_inference \
            --json_path=/root/af_input/\${json_file} \
            --model_dir=/root/models \
            --db_dir=/root/public_databases \
            --output_dir=/root/af_output

        # Move newly computed MSAs to cache
        #mv /root/af_output/msas/\${proteinA}.a3m "\$MSA_A"
        #mv /root/af_output/msas/\${proteinB}.a3m "\$MSA_B"
        #mv /root/af_output/msas/\${proteinA}_\${proteinB}.a3m "\$PAIRED_MSA"
        #mv /root/afm_output/features/\${proteinA}_\${proteinB}.afm "\$AFM_FILE"

    fi

    cp "\$AFM_FILE" ./
    """
}


process AF3_ALIGNMENT {

    label 'af3_alignment'

    publishDir "$params.msa_cache", mode : 'copy', pattern : "*_data.json"

    input :
        val json_file_path
    
    output : 
        path("*_data.json"), emit: af3_alignment

    script:
    """ 
    #!/bin/bash

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