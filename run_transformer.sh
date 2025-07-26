#!/bin/bash

export LD_LIBRARY_PATH=/opt/share/gcc-10.1.0//lib64:/opt/share/gcc-10.1.0//lib:/usr/local/cuda-12.2/lib64

run_dir=./runs

# Run replication for structure prediction from spectra

# Base Transformer
# python ./benchmark/generate_input.py \
#         --analytical_data  ./data/multimodal_spectroscopic_dataset \
#         --out_path ${run_dir}/h_and_c_nmr_v1 \
#         --formula \
#         --h_nmr \
#         --c_nmr

# python ./benchmark/start_training.py \
#         --output_path ${run_dir}/h_and_c_nmr_v1 \
#         --template_path ./benchmark/base_transformer_template.yaml \
#         --seed 3435

# Finetune
# python ./benchmark/start_training.py \
#         --output_path ${run_dir}/finetune_v1 \
#         --template_path ./benchmark/finetune_template.yaml \
#         --seed 3435 \
#         --train_from ${run_dir}/h_and_c_nmr/model_step_250000.pt
