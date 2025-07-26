# multimodal_nmr_transformer

## Installation

Clone the repository by running

```
git clone https://github.com/rxn4chemistry/multimodal-spectroscopic-dataset.git
```
#### Conda install 

In the project directory, create a new conda environment:
```
conda create -n multinmrtransformer python=3.10
```
Activate the created conda enviroment and insall the required dependencies

```
conda activate multinmrtransformer
pip install -r requirements.txt
```

## Data Download

Please run `download_spectra.sh` to download the simulated data for pre-training.

Please run `download_nmrshiftdb2.sh` to download experimental data from NMRShiftDB2 (used for fine-tuning).

## Training and Inference

Please run `run_transformer.sh` to generate simulated pre-training data, pre-train the model, and fine-tune the model.

Please run `python create_nmrshiftdb2_data.py` to generate the fine-tuning data.

All hyperparameters in the yaml and args in the shell script can be customized.

To run inference, use the following commands:

```
onmt_translate 
  -model <model_path> 
  -src <src_path> 
  -output <out_file> 
  -beam_size 10 
  -n_best 10 
  -min_length 3 
  -gpu 0
```

```
python benchmark/analyse_results.py 
  --pred_path <Path to the predictions made with onmt_translate>
  --test_path <Path to the ground truth test data>
```

## Acknowledgements
Would like to acknowledge the authors of "Unraveling Molecular Structure: A Multimodal Spectroscopic Dataset for Chemistry" for the simulated data and the starting point for this repository.

```
@article{alberts2024unraveling,
  title={Unraveling Molecular Structure: A Multimodal Spectroscopic Dataset for Chemistry},
  author={Alberts, Marvin and Schilter, Oliver and Zipoli, Federico and Hartrampf, Nina and Laino, Teodoro},
  year={2024},
  url={https://arxiv.org/abs/2407.17492}, 
}
```
