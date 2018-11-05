# Bayesian Learning of Latent Representations of Language Structures

## Requirements

- Python3
  - numpy
- R (for missing data imputation)
  - missMDA package (required for now but it is not difficult to remove dependency on it)
  - NPBayesImpute (only for comparison)
 
## Preprocessing

### WALS

- Download wals_language.csv.zip from WALS http://wals.info/ to obtain data/wals/language.csv (already in our repository)

- Convert the CSV into two JSON files

```sh
python format_wals.py ../data/wals/language.csv ../data/wals/langs.json ../data/wals/flist.json
```

- Missing data imputation for initialization

```sh
python -mmv.json2tsv ../data/wals/langs.json ../data/wals/flist.json ../data/wals/langs.tsv
R --vanilla -f mv/impute_mca.r --args ../data/wals/langs.tsv ../data/wals/langs.filled.tsv
python -mmv.tsv2json ../data/wals/langs.json ../data/wals/langs.filled.tsv ../data/wals/flist.json ../data/wals/langs.filled.json
```

TODO: Remove the dependency on missMDA as our model is now insensitive to initialization.

### Autotyp

- Suppose we are at ~/download. First download the database.
  
```sh
git clone git@github.com:autotyp/autotyp-data.git
```

- Convert the data into two JSON files

```sh
python format_autotyp.py ~/download/autotyp-data ../data/autotyp/langs.json ../data/autotyp/flist.json
```

- Missing data imputation for initialization

```sh
python -mmv.json2tsv ../data/autotyp/langs.json ../data/autotyp/flist.json ../data/autotyp/langs.tsv
R --vanilla -f mv/impute_mca.r --args ../data/autotyp/langs.tsv ../data/autotyp/langs.filled.tsv
python -mmv.tsv2json ../data/autotyp/langs.json ../data/autotyp/langs.filled.tsv ../data/autotyp/flist.json ../data/autotyp/langs.filled.json
```
 
## Run the model

- The hyperparameter settings must be changed properly.

```sh
python train_mda.py --seed=10 --K=100 --iter=1000 --bias --hmc_epsilon=0.025 --maxanneal=100 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --output ../data/wals/mda_K100.pkl ../data/wals/langs.filled.json ../data/wals/flist.json
```

```sh
python train_mda.py --seed=10 --K=50 --iter=1000 --bias --maxanneal=100 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --output ../data/autotyp/mda_K50.pkl ../data/autotyp/langs.filled.json ../data/autotyp/flist.json
```

- Collect samples

```sh
python sample_auto.py --seed=10 --a_repeat=5 --iter=100 ../data/wals/mda_K100.pkl.final - | bzip2 -c > ../data/wals/mda_K100.xz.json.bz2
python convert_auto_xz.py --burnin=0 --update --input=../data/wals/mda_K100.xz.json.bz2 ../data/wals/langs.filled.json ../data/wals/flist.json > ../data/wals/mda_K100.xz.merged.json
```

```sh
python sample_auto.py --seed=10 --a_repeat=5 --iter=100 ../data/autotyp/mda_K50.pkl.final - | bzip2 -c > ../data/autotyp/mda_K50.xz.json.bz2 &
python convert_auto_xz.py --burnin=0 --update --input=../data/autotyp/mda_K50.xz.json.bz2 ../data/autotyp/langs.filled.json ../data/autotyp/flist.json > ../data/autotyp/mda_K50.xz.merged.json
```

## Evaluation of missing data imputation

```sh
make -j 20 -f eval_mv.make DATATYPE=wals CV=10 MODEL_PREFIX=mda TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --hmc_epsilon=0.025 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if" mda
make -j 20 -f eval_mv.make DATATYPE=wals CV=10 MODEL_PREFIX=mda_dv TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --hmc_epsilon=0.025 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --drop_vs" mda
make -j 20 -f eval_mv.make DATATYPE=wals CV=10 MODEL_PREFIX=mda_dh TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --hmc_epsilon=0.025 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --drop_hs" mda
make -j 20 -f eval_mv.make DATATYPE=wals CV=10 MODEL_PREFIX=mda_oa TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --hmc_epsilon=0.025 --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --only_alphas" mda
make -j 100 -f eval_mv.make al DATATYPE=wals CV=10
```

```sh
make -j 20 -f eval_mv.make DATATYPE=autotyp CV=10 MODEL_PREFIX=mda TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --norm_sigma=10.0 --gamma_scale=1.0 --resume_if" mda
make -j 20 -f eval_mv.make DATATYPE=autotyp CV=10 MODEL_PREFIX=mda_dv TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --drop_vs" mda
make -j 20 -f eval_mv.make DATATYPE=autotyp CV=10 MODEL_PREFIX=mda_dh TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --drop_hs" mda
make -j 20 -f eval_mv.make DATATYPE=autotyp CV=10 MODEL_PREFIX=mda_oa TRAIN_OPTS="--maxanneal=100 --iter=500 --bias --norm_sigma=10.0 --gamma_scale=1.0 --resume_if --only_alphas" mda
make -j 100 -f eval_mv.make al DATATYPE=autotyp CV=10
```


# Bayesian Analysis of Correlated Evolution Involving Multiple Discrete Features

## About

Yugo Murawaki. [Analyzing Correlated Evolution of Multiple Features Using Latent Representations](http://aclweb.org/anthology/D18-1468). In Proceedings of the 2018 Conference on Empirical Methods in Natural Language Processing (EMNLP2018), pp. 4371-4382, Brussels, Belgium, 2018.

## Preprocessing

- convert Glottolog trees

```sh
python newick_tree.py ../data/glottolog/tree_glottolog_newick.txt ../data/glottolog/trees_all.pkl
```

- combine WALS languages and Glottolog trees

```sh
python merge_glottolog.py --npriors ../data/node_priors.json  ../data/wals/langs.json ../data/glottolog/trees_all.pkl ../data/wals/trees_attached.pkl
```

## Train the model

- the main inference

```sh
nice -19 python train_bin_ctmc.py --has_bias --resume_if --seed=0 --npriors ../data/node_priors.json ../data/wals/trees_attached.pkl ../data/wals/mda_K100.0.xz.merged.json ../data/wals/paramevo_K100.0.tree.pkl 2>&1 | tee -a ../data/wals/paramevo_K100.0.tree.log
```

- collect samples

```sh
nice -19 python train_bin_ctmc.py --iter=1100 --save_interval=10 --has_bias --resume ../data/wals/paramevo_K100.0.tree.pkl.final --seed=0 --npriors ../data/node_priors.json ../data/wals/trees_attached.pkl ../data/wals/mda_K100.0.xz.merged.json ../data/wals/paramevo_K100.0.tree_plus.pkl 2>&1 | tee -a ../data/wals/paramevo_K100.0.tree_plus.log
```

- estimate CTMC parameters for the surface feature

```sh
nice -19 python train_surface_ctmc.py --seed=0 ../data/wals/paramevo_K100.0.tree.pkl.final ../data/wals/flist.json ../data/wals/mda_K100.0.xz.merged.json ../data/wals/paramevo_K100.0.surface_tree.pkl 2>&1 | tee ../data/wals/paramevo_K100.0.surface_tree.log
```

(TODO) clean up a Jupyter Notebook (used for further analysis) and add it to the repository
