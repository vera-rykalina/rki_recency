# Recency Pipeline (HIV-phyloTSI)
This is a seamless pipeline developted to automate HIV-phyloTSI method. The [method](https://www.medrxiv.org/content/10.1101/2022.05.15.22275117v1) is meant to estimate time since infection (TSI) from HIV deep-sequencing data. Here is a link for the HIV-phyloTSI GitHub [page](https://github.com/BDI-pathogens/HIV-phyloTSI/tree/main).


## Tools
The pipeline includes the following tools:
- KRAKEN2
- SHIVER (C. Wymant)
- KALLISTO
- IVA
- PHYLOSCANNER (C. Wymant & M. Hall)
- IQTREE
- HIV-phyloTSI (T. Golubchik)

## Installation
Make sure that R with the phyloscannerR package and pyfastaq installed globally. All other software tools with their dependencies get installed automatically via conda directives. 

Download Kraken2 database and change the path accordingly in the shiver_phyloscanner_tsi_pipeline.nf if Kraken2 reports are desired.

>params.krakendb = "/scratch/databases/kraken2_20230314/"

The database used by the pipeline is PlusPFP (3 March, 2023) which is stored within HPC resources. 

Alternativelly, one can install all the tools manually, using .yml recipes (see Environments folder). If so, pyfastaq can be installed within a shiver's environment:

```sh
conda activate shiver
```

```sh
pip3 install --ignore-installed pyfastaq
```

## Usage
Populate RawData/ folder with *.fastq.gz files.

Change directory to the Pipeline/: 
```sh
cd Pipeline/
```

Create a nextflow environment, using nextflow.yml file and activate it:

```sh
conda env create -n nextflow --file nextflow.yml
conda activate nextflow
```

Run the pipeline: 

```sh
nextflow /Scripts/shiver_phyloscanner_tsi_pipeline.nf --outdir Results -c Script/rki_profile.config -profile rki_slurm,rki_mamba
```

## Flowchart
Below is a flowchart of the pipeline:
[![](https://mermaid.ink/img/pako:eNqNVm1vm0gQ_isRn1ypjWAXv-XDSTTGLRcHO4ZWl8PVagOLvXcYKKx7qZr-99vF7NgQX3X-5HnmfXbmET-MuEiYcWOkWfFPvKOVuApnm_xK_kpzMLjd0Txn2XVaFfsVFbs3b1qd1dXNecZWlFc1GKDB4AoEHN2tnTvXR-R24QSBN3_80mrsc7PhuTA6F8aR53uh5yy8wAm9pa_dJ9Gds5BguCSeP3P_0Ph0EMXF_onn7AtUbJ5sHz45fqhtLSt67wYhkcE_-PfumQJFcycIH4g3Ix9dZ-auQYMj77NDbpd-6H0IALWjJsYrfHihmtEFbHwBm0T3zgoiTaU01xIylX2WsViAPbJex0Ao-n3p-US6Qk0IR--de7J256q72-AzKGwIqh71FEQ2saflSZYNpBkVguUnbHyhINXAnUvCteuesk8722F2JKsjoddpMJYYzwSrTpAdeQ8qhc6AO7uEO8uEL5SJJ5FyJ47vLB4DDyrFnUrtTqV2p1K7s_E2jlYfHxdLEgYe7LodrdZuGMr1l3lmZLVYwq7ZnXrtTr3m1bt3v73EO1KxlOyen9CLeme9ElqZ0lp8JaU6QqnHetKNOnk6g7CC5H11xWFXtPXhHUWrI6t0POci4dW1-8xrwfPtmqX1p3xLy5IlMtnkv-y9POHfeHKgmfJ4Ueenb_mYavrrxrR6qtV_0yyTJRQyQ8Ke5Yxo0oQ14e7bFmBiVgsMf53JQsAD_9MBAz30pmbZLTAGSmiBUd8FgFHfBUYpz5DQatuk7M65HfK5wmrHasFcp0cAmUAjLWABh-hMatlkulQtHAI20dq_CskySau3cS8hghVsW0M2UIyO8ET3zVLzhMhXVFn6ZaMh8M9Z-3LHSHEQJJevrrYNjYCVWrcxUFILTICPjgBGfcDqA2YPQFMgpdYCTgrr4mjGt2oozRqSsqi54EVeE_Yc87q5DGwDZ8FtfBUVY0o37ugUcfUBuNO2MTwBBmuApoicZt9lNqLC1tclFbWggsT1t_OHah0kcfUBqw-YPQBPgeZaCxs47liFqLnKR-q4OOTiRTJaz6DMCqG1iv6Mt8aeVXvKE_kx8kMZbwyxY3u2MW7k34Sl9JCJjbHJf0pTehBF8D2PjRtRHdhb41AmVLAZp9uK7o2blGa1REua_1kUWv75LyBDf6Q?type=png)](https://mermaid.live/edit#pako:eNqNVm1vm0gQ_isRn1ypjWAXv-XDSTTGLRcHO4ZWl8PVagOLvXcYKKx7qZr-99vF7NgQX3X-5HnmfXbmET-MuEiYcWOkWfFPvKOVuApnm_xK_kpzMLjd0Txn2XVaFfsVFbs3b1qd1dXNecZWlFc1GKDB4AoEHN2tnTvXR-R24QSBN3_80mrsc7PhuTA6F8aR53uh5yy8wAm9pa_dJ9Gds5BguCSeP3P_0Ph0EMXF_onn7AtUbJ5sHz45fqhtLSt67wYhkcE_-PfumQJFcycIH4g3Ix9dZ-auQYMj77NDbpd-6H0IALWjJsYrfHihmtEFbHwBm0T3zgoiTaU01xIylX2WsViAPbJex0Ao-n3p-US6Qk0IR--de7J256q72-AzKGwIqh71FEQ2saflSZYNpBkVguUnbHyhINXAnUvCteuesk8722F2JKsjoddpMJYYzwSrTpAdeQ8qhc6AO7uEO8uEL5SJJ5FyJ47vLB4DDyrFnUrtTqV2p1K7s_E2jlYfHxdLEgYe7LodrdZuGMr1l3lmZLVYwq7ZnXrtTr3m1bt3v73EO1KxlOyen9CLeme9ElqZ0lp8JaU6QqnHetKNOnk6g7CC5H11xWFXtPXhHUWrI6t0POci4dW1-8xrwfPtmqX1p3xLy5IlMtnkv-y9POHfeHKgmfJ4Ueenb_mYavrrxrR6qtV_0yyTJRQyQ8Ke5Yxo0oQ14e7bFmBiVgsMf53JQsAD_9MBAz30pmbZLTAGSmiBUd8FgFHfBUYpz5DQatuk7M65HfK5wmrHasFcp0cAmUAjLWABh-hMatlkulQtHAI20dq_CskySau3cS8hghVsW0M2UIyO8ET3zVLzhMhXVFn6ZaMh8M9Z-3LHSHEQJJevrrYNjYCVWrcxUFILTICPjgBGfcDqA2YPQFMgpdYCTgrr4mjGt2oozRqSsqi54EVeE_Yc87q5DGwDZ8FtfBUVY0o37ugUcfUBuNO2MTwBBmuApoicZt9lNqLC1tclFbWggsT1t_OHah0kcfUBqw-YPQBPgeZaCxs47liFqLnKR-q4OOTiRTJaz6DMCqG1iv6Mt8aeVXvKE_kx8kMZbwyxY3u2MW7k34Sl9JCJjbHJf0pTehBF8D2PjRtRHdhb41AmVLAZp9uK7o2blGa1REua_1kUWv75LyBDf6Q)