### Data
- raw_data/lung_phyloseq.Rdata: lotus2(v2.23, uparse)
- raw_data/gut_phyloseq.Rdata: lotus2(v2.23, uparse)
- raw_data/lung_mebo_pre.csv: lung metabolome profiles
- raw_data/lung_mebo_pre.csv: plasma metabolome profiles


### Command
### set cpu number
edit  0.make_targets_config.R at first, then
```shell
Rscript make_targets_config.R
```
#### run pipeline
```shell
Rscript targets.R
```