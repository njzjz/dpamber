# dpamber

Some useful tools related to Amber and DP.

## Installation

```sh
pip install dpamber
```

## Tools
### corr: generating data for DPRc models

[![DOI:10.1021/acs.jctc.1c00201](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.1c00201-blue)](https://doi.org/10.1021/acs.jctc.1c00201)
[![Citations](https://citations.njzjz.win/10.1021/acs.jctc.1c00201)](https://doi.org/10.1021/acs.jctc.1c00201)

`corr` tool generates [DeePMD-kit](https://github.com/deepmodeling/deepmd-kit) training data for DPRc from AMBER sander low-level QM/MM data and high-level data. For details of DPRc, read the [DPRc paper](https://doi.org/10.1021/acs.jctc.1c00201).

Before using this tool, one need to prepare low-level and high-level QM/MM data:

$$
E_\text{hl}(\mathbf R)=E_\text{hl,QM}(\mathbf R)+E_\text{hl,QM/MM}(\mathbf R)+E_\text{MM}(\mathbf R)
$$

$$
E_\text{ll}(\mathbf R)=E_\text{ll,QM}(\mathbf R)+E_\text{ll,QM/MM}(\mathbf R)+E_\text{MM}(\mathbf R)
$$

Low-level and high-level data should use the same coordinate and the same MM method, but different QM methods. So, the correction energy for training will be

$$
\Delta E (\mathbf R) = E_\text{hl}(\mathbf R) - E_\text{ll}(\mathbf R) = (E_\text{hl,QM}(\mathbf R) - E_\text{ll,QM}(\mathbf R)) + (E_\text{hl,QM/MM}(\mathbf R) - E_\text{ll,QM/MM}(\mathbf R))
$$

An example of the command is
```sh
dpamber corr --cutoff 6. --qm_region ":1" --parm7_file some_param.param7 --nc some_coord.nc --hl high_level --ll low_level --out dataset
```
where `--cutoff` takes cutoff radius of the QM/MM interaction for training. `--qm_region` takes AMBER mask format for the QM region. `--parm7_file` and `--nc` take the PARM7 file and the trajectory (NetCDF) file, respectively. `--ll` and `--hl` are the prefixes of low-level and high-level files, including the mdout file (`.mdout`), the mden file (`.mden`) and the mdfrc file (`.mdfrc`). The output dataset directory should be put in `--out`.

See details from `dpamber corr -h`.

### devi: calculate model deviation

`devi` can be used to calculate the model deviation of a given trajectory.
You need to install DeePMD-kit using
```sh
pip install dpamber[dpgpu]
```

See `dpamber devi -h` for details.
