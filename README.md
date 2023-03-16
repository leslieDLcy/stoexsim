# stoexsim

A stochastic implementation of EXSIM

## Introduction

This repo is a wrapper for the implementation of the EXSIM (stochasict finite-fault model). Particularly, we provide extra capabilities in modeling some key parameters as random variables to account for variability.

> A Python implementation can be found in another repo.

## Usage

```shell
$ "mechanism=N;depth=10;Mw=6.5;Repi=10; csh do_exsim.csh $mechanism $depth $Mw $Repi"
```

- [x] A CLI interface facilitating the use of stochastic finite fault model;
- [x] Aleatoric uncertainty on the region-specific parameters

## References

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`stoexsim` was created by Y. Chen. It is licensed under the terms of the MIT license.

<!-- TODO: Tweak the standart out look of the CLI -->
