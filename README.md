# Optimizing Linear Correctors: A Tight Output Min-Entropy Bound and Selection Technique

This repository contains Python code and data for the paper: *Optimizing Linear Correctors: A Tight Output Min-Entropy Bound and Selection Technique*.

## Code

The directory `Code/` contains necessary Python code to reproduce figures from the paper. All figures can be generated directly from [main.py](Code/main.py). The content of [parameters.py](Code/parameters.py) can be modified to speed up the code execution by lowering `Decimal` context precision, setting another output min-entropy lower bound or to modify the parameters of the bisection method. 

## Files and data description

1. [NBC_list_and_weight_distributions](NBC_list_and_weight_distributions.zip): List with parameters and weight distributions of binary linear codes that are used to calculate the new bound of the corresponding correctors (NBC set)
2. [NBCCYC_list_and_weight_distributions](NBCCYC_list_and_weight_distributions.txt): List with parameters and weight distributions of cyclic binary linear codes that are used to calculate the new bound of the corresponding correctors (NBCCYC set)
3. [OBC_list_and_parameters](OBC_list_and_parameters.txt): List with parameters of binary linear codes that are used to calculate the old bound of the corresponding correctors (OBC set)
4. [OBCCYC_list_and_parameters](OBCCYC_list_and_parameters.txt): List with parameters of cyclic binary linear codes that are used to calculate the old bound of the corresponding correctors (OBCCYC set)
5. [NBC_H_out_0.999_PF_correctors](NBC_H_out_0.999_PF_correctors.txt): Optimal correctors for output min-entropy rate of 0.999 according to the new bound
6. [NBCCYC_H_out_0.999_PF_correctors](NBCCYC_H_out_0.999_PF_correctors.txt): Optimal correctors for output min-entropy rate of 0.999 based only on the cyclic codes according to the new bound
7. [OBC_H_out_0.999_PF_correctors](OBC_H_out_0.999_PF_correctors.txt): Optimal correctors for output min-entropy rate of 0.999 according to the old bound
8. [OBCCYC_H_out_0.999_PF_correctors](OBCCYC_H_out_0.999_PF_correctors.txt): Optimal correctors for output min-entropy rate of 0.999 based only on the cyclic codes according to the old bound
9. [Modified_generator_matrices_OBC_NBC](Modified_generator_matrices_OBC_NBC.zip): Modified generator matrices of the best known linear codes that originally contained one or multiple all-zero columns (BKLCs from  M. Grassl, *''Bounds on the minimum distance of linear codes and quantum codes,''* Online available at: http://www.codetables.de) 

	
Files [1](NBC_list_and_weight_distributions.zip) and [2](NBCCYC_list_and_weight_distributions.txt) contain in each line:
- parameters of the code on which the corresponding corrector is based (`n`, `k`, `d`), 
- reference where the complete code description and/or its weight distribution can be found (`source`),
- sequence of tuples which represents code's weight distribution (`Weight Distribution`), where the *i*-th tuple `<w_i, a_i>` represents number of codewords `a_i` with weight `w_i`,
- in cases where MAGMA uses PRNG we additionally provide `MAGMA Seed` and `MAGMA Seed iteration` values in brackets that were used for our constructions.

Files [3](OBC_list_and_parameters.txt) and [4](OBCCYC_list_and_parameters.txt) contain in each line:
- parameters of the code on which the corresponding corrector is based (`n`, `k`, `d`), 
- reference where the complete code description can be found (`source`).

Files [5](NBC_H_out_0.999_PF_correctors.txt), [6](NBCCYC_H_out_0.999_PF_correctors.txt), [7](OBC_H_out_0.999_PF_correctors.txt) and [8](OBCCYC_H_out_0.999_PF_correctors.txt) contain in each line:
- parameters of the code on which the corresponding corrector is based (`n`, `k`, `d`), 
- minimum required min-entropy of the input raw bits to achieve the output min-entropy rate of 0.999 (`H_in_req`),
- code rate = throughput reduction (`CR`),
- efficiency of the corrector at H_in_req (`efficiency`).

File [9](Modified_generator_matrices_OBC_NBC.zip) contains:
- parameters of the code on which the corresponding corrector is based (`n`, `k`, `d`), 
- in cases where MAGMA uses PRNG we additionally provide `MAGMA Seed` and `MAGMA Seed iteration` values in brackets that were used for our constructions,
- followed by `k` rows with `n` entries of the code's generator matrix.

## References

All constructions are based on binary linear codes whose descriptions can be found here:
- M. Grassl, *''Bounds on the minimum distance of linear codes and quantum codes,''* Online available at: http://www.codetables.de.
- M. Terada, J. Asatani, and T. Koumoto, *''Weight Distribution,''* Online available at: https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/index.html.
- N. J. Sloane, *''List of weight distributions in the on-line encyclopedia of integer sequences,''* Online available at: https://oeis.org/wiki/List_of_weight_distributions.
- S. Lin and D. J. Costello, *''Error Control Coding: Fundamentals and Applications,''* Pearson-Prentice Hall, 2004.
- D. Schomaker and M. Wirtz, *''On binary cyclic codes of odd lengths from 101 to 127,''* IEEE Transactions on Information Theory, vol. 38, no. 2, pp. 516–518, 1992.
- T. Sugita, T. Kasami, and T. Fujiwara, *''The weight distribution of the third-order Reed-Muller code of length 512,''* IEEE Transactions on Information Theory, vol. 42, no. 5, pp. 1622–1625, Sep. 1996.
- Y. Desaki, T. Fujiwara, and T. Kasami, *''The weight distributions of extended binary primitive BCH codes of length 128,''* IEEE Transactions on Information Theory, vol. 43, no. 4, pp. 1364–1371, 1997.
- T. Fujiwara and T. Kasami, *''The weight distribution of (256, k) extended binary primitive BCH code with k<= 63, k>= 207,''* IEICE, IT97, Technical Report, 1993.
- T.-K. Truong, Y. Chang, and C.-D. Lee, *''The weight distributions of some binary quadratic residue codes,''* IEEE Transactions on Information Theory, vol. 51, no. 5, pp. 1776–1782, 2005.
- M. Tomlinson, C. J. Tjhai, M. A. Ambroze, M. Ahmed, and M. Jibril, *''Error-Correction Coding and Decoding: Bounds, Codes, Decoders, Analysis and Applications,''* Springer Nature, 2017.


## Further Information

For more details on the new bound and selection procedure, please consult the paper or contact the authors: milos.grujic@esat.kuleuven.be, ingrid.verbauwhede@esat.kuleuven.be
