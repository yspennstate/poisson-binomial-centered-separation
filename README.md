# A lower bound for the Kolmogorov distance between centered Poisson binomial pairs

Yitzchak Shmalo, 16 November 2026.

Let X and Y be sums of n independent Bernoulli variables with parameters p_i >= q_i, and let g be the mean of the gaps p_i - q_i. Aryeh Kontorovich asked on MathOverflow (question 477402) whether, for even n, the Kolmogorov distance between X and Y is at least the distance between Bin(n, (1+g)/2) and Bin(n, (1-g)/2). The paper proves this in two centered cases, where sum_i (p_i + q_i) = n:

- all gaps equal, any even n >= 4 (Theorem 1.1): both sums are built from one sum of independent Bernoulli variables with integer mean, and Hoeffding's interval comparison for that sum writes the excess over the binomial value as a sum of nonnegative terms, bounded below by a positive multiple of sum_i (p_i + q_i - 1)^2;
- four summands with arbitrary gaps (Theorem 1.2), with an exact formula for the excess of max(D_2, D_3) over b_4(g) and a sharp constant 8/3.

When the gaps are not all zero, equality holds in both cases only for the symmetric binomial pair. The inequality itself, for every pair and without centering, and its equality cases are due to C. D. Long (github.com/long-mathematics/poisson-binomial-comparison); the paper gives the excess over the binomial value.

## Files

- `main.tex`, `main.pdf`: the paper.
- `CREDITS.md`: which statements are known and which are new, and whose work they rest on.
- `code/certificate_checks.py`: verifies in exact rational arithmetic the polynomial identities of Section 3 (the four-summand averages and differences, the coordinates of Q, the affine form of L, the regrouping, the gap identity and the vertex expansion in Lemma 3.1). Output: `code/certificate_checks.json`.
- `code/general_uniform.py`: checks the identities of Section 2 (the average, the jumps, the exact excess, the variance identity) by exact convolution on 1,335 rational parameter vectors with n <= 20, and the bound of Theorem 1.1 on the 1,280 of them with n >= 4. Output: `code/uniform_checks.json`.
- `code/print_formulas.py`: prints the low-dimensional coefficients displayed in Section 2 as exact rationals.
- `code/verify_paper_identities.py`: a second check of the same identities written in the notation of the paper, with SymPy, together with the equal-gap formulas for m = 2, ..., 6 in exact rationals.

`certificate_checks.py`, `general_uniform.py` and `print_formulas.py` use only the Python 3 standard library; `verify_paper_identities.py` needs SymPy. Each program stops with an error if a check fails.

## License

The paper is licensed CC BY 4.0 and the code MIT.
