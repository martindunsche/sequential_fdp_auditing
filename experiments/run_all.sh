#!/bin/bash

for mu in 5 7 10; do
    ./dpsgd.R "$mu"
done

for mu in 0.5 0.8 1.0; do
    ./gaussian.R "$mu"
done

for mu in 0.01 0.1 0.5; do
    ./gaussian_comp_other_paper.R "$mu"
done

for mu in 0.5 0.8 1.0; do
    ./laplace.R "$mu"
done

for mu in 0.01 0.1 0.5; do
    ./laplace_comp_other_paper.R "$mu"
done
