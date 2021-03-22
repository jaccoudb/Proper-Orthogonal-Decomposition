## About

Proper Orthogonal Decomposition (POD) is a powerful method for low-order approximation of some high dimensional processes. It is widely used in the situations where model reduction is required. The most favorable feature of the method is its optimality: it provides the most efficient way of capturing the dominant components of high-dimensional processes with, sometimes surprisingly small number of "modes".

In this source code was implemented a derivation of the POD basis that will be presented here is the Singular Value Decomposition (SVD). SVD can be viewed as the extension of the eigenvalue decomposition for the case of non-square matrices. In general, this decomposition states, that for any real rectangular $N\times{M}$ matrix $\Sigma$, there exists orthogonal matrices, $N\times{N}$ matrix $V$ and $M\times{M}$ $Z$ such that
