# A Conic Quadratic Programming Approach to Physical Layer Multicasting for Large-Scale Antenna Arrays

This repo contains the code for implementing the iterative SOCP-based method presented in the following scientific paper:

Le-Nam Tran, Muhammad-Fainan Hanif, and Markku Juntti, "A Conic Quadratic Programming Approach to Physical Layer Multicasting for Large-Scale Antenna Arrays," IEEE Signal Process. Lett., vol. 21, no. 1, pp. 114-117, Jan. 2014. 

## Instructions
The Matlab code makes use of [Yalmip](https://yalmip.github.io/) as a parser and [MOSEK](https://www.mosek.com/) as the internal convex conic solver for speed.

Recently we have found that [MOSEK Fusion API](https://docs.mosek.com/latest/pythonfusion/index.html) is particularly useful if one wishes to solve a same convex problem for different problem data. In this regard, we also provide a Python script which is much more efficient if you want to run the algorithm for many channel realizations, compared to the Matlab implementation.

We also publish the codes on [Code Ocean](https://codeocean.com/capsule/9830614/tree/v1) where you can run them online without the need to install the required packages.
