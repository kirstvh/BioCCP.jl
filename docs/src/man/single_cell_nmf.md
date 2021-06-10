# SingleCellNMF

This project implements coupled non-negative matrix factorization(NMF) for multiomics single-cell data analysis.

## Motivation

Single-cell technologies allow studying cellular heterogeneity at unprecedented resolution thousands of cells. It is possible to sequence transcriptome, epigenome, proteome, and other -omes of single cells, and the price is constantly dropping. As it is possible to profile both epigenetic features and the transcriptome, we can try to connect these features to uncover the mechanisms of regulatory heterogeneity. However, integrating the data from multiple experiments and modalities can be quite challenging due to technical and biological variability. Single-cell data analysis is full of challenges and most tools are relatively new. Notably, chromatin accessibility measurements at the single-cell level are extremely sparse, and it is not always possible to reliably recover all the cell types from such measurements.

## What is matrix factorization?

Matrix factorization (MF) has been been used to solve the challenges of single-cell data intergration. The main idea is to jointly factorize data matrices from multiple modalities to infer a "latent space" - shared low-dimensional representation of all modalities. Several authors proposed methods based on MF:

- [Coupled NMF](https://doi.org/10.1073/pnas.1805681115)
- [LIGER](https://doi.org/10.1016/j.cell.2019.05.006)
- [MOFA and MOFA+](https://doi.org/10.1186/s13059-020-02015-1)
- [scAI](https://doi.org/10.1186/s13059-020-1932-8)
- [DC3](https://doi.org/10.1038/s41467-019-12547-1)

Let us define the problem for one input matrix first. Given an input matrix $X$, we seek such matrices $W$ and $H$ that minimizes $|X - WH|^2_F$. The product $WH$ approximates $X$ in a low-dimensional space and captures relevant features of the data. NMF introduces additional constraint such what W and H have to be component-wise non-negative:

$\underset{W\geq0, H\geq0}{min} \Vert{X-WH}\rVert_{F}^{2}$

Additionally, NMF has an inherent clustering property, which is useful since clustering is part of typical single-cell analysis.

In case there are two data sets (e.g. transcriptome and chromatin accessibility), the input matrices can be simultaneously factorized to infer a shared latent space. Such a procedure allows capturing shared variability in both datasets. scAI authors proposed including additional terms into the objective to aggregate sparse epigenetic signal over the cells, which improves the clustering and infers more meaningful latent factors. The optimization objective is formulated as:

${\min}_{W_1,{W}_2,H,Z\ge 0}\alpha {\left\Vert {X}_1-{W}_1H\right\Vert}_F^2+{\left\Vert {X}_2\left(Z\circ R\right)-{W}_2H\right\Vert}_F^2+\lambda {\left\Vert Z-{H}^TH\right\Vert}_F^2+\gamma \sum \limits_j{\left\Vert {H}_{.j}\right\Vert}_1^2,$

Here $X_{1}$ is a gene expression matrix with cells in the columns and features in the rows. $X_{2}$ is ATAC matrix with cells in the columns and features in the rows. $Z$ is an aggregation matrix, $R$ is a matrix with elements drawn from Bernoulli distribution. Element-wise multiplication of $Z$ and $R$ prevents over-aggregation of epigenomic signal. Multiplying $X_{2}$ by the product of $Z$ and $R$ effectively aggregates sparse signal over similar cells. $\alpha$, $\lambda$, and $\gamma$ are regularization terms.

## Solving NMF

The coordinate descent approach is commonly used to solve the optimization problem. Therein, one of the factors is updated, while others are fixed. Below is the pseudocode for factorizing two input matrices:
```
Given matrices X1(p by n) and X2(q by n):
Randomly initialize W1(p by k), W2(q by k), H(k by n), Z(n by n)
for i = 1...number of iterations:
  W1(i) = update(X1, X2, H(i-1), W1(i-1), W2(i-1), Z(i-1))
  W2(i) = update(X1, X2, H(i-1), W1(i), W2(i-1), Z(i-1))
  H(i) = update(X1, X2, H(i-1), W1(i), W2(i), Z(i-1))
  Z(i) = update(X1, X2, H(i), W1(i), W2(i), Z(i-1))
end
```
The terms can be updated using multiplicative update rule, which is trivial to implement:

${W}_1^{ij}\leftarrow {W}_1^{ij}\frac{{\left({X}_1{H}^T\right)}^{ij}}{{\left({W}_1H{H}^T\right)}^{ij}}$
${W}_2^{ij}\leftarrow {W}_2^{ij}\frac{{\left({X}_2\left(Z\circ R\right){H}^T\right)}^{ij}}{{\left({W}_2H{H}^T\right)}^{ij}}$
${H}^{ij}\leftarrow {H}^{ij}\frac{{\left(\alpha {W}_1^T{X}_1+{W}_2^T{X}_2\left(Z\circ R\right)+\lambda H\left(Z+{Z}^T\right)\right)}^{ij}}{{\left(\left(\alpha {W}_1^T{W}_1+{W}_2^T{W}_2+2\lambda H{H}^T+\gamma e{e}^T\right)H\right)}^{ij}}$
${Z}^{ij}\leftarrow {Z}^{ij}\frac{{\left(\left({X}_2^T{W}_2H\right)\circ R+\lambda {H}^TH\right)}^{ij}}{{\left(\left({X}_2^T{X}_2\left(Z\circ R\right)\right)\circ R+\lambda Z\right)}^{ij}}$

## Choice of dimensionality

Gene and locus loading matrices $W_{1}$ and $W_{2}$ have dimensionality of $(N_{features}, k)$, and $H$ is $(N_{cells}, k)$. As we seek a low-dimensional representation of data, $k$ should be lower than number of features (genes and loci) and cells. We can use our knowledge about biological system at hand to approximately choose $k$, and to further improve the choice of $k$ by trying out different values.

```@docs
perform_nmf
reduce_dims_atac
```
