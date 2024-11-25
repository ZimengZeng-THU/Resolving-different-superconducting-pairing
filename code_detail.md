The paper[1] distills the core idea of our idea, but for the convenience of numerical calculations and the beauty of the results, our program has some additional steps. Reading this note can help you understand every detail of our program. I will introduce in the following order: 
1. How the interaction matrix is calculated in our code.
2. How do we diagonalize this matrix and the meaning of its eigenvalues. 
3. Our agreement on the gap function $\Delta_{\mathbf{n}k}$.
4. How to eliminate solutions that have no physical meaning. 
5. Introduce how can be determined what irrep the gap function belong . 
6. Show how we drew an image of the gap function.

## interaction matrix calculation(cal_Delta.py CalVkk)
In this section, I will show you how to calculate the interaction matrix by our code.
It is mentioned in the paper that we can calculate the eigenvectors of the interaction matrix $V^S$ and the triplet matrix $V^T$ of singlet respectively to obtain the superconducting energy gap function. The basis of $V^S$ and $V^T$ in the text is
$$
\begin{aligned}
\hat\Delta^S(n\mathbf{k})=&\frac{1}{2}[\hat\Delta(n\mathbf{k}+-)-\hat\Delta(n\mathbf{k}-+)],\\
\hat\Delta^{T_0}(n\mathbf{k})=&\frac{1}{2}[\hat\Delta(n\mathbf{k}+-)+\hat\Delta(n\mathbf{k}-+)] \nonumber \\
\hat\Delta^{T_{\pm}}(n\mathbf{k})=&\hat\Delta(n\mathbf{k}\pm\pm)
\end{aligned}
$$
here
$$
\begin{aligned}
   \hat \Delta(n\mathbf{k}s_1s_2)=\tilde{c}_{n,\mathbf{-k},s_2}c_{n\mathbf{k}s_1},  \\
   \mathcal{T}c_{n\mathbf{k}s}\equiv s\tilde{c}_{n,-\mathbf{k},-s},    
\end{aligned}
$$
Defining spin metrics in Kohn Sham orbital is difficult when considering SOC effects.


For the convenience of numerical implementation, I define
$$
\begin{aligned}
\mathcal{T}c_{n\mathbf{k}s}\equiv \tilde{c}_{n,-\mathbf{k},-s},    
\end{aligned}
$$
here
$$
\begin{aligned}
\hat \Delta(n\mathbf{k}s_1s_2)=\tilde{c}_{n,\mathbf{-k},s_2}c_{n\mathbf{k}s_1} \\
\end{aligned}
$$
write the interaction matrix by this basis
$$
\begin{aligned}
   V_{s_2,s_1;s_3,s_4}(n\mathbf{k},m\mathbf{k}')&=&\sum_\upsilon \frac{g_{n\mathbf{k}s_1,m\mathbf{k'}s_3}^{\upsilon}(g_{n\mathbf{k},-s_2,m\mathbf{k'},-s_4}^{\upsilon})^*}{\omega_{\mathbf{k'}-\mathbf{k},\upsilon}}  \nonumber
N_{\mathbf{n}k}
\end{aligned}
$$
Here, k and k' can be valued in the range of those that satisfy the $abs(E_{\mathbf{n}k}-E_f)<fthick$ so $delta_{\mathbf{n}k}$[[(https://docs.epw-code.org/doc/Electron-phononCoupling.html)]] in the body does not appear here.

In the interaction matrix element, the elestron phonon coupling matrix element $g_{\mathbf{n}k,\mathbf{m}k'}^{\nu}$is the eletron phonon coupling matrix element calculated by EPW.According to the EPW Code's method[(https://docs.epw-code.org/doc/Electron-phononCoupling.html)], the density of states $N_{\mathbf{n}k}$ of $|psi_{\mathbf{n}k}$ can be simulated by the Gaussian distribution function. It is worth mentioning that the elestron phonon coupling of EPW is based on the formula in the link, where the improved Gaussian distribution function is used to simulate the gap function (equivalent to the density of states at k points[]), and my calculation uses the Gaussian distribution function output by EPW.
## Diagonalizing Interaction Matrix(cal_Delta.py CalEigen NormVal)

As mentioned in the text, the eigenvalue obtained by the diagonal interaction matrix is the superconducting intensity of the corresponding eigenvector (pairing channles). It should be noted that different Gaussian distributions ($w0gauss(\sigma(\epsilon_{\mathbf{n}k}),0)$, $w0gauss(\sigma(\epsilon_{\mathbf{n}k}),1)$) are applied to $Delta_{mathbf{n}k,dos}$ in the density of states calculation and $Delta_{mathbf{n}k,epc}$ in the elestron phonon coupling calculation. To align our calculations, I'll do the following
$$
\begin{aligned}
\lambda_i -> \lambda_i \frac{N_{epc}}{N_{dos}} \\
N_{epc} = \sum\frac{w0g(\mathbf{n}k)}{N_k} \\
N_{epc} = \sum\frac{w1g(\mathbf{n}k)}{N_k}
\end{aligned}
$$
With this processing, I ensure that when the electron phonon coupling matrix element is constant, the intensity of our leading channel is equal to the total electron phonon coupling strength calculated by EPW.


## notation of gap function(irrep_Delta.py CheckSolution)
We write the eigenvectors of the interaction matrix in the form of a 2 * 2 matrix at each k point

For the solution with even parity，$\Delta_{\mathbf{n}k}$proportional to $\sigma_0$;for the solution with odd parity，$\Delta_{\mathbf{n}k}$proportional to$\alpha_1\sigma_1+\alpha_2\sigma_2+\alpha_3\sigma_3$。For ease of presentation, we would like to define $\ Delta_ {\ mathbf {n}k }$ as a single valued function about nk, and the following convention is made for this purpose
For even parity:
$$
\begin{aligned}
\Delta_{\mathbf{n}k} = \Delta_{\mathbf{n}k,00}
\end{aligned}
$$
For odd parity:
$\Delta_{\mathbf{n}k}$is Hermitian，We can always find the unitary matrix $G_k $such that $G_k\Delta_{\mathbf{n}k}G_k^{+}$is a diagonal matrix
$$
\begin{aligned}
\Delta_{\mathbf{n}k} -> (G_k\Delta_{\mathbf{n}k}G_k^{+})_{00}
\end{aligned}
$$
After the transformation of $G_k$,  basis of $G_k\Delta_{\mathbf{n}k}G_k^{+}$ is different from $\Delta_{\mathbf{n}k}$. To record the basis of the gap function, we also perform the same transformation on the electron wave function and record its spin expectation value.
$$
\begin{aligned}
\Psi_{\mathbf{n}k} = G_k\phi_{\mathbf{n}k}
\end{aligned}
$$
here $\phi_{\mathbf{n}k}$is a matrix of 2 * n, with n numbers in each row representing the wave function distribution under the Wannier basis, and two columns corresponding to two Kramers degenerate states.
Please note that at this time, our basis $\Psi_ {\mathbf{n}k}$and $\Psi_{\mathbf{n}-k}$ .The relationship between is unknown, and due to inversion symmetry, their relationship have only two possibilities:
(1)
$$
\begin{aligned}
\Psi_{\mathbf{n}k} = \hat{I}\Psi_{\mathbf{n}-k}
\end{aligned}
$$
in this condition, $\Delta_{\mathbf{n}k} = -\Delta_{\mathbf{n}-k}$


(2)
$$
\begin{aligned}
\Psi_{\mathbf{n}k} = \hat{T}\Psi_{\mathbf{n}-k}
\end{aligned}
$$
in this condition, $\Delta_{\mathbf{n}k} = \Delta_{\mathbf{n}-k}$


For the aesthetic presentation of numerical values, we hope that regardless of the basis, the gap function satisfies
$$
\Delta_{\mathbf{n}k} = \Delta_{\mathbf{n}-k}
$$
So we add the following agreement: for $\Delta_{\mathbf{n}k }$, if the x-component of spin of its basis ($\Psi_{\mathbf{n}k0}$) is less than zero, let $\Delta_{\mathbf{n}k }-> -\Delta_{\mathbf{n}k}$. Even in the case of (2), the gap function we agreed upon is still odd parity.


## Eliminate redundant solutions(irrep_Delta.py CheckSolution)
The eigenvector of the interaction matrix $\hat{V}$ can be considered as representing the ratio between the order parameter  $<\tilde c_k{c_{-k}}>$ for different energy bands and momenta. In this note, I denote the element of the eigrnvector by $<\tilde c_k {c_{-k}}>$


It must be emphasized that the eigenvector $\Delta$ contains order parameters linked by commutative relationships($<\tilde c_k {c_{-k}}>,<\tilde c_{-k} {c_{k}}>$). Obviously, they should be equal(These two order parameters undergo time difference inversion and fermion pair inversion, with each operation providing a negative sign). But our eigenvectors do not necessarily have this feature. We need to manually exclude those eigenvectors that do not satisfy this feature, as they have no physical meaning. Therefore, according to the agreement in the previous section, we can exclude solutions that have no physical meaning
对于even parity:
$$
\Delta_{\mathbf{n},-k} = -\Delta_{\mathbf{n},k}
$$
When $\ Delta_{\ mathbf {n},-k} $in matrix format is proportional to $\sigma_0 $, the pairing should be even parity, but the calculated result is an odd parity, which is obviously unreasonable.


For odd parity, in last section notation:
$$
\begin{aligned}
\Delta_{\mathbf{n},k}  = <\tilde c_{-k,0} {c_{k,0}}>
\Delta_{\mathbf{n},-k} = -<\tilde c_{k,0} {c_{-k,0}}>
\end{aligned}
$$
if the results shows =
$$
\Delta_{\mathbf{n},-k} = \Delta_{\mathbf{n},k}
$$
It has no physical meaning

## Confirm which irreducible representation the gap function belongs to(irrep_Delta.py CheckRepC2x)
In order to determine the irreducible representation of a certain gap function, we need to apply a symmetry operation to this gap function, check the transformation of this gap function under the symmetry operation, and refer to the characteristic label table to determine which irreducible representation basis function this gap function belongs to.

Specifically, for the gap function $\Delta_{\mathbf{n},k} $, its basis is $c_{\Psi_{\mathbf{n}k0}}\tilde{c}_{\Psi_{\mathbf{n}k0}}$. To determine its features, it is necessary to compare $G(\Delta_{\mathbf{n}k})$and $\Delta_{\mathbf{n}k}$.Consider a simple scenario: one-dimensional representation. At this point, $\alpha=\pm1 $.Known:

$$
\begin{aligned}
\Delta_{\mathbf{n}G(k)}=\beta\Delta_{\mathbf{n}k} \\
G(\Delta_{\mathbf{n}k}) =& G(c_{\Psi_{\mathbf{n}k0}})G(\tilde{c}_{\Psi_{\mathbf{n}k0}}) \\
=& c_{G(\Psi_{\mathbf{n}k0})}\tilde{c}_{G(\Psi_{\mathbf{n}k0})} \\
=& \alpha\Delta_{\mathbf{n}G(k)} =& \alpha\beta\Delta_{\mathbf{n}k}
\end{aligned}
$$
考虑简单的情况：一维表示。此时$\alpha=\pm1$。

There are only two possibilities for G to transform the basis
(1)
$$
G(\Psi_{\mathbf{n}k0})=\Psi_{\mathbf{n}G(k)0}
$$
In this condition $\alpha=1$

(2)
$$
G(\Psi_{\mathbf{n}k0})=\Psi_{\mathbf{n}G(k)1}
$$
In this condition $alpha=-1$

$\alpha\beta$ is charater of this opration.


In the CheckSym function, we studied the characteristics of Pb's symmetric operation C2x in the leading odd parity channel. The transformation matrix is obtained from the output file scf.out of QE. By running this function and observing its output, it can be observed that for any k, $\ alpha \ beta=1 $. This indicates that the characteristic label of C2x is -1, and this channel is the basis function represented by A1u

## interpolation and plot

The drawing software we use is fermisuffer, which takes electronic structure data E(ibnd,kx,ky,kx) and another set of functions c(ibnd,kx,ky,kz) as inputs. We use the electronic structure data to draw the Fermi surface, and then use c(ibnd,kx,ky,kz) as its color.

In fact, due to the limited grid density of electronic structural data, this software automatically interpolates on denser k-grids during drawing to achieve smooth patterns (the specific interpolation method is uncertain). So this has led to difficulties in our plotting, as the gap function we calculated is only non-zero when $|E(ibnd, kx, ky, kz) - E_f|<fthick$, and the software's automatic interpolation results in many gap functions at k points being 0, which is clearly unreasonable.

To solve the problem of numerical presentation, I manually interpolated the gap function and extrapolated it to the entire energy band to ensure that there were no issues when plotting. It should be emphasized here that our manual extrapolation does not change the energy gap function near the Fermi surface. This is just a method to make the pattern more beautiful and does not affect the accuracy of the data.

After completing the above steps, I will write our data to the. frmsf file in the format of the input file from fermisuffer, which needs to be opened by fermisuffer.