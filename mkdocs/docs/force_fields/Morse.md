# Morse potential


Potential form of the Morse potential is,

$$\begin{equation*}
E_{ij}(R_{ij}) = D_{ij} \left\{ [\exp (\alpha_{ij} (R_{0,ij} -R_{ij}))-1]^2 -1\right\}
\end{equation*}$$

The potential requires a parameter file `in.params.Morse` that has the
following format.

    #  isp, jsp, D_ij, alpha_ij, R0_ij
    1    1    3.7    2.0    1.68
    1    2    3.2    1.9    2.5
    2    2    2.5    2.1    2.3

In the above case, there are three interactions between 1-1, 1-2 and
2-2. Each interaction has three parameters:

-   *D_ij*: depth of the potential curve (eV)
-   *alpha_ij*: related to the width of the potential well
    (Ang.^{-1}).
-   *R0_ij*: position of the potential minimum (Ang.)
