# Density Functional Theory Library

Description

## Intallation

Download stuff

```bash
cd src && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../DFT-install
cmake --build .
cmake --install .
```

## Usage

Make C++ file and CMakeList as in example and run
```bash
mkdir build && cd build
cmake ..
cmake --build .
```

LSDA

Let the collected spin density be,\
$$ n = n_\uparrow + n_\downarrow, $$\
and the polarization denisty,\
$$ \zeta = \frac{n_\uparrow - n_\downarrow}{n}.$$\

Then the exchange energy will be\
$$E_{X}[n]=\frac{1}{2}\left[E^H[2n_\uparrow] + E^H[2n_downarrow]\right],$$\
with potential\
$$v^H[n]=(6\pi n_uparrow)^{1/3} + (6\pi n_downarrow)^{1/3}.$$\

The correlation energy density given by Vosko, Wilk, Nusair ([VWN](https://doi.org/10.1139/p80-159)),\
$$e_c^{VWN}(r_s,\zeta)=e_c(r_s,0)+\Delta e_c(r_s,\zeta).$$\
$$e_c(n,\{0,1\)&=n\frac{1-\text{ln}2}{\pi^2}\Bigg\{\text{ln}\frac{x^2}{X(x)}+\frac{2b}{Q}\text{arctan}\left(\frac{Q}{2x+b}\right)-\frac{bx_0}{X(x_0)}\left[\text{ln}\frac{(x-x_0)^2}{X(x)}+\frac{2(2x_0+b)}{Q}\text{arctan}\left(\frac{Q}{2x+b}\right)\right]\Bigg\}$$,\
Here $X(x)=x^2+bx+c$, and $Q=\sqrt{4c-b^2}$, with the Wigner-Seitz radius in 3D $\frac{4}{3}\pi r_s^3 n=1$. $x_0$, $b$, and $c$ are fit parameters, while $x=\sqrt{r_s}$.\
Then define $\Delta e_c(r_s,1)=e_c(r_s,0)-e_c(r_s,1)$,\
$$f(\zeta)=\frac{(1+\zeta)^{4/3}+(1-\zeta)^{4/3}-2}{2(2^{1/3}-1)},$$\
and\
$$\alpha_c=(r_s)=-\frac{\ln r_s}{3\pi^2}+C_\alpha$$.\

