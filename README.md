# BeamHardening-Correction
This repository allows to correct the beam-hardening effect in computerized tomography assuming to know the spectrum of the source and the type of materials.


## Beam-hardening

### Monochromatic source 

The variation in intensity for an incident beam of intensity $\Phi_{0,E}$ after crossing the material on length $dx$ is given by 
$$d\Phi_E = - f_{E} \Phi_{0,E} dx \% $$
with $f_{E}$ the energy dependant linear attenuation coefficient. For monochromatic radiation (energy $E_0$), this leads after integration to the well-known Beer-Lambert law which describes the variation of the intensity of a beam crossing a medium of attenuation $f_{E_0}$ along the straight lines $L$
$$\Phi(L) = \Phi_0 e^{-\int_L f_{E_0}(x)dx}. \%$$ 
Denoting by $L(\mathbf{x},\theta) = \{\mathbf{x}+t \theta, t\in\mathbb{R}^+\}$ the straight line starting at $\mathbf{x}\in \mathbb{R}^3$ with direction $\theta \in S^2$, the measure of the attenuated beam for monochromatic radiation delivers the well-known Ray transform $\mathcal{P}$
$$\ln\frac{\Phi_0}{\Phi(\mathbf{x},\theta)} = \int_{0}^\infty f_{E_0}(\mathbf{x}+t \theta)dt =: \mathcal{P}f_{E_0}(\mathbf{x},\theta).\%$$

### Polychromatic source

When the radiation is polychromatic, the variation in intensity for an incident beam of intensity $\Phi_{0,E}$ after crossing the material on length $dx$ is given by 
$$d\Phi_E = - f_{E} \Phi_{0,E} dx \%$$
with $f_{E}$ the energy dependent linear attenuation coefficient. For polychromatic radiation (energy domain $\mathbb{E}=[E_m,E_M]$), the linearity of the model does not hold anymore due to the energy dependency of the attenuation. Integrating over the energy, we obtain the beam-hardening formula,
$$\frac{\Phi(\mathbf{x},\theta)}{\Phi_0} = \int_\mathbb{E} T(E) e^{-\int_{0}^\infty f_{E_0}(\mathbf{x}+t \theta)dt} dE = \int_\mathbb{E} T(E) e^{-\mathcal{P} f_{E}(\mathbf{x},\theta)} dE \%$$
with $T(E)$ the normalized spectrum taking into account the response of the detector. As a consequence, the modeling based on the Ray/Radon transform stops to hold for a polychromatic source. 


From Stonestrom et al., the dependency regarding the energy of the attenuation coefficient can be approximated by
$$ f_E(\mathbf{x}) = E^{-3} f_1(\mathbf{x}) + C(E) f_2(\mathbf{x}) \% $$
with $C(E)$ the Klein-Nishina function and $f_1,f_2$ two characteristic maps  of the attenuating medium. The first part describes the effect of photoelectric absorption whereas the second one gives the attenuation due to Compton scattering.


### A model for correction 


In practice, $C(E)$ varies only a little due to the size and position of the X-ray energy range. This is why, to simplify the algebra, we consider in our study $C(E)$ to be a constant $a$. This approximation is crucial because it enables to rewrite the beam-hardening phenomenon as
$$B(\mathcal{P}f_1,\mathcal{P}f_2)(\mathbf{x},\theta) := e^{-a \mathcal{P}f_2(\mathbf{x},\theta)}\int_0^\infty T(E) e^{- \mathcal{P}f_1(\mathbf{x},\theta)/E^{3}} dE =:  e^{-a \mathcal{P}f_2(\mathbf{x},\theta)} P(\mathcal{P}f_1(\mathbf{x},\theta)) =: e^{-g_{BH}}(\mathbf{x},\theta) \%$$
with
$$P(z) := \int_0^\infty T(E) e^{-z/E^3} dE. \%$$
The previous decomposition of the function $B(\cdot)$ shows that the strategy to compute a good approximation of the beam-hardening is to design the suited representation for the polychromatic part $P(z)$. Assuming to know the function $T(E)$, we propose to approximate $P(z)$ by 
$$\tilde{P}(z) = \left(1+bz\right)^{-c}, \ b,c >0. \%$$ 
The estimation of parameters $b$ and $c$ can be performed using a Gauss-Newton algorithm or a Levenberg-Marquardt algorithm for a more robust convergence. 

These approximations on the polychromatic part and on the monochromatic part lead to approximate the whole data by 
$$-\ln (B(\mathcal{P}f_1,\mathcal{P}f_2)) \approx a \mathcal{P}f_2 + c \ln\left(1+b\mathcal{P}f_1 \right).\%$$
Without more prior information about $\mathcal{P}f_1$ and $\mathcal{P}f_2$, this formula does not provide a unique solution to reverse the effect of beam-hardening but enables a simpler interpretation. However, if we consider homogeneous objects, the formula becomes invertible.


#### Inversion for homogeneous objects

In the case of homogeneous objects the functions $f_1$ and $f_2$ are constant with respective values $a_1,a_2>0$, i.e.
$$f_E(\mathbf{x}) = (a_1 E^{-3} + a_2 C(E)) f(\mathbf{x})\%$$
with $f(\mathbf{x}) = 1$ if $\mathbf{x}$ belongs to the support of the object, 0 otherwise. 

Putting $\alpha=a_2 a$ and $\beta = a_1 b$, the approximation for the data altered by beam-hardening becomes in the homogeneous case
$$-\ln ({B}(\mathcal{P}f)) \approx \alpha \mathcal{P}f + c \ln\left((1+\beta \mathcal{P}f ) \right) =: \tilde{g}. \%$$
This analytic model presents the advantages of being invertible thanks to the W-Lambert function, $W(x)$. More precisely, we get  
$$\mathcal{P}f = \frac{\beta c W\left( \frac{\alpha}{\beta c}e^{\frac{\alpha+\beta \tilde{g}}{\beta c}} \right)-\alpha}{\alpha\beta} = \mathcal{C}(\tilde{g},\alpha,\beta,c). \%$$
From this model we can therefore correct the effect of beam-hardening and recover the sought-for data $\mathcal{P}f$ in the homogeneous case. By applying the derived inversion formula on the real data $g_{BH}$, we produce corrected data with $g_{corr}  =  \mathcal{C}(g_{BH},\alpha,\beta,c)$.
The image reconstruction consists then to apply a standard reconstruction technique suited to $\mathcal{P}$ on the corrected data $g_{corr}$. 

This is performed by the function <code>BHcorrection()</code> as illustrated below.


<code>BHcorrection(g_bh,T_E,E,intervalg,absorption_coef,diffusion_coef,tau=0.5,tol=10**(-12))</code> computes the polychromatic function $P()$, the curve fitting and then the correction with the help of the W-Lambert function.

The curve fitting is performed by the Levenberg-Marquardt algorithm via <code>params = LM(P,x,p0,tol)</code> where 
* P denotes the polychromatic part computed from the spectrum $T_E$. 
* $x$ denotes the variable of integration and contains 100 points equally distributed on the interval <code>intervalg</code> given as input of <code>BHcorrection()</code>. For lower energies or extremly dense materials, the range of this variable can be huge. In that case, it can be important to increase the lower bound of <code>intervalg</code> in order to focus on the correction of the higher intensities which are affected the most by beam-hardening. 
* p0 is an initial value for the parameters b and c. 
* tol stands for the tolerance for convergence.
    
The correction is then performed by <code>EvalCorrection(g, alpha, beta, c)</code> where g stands for the CT-data and 
* <code>alpha</code> = $(C_E[0]+$<code>tau</code>$*(C_E[-1]-C_E[0]))$*<code>diffusion_coef</code> with <code>tau</code> given as input of <code>BHcorrection()</code>. Since the approximation of the diffusion coefficient is more difficult for lower energies, this parameter helps to control the contribution of the monochromatic part of the CT-data in the correction. 
* <code>beta</code> = <code>absorption_coef</code>*<code>params[0]</code>
* <code>c</code> = <code>params[1]</code>

It is of course possible to use these functions directly instead of <code>BHcorrection()</code>.

Standard reconstruction approaches such as the FBP can then be applied.
