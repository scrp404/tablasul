# Optical Theory & Mathematical Principles in Tablasul

`tablasul` is founded strictly upon first-principles wave optics and paraxial Fourier diffraction theory. This document provides the mathematical derivations for each component in the optical train.

---

## 1. Complex Wavefield Representation

The electric field of an optical beam at an entrance pupil plane $(x_p, y_p)$ is modeled as a complex 2D scalar field:

$$E_0(x_p, y_p) = A(x_p, y_p) \exp\left[ i \phi(x_p, y_p) \right]$$

where:
* $A(x_p, y_p) \in [0, 1]$ is the real amplitude transmission mask of the telescope aperture.
* $\phi(x_p, y_p)$ is the total optical phase delay in radians, given by:
  $$\phi(x_p, y_p) = \frac{2\pi}{\lambda} \cdot \text{OPD}(x_p, y_p)$$
* $\text{OPD}$ is the physical Optical Path Difference introduced by atmospheric turbulence, telescope mirror figure errors, or static aberrations.

---

## 2. Paraxial Fraunhofer Propagation

Under the Fraunhofer approximation, the complex electric field $E_f(u, v)$ in the focal plane is proportional to the 2D spatial Fourier transform of the pupil field $E_0(x_p, y_p)$:

$$E_f(u, v) = \frac{1}{i \lambda f} \iint_{-\infty}^{\infty} E_0(x_p, y_p) \exp\left[ -i \frac{2\pi}{\lambda f} (x_p u + y_p v) \right] dx_p dy_p$$

### Discrete Implementation & Energy Conservation
Numerically, propagation is computed using a centered 2D Fast Fourier Transform on an $N \times N$ discrete grid:

$$E_f = \frac{1}{N} \text{fftshift}\left( \text{fft2}\left( \text{ifftshift}(E_0) \right) \right)$$

This unitary normalization guarantees exact Parseval energy conservation between pupil and focal planes:

$$\sum_{i, j} |E_f[i, j]|^2 = \sum_{i, j} |E_0[i, j]|^2$$

### Focal Plane Physical Pixel Scale
The focal plane sampling resolution element $\lambda / D$ (radians) corresponds to:

$$\Delta u_{\text{pixel}} = \frac{N_{\text{pupil}}}{N_{\text{grid}}} \cdot \left(\frac{\lambda}{D}\right)$$

where $N_{\text{pupil}}$ is the number of grid pixels spanning the entrance pupil diameter $D$.

---

## 3. Pyramid Wavefront Sensor (PyWFS)

The 4-faceted pyramid prism acts as a focal plane phase mask with apex located on the optical axis $(u=0, v=0)$.

### Pyramid Focal Phase Mask
Each of the four quadrants tilts the transmitted light outward by an angle $(\pm \alpha, \pm \beta)$ corresponding to a spatial shift $\Delta$ in the conjugate pupil plane:

$$\Phi_{\text{pyr}}(u, v) = \frac{2\pi \Delta}{N} \left( |u| + |v| \right)$$

This phase transformation deflects the beam into four non-overlapping pupil images $I_1, I_2, I_3, I_4$ upon propagation to the re-imaging detector:

$$E_{\text{det}} = \text{fftshift}\left( \text{ifft2}\left( \text{ifftshift}\left( E_f \cdot e^{i \Phi_{\text{pyr}}} \right) \right) \right) \cdot N$$

### Dynamic Circular Tip-Tilt Modulation
When modulation with radius $r_{\text{mod}}$ (in $\lambda/D$) is active, the beam traces a circle around the pyramid apex. This is modeled by integrating over $N_{\text{mod}}$ discrete pupil tilts $\theta_k = \frac{2\pi k}{N_{\text{mod}}}$:

$$\phi_{\text{mod}, k}(x_p, y_p) = \pi r_{\text{mod}} \left( \frac{x_p}{R_{\text{pup}}} \cos\theta_k + \frac{y_p}{R_{\text{pup}}} \sin\theta_k \right)$$

The final recorded detector intensity is the incoherently integrated average:

$$I_{\text{det}}(x, y) = \frac{1}{N_{\text{mod}}} \sum_{k=1}^{N_{\text{mod}}} |E_{\text{det}, k}(x, y)|^2$$

### Quad-Cell Slope Estimation
From the extracted 4 pupil quadrants $I_1, I_2, I_3, I_4$, normalized slopes $(S_x, S_y)$ are computed:

$$S_x = \frac{(I_1 + I_4) - (I_2 + I_3)}{I_1 + I_2 + I_3 + I_4 + \epsilon}, \quad S_y = \frac{(I_1 + I_2) - (I_3 + I_4)}{I_1 + I_2 + I_3 + I_4 + \epsilon}$$

---

## 4. Classical Lyot Coronagraph

The Lyot coronagraph operates in three sequential optical planes:

1. **Focal Plane Mask (FPM)**:
   An opaque occulting disk of angular radius $r_{\text{FPM}}$ (typically $2.0\text{--}3.5\,\lambda/D$) suppresses on-axis starlight. With smooth edge apodization width $\Delta r$, transmission is:
   $$M_{\text{FPM}}(r) = \begin{cases} 0 & r \le r_{\text{in}} \\ \frac{1}{2}\left(1 - \cos\left(\pi \frac{r - r_{\text{in}}}{\Delta r}\right)\right) & r_{\text{in}} < r < r_{\text{out}} \\ 1 & r \ge r_{\text{out}} \end{cases}$$
   $$E_{f, \text{masked}} = E_f(u, v) \cdot M_{\text{FPM}}(u, v)$$

2. **Lyot Stop Pupil Plane**:
   Light diffracted by the occulting mask is re-imaged to a conjugate pupil plane:
   $$E_{\text{Lyot}} = \text{IFFT}(E_{f, \text{masked}})$$
   The diffracted on-axis starlight forms bright rings along the outer edge of the pupil and central obscuration. An undersized Lyot diaphragm $M_{\text{Lyot}}(x_p, y_p)$ ($85\%\text{--}90\%$ of pupil diameter) truncates these diffraction edges.

3. **Science Focal Plane**:
   The final science focal plane is obtained by re-focusing:
   $$E_{\text{sci}} = \text{FFT}(E_{\text{Lyot}} \cdot M_{\text{Lyot}}), \quad I_{\text{sci}} = |E_{\text{sci}}|^2$$

---

## 5. Non-Common Path Aberrations (NCPA)

Non-Common Path Aberrations are static or quasi-static optical path differences present in the science arm that are invisible to the WFS arm:

$$\phi_{\text{sci}}(x, y) = \phi_{\text{common}}(x, y) \cdot \left(\frac{\lambda_{\text{wfs}}}{\lambda_{\text{sci}}}\right) + \phi_{\text{ncpa}}(x, y)$$

$\phi_{\text{ncpa}}$ is expanded in orthonormal Zernike polynomials $Z_j(\rho, \theta)$ (Noll indexing):

$$\phi_{\text{ncpa}}(\rho, \theta) = \sum_{j} a_j Z_j(\rho, \theta)$$

---

## 6. Matrix Fourier Transform (MFT) Decoupled Science Camera

In standard discrete 2D FFT propagation, the focal plane field of view and resolution are strictly coupled to the pupil grid dimensions:

$$\text{FOV}_{\text{FFT}} = N_{\text{pupil}} \cdot \left(\frac{\lambda}{D}\right), \quad \Delta u = \frac{N_{\text{pupil}}}{N_{\text{grid}}} \left(\frac{\lambda}{D}\right)$$

For high-contrast imaging systems where the science detector needs an arbitrary field of view (e.g. zooming in on the inner $\pm 15\,\lambda/D$ or zooming out to $\pm 60\,\lambda/D$) and arbitrary detector pixel dimensions ($N_{\text{det}} \times N_{\text{det}}$) without zero-padding the entire simulation grid, `tablasul` incorporates Soummer's Matrix Fourier Transform (MFT):

$$E_f(u, v) = \mathbf{M}_y \cdot E_{\text{pupil}} \cdot \mathbf{M}_x^T$$

where the 1D kernels are:

$$\mathbf{M}_x[j, k] = \exp\left[ -i 2\pi \cdot u_j \cdot x_k \right] \Delta x \sqrt{\frac{\Delta u \Delta v}{\lambda^2}}$$

This decouples the science camera sampling and field of view entirely from the PyWFS detector grid, while preserving exact energy scaling.

---

## 7. Measured Entrance Pupils (`PupilImageAperture`)

Experimental high-contrast testbeds and telescopes (e.g., ground-based telescopes or optical testbeds) often have non-ideal pupil geometries with measured spider distortions, actuator dead zones, or segment gaps. `tablasul` accommodates this via `PupilImageAperture`:
* Reads 2D arrays directly from astronomical FITS files, NumPy `.npy` files, or image files.
* Automatically normalizes transmission peak to unity.
* Resamples to the target simulation grid $N_{\text{grid}} \times N_{\text{grid}}$ with bilinear interpolation.
* Automatically estimates the effective pupil pixel span for angular diffraction scaling $\lambda/D$.
