## Fourier transforms conventions

This document presents a brief overview of Fourier transforms in the tight-binding setting, focusing on two conventions that have been used during the development of this codebase.

Specifically, we distinguish between two common conventions for the Fourier transform. For lack of standard terminology, we will refer to them as:

1. **Convention 1:** This is the convention used by [PythTB](https://www.physics.rutgers.edu/pythtb/index.html), and it is particularly well-suited for the computation of Berry phases and related quantities.

2. **Convention 2:** This convention is adopted in [this review article](https://www.annualreviews.org/content/journals/10.1146/annurev-conmatphys-041720-124134), and is generally more prevalent in the literature.

In the following sections, we will carefully examine both conventions and demonstrate how they are related to each other.

### Convention 1

This convention defines the Fourier transform as:

```math
φ_{I,𝐤}(𝐫) = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ψ_I(𝐫-𝐭),
```
which yields the inverse Fourier transform:

```math
ψ_I(𝐫-𝐭) = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·(𝐭+𝐪_α)} φ_{I,𝐤}(𝐫)
```

In our framework, we are also interested in the quantization of these functions in order to construct the tight-binding Hamiltonian. We can express the same relations in the second-quantized notation as:

```math
\ket{φ_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} \ket{ψ_{I,𝐭}},
```
and:

```math
\ket{ψ_{I,𝐭}} = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·(𝐭+𝐪_α)} \ket{φ_{I,𝐤}}
```

This quantization also determines how the creation and annihilation operators transform under this convention:

```math
â_{I,𝐤}^† \ket{\text{vac}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ĉ_{I,𝐭}^† \ket{\text{vac}} ⇒ â_{I,𝐤}^† = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} ĉ_{I,𝐭}^†,
```
and similarly:

```math
ĉ_{I,𝐭}^† = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·(𝐭+𝐪_α)} â_{I,𝐤}^†
```

### Convention 2

This convention defines the Fourier transform as:

```math
\tilde{φ}_{I,𝐤}(𝐫) = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} ψ_I(𝐫-𝐭),
```
which yields the inverse Fourier transform:

```math
ψ_I(𝐫-𝐭) = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·𝐭} \tilde{φ}_{I,𝐤}(𝐫)
```

In our framework, we are also interested in the quantization of these functions in order to construct the tight-binding Hamiltonian. We can express the same relations in the second-quantized notation as:

```math
\ket{\tilde{φ}_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} \ket{ψ_{I,𝐭}},
```
and:

```math
\ket{ψ_{I,𝐭}} = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·𝐭} \ket{\tilde{φ}_{I,𝐤}}
```

This quantization also determines how the creation and annihilation operators transform under this convention:

```math
\hat{\tilde{a}}_{I,𝐤}^† \ket{\text{vac}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} ĉ_{I,𝐭}^† \ket{\text{vac}} ⇒ \hat{\tilde{a}}_{I,𝐤}^† = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} ĉ_{I,𝐭}^†,
```
and similarly:

```math
ĉ_{I,𝐭}^† = \frac{1}{\sqrt{N}} \sum_𝐤 e^{-i𝐤·𝐭} \hat{\tilde{a}}_{I,𝐤}^†
```

### How they relate to each other

- **Convention 1:** $\ket{φ_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·(𝐭+𝐪_α)} \ket{ψ_{I,𝐭}}$
- **Convention 2:** $\ket{\tilde{φ}_{I,𝐤}} = \frac{1}{\sqrt{N}} \sum_𝐭 e^{i𝐤·𝐭} \ket{ψ_{I,𝐭}}$

Then we have that:

```math
\ket{\tilde{φ}_{I,𝐤}} = e^{-i𝐤·𝐪_α} \ket{φ_{I,𝐤}}
```