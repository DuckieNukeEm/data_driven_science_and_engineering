---
tags: data-driven-science, Fourier, transformations
---

# Chapter 2 - Fourier Transformations
---
## Table of contents
1. [Fourier Series and Fourier Transforms](#intro_fourier)
2. [Discrete Fourier Transfomr (DFT) and Fast Fourier Transform (FFT)](#dftfft)
3. [Transforming Partial Differential Equations](#difeq)
4. [Gabor Transform and the Spectrogram](#Gabor)
5. [Wavelets and Multi-Resolution Analysis](#wavelets)
6. [2D transforms and Image Processing](#image)
----

## 1. Fourier Series and Fourier Transforms <a name="intro_fourier"></a>
So, I'm going to do this a little differently than the book. Generally speaking, in textbook, we have a nice clearn coordinate system: X, Y, and Z are all orthogonal to eachother. Relality...is quiet different. Think about a jet engine. It has a the X,Y,Z coordinates, yes, but also a time dimension, a temperture dimension, a heat dimension, a fueld dimension, a force dimension, and drag dimensions...it as a lot of fucking dimensions. And they may not be orthoginal to eachother. Time and heat are clearly related, as well as force and fuel.

Wouldn't it be nice, however, if we could take the coordinate plane of that jet engine, and move it into an infinte dimensional space, with all dimensions guarneteed to be orthogonial to eachother. Wouldn't that be swell!

That's where Fourier came in. 

$\Biggl\{\begin{align*} 
0 \;for\; x \in [-\pi,-\pi/2 ) \\ 
1+2x/\pi \;for\; x \in [-\pi/2,0) \\
1-sx/\pi \;for\; x\in [0,\pi/2)\\
0 \;for\; x \in [\pi/2,\pi)
\end{align*}$

[Overview Video](https://www.youtube.com/watch?v=jNC0jxb0OxE)
[Fourier Series Part 1](https://www.youtube.com/watch?v=MB6XGQWLV04)
[Fourier Series Part 2](https://www.youtube.com/watch?v=Ud9Xtxsi2HI)
[Inner Product](https://www.youtube.com/watch?v=g-eNeXlZKAQ)
[Complex Fourier Series](https://www.youtube.com/watch?v=4cfctnaHyFM)

## 2. Discrete Fourier Transfomr (DFT) and Fast Fourier Transform (FFT) <a name="dftfft"></a>

## 3. Transforming Partial Differential Equations <a name="difeq"></a>

## 4. Gabor Transform and the Spectrogram <a name="Gabor"></a>

## 5. Wavelets and Multi-Resolution Analysis <a name="wavelets"></a>

## 6. 2D transforms and Image Processing <a name="image"></a>

