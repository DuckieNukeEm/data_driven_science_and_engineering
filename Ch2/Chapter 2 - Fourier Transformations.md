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

So, I'm going to do this a little differently than the book. Generally speaking, in textbook, we have a nice clean coordinate system: X, Y, and Z are all orthogonal to each other. Reality...is quiet different. Think about a jet engine. It has a the X,Y,Z coordinates, yes, but also a time dimension, a temperature dimension, a heat dimension, a fuel dimension, a force dimension, and drag dimensions...it as a lot of fucking dimensions. And they may not be orthogonal to eachother. Time and heat are clearly related, as well as force and fuel.

That's where Fourier came in. Fourier discovered that you can map a complex coordinates, repeating sequence of events onto an infinte series of sine/cosine of increasing frequency. Each step of the series would be at a every increasing frequency, with a given amplitude and association percentage. IE, it maps $\mathbb{R}^{n} -> \Sigma(sin(x) + cosin(x)$. Which is pretty cool. 

This also relates to the Taylor series. You can approximate a function (either absolutely or in infinite) a function by taking the Taylor series. However, some functions are not well suited for a Taylor expansion. However, a Fourier Series *can* approximate these functions. If you want a rough approximation, you could pick, say, the first 10 frequencies. Or first 100...or first 1000.

In order to construct this series of frequencies, we need two things. 

1. Determine the values for each sine/cosine at any given frequency
2. The loading for the frequency

#### 1. Determine the value of sine/cosine at a given frequency

This is pretty easy to solve for: just take the value at the frequency.
$cos(kx)+sin(kx)$

with $k$ being the specific order of frequency. However, it's not clean, lets clean it up a bit by making it whole cycle of $2\pi$.

$cos(2\pi kx)+sin(2\pi kx)$

And what the length of the duration is $L$ long, so like...it's five cycles long, not the standard $2\pi$ we are use to
$cos(\frac{2\pi kx}{L})+sin(\frac{s\pi kx}{L})$

#### 2. The loading factor

So, now we need to figure out how each cosine and sine "load", or how much each of the frequency apply to the decomposition. 

$a_kcos(\frac{2\pi kx}{L})+b_ksin(\frac{s\pi kx}{L})$

where

$a_k = \frac{2}{L}\int_{0}^{L}f(x)cos(\frac{2\pi kx}{L})dx$
$b_k = \frac{2}{L}\int_{0}^{L}f(x)sin(\frac{2\pi kx}{L})dx$

Wait...so what's with this integral of the $f(x)$ against the cosine/sine? THAT! my dear friend, is the [inner product](https://en.wikipedia.org/wiki/Inner_product_space). It's a project of $f(x)$ on to cosine (or sine), scaled by sine. NOw, the issue with the integral is that the large the value $f(x)$, the large the value will be, so we scaled to be between 0 and 1 by multiplying it by $\frac{2}{L}$.

An alternative way to write the inner product is:
$<f(x),cosine(kx)>=\int_{0}^{L}f(x)cos(\frac{2\pi kx}{L})dx$
$\lt f(x),sine(kx)\gt=\int_{0}^{L}f(x)sin(\frac{2\pi kx}{L})dx$

So now we have all the parts. For any given frequencies, we know what the value of the sine/consine is, and how much it loads on it. All we have to do now, is run it across **all** frequencies

$f(x) = \Sigma_0^\infty a_kcos(\frac{2\pi kx}{L})+b_ksin(\frac{s\pi kx}{L})$

with $L$ being the length of the period and $k$ being the increasing frequency.

#### Minimum frequency

If we take a look at the above formula, there is this constant, $\frac{2\pi}{L}$, this defines the *minimum* frequency that can be used (when k = 1). This is the smallest possible frequency we can use, any resolution beyond that can't be measured.

#### Properties within the Fourier Domain:

The Fourier Domain has some really nice properties.

1. A derivative of a function in $\real$ becomes a multiplier of the base frequencies. $F(\frac{d}{dx}f(x)) = i\omega F(f(x))$

2. Fourier Series is a linear operator. $F(\alpha f(x) + \beta f(x)) = \alpha F(f(x)) + \beta F(f(x))$

3. The Fourier transformation preservers the $L_2$ Norm (up to a constant). This is known as Parseval's Theorem. $\int^{\infin}_{-\infin}|\hat f(\omega)|^2d\omega = 2\pi \int^{\infin}_{-\infin}|f(x)|^2$_

4. Convultions are well behaved in the Fourier domain. Multiplying functions in the frequency domain is the same as convolving function in the spatial domain. I'm not going to write out this formula....it's a beast.

##### Links

[Overview Video](https://www.youtube.com/watch?v=jNC0jxb0OxE)
[Fourier Series Part 1](https://www.youtube.com/watch?v=MB6XGQWLV04)
[Fourier Series Part 2](https://www.youtube.com/watch?v=Ud9Xtxsi2HI)
[Inner Product](https://www.youtube.com/watch?v=g-eNeXlZKAQ)
[Complex Fourier Series](https://www.youtube.com/watch?v=4cfctnaHyFM)

[Transforms and derivatives](https://www.youtube.com/watch?v=d5d0ORQHNYs)
[Transforms and Convolutions](https://www.youtube.com/watch?v=mOiY1fOROOg)
[Parsevals Theorum](https://www.youtube.com/watch?v=mOiY1fOROOg)

## 2. Discrete Fourier transform (DFT) and Fast Fourier Transform (FFT) <a name="dftfft"></a>

Generally speaking, in the real world, with us normies. We don't have any process that is continuous. So we have to take it in discrete amounts. This isn't a problem with the Fourier Transformation.

If we have a function, $f$ and we take n number of samples points, each point being in $[f_0, f_1, f_2, ... ,f_{n}]$ then we can compute the DFT equivalent, $\hat f$ as $\hat f_k=\Sigma^n_{j=0}f_je^{jk\frac{2i\pi}{n}} $ and the inverse DFT is very similar: $f_k=\frac{1}{n}\Sigma^{n-1}_{j=0}\hat f_je^{jk\frac{2i\pi}{n}}$

However, this is computationally VERY expensive, taking $O(n^2)$ time. 

The Fast Fourier Transformation, takes $O(nlog(n))$ time and is very efficient. I won't go into how it's done or derived, or even the math behind it (I'm still not sure) but it's used EVERYWHERE.

##### example 2.1 Noise filtering with FFT

```matlab
dt = .001;
t = 0:dt:1;
f = sin(2* pi *50*t) + sin(2* pi *120*t); % Sum of 2 frequencies
f = f + 2.5* randn(size(t)); % Add some noise
%% Compute the Fast Fourier Transform FFT
n = length(t);
fhat = fft(f,n);
% Compute the fast Fourier transform
PSD = fhat.* conj(fhat)/n; % Power spectrum (power per freq)
freq = 1/(dt*n)*(0:n); % Create x-axis of frequencies in Hz
L = 1:floor(n/2);
% Only plot the first half of freqs
%% Use the PSD to filter out noise
indices = PSD>100; % Find all freqs with large power
PSDclean = PSD.*indices; % Zero out all others
fhat = indices.*fhat; % Zero out small Fourier coeffs. in Y
ffilt = ifft(fhat); % Inverse FFT for filtered time signal
%% PLOTS
subplot(3,1,1)
plot(t,f,’r’,’LineWidth’,1.2), hold on
plot(t,f,’k’,’LineWidth’,1.5)
legend(’Noisy’,’Clean’)
subplot(3,1,2)
plot(t,f,’k’,’LineWidth’,1.5), hold on
plot(t,ffilt,’b’,’LineWidth’,1.2)
legend(’Clean’,’Filtered’)
subplot(3,1,3)
plot(freq(L),PSD(L),’r’,’LineWidth’,1.5), hold on
plot(freq(L),PSDclean(L),’-b’,’LineWidth’,1.2)
legend(’Noisy’,’Filtered’)
```

##### example 2.2 spectral derivatives

```matlab
n = 128;
L = 30;
dx = L/(n);
x = -L/2:dx:L/2-dx;
f = cos(x).* exp(-x.ˆ2/25);
df = -(sin(x).* exp(-x.ˆ2/25) + (2/25)*x.*f);
% Function
% Derivative
%% Approximate derivative using finite Difference...
for kappa=1:length(df)-1
dfFD(kappa) = (f(kappa+1)-f(kappa))/dx;
end
dfFD(end+1) = dfFD(end);
%% Derivative using FFT (spectral derivative)
fhat = fft(f);
kappa = (2* pi/L)*[-n/2:n/2-1];
kappa = fftshift(kappa); % Re-order fft frequencies
dfhat = i*kappa.*fhat;
dfFFT = real(ifft(dfhat));
%% Plotting commands
plot(x,df,’k’,’LineWidth’,1.5), hold on
plot(x,dfFD,’b--’,’LineWidth’,1.2)
plot(x,dfFFT,’r--’,’LineWidth’,1.2)
legend(’True Derivative’,’Finite Diff.’,’FFT Derivative’)
```

###### Videos

[Overview of DFT Video](https://www.youtube.com/watch?v=nl9TZanwbBk)
[Computing the DFT](https://www.youtube.com/watch?v=Xw4voABxU5c)
[Overview of FFT](https://www.youtube.com/watch?v=E8HeD-MUrjY)
[algorithm of FFT](https://www.youtube.com/watch?v=toj_IoCQE-4)
[Denoising FFT](https://www.youtube.com/watch?v=c249W6uc7ho)

## 3. Transforming Partial Differential Equations <a name="difeq"></a>

Ignoring for now, man....I need to come back to this and really get better at it

##### Videos

[Solving the Heat Equation](https://www.youtube.com/watch?v=7haZCrQDHpA)
[Solving PDE with Python (1)](https://www.youtube.com/watch?v=hDeARtZdq-U)
[Solving PDE with Python (2)](https://www.youtube.com/watch?v=mMdIxa5qC9Y)

## 4. Gabor Transform and the Spectrogram <a name="Gabor"></a>

A limitation of of the Fourtier transformation is that it doesn't tell you *when* in time a frequency occures. Recall, up to this point, we were using fixed length or repeating sequences in which knowing in time when that occures isn't needed. But it does become important when we start to consider systems that evolve over time. i.e. music. So, we have the Gabor transofmration. 

The theory is very straight forward. We will just have a sliding window, that has a center point $\tau$ , and run the Fourier transform over it, then we'll move $\tau$ over a hair, and do it again, and again, and again until we've covered the entire time $t$.  

<img src="file:///home/asmodi/.config/marktext/images/2023-04-30-10-40-34-Gabor%20window.png" title="" alt="" width="605">

That localized window, $g_{t,\omega}(\tau)$, is called the kernel and is defined as: $g_{t,\omega}=e^{i\omega\tau}g(\tau - t)$, with $g(\tau-t)$ being any function, but is generally set as to $g(t)=e^{-\frac{(t-\tau)^2}{a^2}}$ for simplicity.

Take the kernel $e^{i\omega t}g(\tau-t)$ and then multiple it by Fourieir transform.

now if we do the Fourier transform on it, we can get a function that takes two inputs, time $t$, and frequency, $\omega$.

$G(f)(t,\omega) = \hat f_g(t,\omega)=\int_{-\infin}^\infin f(\tau)e^{-i\omega \tau}\bar g(\tau - t)d\tau = <f,g_{t,\omega}> $

 $g_{t,\omega}$ is called the kernel, and is written as: $g_{t,\omega}=e^{i\omega\tau}g(\tau - t)$, with $g(t)$, called the kernal, usually set as $g(t)=e^{-\frac{(t-\tau)^2}{a^2}}$, which ultimatly gives us 

$G(f)(t,\omega) = \hat f_g(t,\omega)=\int_{-\infin}^\infin f(\tau)e^{-\frac{t^2}{a^2}}d\tau$ 

This isn't bad. You can see we evaluate function at $\tau$, given time step $t$, at frequency $\omega$.

###### Videos

[Spectrogram](https://www.youtube.com/watch?v=EfWnEldTyPA)
[Uncertaintiy PRinciples](https://www.youtube.com/watch?v=EfWnEldTyPA)

## 5. Wavelets and Multi-Resolution Analysis <a name="wavelets"></a>

##### Videos

[Wavelets and Multi-Resolution Analysis](https://www.youtube.com/watch?v=y7KLbd7n75g)

## 6. 2D transforms and Image Processing <a name="image"></a>

##### Videos

[Image compression](https://www.youtube.com/watch?v=gGEBUdM0PVc)
[Laplace transformation](https://www.youtube.com/watch?v=7UvtU75NXTg)
[Differential Equations](https://www.youtube.com/watch?v=iBde8qOW0h0)
