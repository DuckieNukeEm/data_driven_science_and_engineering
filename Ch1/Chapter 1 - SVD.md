---
tags: data-driven-science, SVD, PCA
---

# Chapter 1 - SVD

# Table of contents
1. [Introduction](#introduction)
2. [Matrix Approximation](#matapprox)
3. [Mathematical properties and Manipulation](#props)
4. [PsuedoInverse, Leaset Squares, and Regression](#OLS)
    1. [Least Squares, and Regression](#OLS_R)
5. [Princiapl Component Analysis](#PCA)
6. [Examples]()
7. [Truncation and Alignment](#TRUNC)
8. [Randomized Singular Value Decomposisition](#RAND)
9. [Tensor Decomposisiont](#TENSE)
10. [Appendix](#APEND)  
----

## 1 introduction <a name="introduction"></a>

The SVC provides a system way to detmine a low-diemnstion approximation  of higher dimensional data in terms of dominat patterns. What does that mean, basically a huge data matrix is split into three, way smaller matricies. Furthermore, there is always a solutions 
1. One Matrix holds groupings of how the columns are similar $U$
2. Another matrix holds groups of how the rows are similar $V^t$
3. And the third matrix is an odd stacked zeroish matrix that can easily be truncted $\Sigma$
    
By design, the matrix has the following formula (And below are restatements of the orignal matrix
    
$\bar X=U\Sigma V^T$  
with:
 * $X=mxn$, is our original input matrix  
 * $U=nxn$, is a matrix that holds 'groupings' or 'similariteis' of the columns and is unitary [ie $UU^t=U^tU=U^{-1}U=UU^{-1}=I$]
 * $V=mxm$, is a matrix that holds 'groupings' or 'similarties' of the **rows** of the $\bar X$ and is also, itself, unitary
 * $\Sigma=mxn$ has all zero's off the diagonal and only real, non negative values along the diagonal. What is intresting, if $n \geq m$, then $\Sigma$ has, **at most** $m$ entries along the diagonal. Due to $\Sigma$ only having non-negative values on the diagonal, the inverse exists.
 
Becuase $\Sigma$ has at most $m$ solution, we can write the above formula in a simplfied format, known as **Economy SVD**.
$\bar X=[\hat U \hat U^t][\frac{\hat \Sigma}{0}]V^t$
     
with $\hat U^t$ spanning the column space that is complementaty and orthogonal to that spaned by $\hat U$. Which, basically means, you can use U and be fine. You can also throw away all the zero's of $\Sigma$ as now it only has rank of $m$. This is also a precise solution. So extra savings baby!

[SVD videos](https://www.youtube.com/watch?v=gXbThCXjZFM&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv) 
[SVD Mathematical Overview](https://www.youtube.com/watch?v=nbBvuuNVfco&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
## 2 Matrix Approximation <a name = 'matapprox'></a>

An intresting property about SVD is that you can get a really good matrix approximation (and reduction). The theorm is

***Theorem Eckart-Young***:*The optimal rank-r approximation to $X$ in a least square sense, is given by the rank-r SVD truncation of $X$*

$\argmin \limits_{\tilde X s.t. rank(\tilde X)=r} \|X - \tilde X\|_{F}=\tilde U \tilde \Sigma \tilde V^T$ 

What this means is that if select the first r rows of $U$ $\Sigma$ and first r columns for $V$, we can get a good approcimation: $X\approx \tilde U \tilde \Sigma \tilde V^T$

Which means, compressions. You can see this at work on images. See example 1

[SVD Matrix Approximation video](https://www.youtube.com/watch?v=xy3QyyhiuY4&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
[SVD Frobenius Norm video](https://www.youtube.com/watch?v=Gt56YxMBlVA&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

### 2.1 Examples

```{python}
# loading images
image = Image.open(IMG_DIR + 'lady_and_street.jpg')
npicture = np.asarray(image)[:,:,0]

# performing SVD
nu, ns, nv = np.linalg.svd(npicture, full_matrices = False)

# selecting various lelvels of r
for r in [1,10, 25, 50, 100, 150 , 250, 400, 500, 750, 900, 1000]:
    nuA = nu[:,0:r]
    nsA =np.diag(ns[0:r])
    nvA = nv[0:r,:]
    xApprox =nuA @ nsA @ nvA
    
    # calcualting  what savings we got and the relative error of it all 
    POC = int((nuA.shape[0] * nuA.shape[1] + nsA.shape[0] + nvA.shape[0] * nvA.shape[1])/(1500*1000) * 100)
    error =  np.sum(np.abs(xApprox - npicture))/np.sum(npicture) * 1000 // 1 / 10
    print(f"new dimensions with the reducitons of r {r}: U = {nuA.shape}, S = {nsA.shape}, V = {nvA.shape}, appox = {xApprox.shape}, percent of total {POC}% , error = {error}")
    # saving it down, we need the np.unit8 so that Pillow can deal with the image :/
    im = Image.fromarray(np.uint8(xApprox))
    im.save(IMG_DIR + 'lady_and_street_r' + str(r) + '.jpg')

```

## 3 Mathematic Properties and Manipulation <a name = 'props'></a>

### 3.1 Dominate Correlations
Eh, not mutch here
[SVD Correlations Video](https://www.youtube.com/watch?v=WmDnaoY2Ivs&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
[SVD Methods of Snapshot video](https://www.youtube.com/watch?v=rs63fnUWJkk&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

### 3.2 SVD Invarance

By design, SVD are invarity to unitary transformation (it basically get absorbed into $U$ or $V$):

Let's define $Y=CX$ then $Y=CX=CU_{x}\Sigma_{x} V^{T}_{x}$
This similarly holds for right hand multiplcation: $Y=XP=U_{x}\Sigma_{x} V^{T}_{x}P$

[SVD Invariance video](https://www.youtube.com/watch?v=46Hpy4FiGls&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

## 4 Psuedo-Inverse, Leaset Squares, and Regression <a name = 'OLS'></a>
[Overview Video](https://www.youtube.com/watch?v=PjeOmOz9jSY&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

### 4.2 Least Squares, and Regressoin <a name = 'OLS_R'></a>
So, recall, what OLS is trying to solve:

$\bar xB=y$

To solve for the beta's, we do the following:

$\bar x^{-1}\bar xB=\bar x^{-1}y$  
$IB=\bar x^{-1}y$
$B=\bar x^{-1}y$

Now, $x$ may **not** be intervertable, but don't worry, SVD to the rescue. Recall, by design that $U$ and $V$ are unitary (have an inverse), and $\Sigma$ has a diagonal matrix with non zero elements, thus the ivnerse exists. We can then do the following to get to the inverse

$\bar xB=y$  
$(U\Sigma V^t)B=y$  
$(U\Sigma V^T)^{-1}(U\Sigma V^T)B=(U\Sigma V^T)^{-1}y$   
$V\Sigma ^{-1}U^{-1}U\Sigma V^T B=V\Sigma^{-1}U^{-1}y$  
$V\Sigma ^{-1}I\Sigma V^T B=V\Sigma^{-1}U^{-1}y$  
$V\Sigma ^{-1}\Sigma V^T B=V\Sigma^{-1}U^{-1}y$  
$VIV^T B=V\Sigma^{-1}U^{-1}y$  
$V V^T B=V\Sigma^{-1}U^{-1}y$  
$I B=V\Sigma^{-1}U^{-1}y$  
$B=V\Sigma^{-1}U^{-1}y$  


and viola, you now have $B$. This is far more effencint to compute, taking roughly $O(n^2)$ to compute, while inverting the original matrix, A, will take $O(n^3)$, so an order of magnatude better.

[Least Squares Video](https://www.youtube.com/watch?v=02QCtHM1qb4&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
[Ax=b Video](https://www.youtube.com/watch?v=N9uf0YGDaWk&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
[Linear Regression Video](https://www.youtube.com/watch?v=EDPCsD6BzWE&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

### 5.3 Example
```{python}
# Loading Data set
housing = pd.DataFrame(np.loadtxt('https://archive.ics.uci.edu/ml/machine-learning-databases/housing/housing.data'),
                       columns = [ 'crime_per_capita', # per capita crime rate by town
                                'zone_proportion', # proportion of residential land zoned for lots over 25,000 sq.ft.per capita crime rate by town
                               'nonretail_business', #proportion of non-retail business acres per town
                               'on_river', #Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
                               'nox_ppm',# nitric oxides concentration (parts per 10 million)
                               'rooms', # average number of rooms per dwelling
                               'age', # proportion of owner-occupied units built prior to 1940
                               'distance_empoyment', # weighted distances to five Boston employment centres
                               'highway_index_access', #index of accessibility to radial highways
                               'tax', # full-value property-tax rate per $10,00.60
                               'pupil-teacher-ratio', #pupil-teacher ratio by town
                               'blank index', # 1000(Bk - 0.63)^2 where Bk is the proportion of blacks  by town
                               'lower_status', # % lower status of the population
                               'median_price' # Median value of owner-occupied homes in $1000's
                                 ])
                                 
# Seperating out dependent variable and adding intercept
b_housing = housing.iloc[:,13]
a_housing = housing.iloc[:,0:13]
a_housing['intercept'] = 1

# performing SVD
U, S, V = np.linalg.svd(a_housing.values, full_matrices = False)

# showing Beta's
np.transpose(V) @ np.linalg.inv(np.diag(S)) @ np.transpose(U) @ b_housing.values

```

## 5 PCA <a name = 'PCA'></a>

A direct application of SVD is being able to easily get the Eigenvectors and eignevalues of a matrix.

First, you need to zero center the data. To do that,  you take a row-wise average (IE, sum up a column and then divded by the number of entires). Formally, define $x_j = \frac{1}{n}  \Sigma^{n}_{i=1}x_ij$, thus, the mean matrix is defined as
$\bar X = [1, 1, 1, .... 1] * [x_1,x_2,....,x_m]$
and the centered data matrix, $X$, can be calculated as
$B=X-\bar X$

Taking an SVD yields
$B=USV^T = P_cV^T$

Where $P_c$ is the principal componentents and $V^T$ holds the principal directions/ais.

So how were we able to use SVD to get PCA?

Recall the covariance matrix:
$C = \frac{B^{-1}B}{n-1}$ 

Okay...so what, how are these related? First, becuase $C$ is hermitian (it's it's own transpose), we can diaganlize it.
$C=VDV^{T}$

$V$ will holds the Eigenvectors and $D$ holds the eigenvalues. So, start from the definition of $B$
$B=USV^T$ (becuase it's a linear shift of $X$)
$C=\frac{B^{-1}B}{1-n}$
$C=\frac{1}{n-1}VSU^TUSV^T$ Note: $U$ is unitary 
$C=\frac{1}{n-1}VS^2V^T$ Note: $S$ is a diagonal matrix so $S=S^T$
$C=V\frac{S^2}{n-1}V^T=VDV^T$

Thus, you can see if we have $V$ from SVD, you automatically have D and the diagnoalization of C.

Thus in $X=USV^T$, $US$ is are the the principal components and $V^T$ is the principal directions/Axis with each column being othogonoal to each other column. 

[More facts from stackoverflow](https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca)

[PCA Overview Video](https://www.youtube.com/watch?v=fkf4IBRSeEc&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

## 7 Truncation and Aligment <a name = 'TRUNC'></a>
Picking where to truncate the data has always been a problem. Similar to K-means clusters, you would look for a elbow or a knee to that would indicate you are at a good point to truncate the SVD (include on the first $p$ ranks of the dataset.) However, there are several methods you can use to find the optimal cut off.

First, if we assume that there is some, erm, noise in the data $X=X_{true} + \gamma X_{guassian noise}$ then we have several solutions that we can find the optimal solutions

1. If $X\exists R^{nxn}$ (ti's square)
Then the hard stop is supper easy. The hard stop $tau$ is $tau =\frac{4}{\sqrt(3)}\sqrt{n}\gamma$

2. if $X\exists R{nxm}$ then
$tau=\lambda(\beta)\sqrt{n}\gamma$ with $\lambda(\beta )=(2(\beta + 1)+\frac{8\beta }{(\beta +1)+(\beta ^2+14\beta+1)^(\frac{1}{2})})^(\frac{1}{2})$  
if $m>>n$ then $\beta =\frac{n}{m}$
if $m<<n$ then $\beta =\frac{m}{n}$


3. If we don't know $\gamma $
then...eh, look at the book. You'll have to solve for this numerical :/

**Alignment**
So, it's important to note, that SVD can't handle rotations, translate or scale. The reason is that, at it's heart, SVD is geometric, IE it's aligned to the axiss of the matrix.  A general $mXn$ matrix is algined along the *x* and *y* axis. When rotating, it s breaks that alignment and you suddenly see an explosion in the number of ranks required to generate a good image (also the principlae componets become a whole lot less effects)

This is a key weakness of SVD.

[Optimal Truncation Video](https://www.youtube.com/watch?v=9vJDjkx825k&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)

## 8. Randomized Singular Value Decomposisition <a name = 'RAND'></a>

Due to the size of matricies now day, we need more efficent ways to compute SVD, hence the randomization of the this process. it's a tad cumberson and I'll add a image of it here. But for right now....this isn't a problems.

[Randomized SVD Video](https://www.youtube.com/watch?v=fJ2EyvR85ro&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)
[Power itterations and oversampling Video](https://www.youtube.com/watch?v=UXXMbpLI7AM&list=PLMrJAkhIeNNSVjnsviglFoY2nXildDCcv)


## 9 Tensor Decomposisiont <a name = 'TENSE'></a>

Now, we can easily expanded out SVD to multiple dimensions. this is useful if you have a movie. A x-y frame moving through the Z time dimension. It's not really gone into so when I find the numpy verison of it, I'll provided it.

[Tensor decomposition Video](https://www.youtube.com/watch?v=tm5am60CId4)


## Appendix A <a name = 'APPEND'></a>
[Source](http://www.databookuw.com/page-2/page-4/ )
