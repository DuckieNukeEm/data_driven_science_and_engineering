# Chapter 1 - SVD

# Table of contents
1. [Introduction](#introduction)
2. [Matrix Approximation](#matapprox)
3. [Mathematical properties and Manipulation](#props)
4. [PsuedoInverse, Leaset Squares, and Regression](#OLS)
    1. [Least Squares, and Regression](#OLS_R)
5. [Princiapl Component Analysis]()
6. [Examples]()
7. [Tyuncation and Alignment]()
8. [Randomized Singular Value Decomposisition]()
9. [Tensor Decomposisiont]()
    
3. [Another paragraph](#paragraph2)

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
 * $U=nxn$, is a matrix that holds 'groupings' or 'similariteis' of the columns and is unitary [ie $UU^t=U^tU=I$]
 * $V=mxm$, is a matrix that holds 'groupings' or 'similarties' of the **rows** of the $\bar X$ and is also, itself, unitary
 * $\Sigma=mxn$ has all zero's off the diagonal and only real, non negative values along the diagonal. What is intresting, if $n \geq m$, then $\Sigma$ has, **at most** $m$ entries along the diagonal. Due to $\Sigma$ only having non-negative values on the diagonal, the inverse exists.
 
Becuase $\Sigma$ has at most $m$ solution, we can write the above formula in a simplfied format, known as **Economy SVD**.
$\bar X=[\hat U \hat U^t][\frac{\hat \Sigma}{0}]V^t$
     
with $\hat U^t$ spanning the column space that is complementaty and orthogonal to that spaned by $\hat U$. Which, basically means, you can use U and be fine. You can also throw away all the zero's of $\Sigma$ as now it only has rank of $m$. This is also a precise solution. So extra savings baby!
    
## 2 Matrix Approximation <a name = 'matapprox'></a>

An intresting property about SVD is that you can get a really good matrix approximation (and reduction). The theorm is

***Theorem Eckart-Young***:*The optimal rank-r approximation to $X$ in a least square sense, is given by the rank-r SVD truncation of $X$*

$\argmin \limits_{\tilde X s.t. rank(\tilde X)=r} \|X - \tilde X\|_{F}=\tilde U \tilde \Sigma \tilde V^T$ 

What this means is that if select the first r rows of $U$ $\Sigma$ and first r columns for $V$, we can get a good approcimation: $X\approx \tilde U \tilde \Sigma \tilde V^T$

Which means, compressions. You can see this at work on images. See example 1

## 3 Mathematic Properties and Manipulation <a name = 'props'></a>

### 3.2 SVD Invarance

By design, SVD are invarity to unitary transformation (it basically get absorbed into $U$ or $V$):

Let's define $Y=CX$ then $Y=CX=CU_{x}\Sigma_{x} V^{T}_{x}$
This similarly holds for right hand multiplcation: $Y=XP=U_{x}\Sigma_{x} V^{T}_{x}P$

## 4 Psuedo-Inverse, Leaset Squares, and Regression <a name = 'OLS'></a>


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