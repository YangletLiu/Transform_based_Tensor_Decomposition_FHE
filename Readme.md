###<center>Secure Tensor Decomposition for Heterogeneous Multimedia Data in Cloud Computing</center>

####Phase I:  The fully homomorphic encryption algorithm coding in C++

1. The fully homomorphic encryption coding based on the Helib library.

   Download this repository via git

```
git clone https://github.com/shaih/HElib.
```

2. Create a new code project called TestMymain.cpp and add corresponding include file.

   This file will encrypt sub-tensors using the software library .

####Phase II: Algorthm for Computing S-tSVD

We use aerial video data as evidence of our experimental testing, which is of MPEG4 format.

1. Load aerial data from YouTube or the video folder address:

   Transform_based_Tensor_Decomposition_FHE-master/Aerial_data.mp4

2.  Run PaperCode.m to caculate secure tensor singular value decomposition

3. This code is also calculated to include: compression degree, dimensionality reduction ratio comparison with traditional svd and tSVD-slice, running time compare traditional svd.

In particular, the new tensor product (t-product) coding is load by t_product.m

####Attachment:  Presentation of experimental results coding in Matlab

on the one hand, the feature matrix and eigenvalues obtained after tensor decomposition are obtained by the last part of code PaperCode.m.

on the other hand, the code of the experimental experiment charts is in AerialData_Experiment.m



