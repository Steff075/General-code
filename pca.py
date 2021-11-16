# see http://sebastianraschka.com/Articles/2015_pca_in_3_steps.html


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import math
from sklearn.preprocessing import StandardScaler


df = pd.read_csv(
    filepath_or_buffer='/homeb/steff075/phd_ecm/data-pca/ecm-model-pars.csv',
    header=1,
    sep=',')

df.columns=['height', 'gamma', 'Fe', 'xi', 'b', 'q2', 'tmax','tshift','class']
df.dropna(how="all", inplace=True) # drops the empty line at file-end

df.tail()


# split data table into data X and class labels y

X = df.iloc[:,0:8].values
y = df.iloc[:,8].values

label_dict = {1: 'Lo-source',
              2: 'Med-source',
              3: 'High-source'}

feature_dict = {0: 'Height',
                1: 'Gamma',
                2: 'Fe',
                3: 'xi',
                4: 'b',
                5: 'q2',
                6: 'tmax',
                7: 'tshift'}

with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(14, 20))
    for cnt in range(8):
        plt.subplot(6, 2, cnt+1)
        for lab in ('Lo-source', 'Hi-source'):
            plt.hist(X[y==lab, cnt],
                     label=lab,
                     bins=10,
                     alpha=0.3,)
        plt.xlabel(feature_dict[cnt])
    plt.legend(loc='upper right', fancybox=True, fontsize=9)

    plt.tight_layout()
    plt.show()


from sklearn.preprocessing import StandardScaler

X_std = StandardScaler().fit_transform(X)

#The mean vector is a d-dimensional vector where each value in this vector represents the sample mean
#of a feature column in the dataset.

mean_vec = np.mean(X_std, axis=0)
cov_mat = (X_std - mean_vec).T.dot((X_std - mean_vec)) / (X_std.shape[0]-1)
print('Covariance matrix \n%s' %cov_mat)

# Covariance Matrix

cov_mat = np.cov(X_std.T)


#Next, we perform an eigendecomposition on the covariance matrix:
eig_vals, eig_vecs = np.linalg.eig(cov_mat)

print('Eigenvectors \n%s' %eig_vecs)
print('\nEigenvalues \n%s' %eig_vals)


# Eigendecomposition of the standardized data based on the correlation matrix:

cor_mat1 = np.corrcoef(X_std.T)

eig_vals, eig_vecs = np.linalg.eig(cor_mat1)

print('Eigenvectors \n%s' %eig_vecs)
print('\nEigenvalues \n%s' %eig_vals)



# Eigendecomposition of the raw data based on the correlation matrix:

cor_mat2 = np.corrcoef(X.T)

eig_vals, eig_vecs = np.linalg.eig(cor_mat2)

print('Eigenvectors \n%s' %eig_vecs)
print('\nEigenvalues \n%s' %eig_vals)


# Singular Value Decomposition

u,s,v = np.linalg.svd(X_std.T)
u


# Selecting Principal Components


#Sorting Eigenpairs

#However, the eigenvectors only define the directions of the new axis, since they have all the same
#unit length 1, which can confirmed by the following two lines of code:

for ev in eig_vecs.T:
    np.testing.assert_array_almost_equal(1.0, np.linalg.norm(ev))
print('Everything ok!')


# rank the eigenvalues from highest to lowest in order choose the top k eigenvectors.

# Make a list of (eigenvalue, eigenvector) tuples
eig_pairs = [(np.abs(eig_vals[i]), eig_vecs[:,i]) for i in range(len(eig_vals))]

# Sort the (eigenvalue, eigenvector) tuples from high to low
eig_pairs.sort(key=lambda x: x[0], reverse=True)

# Visually confirm that the list is correctly sorted by decreasing eigenvalues
print('Eigenvalues in descending order:')
for i in eig_pairs:
    print(i[0])


#The explained variance tells us how much information (variance) can be attributed to each of
#the principal components.

tot = sum(eig_vals)
var_exp = [(i / tot)*100 for i in sorted(eig_vals, reverse=True)]
cum_var_exp = np.cumsum(var_exp)

with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(7, 5))

    plt.bar(range(8), var_exp, alpha=0.6, align='center',
            label='individual explained variance')
    plt.step(range(8), cum_var_exp, where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.legend(loc='best')
    plt.tight_layout()




#Projection Matrix

matrix_w = np.hstack((eig_pairs[0][1].reshape(8,1),
                      eig_pairs[1][1].reshape(8,1)))

print('Matrix W:\n', matrix_w)



#Projection onto the New Feature Space

Y = X_std.dot(matrix_w)

with plt.style.context('seaborn-whitegrid'):
    plt.figure(figsize=(8, 6))
    for lab, col in zip(('Hi-source', 'Med-source','Lo-source'),
                        ('red', 'green', 'blue')):
        plt.scatter(Y[y==lab, 0],
                    Y[y==lab, 1],
                    label=lab,
                    c=col)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/steff075/phd_ecm/data-pca/ecm-pca.png')
