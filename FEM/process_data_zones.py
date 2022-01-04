# script for reading in soln.dat files, computing various stress measures,
# and generating a new filtered.dat file containing these different measures.
# The most important thing about this script is that we preserve the position
# of all the "ZONE I=5,..." lines for plotting in paraview. 
import numpy as np
import pandas as pd
import sys

# set display threshold
pd.set_option('chop_threshold', 0.001)

# read in file location
filename = str(sys.argv[1])

# number of files
nfiles = int(sys.argv[2])

# define yield stresses
Yfibre = 4.0
Ymatrix = Yfibre*0.02

# function for computing the cauchy stress from the second piola kirchoff stress
# and the deformation gradient. Can be to dataframe columns - more efficient than
# looping over rows
def compute_cauchy_stress(x):
    # initialise
    sigma = np.zeros((3,3))
    # read in elements of 2pk stress and deformation gradient in specific order
    S = np.array([[x[0], x[3], x[4]], [x[3], x[1], x[5]], [x[4], x[5], x[2]]])
    F = np.array([[x[6], x[9], x[11]], [x[10], x[7], x[13]], [x[12], x[14], x[8]]])

    # loop over indices
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    sigma[i, j] +=  F[i, k]*S[k, l]*F[j, l]

    return sigma


# same as abover but to rotate the cauchy stress tensor so that it's axis aligns
# with the fibre direction m
def transform_cauchy_stress(x):
    # one column contains the cauchy stress as an entire numpy array, pass in first
    sigma = x[0]
    # ... then the fibre direction
    m = np.array([x[1], x[2], x[3]])
    # define the skew-symmetric matrix required for computing R
    skewSym = np.array([[0, 0, m[0]], [0, 0, m[1]], [-m[0], -m[1], 0]])
    # define R using Rodriguezs formula
    rot = np.identity(3) + skewSym + np.dot(skewSym,skewSym)/(1+m[2])
    # initialise
    sigma_prime = np.zeros((3,3))

    # loop over indices
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    sigma_prime[i,j] += rot[i,k]*sigma[k, l]*rot[j, l]

    return sigma_prime


# compute the von-mises stress from the cauchy stress
# x is cauchy stress, don't need any other inputs
def compute_von_mises_stress(x, Y):
    return np.sqrt(0.5*((x[0,0]-x[1,1])**2 + (x[1,1]-x[2,2])**2 + (x[2,2]-x[0,0])**2
    + 6*(x[0,1]**2 + x[1,2]**2 + x[2,0]**2)))/Y


# function to compute the torsional shear stress (i.e. convert 2pk into cylindrical)
# and take the rtheta component
def compute_rtheta_stress(x):
    # x[0] is x-position, x[1] is y-position, x[2] is S00, then S11, then S01
    return x[0]*x[1]*(x[3] - x[2]) + x[4]*(x[0]**2 - x[1]**2)/(x[0]**2 + x[1]**2)


# compute Hill stress, given cauchy OR transformed cauchy stress + yield stresses
def compute_hill_stress(x, Yf, Ym, tau32, tau12):
        # constants for transverse isotropy
        F = 0.5*(1/Yf**2)
        G = 0.5*(1/Yf**2)
        H = 0.5*(2/Ym**2 - 1/Yf**2)
        L = 0.5*(1/Ym**2)
        M = 0.5*(1/Ym**2)
        N = 0.5*(1/Ym**2)
        # compute hill stress (=1 when the material yields)
        return F*(x[1,1]-x[2,2])**2 + G*(x[2,2]-x[0,0])**2+ H*(x[0,0]-x[1,1])**2
        + 2*L1*x[1,2]**2 + 2*M1*x[2,0]**2 + 2*N1*x[0,1]**2


# loop over all files in the folder
for testno in range(6):
    for solnno in range(nfiles):
        print(solnno)
        # specify the order of the columns
        cols = ["x0", "x1", "x2", "xi0", "xi1", "xi2", "S00", "S11", "S22", "S01",
        "S02", "S12", "F00", "F11", "F22", "F01", "F10", "F02", "F20", "F12", "F21",
        "m0", "m1", "m2", "I4"]
        # read in the .dat file
        df = pd.read_csv(filename + "/test" + str(testno) + "/data/soln" + str(solnno) + ".dat", names=cols,
        engine='python', delimiter=" ", header=None, index_col=False)

        # compute the cauchy stress
        df["sigma"] = df[[ "S00", "S11", "S22", "S01", "S02", "S12", "F00", "F11",
        "F22", "F01", "F10", "F02", "F20", "F12", "F21"]].apply(lambda x: compute_cauchy_stress(x), axis=1)
        # extract the elements we want
        df["sigma_22"] = df["sigma"].apply(lambda x: x[2,2])

        # compute the transformed cauchy stress
        df["sigma_prime"] = df[["sigma", "m0", "m1", "m2"]].apply(lambda x: transform_cauchy_stress(x), axis=1)
        # extract the elements we want
        df["sigma_prime22"] = df["sigma_prime"].apply(lambda x: x[2,2])

        # compute the von mises stress. Change syntax because function has multiple arguments
        df["sigma_vm"] = df.apply(lambda x: compute_von_mises_stress(x["sigma"], Yfibre), axis=1)

        # compute the Hill stress without transforming cauchy stress
        df["hill"] = df.apply(lambda x: compute_hill_stress(x["sigma"], Yfibre, Ymatrix, Ymatrix, Ymatrix), axis=1)
        # compute the Hill stress with the transformed cauchy stress
        df["hill_transformed"] = df.apply(lambda x: compute_hill_stress(x["sigma_prime"], Yfibre, Ymatrix, Ymatrix, Ymatrix), axis=1)
        # compute the torsional shear stress
        # the first two columns can have non-numeric values here because of ZONE.. etc.
        # first generate new columns for these values. DO NOT OVERWRITE THEM, we
        # need to preserve zone positions
        # df["xi0_numeric"] = pd.to_numeric(df["xi0"], errors='coerce')
        # df["xi1_numeric"] = pd.to_numeric(df["xi1"], errors='coerce')
        # df["S_rtheta"] = df[["xi0_numeric", "xi1_numeric", "S00", "S11", "S01"]].apply(lambda x: compute_rtheta_stress(x), axis=1)

        # create a new dataframe to output. specify which columns of df we want
        # output_df = df[["x0", "x1", "x2", "xi0", "xi1", "xi2", "S22", "S_rtheta",
        # "sigma_22", "sigma_prime22", "sigma_vm", "hill", "hill_transformed", "I4"]]
        output_df = df[["x0", "x1", "x2", "xi0", "xi1", "xi2", "S22","sigma_22",
        "sigma_prime22", "sigma_vm", "hill", "hill_transformed", "I4"]]
        # output to csv, should ignore nan by default.
        output_df.to_csv(filename + "/test" + str(testno) + "/data/filtered" + str(solnno)+ ".dat", sep=' ', index=False)
