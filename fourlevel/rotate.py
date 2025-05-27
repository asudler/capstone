import numpy as np
import argparse
import re

# Get filename information through command-line parsing
parser = argparse.ArgumentParser(description="Read and display the contents"
                                 + " of a file.")
parser.add_argument("-r", "--read", type=str, default="out.physical_grid.dat",
                    help="Data file name (default: out.physical_grid.dat)")
parser.add_argument("-w", "--write", type=str, default="out.polaritons.dat",
                    help="Output file (default: out.polaritons.dat)")
parser.add_argument("-i", "--input", type=str, default="in.parameters.dat",
                    help="Simulation inputfile for generating mixing angles "
                         + "(default: in.parameters.dat)")
args = parser.parse_args()

readfile = args.read
writefile = args.write
inputfile = args.input

print("Loading simulation data.")
data = np.genfromtxt(readfile, skip_header=2)

print("Extracting relevant quantities from the input file.")
with open(inputfile, 'r') as f:
    content = f.read()

# Pattern to find lines like: <parameter> : <value>
pattern = r'^(\w+)\s*:\s*([-+]?\d*\.?\d+([eE][-+]?\d+)?)'
matches = re.findall(pattern, content, re.MULTILINE)

# Store values dynamically in global variables
for match in matches:
    var_name = match[0]
    value = float(match[1])
    globals()[var_name] = value  # DANGEROUS: use with care!

print("Initializing mixing angles \\theta and \\phi, rotation matrix.")
nn = 1.6*10**7 # from Shan
print(nn)
def control_p(t):
    return (cap_omega_plus/2.*(np.tanh(rise1*(t - t_on1_pm)) 
                               - np.tanh(fall1*(t - t_off1_pm))
                               + np.tanh(rise2*(t - t_on2_pm))
                               - np.tanh(fall2*(t - t_off2_pm))))
def control_m(t):
    return (cap_omega_minus/2.*(np.tanh(rise1*(t - t_on1_pm)) 
                               - np.tanh(fall1*(t - t_off1_pm))
                               + np.tanh(rise2*(t - t_on2_pm))
                               - np.tanh(fall2*(t - t_off2_pm))))
def f_chi_p(t):
    if (t < t_phase): return 0.
    else: return chi_p
def f_chi_m(t):
    if (t < t_phase): return 0.
    else: return chi_m
def theta(t):
    return np.arctan2(xi_max*mu_alpha/np.sqrt(nn),
                      np.sqrt(control_p(t)**2 + control_m(t)**2))
def phi(t):
    return np.arctan2(control_m(t), control_p(t))
def rotate_to_physical(t):
    return np.array(
        [[np.cos(theta(t)), np.sin(theta(t)), 0],
        [-np.sqrt(nn)*np.exp(-1j*f_chi_m(t))*np.sin(theta(t))*np.sin(phi(t)),
         np.sqrt(nn)*np.exp(-1j*f_chi_m(t))*np.cos(theta(t))*np.sin(phi(t)),
         np.sqrt(nn)*np.exp(1j*f_chi_p(t))*np.cos(phi(t))], 
        [-np.sqrt(nn)*np.exp(-1j*f_chi_p(t))*np.sin(theta(t))*np.cos(phi(t)),
         np.sqrt(nn)*np.exp(-1j*f_chi_p(t))*np.cos(theta(t))*np.cos(phi(t)),
         -np.sqrt(nn)*np.exp(1j*f_chi_m(t))*np.sin(phi(t))]])
def rotate_to_polariton(t):
    return np.array(
        [[np.cos(theta(t)), 
          -1./np.sqrt(nn)*np.exp(1j*f_chi_m(t))
          *np.sin(theta(t))*np.sin(phi(t)),
          -1./np.sqrt(nn)*np.exp(1j*f_chi_p(t))
          *np.sin(theta(t))*np.cos(phi(t))],
         [np.sin(theta(t)),
          1./np.sqrt(nn)*np.exp(1j*f_chi_m(t))
          *np.cos(theta(t))*np.sin(phi(t)),
          1./np.sqrt(nn)*np.exp(1j*f_chi_p(t))
          *np.cos(theta(t))*np.cos(phi(t))],
         [0,
          1./np.sqrt(nn)*np.exp(-1j*f_chi_p(t))*np.cos(phi(t)),
          -1./np.sqrt(nn)*np.exp(-1j*f_chi_m(t))*np.sin(phi(t))]])

#print("Rotate to physical: ")
#print(rotate_to_physical(0))
#print("NumPy inverse: ")
#print(np.linalg.inv(rotate_to_physical(0)))
#print("My inverse: ")
#print(rotate_to_polariton(0))
#print("Check NumPy inverse: ")
#print(rotate_to_physical(0) @ np.linalg.inv(rotate_to_physical(0)))
#print("Check my inverse: ")
#print(rotate_to_physical(0) @ rotate_to_polariton(0))

print("Doing the rotation.")
probe_grid = (np.stack(np.split(data[:,2], nxi)) 
              + 1j*np.stack(np.split(data[:,3], nxi)))
sigma21_grid = (np.stack(np.split(data[:,12], nxi)) 
                + 1j*np.stack(np.split(data[:,13], nxi)))
sigma23_grid = (np.stack(np.split(data[:,16], nxi)) 
                + 1j*np.stack(np.split(data[:,17], nxi)))
xis = np.unique(data[:,0])
taus = np.unique(data[:,1])

output = []
for i in range(len(xis)):
    for j in range(len(taus)):
        check = rotate_to_polariton(taus[j]) @ rotate_to_physical(taus[j])
        assert (check.shape[0] == check.shape[1]) # make sure it's the identity
        assert (np.around(check, 6) == np.eye(check.shape[0])).all()
        lhs = rotate_to_polariton(taus[j]) @ np.array([probe_grid[i,j],
                                                       sigma21_grid[i,j],
                                                       sigma23_grid[i,j]])
        output.append([xis[i], taus[j], np.real(lhs[0]), np.imag(lhs[0]),
                       np.real(lhs[1]), np.imag(lhs[1]), np.real(lhs[2]),
                       np.imag(lhs[2]), theta(taus[j]), phi(taus[j])])
    output.append([])

print("Writing output to file.")
header="xi, tau, Re/Im Psi, Re/Im Phi, Re/Im Z, theta, phi"
with open(writefile, 'w') as f:
    f.write(header + '\n\n')
    for row in output:
        if not row:
            f.write("\n")
        else:
            f.write("\t".join(map(str, row)) + '\n')
