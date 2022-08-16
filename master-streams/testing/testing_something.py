import numpy as np
import matplotlib.pyplot as plt

# Sample gaussian
import ascii_info
import hdfutils
import windows_directories

rng = np.random.default_rng()

table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)

mean = table[33]['dist']
error = table[33]['edist']
with_normal = rng.normal(mean, error, 10000000)
print(np.mean(with_normal), np.std(with_normal))

# Convert to distance modulus
mean_pc = mean*1000
mu = 5*(np.log10(mean_pc) - 1)
err_pc = error*1000
err_mu = err_pc/(0.461*mean_pc)
with_mu = rng.normal(mu, err_mu, 10000000)
with_mu = 1 + (with_mu/5)
with_mu = 10**with_mu
with_mu /= 1000
print(np.mean(with_mu), np.std(with_mu))

fract_above = 0
for i in with_mu:
    if i > mean:
        fract_above += 1
fract_above /= len(with_mu)
print(fract_above)
plt.vlines(mean, 0, 0.18, colors='black')

plt.hist(with_normal, color='red', bins=50, density=True)
plt.hist(with_mu, color='green', alpha=0.5,bins=200, density=True)
plt.show()