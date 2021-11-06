import numpy as np
import pandas as pd
import os
np.set_printoptions(precision=9) # for that 1.000000167

# Define grid spacings and a placeholder for calculated values
spacings = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
spacingvalues = []

# Calculate for each spacing the three different approximations (forward, backward, centered/centred)
for spacing in spacings:
    # Set up the grid (for ease I'm using odd number of cells- centred on 0).
    # Note: Minor Rounding Error introduced.
    n_cells = 9
    n_middle = np.rint(n_cells/2 - 1e-5, out=np.zeros(1, int), casting='unsafe')[0] # middle element index (0)
    max_cell = n_cells - 1 # pythonic maximum cell index
    n_right = np.rint(n_cells/2 - 1e-5, out=np.zeros(1, int), casting='unsafe')[0] # number of spacings to the right
    grid = np.linspace(-1*n_right*spacing, +1*n_right*spacing, n_cells)
    grid = np.exp(grid) # set the values on the grid

    # Forward Method (periodic boundary conditions)
    forwrd_gradient = np.zeros(shape=np.shape(grid))
    for num in range(n_cells):
        try:
            gradient = (grid[num+1] - grid[num])/spacing
            forwrd_gradient[num] = gradient
        # Hit a boundary (the upper, since this is a forward method.)
        except:
            gradient = (grid[0] - grid[num]) / spacing
            forwrd_gradient[num] = gradient

    # Backward Method (periodic boundary conditions)
    backwd_gradient = np.zeros(shape=np.shape(grid))
    for num in range(n_cells):
        try:
            gradient = (grid[num] - grid[num-1])/spacing
            backwd_gradient[num] = gradient
        # Hit a boundary (the lower, since this is a backward method.)
        except:
            gradient = (grid[num] - grid[max_cell]) / spacing
            backwd_gradient[num] = gradient

    # Centred Method (periodic boundary conditions)
    """
    Two options for periodicity: numpy rolling, or numpy tiling. 
    - roll is better for memory/etc- try/except loops would be involved though- no lengthy arrays
    - tiling is what I'll use, since there's no need for loops at all + faster to code 
    """
    centred_gradient = np.zeros(shape=np.shape(grid))
    tiled_grid = np.tile(grid, 3)
    for num in range(n_cells):
        gradient = (tiled_grid[num+n_cells+1] - tiled_grid[num+n_cells-1])/(2*spacing)
        centred_gradient[num] = gradient

    # Select the middle element (i.e.
    forwrd_val = forwrd_gradient[n_middle]
    backwd_val = backwd_gradient[n_middle]
    centre_val = centred_gradient[n_middle]
    vals = [forwrd_val, backwd_val, centre_val]

    # Append
    spacingvalues.append(vals)

# Formatting stuff for export to a table
spacingvalues = np.array(spacingvalues)
df = pd.DataFrame(data=spacingvalues,
                  index=spacings,
                  columns=["Forward","Backward","Centred"])
dir = os.getcwd()
df.to_latex(buf=dir + "\\boden2.1.tex",
            label="Table 2.1, Bodenheimer, Repeated",
            float_format='%11.9f')

# Create and format out a table for the "distance" to the true value
distancevalues = spacingvalues
for i in range(np.shape(spacingvalues)[0]):
    for j in range(np.shape(spacingvalues)[1]):
        distancevalues[i,j] = np.abs(distancevalues[i,j] - 1)
df = pd.DataFrame(data=distancevalues,
                  index=spacings,
                  columns=["Forward","Backward","Centred"])
df.to_latex(buf=dir + "\\boden2.1_Difference.tex",
            label="Table 2.1, Bodenheimer, Repeated, Error Size",
            float_format='%11.9f')



