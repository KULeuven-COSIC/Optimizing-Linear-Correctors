import numpy as np
from decimal import Decimal, getcontext, ROUND_UP, ROUND_HALF_EVEN
import re
# Plotting
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib.ticker import StrMethodFormatter
import matplotlib.legend_handler
import scienceplots

# Setting the precision for Decimal
getcontext().prec = 35 # 1766

# Entropy parameter
minimum_output_entropy_per_bit_bound = Decimal('0.999')

# Bisection parameters
tol = Decimal('1e-8')
max_bisection_iter = 110000

# Plotting parameters
ieee_one_column_inch = 3.5
ieee_two_column_inch = 7.16