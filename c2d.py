import sys
import pandas as pd
import numpy as np
sys.path.append('/home/yao/scripts/toolkit')
import core

structure = core.Writefile(filename=sys.argv[1])
structure.tovasp(outfile_name=sys.argv[1]+'_d',direct=True)