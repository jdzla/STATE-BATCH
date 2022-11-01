from distutils.log import debug
import sys
sys.path.append('../../statebatch/')
from statebatch.QEbatch import Batch

Batch_obj = Batch('input.yaml', debug=True)

Batch_obj.prerun(make_jobscript='parallel')

