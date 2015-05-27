__author__ = 'Anthony G'

import re

text = 'G1 X159.740 Y80.260 E9.31481 F1200.000'

m = re.search('E(.+?)\\s', text)
if m:
    found = m.group(1)
    E = float(found)
    print E
