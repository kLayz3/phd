#!/usr/bin/env python3

import numpy as np

# FOOT box has equidistance two pairs, always ~109.7 mm
# So it is always:
# |FOOT12| - 109.7mm - |FOOT23|  -????-  |FOOT45| - 109.7mm - |FOOT67|

width_outer = 435.5
pair_width = 10.15
foot_width = 0.9
dfixed = 129.9

dx = 5.5 # dx = 5.5 best found now. 
d0 = 48.0 - dx
dist12 = 136.85 + dx

dtarget = -10.1

# Distances always measured with respect to the outer edge:
# ====>>> BEAM DIRECTION ====>>>
#                   (0)    (1)                   (2)    (3)                    (4)    (5)                   (6)    (7)
# | W  |             | FOOT |                     | FOOT |                      | FOOT |                     | FOOT |
# | A  |             | PAIR |                     | PAIR |                      | PAIR |                     | PAIR |
# | L  |             |  0   |                     |  1   |                      |  2   |                     |  3   |
# | L  |             |      |                     |      |                      |      |                     |      |
#                    <------>  `pair_width`       <------> `pair_width`         <------> `pair_width`        <------> `pair_width                         
# <------------------------->  `d0`; CHANGING     <------------------------------------> `dist12`; CHANGING                       
#                    <-----------------------------------> `dfixed`             <-----------------------------------> `dfixed` 
#                                                 
#                                                                               

print('Outer box-edge: 0.0')

dist = []
dist.append( d0 - pair_width - foot_width ); # dist[0]
dist.append( d0              + foot_width ); # dist[1]
dist.append( dist[0] + dfixed - pair_width - foot_width); # dist[2]
dist.append( dist[0] + dfixed              + foot_width); # dist[3]
dist.append( dist[2] + dist12 - pair_width - foot_width); # dist[4]
dist.append( dist[2] + dist12              + foot_width); # dist[5]
dist.append( dist[4] + dfixed - pair_width - foot_width); # dist[6]
dist.append( dist[4] + dfixed              + foot_width); # dist[7]


print(f"\"det_pos\": [{', '.join(f'{x:.2f}' for x in dist)}],")
print(f'Target: {dtarget:.2f}')
