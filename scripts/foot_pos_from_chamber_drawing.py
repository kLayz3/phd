#!/usr/bin/env python3


width_outer = 440
width_inner = 412

target='right' # or left
dist = [89,90,136,100]
dist_from = ['left', 'middle'] # right, outer, inner

# Thickness of box edge, along x/y axis, on only one end.
z_thickness = (width_outer - width_inner)/2

print('Target sitting at nominally at the right-side outer edge: 0.0')

dist_cumulative = 0
for dist_raw in dist:
    dist_cumulative = dist_cumulative + dist_raw

    # Distance from left inner:
    zLI = dist_cumulative - (z_thickness / 2)

    # Distance from right inner:
    zRI = width_inner - zLI

    # Distance from right outer:
    zRO = zRI + z_thickness

    # Correct it by ~1 mm accounting for FOOT non-exact placements
    zROC = zRO + 1.0
    print(f'Right outer: {zRO} | right inner: {zRI} ; nominal outer: {zROC}')


print(f'Box thickness along 1 wall: {z_thickness}')
