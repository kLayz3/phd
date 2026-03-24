#### `tpc_draw_track.C`

- File Name
- `OptionE o = 0`
  - if 1 then will do extra fitting to the difference spectrum. **Mostly not working**
- `int I = -1`
  - specify TPC ID to draw only for specific ID
- `uint64_t max_events = -1`
Will draw TPC tracking through S2, tracked by `_take[N][2]` array specifying which measurements
go into the consideration.

#### `foot_eta_corr.C`
- File Name
- Foot Index
- `{a,b}` cut -> can be left open, will prompt for cut

#### `foot_eta_corr_fit.C`
- File Name
- Foot Index
- `{a,b}` SCI21 cut -> can be left open, will prompt for cut
  - Always taken from SCI21 !

Unlike first macro, this one won't take `{E_Sum, CoG}` but rather the fitted value
for the cluster.



