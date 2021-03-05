E_bv = Av/3.1 (correction for BV axis)
E_ub = slope * E_bv (correction for UB axis)

Av_val= 1.9655172413793103
Av_ngc= 0.9689655172413794
slope =  0.7432432432432426
Syaza sent Yesterday at 22:37
i used this one in my code
Syaza sent Yesterday at 22:37
v_m52 = np.array(V(phottable_plot_dist)) - Av_val   
ub_m52 =  (np.array(U(phottable_plot_dist))-np.array(B(phottable_plot_dist)))-(slope*Av_val/3.1)
bv_m52 = (np.array(B(phottable_plot_dist))-np.array(V(phottable_plot_dist)))-Av_val/3.1