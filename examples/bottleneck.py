import numpy as np
import EWF_pybind as EWF

def build_demography_map(demography, gen_time=10.0):
    """
    demography: array (M,2) with rows [Ne, t_du_in_past]
    Interpreted as: Ne is constant from the previous breakpoint up to time t_du.
    Returns:
      T_du : sorted breakpoints in diffusion units, including 0 as present
      N_seg: Ne per interval [T_du[i], T_du[i+1])  (length = len(T_du)-1)
      Y_yrs: cumulative years at each T_du breakpoint (same len as T_du)
    """
    dem = np.array(demography, dtype=float)
    N_raw = dem[:, 0]
    T_raw = np.abs(dem[:, 1])  # make positive in the past

    # ensure present (0) exists as first breakpoint
    if np.all(T_raw > 0):
        T_raw = np.concatenate(([0.0], T_raw))
        # the first interval [0, T_raw[1]) uses the N attached to T_raw[1],
        # so the N we attach to T=0 is never used; set to a placeholder
        N_raw = np.concatenate(([N_raw[0]], N_raw))

    # sort by time
    order = np.argsort(T_raw)
    T_sorted = T_raw[order]
    N_sorted = N_raw[order]

    # unique times
    T_du, inv_idx = np.unique(T_sorted, return_index=True)
    N_at_T = N_sorted[inv_idx]

    if len(T_du) < 2:
        raise ValueError("Demography needs at least one breakpoint beyond 0.")

    # Interval i is [T_du[i], T_du[i+1]) with Ne = N_at_T[i+1]
    N_seg = N_at_T[1:]                       # length K-1
    delta_t = np.diff(T_du)                  # length K-1
    delta_y = 2.0 * N_seg * gen_time * delta_t
    Y_yrs = np.concatenate(([0.0], np.cumsum(delta_y)))

    return T_du, N_seg, Y_yrs

def years_to_du(years_bp, T_du, N_seg, Y_yrs, gen_time=10.0):
    """
    Invert years->diffusion units using the same semantics as above.
    On interval i (where Y_yrs[i] <= y < Y_yrs[i+1]):
       y = Y_yrs[i] + 2*N_seg[i]*gen_time * (t - T_du[i])
       => t = T_du[i] + (y - Y_yrs[i]) / (2*N_seg[i]*gen_time)
    """
    y = np.asarray(years_bp, dtype=float)
    # find interval
    i = np.searchsorted(Y_yrs, y, side='right') - 1
    i = np.clip(i, 0, len(N_seg) - 1)

    denom = 2.0 * N_seg[i] * gen_time
    t = T_du[i] + (y - Y_yrs[i]) / denom
    return t

def du_to_years(t_du, T_du, N_seg, Y_yrs, gen_time=10.0):
    """
    Convert diffusion time(s) -> years before present using the same
    piecewise-constant Ne mapping defined by build_demography_map.

    On interval i (T_du[i] <= t < T_du[i+1]):
        y = Y_yrs[i] + 2 * N_seg[i] * gen_time * (t - T_du[i])

    Parameters
    ----------
    t_du : float or array-like
        Diffusion time(s) measured from present (>= 0).
    T_du : 1D array
        Breakpoint times in diffusion units (increasing, with T_du[0] = 0).
    N_seg : 1D array
        Effective sizes per interval [T_du[i], T_du[i+1]) (len = len(T_du)-1).
    Y_yrs : 1D array
        Cumulative years at each breakpoint (same length as T_du).
    gen_time : float
        Generation time in years.

    Returns
    -------
    y : ndarray
        Years before present corresponding to t_du (same shape as t_du).
    """
    t = np.asarray(t_du, dtype=float)

    # find interval index i with T_du[i] <= t < T_du[i+1]
    i = np.searchsorted(T_du, t, side='right') - 1
    i = np.clip(i, 0, len(N_seg) - 1)

    y = Y_yrs[i] + 2.0 * N_seg[i] * gen_time * (t - T_du[i])
    return y

if __name__ == "__main__":
    mut_rate = 2.5e-8 
    sel_rate = 3e-4
    selectionCoefficients = np.array([-0.09375, 0.6875, -1.5, 1])
    gen_time = 10
    selectionSetup = 2
    dominance_parameter = 0.0
    selectionPolynomialDegree = 3
    demography = np.loadtxt("bottleneck.demo")

    observations = np.loadtxt("bottleneck.txt")
    
    nSim = 100
    ntimes = 100
    run_bottleneck = True
    if run_bottleneck:
        T_du, N_seg, Y_yrs = build_demography_map(demography, gen_time=10.0)
        t_du_abs = years_to_du(observations[:,0], T_du, N_seg, Y_yrs, gen_time=10.0)
        t0 = t_du_abs[0]
        change_times_since_first = T_du[T_du >= t0] - t0

        obs_du = np.column_stack([t_du_abs, observations[:,1]])
        obs_du_from_first = np.column_stack([np.abs(t_du_abs - t0), observations[:,1]])

        changepoints = (np.flip((demography[:, 1])) + t0)
        changepoints[0] = 0.

        mut_vec = np.array([np.array([2.0 * Ne * mut_rate, 2.0 * Ne * mut_rate]) for Ne in np.flip(N_seg)])
        sel_vec = 2.0 * np.flip(N_seg) * sel_rate

        WF = EWF.WrightFisher(changepoints, mut_vec, True, sel_vec, selectionSetup, dominance_parameter, selectionPolynomialDegree, selectionCoefficients)
        
        simulate_bottleneck = True
        if simulate_bottleneck:
            start_index, end_index = 0, 1
            paths = np.zeros((nSim, int((ntimes-1)*(len(t_du_abs)-1))))
            all_times = []
            time_counter = 0
            chngpt_counter = 1
            s_ind = 0
            for (tstart, tend, x, z) in zip(obs_du_from_first[:-1, 0], obs_du_from_first[1:, 0], obs_du[:-1, 1], obs_du[1:, 1]):
                print("Time counter = " + str(time_counter))
                tvals = np.linspace(tstart, tend, ntimes)
                if tstart == t_du_abs[0]:
                    all_times = np.concatenate((all_times, tvals))
                else:
                    all_times = np.concatenate((all_times, tvals[1:]))
                for s in tvals[1:-1]:
                    if s == tvals[1]:
                        WF.BridgeDiffusionRunner(nSim, x, z, tvals[0], tvals[-1], s, False, "tmp.txt", False)
                        paths[:, s_ind] = np.loadtxt("tmp.txt")
                    else:
                        x_ind = 0
                        for xv in paths[:, s_ind-1]:
                            WF.BridgeDiffusionRunner(1, xv, z, prev_s, tvals[-1], s, False, "tmp.txt", False)
                            val = np.loadtxt("tmp.txt")
                            paths[x_ind, s_ind] = val
                            x_ind += 1
                    s_ind += 1
                    prev_s = s
                paths[:, s_ind] = np.ones(np.shape(paths[:, s_ind]))*z
                s_ind += 1
                start_index += 1
                end_index += 1
                time_counter += 1
                
            np.savetxt("bottleneck_sim_paths.txt", paths)
            np.savetxt("bottleneck_sim_times.txt", all_times)
        else:
            paths = np.loadtxt("bottleneck_sim_paths.txt")
            all_times = np.loadtxt("bottleneck_sim_times.txt")
        all_times = np.insert(all_times, 0, 0.)
        yr_times = -du_to_years(np.abs(all_times - 1.5), T_du, N_seg, Y_yrs, gen_time)
        paths = np.hstack((np.zeros((nSim,1)), paths))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        for i in np.arange(np.shape(paths)[0]):
            ax.plot(yr_times, paths[i, :], linewidth=1.5, alpha=0.075)
        ax.plot(-observations[:, 0], observations[:, 1], 'x', color='black', markersize=8, markeredgewidth=2)
        ax.set_ylabel("Frequency", fontsize=16)
        ax.set_xlabel("Years before present", fontsize=16)
        ax2 = ax.twinx() 
        hline1 = 5000
        hline2 = 40000
        chngpt_yr = -du_to_years(np.abs(changepoints - 1.5), T_du, N_seg, Y_yrs, gen_time) 
        ax2.hlines([hline1, hline1, hline1], [chngpt_yr[0], chngpt_yr[2], chngpt_yr[4]], [chngpt_yr[1], chngpt_yr[3], chngpt_yr[5]],
                colors='tab:red',
                linewidths=2, alpha=0.9)
        ax2.hlines([hline2, hline2, hline2], [chngpt_yr[1], chngpt_yr[3], chngpt_yr[5]], [chngpt_yr[2], chngpt_yr[4], -observations[-1, 0]],
                colors='tab:red',
                linewidths=2, alpha=0.9)
        ax2.set_ylim(-1000, 55000)
        ax2.yaxis.set_ticks([5000,40000])
        ax.set_ylim(0., 1.)
        ymin, ymax = ax2.get_ylim()
        ax2.set_ylabel(r"$N_e(t)$", fontsize=16)
        ax.tick_params(axis='x', which='major', labelsize=12)
        ax2.tick_params(axis='y',   which='major', labelsize=12)
        ax.tick_params(axis='y', which='major', labelsize=12)
        plt.tight_layout()
        plt.savefig("bottleneck.pdf")
        plt.close()

    run_const40k = True
    if run_const40k:
        t_du_abs = observations[:, 0] / (2. * gen_time * 40000.)
        t0 = t_du_abs[0]
        change_times_since_first = t_du_abs - t0

        obs_du = np.column_stack([t_du_abs, observations[:,1]])
        obs_du_from_first = np.column_stack([np.abs(t_du_abs - t0), observations[:,1]])

        changepoints = np.array([0.])

        mut_vec_const40k = np.array([[2.0 * 40000. * mut_rate, 2.0 * 40000. * mut_rate]])
        sel_vec_const40k = np.array([2. * 40000. * sel_rate])
        
        WF_const40k = EWF.WrightFisher(changepoints, mut_vec_const40k, True, sel_vec_const40k, selectionSetup, dominance_parameter, selectionPolynomialDegree, selectionCoefficients)

        simulate_const = True
        if simulate_const:
            start_index, end_index = 0, 1
            paths = np.zeros((nSim, int((ntimes-1)*(len(t_du_abs)-1))))
            all_times = []
            time_counter = 0
            chngpt_counter = 1
            s_ind = 0
            for (tstart, tend, x, z) in zip(obs_du_from_first[:-1, 0], obs_du_from_first[1:, 0], obs_du[:-1, 1], obs_du[1:, 1]):
                print("Time counter = " + str(time_counter))
                tvals = np.linspace(tstart, tend, ntimes)
                if tstart == t_du_abs[0]:
                    all_times = np.concatenate((all_times, tvals))
                else:
                    all_times = np.concatenate((all_times, tvals[1:]))
                for s in tvals[1:-1]:
                    if s == tvals[1]:
                        WF_const40k.BridgeDiffusionRunner(nSim, x, z, tvals[0], tvals[-1], s, False, "tmp.txt", False)
                        paths[:, s_ind] = np.loadtxt("tmp.txt")
                    else:
                        x_ind = 0
                        for xv in paths[:, s_ind-1]:
                            WF_const40k.BridgeDiffusionRunner(1, xv, z, prev_s, tvals[-1], s, False, "tmp.txt", False)
                            val = np.loadtxt("tmp.txt")
                            paths[x_ind, s_ind] = val
                            x_ind += 1
                    s_ind += 1
                    prev_s = s
                paths[:, s_ind] = np.ones(np.shape(paths[:, s_ind]))*z
                s_ind += 1
                start_index += 1
                end_index += 1
                time_counter += 1
                
            np.savetxt("const40k_sim_paths.txt", paths)
            np.savetxt("const40k_sim_times.txt", all_times)
        else:
            paths = np.loadtxt("const40k_sim_paths.txt")
            all_times = np.loadtxt("const40k_sim_times.txt")
        all_times = np.insert(all_times, 0, 0.)
        yr_times = 2. * gen_time * 40000. * (all_times - all_times[0])
        yr_times -= observations[0, 0]
        paths = np.hstack((np.zeros((nSim,1)), paths))
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        for i in np.arange(np.shape(paths)[0]):
            ax.plot(yr_times, paths[i, :], linewidth=1.5, alpha=0.075)
        ax.plot(-observations[:, 0], observations[:, 1], 'x', color='black', markersize=8, markeredgewidth=2)
        ax.set_ylabel("Frequency", fontsize=16)
        ax.set_xlabel("Years before present", fontsize=16)
        ax2 = ax.twinx() 
        hline1 = 40000
        ax2.hlines(hline1, yr_times[0], yr_times[-1],
                colors='tab:red',
                linewidths=2, alpha=0.9)
        ax2.yaxis.set_ticks([40000])
        ax.set_ylim(0., 1.)
        ax2.set_ylim(-1000, 55000)
        ymin, ymax = ax2.get_ylim()
        ax2.set_ylabel(r"$N_e(t)$", fontsize=16)
        ax.tick_params(axis='x', which='major', labelsize=12)
        ax2.tick_params(axis='y',   which='major', labelsize=12)
        ax.tick_params(axis='y', which='major', labelsize=12)
        plt.tight_layout()
        plt.savefig("const40k.pdf")
        plt.close()

    