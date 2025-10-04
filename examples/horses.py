import EWF_pybind as EWF
import numpy as np
import matplotlib.pyplot as plt

def read_demo(demo_file, gen_gap):
    demo = np.loadtxt(demo_file)
    eff_pop_size = demo[:, 0]
    times_diff_pop_size = demo[:, 1]
    times_years_pop_size = np.cumsum(np.diff(np.concatenate(([0.0], times_diff_pop_size))) * 2.0 * eff_pop_size * gen_gap)
    
    return eff_pop_size, times_diff_pop_size, times_years_pop_size

def adjust_demo(required_demo, year_times, gen_gap):
    rebaser = year_times[0]
    year_times -= rebaser
    demo_yrs = np.array([required_demo[0][i] - rebaser for i in range(np.shape(required_demo[0])[0])])
    demo_yrs = np.concatenate(([0.0], demo_yrs[:-1]))
    demo_pop = required_demo[1]
    demo_diff = np.concatenate(([0.0], np.array([demo_yrs[i] / (2.0 * demo_pop[-i] * gen_gap) for i in range(1,np.shape(demo_yrs)[0])])))
    demo_pop = np.flip(demo_pop)

    return year_times, np.vstack((demo_yrs, demo_pop, demo_diff))

def convert_years_to_diffunits(ref_times, times, eff_pop_size, gen_gap):
    diff_times = [0.0]
    max_ind = 0
    for i in range(len(times) - 1):
        t1, t2 = times[i], times[i+1]
        # Find all ref_times strictly inside (t1, t2)
        internal_changepoints = ref_times[(ref_times > t1) & (ref_times < t2)]
        # Partition the interval [t1, t2] at these changepoints
        subintervals = np.concatenate(([t1], internal_changepoints, [t2]))
        sum_contrib = 0.0
        for j in range(len(subintervals) - 1):
            a, b = subintervals[j], subintervals[j+1]
            idx = np.searchsorted(ref_times, a, side='right') - 1
            if idx < 0:
                idx = 0  # before all changepoints, use earliest epoch
            elif idx >= len(eff_pop_size):
                idx = len(eff_pop_size) - 1
            sum_contrib += (b - a) / (2.0 * eff_pop_size[idx] * gen_gap)
            max_ind = max(max_ind, idx)
        diff_times.append(diff_times[-1] + sum_contrib)

    diff_times = np.array(diff_times)

    return diff_times, max_ind

def convert_diffunits_to_years(ref_times, times, eff_pop_size, gen_gap, start_yr):
    years = [start_yr]
    max_ind = 0
    counter = 0
    ischngpt = np.zeros_like(times)
    for i in range(len(times) - 1):
        t1, t2 = times[i], times[i+1]
        internal_changepoints = ref_times[(ref_times > t1) & (ref_times < t2)]
        subintervals = np.concatenate(([t1], internal_changepoints, [t2]))
        sum_contrib = 0.0
        for j in range(len(subintervals) - 1):
            a, b = subintervals[j], subintervals[j+1]
            idx = np.searchsorted(ref_times, a, side='right') - 1
            if idx < 0:
                idx = 0  # before all changepoints, use earliest epoch
            elif idx >= len(eff_pop_size):
                idx = len(eff_pop_size) - 1
            
            sum_contrib += (b - a) * (2.0 * eff_pop_size[idx] * gen_gap)
            max_ind = max(max_ind, idx)
            if a != t1 or b != t2:
                ischngpt[i+1] = 1
        years.append(years[-1] + sum_contrib)
        counter += 1

    years = np.array(years)

    return years, max_ind, np.cumsum(ischngpt)

def read_data(data_file, eff_pop_size, times_years_pop_size, times_diff_pop_size, mut_rate, sel_rate, gen_gap):
    data = np.loadtxt(data_file)
    times = -data[:, 0]
    year_times = times
    ref_times = times_years_pop_size

    index = 0
    while year_times[0] < ref_times[index]:
        index += 1
    required_demo = np.array([ref_times[:index+1], eff_pop_size[:index+1], times_diff_pop_size[:index+1]])
    year_times, adj_demo = adjust_demo(required_demo, year_times, gen_gap)
    
    diff_times, max_ind = convert_years_to_diffunits(adj_demo[0], year_times, adj_demo[1], gen_gap)
    
    changepoints = adj_demo[2]
    mut_vecs = np.array([[mut_rate * 2.0 * eff_pop_size[i], mut_rate * 2.0 * eff_pop_size[i]] for i in range(max_ind,-1,-1)])
    sel_vec = np.array([sel_rate * 2.0 * eff_pop_size[i] for i in range(max_ind,-1,-1)])
    freqs_ASIP = data[:, 2].astype(float) / data[:, 1].astype(float)
    freqs_MC1R = data[:, 3].astype(float) / data[:, 1].astype(float)

    return year_times, diff_times, adj_demo, freqs_ASIP, freqs_MC1R, changepoints, mut_vecs, sel_vec

mut_rate = 1.24e-9 * 5
sel_rate = 7e-4
gen_gap = 5
init_pop_size = 16e3
eff_pop_size, times_diff_pop_size, times_years_pop_size = read_demo("horses.demo", gen_gap)
year_times, diff_times, adj_demo, ASIP, MC1R, changepoints, mut_vecs, sel_vec = read_data("horse_data.data", eff_pop_size, times_years_pop_size, times_diff_pop_size, mut_rate, sel_rate, gen_gap)
print("year_times")
print(year_times)
print("diff_times")
print(diff_times)
print("ASIP")
print(ASIP)
print("MC1R")
print(MC1R)
print("changepoints")
print(changepoints)
print("mut_vecs")
print(mut_vecs)
print("sel_vec")
print(sel_vec)

times_years_pop_size = adj_demo[0]
eff_pop_size = adj_demo[1]
times_diff_pop_size = adj_demo[2]

selectionSetup = 0
dominance_parameter = np.array([0.0, 0.0])
selectionPolynomialDegree = 1
selectionCoefficients = np.array([[], []])
WF = EWF.WrightFisher(changepoints, mut_vecs, True, sel_vec, selectionSetup, dominance_parameter, selectionPolynomialDegree, selectionCoefficients)

nSim = 100
ntimes = 1000
simulate = False
if simulate:
    start_index, end_index = 0, 1
    paths = np.zeros((nSim, int((ntimes-1)*(len(diff_times)-1))))
    all_times = []
    time_counter = 0
    chngpt_counter = 1
    s_ind = 0
    for (tstart, tend) in zip(year_times[:-1], year_times[1:]):
        print("Time counter = " + str(time_counter))
        x = ASIP[start_index]
        z = ASIP[end_index]
        tvals = np.linspace(tstart, tend, ntimes)
        tvals_yrs = tvals
        tvals, _ = convert_years_to_diffunits(times_years_pop_size, tvals, eff_pop_size, gen_gap)
        if tstart == diff_times[0]:
            all_times = np.concatenate((all_times, tvals_yrs))
        else:
            all_times = np.concatenate((all_times, tvals_yrs[1:]))
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
        if tstart <= changepoints[chngpt_counter] and changepoints[chngpt_counter] < tend and chngpt_counter < len(changepoints)-1:
            print([tstart, changepoints[chngpt_counter], tend])
            chngpt_counter += 1
        
    np.savetxt("ASIP_new.txt", paths)
    np.savetxt("ASIP_time_new.txt", all_times)
else:
    paths = np.loadtxt("ASIP_new.txt")
    all_times = np.loadtxt("ASIP_time_new.txt")

paths = np.hstack((np.zeros((nSim,1)), paths))
fig, ax = plt.subplots(1,1)
for i in np.arange(np.shape(paths)[0]):
    ax.plot(all_times-20000, paths[i, :], linewidth=1.5, alpha=0.075)
ax.plot(year_times-20000, ASIP, 'x', color='black', markersize=8, markeredgewidth=2)
ax.set_ylabel("Frequency", fontsize=16)
ax.set_xlabel("Years before present", fontsize=16)
ax2 = ax.twinx()  # share x-axis, independent y-axis
# set your two line levels here:
hline1 = 16000
hline2 = 31423

# draw horizontal lines spanning the current x-limits
xmin, xmax = ax.get_xlim()
ax2.hlines(hline1, -6250, 0,
           colors='tab:red',
           linewidths=2, alpha=0.9)
ax2.hlines(hline2, xmin, -6250,
           colors='tab:red',
           linewidths=2, alpha=0.9)
ax2.set_ylim(-1000, 55000)
ymin, ymax = ax2.get_ylim()
ax2.vlines(-6250, hline1, hline2, colors='tab:red', linewidths=2, alpha=0.9)
ax2.set_ylabel(r"$N_e(t)$", fontsize=16)
ax.tick_params(axis='x', which='major', labelsize=12)
ax2.tick_params(axis='y',   which='major', labelsize=12)
ax.tick_params(axis='y', which='major', labelsize=12)
ax.set_ylim([0., .8])
plt.tight_layout()
plt.savefig("horseTrajectories.png")
plt.savefig("horseTrajectories.eps")
plt.savefig("horseTrajectories.pdf")
plt.close()