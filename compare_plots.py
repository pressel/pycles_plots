import netCDF4 as nc
import pylab as plt
import numpy as np
import collections

out_path = './figs/'

def main():

    data_paths = '/Users/kpressel/GCMFixedData/'

    data_directory = collections.OrderedDict()

    data_directory['0.40x'] = data_paths + 'Stats.GCMVarying.nc'

    parameters = {}
    parameters['t_min'] = 0.0
    parameters['t_max'] = 24 * 60 * 60 * 40

    plot_ref_profiles(data_directory, parameters)
    plot_timeseries(data_directory, parameters)
    plot_profiles(data_directory, parameters)

    return

def plot_ref_profiles(data_directory, parameters):

    list_of_variables = []

    for key in data_directory.keys():

        rt_grp = nc.Dataset(data_directory[key], 'r')

        ref_grp = rt_grp['reference']

        list_of_variables += ref_grp.variables.keys()

        rt_grp.close()

    #now remove dupliate varaibles
    list_of_variables = list(set(list_of_variables))

    #now make plots
    for var_name in list_of_variables:
        print var_name
        if var_name == 't' or var_name == 'z_half' or var_name == 'z':
                pass
        else:
            plt.figure(1)
            for key in data_directory.keys():
                rt_grp = nc.Dataset(data_directory[key], 'r')
                ref_grp = rt_grp['reference']
                plt.plot(ref_grp[var_name][:], ref_grp['zp_half'][:], label=key)
                rt_grp.close()

            plt.title(var_name)
            plt.legend()
            plt.savefig(out_path + var_name + '_ref.pdf')
            plt.close()

    return


def plot_timeseries(data_directory, parameters):


    list_of_variables = []
    for key in data_directory.keys():

        rt_grp = nc.Dataset(data_directory[key], 'r')

        ts_grp = rt_grp['timeseries']

        list_of_variables += (ts_grp.variables.keys())

        rt_grp.close()

    #now remove dupliate varaibles
    list_of_variables = list(set(list_of_variables))

    #now make plots
    for var_name in list_of_variables:
        print var_name
        plt.figure(1)
        for key in data_directory.keys():
            rt_grp = nc.Dataset(data_directory[key], 'r')
            ts_grp = rt_grp['timeseries']

            if var_name == 's_int':
                plt.plot(ts_grp['t'][:]/3600.0,  ts_grp[var_name][:]/ts_grp[var_name][0], label=key)
            else:
                plt.plot(ts_grp['t'][:] / 3600.0, ts_grp[var_name][:], label=key)

            rt_grp.close()

        plt.grid()
        plt.legend()
        plt.title(var_name)
        plt.savefig(out_path + var_name + '_ts.pdf')
        plt.close()


def plot_profiles(data_directory, parameters):
    list_of_variables = []

    for key in data_directory.keys():

        rt_grp = nc.Dataset(data_directory[key], 'r')

        ps_grp = rt_grp['profiles']

        list_of_variables += (ps_grp.variables.keys())

        rt_grp.close()

    #now remove dupliate varaibles
    list_of_variables = list(set(list_of_variables))

    #now make plots
    for var_name in list_of_variables:
        print var_name
        if var_name == 't' or var_name == 'z_half' or var_name == 'z':
                pass
        else:
            plt.figure(1)
            for key in data_directory.keys():
                rt_grp = nc.Dataset(data_directory[key], 'r')
                ps_grp = rt_grp['profiles']
                t = ps_grp['t'][:]
                t_max_indx = np.where(np.abs(t - parameters['t_max']) == np.min(np.abs(t - parameters['t_max'])))[0]
                t_min_indx = np.where(np.abs(t - parameters['t_min']) == np.min(np.abs(t - parameters['t_min'])))[0]

                print t_min_indx, t_max_indx

                print np.shape(ps_grp[var_name][t_min_indx:t_max_indx,:])

                plt.plot(np.mean(ps_grp[var_name][t_min_indx:t_max_indx,:],axis=0), ps_grp['z'][:], label=key)


                rt_grp.close()
            plt.title(var_name)
            plt.legend()
            plt.savefig(out_path + var_name + '_ps.pdf')
            plt.close()

    return


if __name__ == "__main__":
    main()
