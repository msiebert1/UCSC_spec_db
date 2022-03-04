#!/usr/bin/env python3
import copy
import sqlite3 as sq3
import matplotlib
# matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
import glob
import numpy as np
import scipy.signal
from optparse import OptionParser


def basic_format(size=[8,8]):
    plt.rc('font', family='serif')
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(size[0], size[1], forward = True)
    plt.minorticks_on()
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.tick_params(
        which='major', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=20)
    plt.tick_params(
        which='minor', 
        bottom='on', 
        top='on',
        left='on',
        right='on',
        direction='in',
        length=10)

class spectrum(object):
    """A generic class to represent a spectrum and its associated metadata
    """
    def __init__(self, wavelength = None, flux = None, err = None, meta_dict=None):
        self.wavelength = wavelength
        self.flux = flux
        self.err = err
        self.meta_dict = meta_dict


def grab_all_spec_data(sql_input, _local):
    if not _local:
        db_path = '/data2/UCSC_Spectral_Database/'
    else:
        db_path = '/Users/msiebert/Documents/UCSC/Research/UCSC_spec_database/'

    # if db_file is None:
    #     db_file = glob.glob('../*.db')[0]

    db_file = db_path+'UCSC_SPEC_DATA_DEV.db'
    con = sq3.connect(db_file)
    print ("Collecting data from", db_file, "...")
    cur = con.cursor()

    spec_Array = []

    spec_table = cur.execute('PRAGMA TABLE_INFO({})'.format("SPECTRA"))
    spec_cols = [tup[1] for tup in cur.fetchall()]
    cur.execute(sql_input)
 
    spec_metadata = {}
    spec_list = []
    for row in cur:
        """retrieve spectral metadata.
        """
        for i, value in enumerate(row):
            if spec_cols[i] == 'RAW_FLUX':
                flux = np.frombuffer(value, dtype='>f8')
            elif spec_cols[i] == 'RAW_ERR':
                err = np.frombuffer(value, dtype='>f8')
            else:
                spec_metadata[spec_cols[i]] = value
        wave = np.arange(len(flux))*spec_metadata['WAVEDELT']+spec_metadata['MINWAVE']
        spec =  spectrum(wavelength=wave, flux=flux, err=err, meta_dict=spec_metadata)
        spec_list.append(copy.deepcopy(spec))
        
    for spec in spec_list:
        print(spec.meta_dict['FILENAME'])
    print (len(spec_list), 'Total Spectra found')
    return spec_list

def query_and_plot(sql_input):
    spec_list = grab_all_spec_data(sql_input,_local)
    plt.figure(figsize=[15,8])
    buff=0.3
    for i, nspec in enumerate(np.asarray(spec_list)):
        plt.plot(nspec.wavelength, nspec.flux - i*buff, drawstyle='steps-mid', label=nspec.meta_dict['FILENAME'])
        plt.fill_between(nspec.wavelength, nspec.flux - i*buff - nspec.err, nspec.flux - i*buff + nspec.err, color = 'gray')
        plt.legend(fontsize=15, loc=1)
    plt.show() 

    return spec_list


def check_date(date, _local, all_files=False):

    if all_files:
        sql_input = "SELECT * from SPECTRA where DATE_OBS like '%{date}%'".format(date=date)
    else:
        sql_input = "SELECT * from SPECTRA where DATE_OBS like '%{date}%' and FILENAME like '%combined%'".format(date=date)

    spec_list = grab_all_spec_data(sql_input, _local)
    # plt.figure(figsize=[15,8])
    # buff=0.3
    # for i, nspec in enumerate(np.asarray(spec_list)):
    #     plt.plot(nspec.wavelength, nspec.flux - i*buff, drawstyle='steps-mid', label=nspec.meta_dict['FILENAME'])
    #     plt.fill_between(nspec.wavelength, nspec.flux - i*buff - nspec.err, nspec.flux - i*buff + nspec.err, color = 'gray')
    #     plt.legend(fontsize=15, loc=1)
    # plt.show() 

    return spec_list

def plot_sn(sn_name, _local, all_files=False):

    if all_files:
        sql_input = "SELECT * from SPECTRA where OBJECT like '%{sn_name}%'".format(sn_name=sn_name)
    else:
        sql_input = "SELECT * from SPECTRA where OBJECT like '%{sn_name}%' and FILENAME like '%combined%'".format(sn_name=sn_name)

    spec_list = grab_all_spec_data(sql_input,_local)
    plt.figure(figsize=[15,8])
    buff=0.3
    for i, nspec in enumerate(np.asarray(spec_list)):
        print (nspec.meta_dict['DATE_OBS'])
        plt.plot(nspec.wavelength, nspec.flux - i*buff, drawstyle='steps-mid', label=nspec.meta_dict['FILENAME'])
        plt.fill_between(nspec.wavelength, nspec.flux - i*buff - nspec.err, nspec.flux - i*buff + nspec.err, color = 'gray')
        plt.legend(fontsize=15, loc=1)
    plt.show() 

    return spec_list

def host_line_plots(spectra, z, region = [6500, 6800], plot_region = None, legend = False, rescale=False, cont_subtract=True):
    w1 = region[0]
    w2 = region[1]
    plt.figure(figsize=[15,8])
    for spec in spectra:
        wavelength = spec.wavelength/(1.+z)
        roi = (wavelength > w1) & (wavelength < w2)
        # flux_roi_norm = spec.flux[roi]/(np.median(spec.flux[roi]))
        # err_roi_norm = spec.err[roi]/(np.median(spec.flux[roi]))

        scale=1.
        if rescale:
            # flux_roi_norm = spec.flux[roi] - (np.median(spec.flux[roi]))
            flux_norm = spec.flux - (np.median(spec.flux[roi]))
            scale = 1./np.amax(flux_norm[roi])
            flux_norm = scale*flux_norm
        elif cont_subtract:
            flux_norm = spec.flux - (np.median(spec.flux[roi]))
        else:
            flux_norm = spec.flux

        if spec.err is not None:
            err_roi_norm = scale*spec.err[roi]

        if plot_region:
            new_roi = (wavelength > plot_region[0]) & (wavelength < plot_region[1])
            plt.plot(wavelength[new_roi], flux_norm[new_roi], drawstyle='steps-mid', label=spec.meta_dict['FILENAME'])
        else:
            plt.plot(wavelength[roi], flux_norm[roi], drawstyle='steps-mid', label=spec.meta_dict['FILENAME'])
        # plt.plot(wavelength, flux_norm, drawstyle='steps-mid', label=spec.meta_dict['FILENAME'])

        # plt.fill_between(wave_roi, flux_roi_norm - err_roi_norm, flux_roi_norm + err_roi_norm, color = 'gray')


    line_dict = {'H':       ([6562.79, 4861.35, 4340.472, 4101.734], 'mediumblue'),
                 '[O III]': ([4958.911, 5006.843, 4363.210], 'magenta'),
                 '[O II]':  ([3726.032, 3728.815], 'magenta'),
                 '[N II]':  ([6548.050, 6583.460], 'darkorange'),
                 '[S II]':  ([6716.440, 6730.810], 'darkgreen'),
                 'Ca H':    ([3968.5], 'red'),
                 'Ca K':    ([3933.7], 'red'),
                 'G-band':  ([4304.4], 'gold'),
                 'Mg':      ([5175.3], 'purple'),
                 'Na ID':   ([5894.0], 'lime'),
                 '[Fe II]': ([7155.1742], 'slategray'),
                 '[Ca II]': ([7291.47, 7323.89], 'indigo'),
                 '[Ni II]': ([7377.83], 'coral')
                 }

    if plot_region:
        w1 = plot_region[0]
        w2 = plot_region[1]
    for i, line in enumerate(line_dict.keys()):
        waves = line_dict[line][0]
        color = line_dict[line][1]
        for w in waves:
            plt.axvline(w, color=color, linestyle = ':', linewidth=3)

            if w > w1 and w < w2:
                plt.text(w+5, np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), line, fontsize = 15, horizontalalignment='center', rotation=90)

    tell1 = 7600/(1.+z)
    tell2 = 7630/(1.+z)
    if tell1 > w1 and tell1 < w2:
        plt.gca().axvspan(tell1, tell2, alpha=0.2, color='red')
        plt.text((tell1+tell2)/2., np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), '$\oplus$', fontsize = 15, horizontalalignment='center')

    tell1 = 6860/(1.+z)
    tell2 = 6890/(1.+z)
    if tell1 > w1 and tell1 < w2:
        plt.gca().axvspan(tell1, tell2, alpha=0.2, color='red')
        plt.text((tell1+tell2)/2., np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), '$\oplus$', fontsize = 15, horizontalalignment='center')

    if legend:
        plt.legend(fontsize=15, loc=1)
    if plot_region:
        plt.xlim([plot_region[0],plot_region[1]])
    else:
        plt.xlim([region[0],region[1]])
    
    plt.show()

def host_line_plots_fancy(spectra, z, neb = None, region = [6500, 6800], plot_region = None, legend = False, rescale=False, cont_subtract=True):
    w1 = region[0]
    w2 = region[1]
    # plt.figure(figsize=[15,8])
    basic_format(size=[15,8])

        # plt.plot(neb.wavelength[new_roi], flux_filt, drawstyle='steps-mid', color='r')

    labels = ['Galaxy','Firefly Continuum Fit']
    for i, spec in enumerate(spectra):
        wavelength = spec.wavelength/(1.+z)
        roi = (wavelength > w1) & (wavelength < w2)
        # flux_roi_norm = spec.flux[roi]/(np.median(spec.flux[roi]))
        # err_roi_norm = spec.err[roi]/(np.median(spec.flux[roi]))

        scale=1.
        if rescale:
            # flux_roi_norm = spec.flux[roi] - (np.median(spec.flux[roi]))
            flux_norm = spec.flux - (np.median(spec.flux[roi]))
            scale = 1./np.amax(flux_norm[roi])
            flux_norm = scale*flux_norm
        elif cont_subtract:
            flux_norm = spec.flux - (np.median(spec.flux[roi]))
        else:
            flux_norm = spec.flux

        if spec.err is not None:
            err_roi_norm = scale*spec.err[roi]

        if plot_region:
            new_roi = (wavelength > plot_region[0]) & (wavelength < plot_region[1])
            plt.plot(wavelength[new_roi], flux_norm[new_roi], drawstyle='steps-mid', label=labels[i])
        else:
            plt.plot(wavelength[roi], flux_norm[roi], drawstyle='steps-mid', label=spec.meta_dict['FILENAME'])
        # plt.plot(wavelength, flux_norm, drawstyle='steps-mid', label=spec.meta_dict['FILENAME'])

        # plt.fill_between(wave_roi, flux_roi_norm - err_roi_norm, flux_roi_norm + err_roi_norm, color = 'gray')

    if neb:
        new_roi = (neb.wavelength > plot_region[0]) & (neb.wavelength < plot_region[1])
        flux_filt = scipy.signal.medfilt(neb.flux, kernel_size=201)[new_roi]
        plt.plot(neb.wavelength[new_roi], neb.flux[new_roi]-flux_filt, drawstyle='steps-mid', color='g', label='Nebular Emission')

    line_dict = {'H':       ([6562.79, 4861.35, 4340.472, 4101.734], 'mediumblue'),
                 '[O III]': ([4958.911, 5006.843, 4363.210], 'magenta'),
                 '[O II]':  ([3726.032, 3728.815], 'magenta'),
                 '[N II]':  ([6548.050, 6583.460], 'darkorange'),
                 '[S II]':  ([6716.440, 6730.810], 'darkgreen'),
                 'Ca H':    ([3968.5], 'red'),
                 'Ca K':    ([3933.7], 'red'),
                 'G-band':  ([4304.4], 'gold'),
                 'Mg':      ([5175.3], 'purple'),
                 'Na ID':   ([5894.0], 'lime'),
                 '[Fe II]': ([7155.1742], 'slategray'),
                 '[Ca II]': ([7291.47, 7323.89], 'indigo'),
                 '[Ni II]': ([7377.83], 'coral')
                 }

    if plot_region:
        w1 = plot_region[0]
        w2 = plot_region[1]
    # for i, line in enumerate(line_dict.keys()):
    #     waves = line_dict[line][0]
    #     color = line_dict[line][1]
    #     for w in waves:
    #         plt.axvline(w, color=color, linestyle = ':', linewidth=3)

    #         if w > w1 and w < w2:
    #             plt.text(w+5, np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), line, fontsize = 15, horizontalalignment='center', rotation=90)

    # tell1 = 7600/(1.+z)
    # tell2 = 7630/(1.+z)
    # if tell1 > w1 and tell1 < w2:
    #     plt.gca().axvspan(tell1, tell2, alpha=0.2, color='red')
    #     plt.text((tell1+tell2)/2., np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), '$\oplus$', fontsize = 15, horizontalalignment='center')

    # tell1 = 6860/(1.+z)
    # tell2 = 6890/(1.+z)
    # if tell1 > w1 and tell1 < w2:
    #     plt.gca().axvspan(tell1, tell2, alpha=0.2, color='red')
    #     plt.text((tell1+tell2)/2., np.amax(flux_norm[roi]) - .05*np.amax(flux_norm[roi]), '$\oplus$', fontsize = 15, horizontalalignment='center')

    # if legend:
    #     plt.legend(fontsize=15, loc=1)
    if plot_region:
        plt.xlim([plot_region[0],plot_region[1]])
    else:
        plt.xlim([region[0],region[1]])
    
    plt.ylabel(r'$F_{\lambda}$ ($10^{-15}$ ergs cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)', fontsize = 25)
    plt.xlabel(r'Rest Wavelength ($\mathrm{\AA}$)', fontsize = 25)
    plt.legend(loc=2, fontsize=20, frameon=False)
    # plt.savefig('/Users/msiebert/Documents/UCSC/Research/Foundation_Hosts/plots/gal_fit.png', dpi = 300, bbox_inches = 'tight')
    # plt.savefig('/Users/msiebert/Documents/UCSC/Research/Foundation_Hosts/plots/gal_fit.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()

if __name__ == "__main__":

    description = "For quick investigation of database contents"
    usage = "%prog    \t [option] \n Recommended syntax: %prog"
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option("-l", "--local", dest="local", action="store_true",
                      help='Add data to local (msiebert) database')
    parser.add_option("--sn", dest="sn", action="store_true",
                      help='Plot the spectra of a single target')
    parser.add_option("--date", dest="date", action="store_true",
                      help='List the spectra from a single date')

    option, args = parser.parse_args()
    _local= option.local
    _sn= option.sn
    _date= option.date


    if _sn:
        sn_name = args[0]
        spec_list = plot_sn(sn_name, _local)
    if _date:
        date = args[0]
        spec_list = check_date(date, _local)



















