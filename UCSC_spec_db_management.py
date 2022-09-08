#!/usr/bin/env python3
import sys
import os
import numpy as np 
import glob
import matplotlib.pyplot as plt
import sqlite3 as sq3
from astropy.io import fits
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
import copy
import UCSC_spec_db_interaction as spec_db_interact
from optparse import OptionParser


###################
# example command to add pre_reduced data to ziggy
# rsync -aPvz ./pre_reduced/* msiebert@ziggy.ucolick.org:/data2/Lick/Kast/ut190831/pre_reduced/
###################


def process_keck_file(spec):
    meta_dict = {}
    spec_fits=fits.open(spec)
    
    filename = spec.split('/')[-1]
    head=spec_fits[0].header
    meta_dict['FULL_PATH']  = spec
    meta_dict['FILENAME']  = filename
    meta_dict['OBJECT']  = head['OBJECT']
    
    flux, err = np.transpose(spec_fits[0].data)
    snr = np.round(np.median(flux/err),3)
    meta_dict['SNR'] = snr
    
    wavezero=float(head['CRVAL1'])
    wavedelt=float(head['CDELT1'])
    wave=np.arange(len(flux))*wavedelt+wavezero
    airmass=float(head['AIRMASS'])
    
    if head.get('EXPTIME', None) is None:
        if head.get('TTIME', None) is None:
            meta_dict['EXPTIME']  = float(head['ELAPTIME'])
        else:
            meta_dict['EXPTIME']  = float(head['TTIME'])
    else:
        meta_dict['EXPTIME']  = float(head['EXPTIME'])

    meta_dict['AIRMASS'] = np.round(float(head['AIRMASS']),3)
    meta_dict['MINWAVE'] = np.round(float(head['W_RANGE'].split()[0]),4)
    meta_dict['MAXWAVE'] = np.round(float(head['W_RANGE'].split()[1]),4)
    meta_dict['WAVEDELT'] = np.round(wavedelt,4)
    meta_dict['CRVAL'] = np.round(wavezero,4)
    meta_dict['DATE_OBS'] = head['DATE-OBS'].strip().split('T')[0]

    if head.get('MJD-OBS', None) is None:
        ut_date = head['DATE-OBS'].strip()
        t = Time(ut_date, format='isot')
        meta_dict['MJD'] = float(t.mjd)
    else:
        meta_dict['MJD'] = float(head['MJD-OBS'])
    meta_dict['UTC']  = head['UTC']
    meta_dict['POS_ANG'] = np.round(float(head['ROTPOSN'])+90, 2)
    meta_dict['FLUX_OBJ'] = head['FLUX_OBJ']
    
    radec = SkyCoord(head['RA'], head['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
    # meta_dict['RA_OBS'] = np.round(radec.ra.deg, 4)
    # meta_dict['DEC_OBS'] = np.round(radec.dec.deg, 4)
    meta_dict['RA'] = np.round(radec.ra.deg, 4)
    meta_dict['DEC'] = np.round(radec.dec.deg, 4)
#     apo = Observer.at_site("Keck")
#     date = head['DATE_BEG']
#     par_ang = apo.parallactic_angle(date.split('T')[0] + ' ' + date.split('T')[1], radec)
#     print (par_ang)
    
    meta_dict['OBSERVER'] = head['OBSERVER'].strip()
    meta_dict['REDUCER'] = head['REDUCER']
    meta_dict['RED_DATE'] = head['RED_DATE'].strip()
    
    meta_dict['INSTRUMENT'] = head['INSTRUME'].strip()
    meta_dict['GRATING'] = head['GRANAME'].strip()
    meta_dict['GRISM'] = head['GRISNAME'].strip()
    meta_dict['DICHROIC'] = head['DICHNAME'].strip()
    meta_dict['SLIT'] = head['SLITNAME'].strip()
    meta_dict['TELESCOPE'] = head['TELESCOP'].strip()
    
    if 'combined' in filename:
        meta_dict['COMBINED'] = 1
    else:
        meta_dict['COMBINED'] = 0
        
    if 'BAD' in filename:
        meta_dict['AP_GREATER_SEEING'] = 0
    else:
        meta_dict['AP_GREATER_SEEING'] = 1
        
    if 'arcsec' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('arcsec')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        if 'SN' in filename.split('arcsec')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'kpc' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('kpc')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'kpc'
        if 'SN' in filename.split('kpc')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'rkron' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('rkron')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 1
    else:
        meta_dict['AP_SIZE'] = -99
        meta_dict['AP_UNIT'] = 'None'
        meta_dict['AP_LOC'] = 'None'
        meta_dict['IS_KRON_RAD'] = 0
        
    return wave, flux, err, meta_dict
    

def process_lick_file(spec):
    meta_dict = {}
    spec_fits=fits.open(spec)
    
    filename = spec.split('/')[-1]
    head=spec_fits[0].header
    meta_dict['FULL_PATH']  = spec
    meta_dict['FILENAME']  = filename
    meta_dict['OBJECT']  = head['OBJECT']
    
    flux, err = np.transpose(spec_fits[0].data)
    snr = np.round(np.median(flux/err),3)
    meta_dict['SNR'] = snr
    
    wavezero=float(head['CRVAL1'])
    wavedelt=float(head['CDELT1'])
    wave=np.arange(len(flux))*wavedelt+wavezero
    airmass=float(head['AIRMASS'])
    
    if head.get('EXPTIME', None) is None:
        if head.get('TTIME', None) is None:
            meta_dict['EXPTIME']  = float(head['ELAPTIME'])
        else:
            meta_dict['EXPTIME']  = float(head['TTIME'])
    else:
        meta_dict['EXPTIME']  = float(head['EXPTIME'])
    meta_dict['AIRMASS'] = np.round(float(head['AIRMASS']),3)
    meta_dict['MINWAVE'] = np.round(float(head['W_RANGE'].split()[0]),4)
    meta_dict['MAXWAVE'] = np.round(float(head['W_RANGE'].split()[1]),4)
    meta_dict['WAVEDELT'] = np.round(wavedelt,4)
    meta_dict['CRVAL'] = np.round(wavezero,4)
    meta_dict['DATE_OBS'] = head['DATE-OBS'].strip().split('T')[0]
    ut_date = head['DATE-OBS'].strip()
    t = Time(ut_date, format='isot')
    meta_dict['MJD'] = float(t.mjd)
    meta_dict['UTC']  = ut_date.split('T')[1]
    meta_dict['POS_ANG'] = np.round(float(head['TUB']), 2) #TODO:check that tub is correct
    meta_dict['FLUX_OBJ'] = head['FLUX_OBJ']
    
    radec = SkyCoord(head['RA'], head['DEC'], frame='icrs', unit=(u.hourangle, u.deg))
    meta_dict['RA'] = np.round(radec.ra.deg, 4)
    meta_dict['DEC'] = np.round(radec.dec.deg, 4)
#     apo = Observer.at_site("Keck")
#     date = head['DATE_BEG']
#     par_ang = apo.parallactic_angle(date.split('T')[0] + ' ' + date.split('T')[1], radec)
#     print (par_ang)
    
    meta_dict['OBSERVER'] = head['OBSERVER'].strip()
    meta_dict['REDUCER'] = head['REDUCER']
    meta_dict['RED_DATE'] = head['RED_DATE'].strip()
    
    meta_dict['INSTRUMENT'] = head['VERSION'].strip()
    meta_dict['GRATING'] = head['GRATNG_N'].strip()
    meta_dict['GRISM'] = head['GRISM_N'].strip()
    meta_dict['DICHROIC'] = head['BSPLIT_N'].strip()
    meta_dict['SLIT'] = head['SLIT_N'].strip()
    meta_dict['TELESCOPE'] = head['OBSERVAT'].strip()
    
    if 'combined' in filename:
        meta_dict['COMBINED'] = 1
    else:
        meta_dict['COMBINED'] = 0
        
    if 'BAD' in filename:
        meta_dict['AP_GREATER_SEEING'] = 0
    else:
        meta_dict['AP_GREATER_SEEING'] = 1
        
    if 'arcsec' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('arcsec')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        if 'SN' in filename.split('arcsec')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'kpc' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('kpc')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'kpc'
        if 'SN' in filename.split('kpc')[1]:
            meta_dict['AP_LOC'] = 'SN'
        else:
            meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 0
    elif 'rkron' in filename:
        meta_dict['AP_SIZE'] = float(filename.split('rkron')[0].split('_')[-2])
        meta_dict['AP_UNIT'] = 'arcsec'
        meta_dict['AP_LOC'] = 'NUC'
        meta_dict['IS_KRON_RAD'] = 1
    else:
        meta_dict['AP_SIZE'] = -99
        meta_dict['AP_UNIT'] = 'None'
        meta_dict['AP_LOC'] = 'None'
        meta_dict['IS_KRON_RAD'] = 0
        
    return wave, flux, err, meta_dict


def process_soar_file(head):
    return


def setup_db(spec, new_db_name, telescope):
    
    print ('Wiping database and creating new table...')
    con = sq3.connect(new_db_name)
    cur = con.cursor()
    con.execute("""DROP TABLE IF EXISTS SPECTRA""")

    db_setup_fits = spec
    if telescope == 'keck':
        wave, flux, err, meta_dict = process_keck_file(spec)
    elif telescope == 'lick':
        wave, flux, err, meta_dict = process_lick_file(spec)
    else:
        print ('Only lick and keck implemented so far')
        return 
    print (meta_dict['FILENAME'])
    db_setup_str = """CREATE TABLE IF NOT EXISTS SPECTRA (FILENAME TEXT PRIMARY KEY, """
    type_dict = {'str': 'TEXT', 'float64': 'REAL', 'float': 'REAL', 'int': 'BIT'}
    for m in meta_dict:
        if m != 'FILENAME' and type(meta_dict[m]) != None:
            db_setup_str = db_setup_str + m + ' ' + type_dict[type(meta_dict[m]).__name__] + ', '
    db_setup_str = db_setup_str + 'RAW_FLUX BLOB, '
    db_setup_str = db_setup_str + 'RAW_ERR BLOB'
    db_setup_str = db_setup_str + ')'
    
    con.execute(db_setup_str)
    return db_setup_str

def add_spectrum_to_db(new_db_name, spec, telescope):
    con = sq3.connect(new_db_name)

    if telescope == 'keck':
        wave, flux, err, meta_dict = process_keck_file(spec)
    elif telescope == 'lick':
        wave, flux, err, meta_dict = process_lick_file(spec)
    else:
        print ('Only lick and keck implemented so far')
        return 
    flux_binary = flux.tobytes()
    err_binary = err.tobytes()
    print ('Adding:', meta_dict['FILENAME'])
    
    add_str = """INSERT OR REPLACE INTO SPECTRA( """
    meta_data = ()
    for m in meta_dict:
        add_str = add_str + m + ', '
        meta_data = meta_data + (meta_dict[m],)

    add_str = add_str + 'RAW_FLUX, '
    meta_data = meta_data + (flux_binary,)
    add_str = add_str + 'RAW_ERR)'
    meta_data = meta_data + (err_binary,)
    add_str = add_str + ' VALUES('
    for i in range(len(meta_dict)+2):
        add_str = add_str + '?, '
    add_str = add_str[0:-2] + ')'
    con.execute(add_str, meta_data)
    con.commit()


def delete_date(local, date):
    if not local:
        db_path = '/data2/UCSC_Spectral_Database/'
    else:
        db_path = '/Users/msiebert/Documents/UCSC/Research/UCSC_spec_database/'
    con = sq3.connect(db_path+'UCSC_SPEC_DATA_DEV.db')
    cur = con.cursor()
    cur.execute("DELETE FROM Spectra where FILENAME like '%{date}%'".format(date=date))
    con.commit()

def add_final_reductions(local):
    # if local this should be run in the pre_reduced directory

    path = os.getcwd()
    if 'keck' in path.lower():
        telescope = 'keck'
    elif 'lick' in path.lower():
        telescope = 'lick'
    else:
        print ('Telescope not implemented')
        return

    if not local:
        db_path = '/data2/UCSC_Spectral_Database/'
    else:
        db_path = '/Users/msiebert/Documents/UCSC/Research/UCSC_spec_database/'

    spec_files = glob.glob(path+"/final_reductions/*.fits")
    print ('Adding spectra to: ', db_path+'UCSC_SPEC_DATA_DEV.db')

    for spec in spec_files:
        add_spectrum_to_db(db_path+'UCSC_SPEC_DATA_DEV.db', spec, telescope)
    return

def update_spec_database(overwrite=False):
    db_path = '/data2/UCSC_Spectral_Database/'
    kast_path = '/data2/Lick/Kast/'
    lris_path = '/data2/Keck/Keck/LRIS'

    final_reductions = []
    for root, dirs, files in os.walk(lris_path, topdown=False):
        for name in files:
            if 'final_reductions' in os.path.join(root, name):
                # print(os.path.join(root, name))
                final_reductions.append(os.path.join(root, name))
    for root, dirs, files in os.walk(kast_path, topdown=False):
        for name in files:
            if 'final_reductions' in os.path.join(root, name):
                # print(os.path.join(root, name))
                final_reductions.append(os.path.join(root, name))
    print ('Found', len(final_reductions), 'final reductions')

    db_specs = spec_db_interact.grab_all_spec_data("SELECT * from SPECTRA", False)
    db_files = []
    for spec in db_specs:
        # print(spec.meta_dict['FILENAME'])
        db_files.append(spec.meta_dict['FILENAME'])
    if overwrite == True:
        db_files = []
        
    count_new_spec = 0
    for f in final_reductions:
        in_db = False
        for db_f in db_files:
            if db_f in f:
                in_db = True
                break
        if not in_db:
            if 'Kast' in f:
                add_spectrum_to_db(db_path+'UCSC_SPEC_DATA_DEV.db', f, 'lick')
                count_new_spec = count_new_spec + 1
            elif 'LRIS' in f:
                add_spectrum_to_db(db_path+'UCSC_SPEC_DATA_DEV.db', f, 'keck')
                count_new_spec = count_new_spec + 1
    print ('Added', count_new_spec, 'files to database')


if __name__ == "__main__":

    description = "For data model database structure changes "
    usage = "%prog    \t [option] \n Recommended syntax: %prog"
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option("-s", "--setup", dest="setup", action="store_true",
                      help='Run database setup again (will rewrite database)')

    option, args = parser.parse_args()
    _setup= option.setup

    db_path = '/data2/UCSC_Spectral_Database/'
    if _setup:
        #example file, will use header to set database schema
        example_file = '/data2/Keck/Keck/LRIS/20171117/pre_reduced/final_reductions/SN2017coa-blue-20171117_ap1.fits'
        setup_db(example_file, db_path+'UCSC_SPEC_DATA_DEV.db', 'keck')










