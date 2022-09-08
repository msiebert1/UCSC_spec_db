#!/usr/bin/env python3
import UCSC_spec_db_management as spec_db 
from optparse import OptionParser

if __name__ == "__main__":
    description = "For adding spectra to UCSC database"
    usage = "%prog    \t [option] \n Recommended syntax: %prog"
    parser = OptionParser(usage=usage, description=description, version="0.1" )
    parser.add_option("-l", "--local", dest="local", action="store_true",
                      help='Add data to local (msiebert) database')
    parser.add_option("--wipe-date", dest="wipedate", action="store_true",
                      help='Delete reductions from specified date')
    parser.add_option("--update-all", dest="updateall", action="store_true",
                      help='Update database with all final reductions on ziggy')
    parser.add_option("--overwrite", dest="overwrite", action="store_true",
                      help='Overwrite existing spectra in database (meant to be a complete database refresh)')

    option, args = parser.parse_args()
    _local= option.local
    _wipedate= option.wipedate
    _updateall= option.updateall
    _overwrite= option.overwrite

    if _updateall:
        if _overwrite:
            spec_db.update_spec_database(overwrite = True)
        else:
            spec_db.update_spec_database()
    else:
        if _wipedate:
            spec_db.delete_date(_local, args[0])
        spec_db.add_final_reductions(_local)