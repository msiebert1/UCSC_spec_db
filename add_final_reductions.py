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

    option, args = parser.parse_args()
    _local= option.local
    _wipedate= option.wipedate
    
    if _wipedate:
        spec_db.delete_date(_local, args[0])
    spec_db.add_final_reductions(_local)