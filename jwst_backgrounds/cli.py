"""
Usage: 
    jwst_backgrounds [options] <ra> <dec> <wavelength> 

Options:
    --help                        show this help
    --thresh <float>              threshold factor relative to the minimum background [default: 1.1]
    --day <integer>               which day in the year for which to extract the background 
    --showsubbkgs                 show background components in the bathtub plot
    --background_file <string>    output file name for the background [default: background.txt]
    --bathtub_file <string>       output file name for the bathtub curve [default: background_versus_day.txt]

Help:
    For help using this tool, please contact the jwst help desk at jwsthelp.stsci.edu

"""

from jwst_backgrounds import jbt
from jwst_backgrounds.docopt import docopt

def main():
    """Main CLI entrypoint."""
    opt = docopt(__doc__, options_first=True)
    
    if opt['--day'] is None:
        thisday = None
    else:
        thisday = int(opt['--day'])

    jbt.get_background(float(opt['<ra>']),float(opt['<dec>']),float(opt['<wavelength>']), \
                            thresh=float(opt['--thresh']),thisday=thisday,showsubbkgs=opt['--showsubbkgs'], \
                            background_file=opt['--background_file'],bathtub_file=opt['--bathtub_file'])
