"""
Usage: jwst_backgrounds <ra> <dec> <wavelength> [options]

    -h --help              show this
    -t --thresh <float>    threshold factor relative to the minimum background (default is 1.1)

Help:
    For help using this tool, please contact the jwst help desk at jwsthelp.stsci.edu

"""

from jwst_backgrounds import bg_tools
from jwst_backgrounds.docopt import docopt

def main():
    """Main CLI entrypoint."""
    opt = docopt(__doc__)
    bg_tools.get_background(float(opt['<ra>']),float(opt['<dec>']),float(opt['<wavelength>']))
