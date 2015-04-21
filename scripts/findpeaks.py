#!/usr/bin/env python

import argparse
import os
import shutil
import sys
import textwrap

from ImageD11 import columnfile as ic


def process_layer(args, layer):
    ## Checks for existance of directory for saving files
    pdir = '%s/%02d' % (args.output, layer)
    if os.path.exists(pdir):
        if args.force:
            shutil.rmtree(pdir)
        else:
            raise IOError('%s already exists' % pdir)
    ## Creats directory for saving files
    os.mkdir(pdir)

    ## Defines first image and last image number in layer
    first_im =  args.nstart + args.nimages * layer###
    last_im = args.nimages + first_im -1

    ## Creats median file using fable 'median.py'
    ndigits = None
    for i in range(7):
        fmtstr = '%s%0' + str(i) + 'd%s'
        first_im_name = fmtstr % (args.stem, args.nstart, args.image_format)
        if os.path.exists(first_im_name):
            ndigits = i
            break
    if ndigits is None:
        raise RuntimeError('format string problem, revisit')
    ## defines directory path for raw data/images
    ## first image, last image, step size between images
    command = ('median.py -i %s -f %i -l %i -s %i -o %s/median' %
               (first_im_name, first_im, int(args.nimages/args.medianstep),
                args.medianstep, pdir)
              )
    if args.verbose:
        print(command)
    os.system(command)

    ## Runs peaksearch.py on images using specified thresholds
    command = (
        'peaksearch.py -n %s -F %s -f %i -l %i -o %s/peaks '
        '-d %s/median.edf -p Y --ndigits=%i --OmegaOverRide '
        '-S %.2f -T %.3f' %
        (args.stem, args.image_format, first_im, last_im, pdir, pdir, ndigits,
         args.step_ome, args.start_ome)
        )

    ## Adds threshold values onto command
    for t in args.thresholds:
        command +=' -t %s' % t
    if args.verbose:
        print command
    os.system(command)

    ## removes any peak entries arising from spots containing less than N pixels
    for t in args.thresholds:
         d=ic.columnfile('%s/peaks_t%i.flt' % (pdir, t))
         d.removerows('Number_of_pixels', range(args.pixels))
         d.writefile('%s/peaks_min%d_t%i.flt' % (pdir, args.pixels, t))


def main(args):
    ## Define and create directory for saving processed data
    if args.output is None:
        if args.experiment:
            args.output = '%s_foundpeaks' % args.experiment
        else:
            args.output = 'foundpeaks'
    ## Create directory
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    ## Provide default thresholds if not given by user
    if args.thresholds is None:
        args.thresholds=[100,300,500,1000]

    if args.pixels < 1:
        raise RuntimeError('Pixel threshold must be > 0, %s specified' % args.pixels)

    ## Loops through specified number of layers
    for layer in range(args.nlayers):
        if args.verbose:
            print("\nProcessing layer %d of %d" % (layer+1, args.nlayers))
        process_layer(args, layer)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='findpeaks',
        description='Finds peaks in area detector image series.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
            examples:

            $ findpeaks FF_data/slice/Ti-7_14_img 1 1440 -179.875 0.25 -n 11 -x FF
            $ findpeaks FF_data/silicon/raw/a 20 732 171.84 -0.48 -x FF
            $ findpeaks NF_data/a 20 751 171.36 -0.48 -x NF
            ''')
        )
    parser.add_argument(
        'stem', type=str, help='image file name stem (path/basename)'
        )
    parser.add_argument(
        'nstart', type=int, help='file number of first image'
        )
    parser.add_argument(
        'nimages', type=int, help='number of images per layer'
        )
    parser.add_argument(
        'start_ome', type=float, help='starting angle'
        )
    parser.add_argument(
        'step_ome', type=float, help='omega step size'
        )
    parser.add_argument(
        '-n', '--nlayers', type=int, default=1, help='number of layers'
        )
    parser.add_argument(
        '-o', '--output', type=str, help='output directory'
        )
    parser.add_argument(
        '-f', '--force', default=False, action='store_true',
        help='force removal of existing layer directory'
        )
    parser.add_argument(
        '-t', '--thresholds', action='append', type=int,
        help='threshold value(s) for peaksearcher'
        )
    parser.add_argument(
        '-p', '--pixels', type=int, default=10,
        help='minimum number of pixels threshold for a peak to be a peak'
        )
    parser.add_argument(
        '-x', '--experiment', type=str, choices=['nf', 'ff', 'NF', 'FF'],
        default='', help='near-field or far-field experiment'
        )
    parser.add_argument(
        '-F', '--image_format', type=str, choices=['.edf', '.tif'],
        default='.tif', help='format of input files'
        )
    parser.add_argument(
        '-m', '--medianstep', type=int, default=10,
        help='add every mth image to generate a median background image'
        )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report progress in terminal'
        )
    args = parser.parse_args()
    main(args)
