import argparse

from .model_devi import run as model_devi_run
from .corr import run as corr_run


def run():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    parser_model_devi = subparsers.add_parser(
        "devi", help="Calculate model deviations.")
    parser_model_devi.add_argument('--models', type=str, nargs='+',
                                   help="models")
    parser_model_devi.add_argument('--cutoff', type=float,
                                   help="cutoff")
    parser_model_devi.add_argument('--parm7_file', type=str, default="qmmm.parm7",
                                   help="parm7_file")
    parser_model_devi.add_argument('--qm_region', type=str, default=":1",
                                   help="qm_region")
    parser_model_devi.set_defaults(func=model_devi_run)

    parser_corr = subparsers.add_parser(
        "corr", help="Generate systems for DPRc.")
    parser_corr.add_argument('--cutoff', type=float,
                             help="cutoff")
    parser_corr.add_argument('--parm7_file', type=str, default="qmmm.parm7",
                             help="parm7_file")
    parser_corr.add_argument('--qm_region', type=str, default=":1",
                             help="qm_region")
    parser_corr.add_argument('--nc', type=str,
                             help="nc file")
    parser_corr.add_argument('--hl', type=str,
                             help="high level file")
    parser_corr.add_argument('--ll', type=str,
                             help="low level file")
    parser_corr.add_argument('--out', type=str, default="dataset",
                             help="output directory (default is dataset)")
    parser_corr.set_defaults(func=corr_run)
    args = parser.parse_args()
    args.func(args)
