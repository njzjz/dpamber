import argparse

from .corr import run as corr_run
from .model_devi import run as model_devi_run
from .qmbuffer import run as qmwater_run
from .qmwater_corr import run as qmwater_corr_run


def run():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    parser_model_devi = subparsers.add_parser(
        "devi", help="Calculate model deviations."
    )
    parser_model_devi.add_argument("--models", type=str, nargs="+", help="Models")
    parser_model_devi.add_argument("--cutoff", type=float, help="QM/MM cutoff radius")
    parser_model_devi.add_argument(
        "--parm7_file", type=str, default="qmmm.parm7", help="AMBER parm file"
    )
    parser_model_devi.add_argument(
        "--qm_region",
        type=str,
        default=":1",
        help="QM region with AMBER mask. Add quote if containing space",
    )
    parser_model_devi.add_argument("--prefix", type=str, help="prefix to the .nc file")
    parser_model_devi.set_defaults(func=model_devi_run)

    parser_corr = subparsers.add_parser("corr", help="Generate systems for DPRc.")
    parser_corr.add_argument("--cutoff", type=float, help="QM/MM cutoff radius")
    parser_corr.add_argument(
        "--parm7_file", type=str, default="qmmm.parm7", help="AMBER parm file"
    )
    parser_corr.add_argument(
        "--qm_region",
        type=str,
        default=":1",
        help="QM region with AMBER mask. Add quote if containing space",
    )
    parser_corr.add_argument("--nc", type=str, help="AMBER coordinates (nc) file")
    parser_corr.add_argument(
        "--hl", type=str, help="Prefix to high-level files, including .mdout, .mdfrc"
    )
    parser_corr.add_argument(
        "--ll", type=str, help="Prefix to low-level files, including .mdout, .mdfrc"
    )
    parser_corr.add_argument(
        "--out",
        type=str,
        default="dataset",
        help="output directory or hdf5 file for DeePMD-kit data (default is dataset)",
    )
    parser_corr.add_argument("--suffix_mdfrc", type=str, help="suffix of mdfrc file")
    parser_corr.set_defaults(func=corr_run)

    parser_qmwater = subparsers.add_parser(
        "qmwater", help="Generate QM water parm and mdin file."
    )
    parser_qmwater.add_argument("--cutoff", type=float, help="cutoff")
    parser_qmwater.add_argument(
        "--parm7_file", type=str, default="qmmm.parm7", help="parm7_file"
    )
    parser_qmwater.add_argument("--qm_region", type=str, default=":1", help="qm_region")
    parser_qmwater.add_argument("--nc", type=str, help="nc file")
    parser_qmwater.add_argument("--hl", type=str, help="high level mdin file")
    parser_qmwater.add_argument("--ll", type=str, help="low level mdin file")
    parser_qmwater.add_argument(
        "--charge", type=str, default="dataset", help="Charge of QM region"
    )
    parser_qmwater.set_defaults(func=qmwater_run)

    parser_qmwater_corr = subparsers.add_parser(
        "qmwater_corr", help="Generate systems for QM water."
    )
    parser_qmwater_corr.add_argument("--cutoff", type=float, help="cutoff")
    parser_qmwater_corr.add_argument(
        "--qm_region", type=str, default=":1", help="qm_region"
    )
    parser_qmwater_corr.add_argument(
        "--out",
        type=str,
        default="dataset",
        help="output directory (default is dataset)",
    )
    parser_qmwater_corr.set_defaults(func=qmwater_corr_run)

    args = parser.parse_args()
    args.func(args)
