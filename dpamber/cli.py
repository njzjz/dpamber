import argparse

from .model_devi import run as model_devi_run

def run():
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers()
        parser_model_devi = subparsers.add_parser(
                "devi", help="Generating initial data for surface systems.")
        parser_model_devi.add_argument('--models', type=str, 
                                help="models")
        parser_model_devi.add_argument('--cutoff', type=float, 
                                help="cutoff")
        parser_model_devi.add_argument('--parm7_file', type=str,default="qmmm.parm7",
                        help="parm7_file")
        parser_model_devi.add_argument('--qm_region', type=str,default=":1",
                        help="qm_region")
        parser_model_devi.set_defaults(func=model_devi_run)
        args = parser.parse_args()
        args.func(args)
