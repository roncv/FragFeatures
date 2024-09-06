"""
Console Script for FragFeatures - preparing, summarising, and analysing experiments.
"""

import argparse
import time

from FragFeatures.duck.analysis.experiment import DUckExperiment
from FragFeatures.duck.input import DUckInput
from FragFeatures.utils import timefunction


def hello(hello_who=None):
    """Hello World function"""
    if hello_who is not None:
        print(f"Hello, {hello_who}!")
    else:
        print("Hello, World!")


def prepare_duck_experiment(
        compound_selection,
        experiment_name,
        target_dir,
        all_compounds,
        verbosity
        ):
    """Prepare an experiment for DUck.

    Parameters
    ----------
    compound_selection : list | str
            Compound selection for the experiment.
    experiment_name : str
            Name of the experiment.
    target_dir : str
            Path to the Fragalysis target's directory.
    """
    start_time = time.time()

    if verbosity >= 1:
        verbose = True
    else:
        verbose = False
    if verbosity >= 2:
        verbose_l2 = True
    else:
        verbose_l2 = False

    # Prepare the input for DUck simulation
    # TODO: add verbosity options
    duck_input = DUckInput(
        compound_selection=compound_selection,
        experiment_name=experiment_name,
        target_dir=target_dir,
        all_compounds=all_compounds,
        verbose=verbose,
        verbose_l2=verbose_l2,
    )
    # print(duck_input.compound_codes)
    # duck_input.prepare_experiment()
    duck_input.prepare_experiment_prolif()

    end_time = time.time()
    print(f"\n`prepare-duck` executed in {end_time - start_time:.4f} seconds")


def summarize_duck_experiment(experiment_dir, output_dir, wqb_filename):
    """Parse the output of a DUck simulation."""
    print("Parsing the output of a DUck simulation..")

    # Add some timing
    start_time = time.time()

    # Iniitialise the DUck experiment
    duck_output = DUckExperiment(
        experiment_dir=experiment_dir, output_dir=output_dir, wqb_filename=wqb_filename
    )
    duck_output.summarise_experiment()

    # Timing
    end_time = time.time()
    print(
        f"`summarize_duck_experiment` executed in {end_time - start_time:.4f} seconds"
    )


def args_sanitation(parser, modes):
    """Sanitize the parser arguments."""
    args = parser.parse_args()

    ### HELLO ###
    # check if everything is ok
    if args.mode == "hello":
        if args.who is None:
            # This overwrites the function if condition is not met
            modes.choices["hello"].error("You didn't specify who to say hello to.")
            # print("You didn't specify who to say hello to.")
        else:
            pass


    ### PREPARE-DUCK ###
    # TODO: Type check for compound_selection
    elif args.mode == "prepare-duck":
        print(args.compound_selection, args.experiment_name, args.target_dir)
        if (
            (args.compound_selection is None and args.all_compounds is False)
            or (args.experiment_name is None)
            or (args.target_dir is None)
        ):
            # This overwrites the function if condition is not met
            modes.choices["prepare-duck"].error(
                "You didn't specify all the required arguments."
            )
            # print("You didn't specify all the required arguments.")


    ### SUMMARISE-DUCK ###
    elif args.mode == "summarise-duck":
        if args.experiment_dir is None:
            # This overwrites the function if condition is not met
            modes.choices["summarise-duck"].error(
                "You didn't specify the experiment directory."
            )
    else:
        pass

    return args


def parse_input():
    """Main FragFeatures parser, subparsers define action modes"""
    formatter = lambda prog: argparse.HelpFormatter(
        prog, max_help_position=20, indent_increment=1, width=None
    )
    parser = argparse.ArgumentParser(
        formatter_class=formatter,
        description="Open-source toolkit for extracting fragment features from protein-ligand complexes."
    )


    ### HELLO ###
    parser.set_defaults(mode=None)
    modes = parser.add_subparsers(
        title="Subcommands", help=None, metavar="                                  "
    )
    # Arguments for hello function (hello)
    hello = modes.add_parser(
        "hello",
        help="Hello function from FragFeatures (testing).",
        description='A simple "Hello World" function.',
    )
    hello.add_argument(
        "-w", "--who", type=str, default=None, help="Say hello to this person."
    )
    hello.set_defaults(mode="hello")


    ### PREPARE-DUCK ###
    # Arguments for preparing a DUck experiment (prepare-duck)
    prepare_duck = modes.add_parser(
        "prepare-duck", help="Prepare an experiment for DUck."
    )
    prepare_duck.add_argument(
        "-c",
        "--compound-selection",
        type=str,
        nargs="+",
        default=None,
        help="Compound selection for the experiment.",
    )
    prepare_duck.add_argument(
        "-e",
        "--experiment-name",
        type=str,
        default=None,
        help="Name of the experiment.",
    )
    prepare_duck.add_argument(
        "-t",
        "--target-dir",
        type=str,
        default=None,
        help="Path to the Fragalysis target's directory.",
    )
    prepare_duck.add_argument(
        "-a",
        "--all-compounds",
        action="store_true",
        default=False,
        help="Select all available compounds from a target.",
        )

    prepare_duck.add_argument(
        '-v','--verbose',
        action='count',
        default=0,
        help='Increase verbosity. Use -v for verbose and -vv for additional verbosity.'
        )

    prepare_duck.set_defaults(mode="prepare-duck")


    ### SUMMARISE-DUCK ###
    # Arguments for summarising a DUck experiment (summarise-duck)
    summarise_duck = modes.add_parser(
        "summarise-duck", help="Summarise the output of a DUck experiment."
    )
    summarise_duck.add_argument(
        "-e",
        "--experiment-dir",
        type=str,
        default=None,
        help="Path to the experiment directory.",
    )
    summarise_duck.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default="analysis",
        help="Output directory for the analysis.",
    )
    summarise_duck.add_argument(
        "-w",
        "--wqb-filename",
        type=str,
        default="wqb.txt",
        help="Filename for the WQB output.",
    )
    summarise_duck.set_defaults(mode="summarise-duck")


    ### SANITATION ###
    args = args_sanitation(parser, modes)

    return args, parser


@timefunction
def main():
    # Parse and sanitize the inputs
    args, parser = parse_input()

    ### HELLO ###
    if args.mode == "hello":
        hello(hello_who=args.who)

    ### PREPARE-DUCK ###
    elif args.mode == "prepare-duck":
        prepare_duck_experiment(
            compound_selection=args.compound_selection,
            experiment_name=args.experiment_name,
            target_dir=args.target_dir,
            all_compounds=args.all_compounds,
            verbosity=args.verbose
        )

    ### SUMMARISE-DUCK ###
    elif args.mode == "summarise-duck":
        summarize_duck_experiment(
            experiment_dir=args.experiment_dir,
            output_dir=args.output_dir,
            wqb_filename=args.wqb_filename,
        )

    ### HELP ###
    else:
        parser.print_help()



if __name__ == "__main__":
    main()
