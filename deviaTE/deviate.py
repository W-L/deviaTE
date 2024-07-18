import logging

from deviaTE.reference import InputFile
import deviaTE.config




def main() -> None:
    """
    Main entry point for deviate analysis
    :return:
    """
    # define the configuration of the experiment and load resources
    conf = deviaTE.config.Config()
    # analyse all files in a loop
    for i in conf.input:
        inf = InputFile(conf=conf, infile=i)
        inf.analyse_coverage()
        inf.analyse_families()
        if conf.args.no_viz:
            continue
        else:
            inf.visualise()
    logging.info("done")



if __name__ == "__main__":
    main()




