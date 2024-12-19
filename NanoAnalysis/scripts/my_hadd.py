import os
import subprocess

store = "root://eos.grif.fr/"

def get_store_files(direc):
    result_string = subprocess.check_output("xrdfs {} ls {}".format(store, direc), shell=True, text=True)
    result_list = result_string.split("\n")
    return [store+file for file in result_list[:-1]]

def get_files(direc):
    return [os.path.join(direc, file, "ZZ4lAnalysis.root") for file in os.listdir(direc)]

def hadd_files(args):
    out = args["filename"] if args["filename"] != "" else "ZZ4lAnalysis.root"
    direc = args["directory"]

    outfile = os.path.join(direc, out)

    if args["eos_grif"]:
        hadd_comm = "hadd -f {}{} ".format(store, outfile)
        files = get_store_files(direc)
    else:
        hadd_comm = "haddnano.py {} ".format(outfile)
        files = get_files(direc)

    for file in files:
        hadd_comm += " {}".format(file)

    breakpoint()
    os.system(hadd_comm)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Hadd arguments")
    parser.add_argument("--directory", required=True)
    parser.add_argument("--eos_grif", action="store_true")
    parser.add_argument("--filename", default="")
    args = vars(parser.parse_args())

    hadd_files(args)