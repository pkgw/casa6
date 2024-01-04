## Imports
import os, argparse, subprocess, json, sys
from urllib.request import urlopen
from pathlib import Path

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

## Functions
def check_for_git_lfs():
    git_string = subprocess.check_output(["git","lfs","version"], encoding="utf-8")
    if ("is not a git command" in git_string) or ("command not found" in git_string):
        print("ERROR git-lfs is not found. Please Check Configuration")
        print(git_string)
        raise Exception("git-lfs not found")

def run_script(script):
    local_shell = os.environ['SHELL']
    if os.path.exists("/bin/bash"):
        os.environ['SHELL'] = "/bin/bash"
    cmd = ("{} {}".format(os.environ['SHELL'], script)).split()
    subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
    os.environ['SHELL'] = local_shell

def download_data(testfiles: list):
    """
    Usage 
        from casatestutils import sparse_check
        sparse_check.download_data(["ANTEN_sort_hann_for_cvel_reg.ms"])
        sparse_check.download_data(["ngc5921.ms"])
        sparse_check.download_data(["gaincaltest2.ms"])
    """
    check_for_git_lfs()
    paths = []

    sh_filename = "checkout_unit_dir.sh"
    bashFile = open(sh_filename, 'w')
    print("git clone --depth 1 --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
    print("cd casatestdata",file = bashFile)
    print("git ls-tree --name-only --full-tree -r HEAD >> datafile_list.txt",file = bashFile)
    print("mv datafile_list.txt ..",file = bashFile)
    print("cd ..",file = bashFile)
    print("rm -rf casatestdata",file = bashFile)
    bashFile.close()
    run_script(sh_filename)
    os.remove(sh_filename)

    datafile = open("datafile_list.txt","r")
    gitpaths = datafile.readlines()

    fetch_path = []
    failedfind = []

    for testfile in testfiles:
        startlen = len( fetch_path) 
        datapaths = [gitpaths[x] for x in [i for i, x in enumerate(gitpaths) if testfile in x]]
        datapaths = list(set(["/".join(x.split("/")[:-2]) if not x.startswith(tuple(["text","fits"])) else x for x in datapaths]))
        datapaths = [x.strip() for x in datapaths]

        # Check For Path explicitly 1st 
        for datapath in datapaths:
            if datapath.endswith(testfile) and datapath.split("/")[-1] == testfile: 
                val = datapath
                fetch_path.append(val)
                break

        # Check For Path implicitly if not found initially 
        if len(fetch_path) == 0:
            for datapath in datapaths:
                if testfile in datapath:
                    val = datapath
                    fetch_path.append(val)
                    break
        if len(fetch_path) == startlen:
            failedfind.append(testfile)
    datafile.close()
    #sys.exit()
    if len(failedfind) > 0:
        print("Cannot Find Requested Data File(s) {}:.\nPlease Check File Names or contact Verification Team".format(",".join(failedfind)))
        return

    sh_filename = "checkout_unit_dir.sh"
    bashFile = open(sh_filename, 'w')
    print("git clone --depth 1 --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
    print("cd casatestdata",file = bashFile)
    print("git config core.sparseCheckout true",file = bashFile)
    print("git config --global filter.lfs.required true",file = bashFile)
    print('git config --global filter.lfs.clean "git-lfs clean -- %f"',file = bashFile)
    print('git config --global filter.lfs.smudge "git-lfs smudge -- %f"',file = bashFile)
    print('git config --global filter.lfs.process "git-lfs filter-process"',file = bashFile)

    for path in fetch_path:
        substring = ''
        if path.endswith("/"):
            path = path.rstrip(path[-1])
        if path.startswith(tuple(["text","fits"])): # If in text section of casatestdata
            if path.endswith("jyperk_web_api_response"):
              substring  = substring + "{}/*".format(path) 
            else:
              substring  = substring + "{}".format(path) 
        else:
            substring  = substring + "{}/*".format(path) 
        print('echo {} >> .git/info/sparse-checkout'.format(substring),file = bashFile)
    print("git checkout master",file = bashFile)
    for path in fetch_path:
        print("mv {} ..".format(path),file = bashFile)
    print("cd ..",file = bashFile)
    print("rm -rf casatestdata",file = bashFile)
    bashFile.close()
    print("Fetching ", *testfiles)
    run_script(sh_filename)
    os.remove(sh_filename)
    os.remove("datafile_list.txt")

    return

def fetch_data_dir(directory):
    check_for_git_lfs()
    sh_filename = "checkout_unit_dir.sh"
    bashFile = open(sh_filename, 'w')
    print("git clone --depth 1 --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
    print("cd casatestdata",file = bashFile)
    print("git ls-tree --name-only --full-tree -r HEAD >> datafile_list.txt",file = bashFile)
    print("mv datafile_list.txt ..",file = bashFile)
    print("cd ..",file = bashFile)
    print("rm -rf casatestdata",file = bashFile)
    bashFile.close()
    run_script(sh_filename)
    os.remove(sh_filename)
    datafile = open("datafile_list.txt","r")
    gitpaths = datafile.readlines()
    os.remove("datafile_list.txt")
    return gitpaths

def build_checkout(testnames, directory):
    check_for_git_lfs()
    testdata_dir = fetch_data_dir(directory)
    datasets = []
    paths = []
    for test in testnames:
        testname = test.split("_")[-1]
        for datadir in [x for x in testdata_dir if x.startswith('unittest')]:
            if testname in datadir:
                paths.append("/".join(datadir.split("/", 2)[:2]))
                datasets.append(datadir)
    paths = list(set(paths))
    sh_filename = "checkout_unit_dir.sh"
    bashFile = open(sh_filename, 'w')
    print("git clone --depth 1 --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
    print("cd casatestdata",file = bashFile)
    print("git config core.sparseCheckout true",file = bashFile)
    print("git config --global filter.lfs.required true",file = bashFile)
    print('git config --global filter.lfs.clean "git-lfs clean -- %f"',file = bashFile)
    print('git config --global filter.lfs.smudge "git-lfs smudge -- %f"',file = bashFile)
    print('git config --global filter.lfs.process "git-lfs filter-process"',file = bashFile)

    for path in paths:
        substring = ''
        if path.endswith("/"):
            path = path.rstrip(path[-1])
        if path.startswith(tuple(["text","fits"])): # If in text section of casatestdata
            if path.endswith("jyperk_web_api_response"):
              substring  = substring + "{}/*".format(path) 
            else:
              substring  = substring + "{}".format(path) 
        else:
            substring  = substring + "{}/*".format(path) 
        print('echo {} >> .git/info/sparse-checkout'.format(substring),file = bashFile)
    print("git checkout master",file = bashFile)
    print("cd ..",file = bashFile)
    bashFile.close()
    run_script(sh_filename)
    os.remove(sh_filename)

    for x_path in paths:
        for root, dirs, files in os.walk("casatestdata/{}".format(x_path), topdown=False):
            #print(root, dirs, files)
            for name in files:
                try:
                    path = os.readlink("casatestdata/{}/{}".format(x_path,name))
                except OSError:
                    if os.path.isdir("casatestdata/{}/{}".format(x_path,name)):
                        path = "{}/{}".format(x_path,name)
                    else:
                        #print("casatestdata/{}/{}".format(x_path,name))
                        #raise
                        pass
                paths.append(path)
    paths = [x.split("../../")[-1] for x in paths]
    paths = [x.replace("/","",1 ) if x.startswith("/") else x for x in paths]
    paths = list(set(paths))
    paths = sorted(paths)

    headstring = """
## This file will provide the instructions to allow a sparse checkout of data
##
## This file is intended to be used by piping its contents into bash in a
## git clone that has been cloned with --no-checkout, see readme.md at:
##
##  https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git
##

git config core.sparseCheckout true
cat > .git/info/sparse-checkout <<'EOF'
{}/* \n""".format(directory)

    substring = ""
    for path in paths:
        if path.endswith("/"):
            path = path.rstrip(path[-1])
        if path.startswith(tuple(["text","fits"])) or path.endswith(tuple([".txt"])): # If in text section of casatestdata
            if path.endswith("jyperk_web_api_response"):
              substring  = substring + "{}/*".format(path) + "\n"
            else:
              substring  = substring + "{}".format(path) + "\n"
        else:
            substring  = substring + "{}/*".format(path) + "\n"
        
    tailstring = """readme.md
EOF
"""
    string = headstring + substring + tailstring
    return string

## Main
if __name__ == "__main__":
    url = "https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/raw/casatestutils/casatestutils/component_to_test_map.json?at=refs%2Fheads%2Fmaster"
      
    # store the response of URL
    try:
        response = urlopen(url)
    except:
        import ssl
        context =  ssl._create_unverified_context()
        response = urlopen(url, context = context)

    # storing the JSON response
    # from url in data
    component_to_test_map = json.loads(response.read())
    
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('-j','--test_group',  help='Filter tests by a comma separated list of components', required=False)
    parser.add_argument('--bash', help='Generate Full Sparse Checkout Script',  action='store_true')

    args, unknownArgs = parser.parse_known_args()

    directory = "unittest"

    testnames = []
    if args.test_group is not None:
        if args.test_group in ['regression', "performance", "stakeholder"]:
            headstring = """
## This file will provide the instructions to allow a sparse checkout of data
##
## This file is intended to be used by piping its contents into bash in a
## git clone that has been cloned with --no-checkout, see readme.md at:
##
##  https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git
##

git config core.sparseCheckout true
cat > .git/info/sparse-checkout <<'EOF'
{}/* \n
readme.md
EOF""".format(args.test_group)
            filename = "{}-data".format("".join(args.test_group.split()))
            sourceFile = open(filename, 'w')
            print(headstring, file = sourceFile)
            sourceFile.close()
        else:
            components = args.test_group
            components = [x.strip() for x in components.split(",")]
            print("Testing Components" + str(components))
            print("")
            no_test_components = []
            for c in components:
                _isComponent = False
                component = c.strip()
                for myDict in component_to_test_map["testlist"]:
                    #print(component, myDict["testGroup"])
                    if component in myDict["testGroup"] or component in myDict["testType"]:
                        _isComponent = True
                        if (myDict["testScript"] not in testnames):
                            testnames.append(myDict["testScript"])
                if not _isComponent:
                    print("No Tests for Component: {}".format(component))
                    no_test_components.append(component)

            if len(testnames)==0:
                if len(no_test_components) > 0:
                    print("No Test Suite for Component(s): {}".format(no_test_components))
                print("Generating Suite Using Component 'default'")
                component = 'default'
                for myDict in component_to_test_map["testlist"]:
                    if component in myDict["testGroup"]:
                        _isComponent = True
                        testnames.append(myDict["testScript"])
            filename = "-".join(components) + "-data"
            filename = filename.replace(" ", "-")
            sourceFile = open(filename, 'w')
            print(build_checkout(testnames, directory), file = sourceFile)
            sourceFile.close()
            print("File Saved as: {}".format(filename))

    if args.bash:
        sh_filename = "On_demand_sparse_checkout.sh"
        bashFile = open(sh_filename, 'w')
        print("mkdir data",file = bashFile)
        print("cd data",file = bashFile)
        print("git clone --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
        print("cd casatestdata",file = bashFile)
        print("git config core.sparseCheckout true",file = bashFile)
        print("git config --global filter.lfs.required true",file = bashFile)
        print('git config --global filter.lfs.clean "git-lfs clean -- %f"',file = bashFile)
        print('git config --global filter.lfs.smudge "git-lfs smudge -- %f"',file = bashFile)
        print('git config --global filter.lfs.process "git-lfs filter-process"',file = bashFile)
        if args.test_group is not None:
            print("cp ../../{0} ./{0}".format(filename),file = bashFile)
            print("source {}".format(filename),file = bashFile)
        else:
            print('echo {}/* >> .git/info/sparse-checkout'.format(directory),file = bashFile)
        print("git checkout master",file = bashFile)
        print("cd ..",file = bashFile)
        bashFile.close()
        print("File Saved as: {}".format(sh_filename))

    if len(sys.argv) == 1:
        print("No Arguments Given")
        parser.print_help()
        sys.exit(1)

