import os
import argparse
import subprocess
import shutil
from urllib.request import urlopen
from pathlib import Path
import json

def fetch_unittest_dir():
    # store the URL in url as
    # parameter for urlopen
    url = "https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/raw/casatestutils/casatestutils/component_to_test_map.json?at=refs%2Fheads%2Fmaster"
      
    # store the response of URL
    response = urlopen(url)
      
    # storing the JSON response
    # from url in data
    component_to_test_map = json.loads(response.read())
  
    if os.path.exists("unittest"):
        unittest_dir = os.listdir("unittest")
    else:
        sh_filename = "checkout_unit_dir.sh"
        bashFile = open(sh_filename, 'w')
        print("git clone --depth 1 --no-checkout https://open-bitbucket.nrao.edu/scm/casa/casatestdata.git",file = bashFile)
        print("cd casatestdata",file = bashFile)
        print("git config core.sparseCheckout true",file = bashFile)
        print("git config --global filter.lfs.required true",file = bashFile)
        print('git config --global filter.lfs.clean "git-lfs clean -- %f"',file = bashFile)
        print('git config --global filter.lfs.smudge "git-lfs smudge -- %f"',file = bashFile)
        print('git config --global filter.lfs.process "git-lfs filter-process"',file = bashFile)
        print('echo unittest/* >> .git/info/sparse-checkout',file = bashFile)
        print("git checkout master",file = bashFile)
        print("mv unittest ..",file = bashFile)
        print("cd ..",file = bashFile)
        print("rm -rf casatestdata",file = bashFile)
        bashFile.close()

        cmd = ("{} checkout_unit_dir.sh".format(os.environ['SHELL'])).split()
        subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
        os.remove(sh_filename)

        unittest_dir = os.listdir("unittest")
    return unittest_dir

def build_checkout(testnames):

    unittest_dir = fetch_unittest_dir()

    paths = []
    for test in testnames:
        for datadir in unittest_dir:
            if datadir in test:
                datasets = os.listdir("{}/{}".format("unittest",datadir))
                for dataset in datasets:
                    try:
                        path = os.readlink("{}/{}/{}".format("unittest",datadir,dataset))
                        #print(path)
                    except OSError:
                        if os.path.isdir("{}/{}/{}".format("unittest",datadir,dataset)):
                            path = "{}/{}/{}".format("unittest",datadir,dataset)
                            #print(path)
                        else:
                            raise
                    paths.append(path)
    paths = [x.split("../../")[-1] for x in paths]
    paths = [x.replace("/","",1 ) if x.startswith("/") else x for x in paths]
    paths = list(set(paths))
    paths = sorted(paths)
    #print(paths)

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
unittest/* \n"""

    substring = ""
    for path in paths:
        if path.endswith("/"):
            path = path.rstrip(path[-1])
        if path.startswith(tuple(["text","fits"])): # If in text section of casatestdata
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

def download_data(testfiles: list):
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

    cmd = ("{} checkout_unit_dir.sh".format(os.environ['SHELL'])).split()
    subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
    os.remove(sh_filename)
    
    datafile = open("datafile_list.txt","r")
    gitpaths = datafile.readlines()

    fetch_path = []
    for testfile in testfiles:
        datapaths = [gitpaths[x] for x in [i for i, x in enumerate(gitpaths) if testfile in x]]
        datapaths = [x.rstrip() for x in datapaths if not x.startswith("unittest")]
        datapaths = list(set(["/".join(x.split("/")[:-2]) if not x.startswith(tuple(["text","fits"])) else x for x in datapaths]))
        for datapath in datapaths:
            if datapath.endswith(testfile): 
                val = datapath
                fetch_path.append(val)
                break
            elif testfile in datapath:
                val = datapath
                fetch_path.append(val)
                break
    datafile.close()

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
    cmd = ("{} checkout_unit_dir.sh".format(os.environ['SHELL'])).split()
    subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
    os.remove(sh_filename)
    os.remove("datafile_list.txt")

    return

if __name__ == "__main__":
    url = "https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/raw/casatestutils/casatestutils/component_to_test_map.json?at=refs%2Fheads%2Fmaster"
      
    # store the response of URL
    response = urlopen(url)
      
    # storing the JSON response
    # from url in data
    component_to_test_map = json.loads(response.read())
    
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument('-j','--test_group',  help='Filter tests by a comma separated list of components', required=False)
    parser.add_argument('--bash', help='Generate Full Sparse Checkout Script',  action='store_true')

    args, unknownArgs = parser.parse_known_args()
    
    print(args)
    testnames = []
    if args.test_group is not None:
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
        sourceFile = open(filename, 'w')
        print(build_checkout(testnames), file = sourceFile)
        sourceFile.close()

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
            print('echo unittest/* >> .git/info/sparse-checkout',file = bashFile)
        print("git checkout master",file = bashFile)
        print("cd ..",file = bashFile)
        bashFile.close()
        print("File Saved as: {}".format(sh_filename))


    print("")

