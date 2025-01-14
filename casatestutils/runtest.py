########################################################################################################################
############################################            Imports            #############################################
########################################################################################################################

import os, sys, re, json, unittest, shlex
import argparse, subprocess, traceback
import shutil, datetime, platform
import socket
import xml.etree.ElementTree as ET
import signal

default_timeout = 1800
sys.path.insert(0,'')

#### PYTEST IMPORT
HAVE_PYTEST = True
try:
    import pytest
except ImportError:
    HAVE_PYTEST = False

verbose = False

# JIRA BRANCH TO CHECKOUT
JIRA_BRANCH = None

# Dry run of Tests
DRY_RUN = False

########################################################################################################################
###########################################            Functions            ############################################
########################################################################################################################

def write_conftest(filepath):
    """
    Generates and writes a pytest configuration file to the specified filepath.
    
    This function creates a pytest configuration file (`conftest.py`) with
    specific setup based on the operating system. The generated configuration
    includes necessary imports and a pytest configuration function that registers
    a `TestDescriptionPlugin`. The content differs slightly for Linux and Darwin
    (macOS) platforms.

    Args:
        filepath (str): The path where the `conftest.py` file will be written.
    
    Notes:
        - For Darwin (macOS), the configuration attempts to unregister the
          `TestDescriptionPlugin` before registering it again.
        - For Linux, the configuration directly registers the `TestDescriptionPlugin`.

    """
    platform_os='Linux'
    if platform.system() == 'Darwin':
        platform_os = 'Darwin'
    string = """
import pytest
import inspect
import os
    """
    if platform_os == 'Darwin':
        string = string + """
@pytest.mark.trylast
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin('terminalreporter')
    try:
        config.pluginmanager.unregister(TestDescriptionPlugin(terminal_reporter), 'testdescription')
    except:
        pass
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), 'testdescription')
        """
    if platform_os == 'Linux':
        string = string + """
@pytest.mark.trylast
def pytest_configure(config):
    terminal_reporter = config.pluginmanager.getplugin('terminalreporter')
    config.pluginmanager.register(TestDescriptionPlugin(terminal_reporter), 'testdescription')
        """

    string = string + """
class TestDescriptionPlugin:

    def __init__(self, terminal_reporter):
        self.terminal_reporter = terminal_reporter
        self.desc = None
        self.funcn = None

    def pytest_runtest_protocol(self, item):
        #from pprint import pprint
        #d = item.__dict__
        #pprint(d, indent=2)
        self.desc = inspect.getdoc(item.obj)
        #print(item._nodeid)
        self.funcn = item._nodeid

    @pytest.hookimpl(hookwrapper=True, tryfirst=True)
    def pytest_runtest_logstart(self, nodeid, location):
        #print("Verbosity Level: {}".format(self.terminal_reporter.verbosity))
        if self.terminal_reporter.verbosity == 0:
            yield
            self.terminal_reporter.write(f'\\n{self.funcn} \\n')
        else:
            self.terminal_reporter.write('\\n')
            yield
            if self.desc:
                    self.terminal_reporter.write(f'\\n{self.desc} \\n')
            else:
                    self.terminal_reporter.write(f'\\n')
    """
    if platform_os == 'Linux':
        string = string + """
    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        outcome = yield
        report = outcome.get_result()
        #print(dir(report))
        report.start = call.start
        report.stop = call.stop
        if report.when=='teardown':
            filepath = os.path.join(os.getcwd(),'short_summary.log')

            file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
            file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid,))
            file_obj.close()
        """

    string = string + """
    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_makereport(item, call):
        outcome = yield
        report = outcome.get_result()
        if report.when=='call':
            filepath = os.path.join(os.getcwd(),'short_summary.log')
            # write short summary to file
            file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
            file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
            file_obj.close()

            # Write not pass to Textfile
            if report.outcome != 'passed':
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
                filepath = os.path.join(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..')),'summary_of_failed.log')
                file_obj = open(filepath, 'a' if os.path.isfile(filepath) else 'w')
                file_obj.write("{} {}\\n".format(report.outcome.upper(), report.nodeid))
                file_obj.write("\\tDuration: {}s\\n".format(round(report.duration,5)))
                if report.outcome == 'failed':
                    file_obj.write("\\tMessage : {}\\n".format(report.longrepr.reprcrash.message))
                file_obj.close()
    """
    file_obj = open(filepath,'w')
    file_obj.write(string)
    file_obj.close()


def write_pytestini(filepath, testname):
    """
    Generates and writes a pytest.ini configuration file with a specified test suite name.

    This function creates a pytest.ini file containing a `[pytest]` section with
    a `junit_suite_name` entry. The test suite name is specified by the `testname`
    argument.

    Args:
        filepath (str): The path where the `pytest.ini` file will be written.
        testname (str): The name to be used for the `junit_suite_name` in the pytest configuration.

    """
    string = """
[pytest]
junit_suite_name = '{}'
    """.format(testname)

    file_obj = open(filepath,'w')
    file_obj.write(string)
    file_obj.close()

def clean_working_directory(workdir):
    """
    Removes the specified working directory and all its contents if it exists.

    Args:
        workdir (str): The path to the working directory to be cleaned.

    Notes:
        - If the directory does not exist, no action is taken.
        - Use with caution as this function will permanently delete all contents
          of the specified directory.

    """
    print("Cleaning: {}".format(workdir))
    if os.path.exists(workdir):
        shutil.rmtree(workdir)

def list_tests(local_path, test_paths=[]):
    """
    Lists and prints all unit tests from specified test paths, organizing them in a test list directory.
    Args:
        local_path (str): The local path where the 'testlist' directory will be created and managed.
        test_paths (list, optional): A list of test paths to fetch and gather tests from. Defaults to an empty list.

    Notes:
        - The function relies on `fetch_tests` and `gather_all_tests` functions to fetch and gather test files.
        - Only files that start with "test_" will be listed and printed.
        - Existing 'testlist' directory and its contents in the specified `local_path` will be removed before creating a new one.
    """
    print('Full list of unit tests')
    print('-----------------------')
    if os.path.isdir(local_path +"/testlist/"):
        shutil.rmtree(local_path +"/testlist/")
    os.makedirs(local_path +"/testlist/")
    testpaths = fetch_tests(local_path +"/testlist/", 'master', test_paths)
    for path in testpaths:
        gather_all_tests(path, local_path +"/testlist/")
    tests = sorted(os.listdir(local_path +"/testlist/"))
    for test in tests:
        if test.startswith("test_"):
            print(test)

def gather_all_tests(path, workpath):

    if sys.version_info[0] > 2:
        import pathlib
        for filename in pathlib.Path(path).rglob('test_*.py'):
            shutil.copy2(filename, workpath)

def gettests(testfile):
    '''Get the list of specific tests from the command-line
       Ex: from test_clean[test1,test3] returns [test1,test3]'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        temp = testfile[n0+1:n1]
        tests = temp.split(',')
        return tests

def getname(testfile):
    '''Get the test name from the command-line
       Ex: from test_clean[test1], returns test_clean'''
    n0 = testfile.rfind('[')
    n1 = testfile.rfind(']')
    if n0 != -1:
        return testfile[:n0]

def write_xml(name, runtime, testname, classname, fMessage, filename, result):
    """
    Generates and writes an XML file in the JUnit format to report a test case result.

    This function creates an XML structure that conforms to the JUnit report format,
    representing a test suite containing a single test case that has failed. The XML
    is then written to a specified file.

    Args:
        name (str): The name of the test suite and test case.
        runtime (float): The runtime of the test case.
        testname (str): The name of the test case.
        classname (str): The class name of the test case.
        fMessage (str): The failure message to be included in the test case result.
        filename (str): The name of the file to which the XML content will be written.
        result (str): The result of the test case (e.g., 'failure', 'success').

    Notes:
        - The function sets default values for errors, failures, skipped, tests, and time attributes.
        - The current timestamp and hostname are included in the XML.
        - The function creates a failure element with the provided failure message.
        - The generated XML is written to the specified file in binary mode.

    Example:
        write_xml('test_suite', 0.01, 'test_case', 'TestClass', 'An error occurred', 'test_results.xml', 'failure')
    """
    e = datetime.datetime.now()
    timestamp = e.strftime('%Y-%m-%dT%H:%M:%S.%f')

    data = ET.Element('testsuites')

    element1 = ET.SubElement(data, 'testsuite')
    element1.set('name', "'{}'".format(name))
    element1.set('errors', "0")
    element1.set('failures', "1")
    element1.set('skipped', "0")
    element1.set('tests', "1")
    element1.set('time', "0.01")
    element1.set('timestamp', timestamp)
    element1.set('hostname', socket.gethostname())

    s_elem1 = ET.SubElement(element1, 'testcase')
    s_elem1.set('classname', "{}.SomeClass".format(name))
    s_elem1.set('name', "{}".format(name))
    s_elem1.set('time', "0.01")

    ss_elem1 = ET.SubElement(s_elem1, 'failure')
    ss_elem1.set('message', fMessage)
    ss_elem1.text = fMessage

    b_xml = ET.tostring(data)

    with open(filename, "wb") as f:
        f.write(b_xml)

def update_xml(filename, result, name="", runtime="", testname="", classname="", fMessage=""):
    """
    Updates an existing XML file in the JUnit format or generates a new one if it does not exist.

    This function performs the following steps:
    1. Checks if the specified XML file exists.
        - If the file does not exist, it generates a failure message based on the result's return code,
          prints a message indicating the creation of a new XML file, and calls `write_xml` to create it.
    2. Parses the existing XML file.
    3. Iterates through the XML elements and updates the 'testcase' elements:
        - Modifies the 'classname' attribute to contain only the test script name.
        - Modifies the 'name' attribute to contain the test class and test name separated by a period.
    4. Writes the updated XML back to the file with UTF-8 encoding and an XML declaration.

    Args:
        filename (str): The name of the XML file to be updated or created.
        result (obj): The result object containing the return code used to generate a failure message.
        name (str, optional): The name of the test suite and test case. Defaults to an empty string.
        runtime (str, optional): The runtime of the test case. Defaults to an empty string.
        testname (str, optional): The name of the test case. Defaults to an empty string.
        classname (str, optional): The class name of the test case. Defaults to an empty string.
        fMessage (str, optional): The failure message to be included in the test case result. Defaults to an empty string.

    Notes:
        - The function relies on the `write_xml` function to generate a new XML file if it does not exist.
        - The failure message is derived from the return code of the result object if not provided.
        - The function assumes the XML structure conforms to the JUnit report format.
    
    Example:
        update_xml('test_results.xml', result, 'test_suite', '0.01', 'test_case', 'TestClass', 'An error occurred')
    """
    if not os.path.isfile(filename):
        try: fMessage = signal.strsignal(abs(result.returncode))
        except: fMessage = signal.Signals(abs(result.returncode)).name
        print("Nose File Not Generated. Generating: {}".format(filename))
        write_xml(name, runtime, testname, classname, fMessage, filename, result)

    xmlTree = ET.parse(filename)

    rootElement = xmlTree.getroot()
    for element in rootElement.iter():
        if element.tag == 'testcase':
            testname = element.attrib['name']
            testscript = element.attrib['classname'].split(".")[0]
            testclass = element.attrib['classname'].split(".")[1]
            #print(name,testscript,testclass)
            element.set("classname",testscript)
            element.set("name",'.'.join([testclass,testname]))
    xmlTree.write(filename,encoding='UTF-8',xml_declaration=True)


class casa_test:
    def __init__(self,
                 name,
                 path,
                 test_group=None,
                 test_type=None,
                 maintainer=None,
                 email=None,
                 options=None,
                 comment=None,
                 timeout = default_timeout):
        self.name = name
        self.path = path
        self.test_type = test_type
        self.maintainer = maintainer
        self.email = email
        self.options = options
        self.test_group = test_group
        self.comment = comment
        self.timeout = timeout

    def __eq__(self, other):
        return other is not None and \
               self.name == other.name and \
               self.path == other.path and \
               self.options == other.options

    def __hash__(self):
        return hash(('name', self.name,'path', self.path, 'options', self.options))

def read_conf(conf):
    """
    Reads a configuration file and returns its content as a dictionary, with specific formatting for development versions.

    Args:
        conf (str): The path to the configuration file to be read.

    Returns:
        dict: A dictionary containing key-value pairs from the configuration file, with specific modifications for '.dev' versions.

    Notes:
        - The function assumes the configuration file has lines formatted as 'key==value'.
        - If a value contains '.dev', it is modified to follow the 'CAS-' prefix format.

    """
    with open(conf) as f:
        lines = [line.rstrip() for line in f]
    outDict = dict(x.split('==') for x in lines)
    for key in list(outDict.keys()):
        if (".dev" in outDict[key]) and ("casaconfig" not in key):
            tag = re.findall(r"a([\s\S]*)$",outDict[key])[0]
            outDict[key] = "CAS-" + tag.replace(".dev","-")
    return outDict

def run_shell_command(cmd, run_directory):
    """
    Executes a shell command in the specified directory, using ShellRunner if available, or falling back to subprocess.

    This function attempts to execute the given shell command using the ShellRunner utility. If ShellRunner encounters
    an error, the function falls back to using the subprocess module to run the command. The command is executed in the
    specified run directory.

    Args:
        cmd (str): The shell command to be executed.
        run_directory (str): The directory in which to execute the command.

    Notes:
        - The function first tries to use ShellRunner to execute the command with a default timeout.
        - If ShellRunner fails, the function changes the current working directory to the specified run directory,
          executes the command using subprocess, and then changes back to the original working directory.
        - The command's output is suppressed by redirecting stdout and stderr to subprocess.DEVNULL.

    Example:
        run_shell_command('ls -la', '/home/user/directory')

    Raises:
        Exception: If both ShellRunner and subprocess fail to execute the command.
    """
    try:
        r = ShellRunner()
        r.runshell(cmd, default_timeout, run_directory)
    except:
        cwd = os.getcwd()
        os.chdir(run_directory)
        subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr=subprocess.STDOUT)
        os.chdir(cwd)

def is_in_remote(branch,repo_path, repo):
    """
    Checks if a given branch exists in the remote repository.

    This function determines whether a specified branch exists in the remote repository by using the
    `git ls-remote` command. If the branch does not exist, it returns False, indicating that the branch
    is not present in the remote repository. For the 'master' branch, it always returns True.

    Args:
        branch (str): The name of the branch to check.
        repo_path (str): The path to the remote repository.
        repo (str): The name of the repository.

    Returns:
        bool: True if the branch exists in the remote repository, False otherwise.

    Notes:
        - If the branch is not 'master', the function constructs a `git ls-remote` command to check for the branch.
        - If the branch name starts with 'origin', it extracts the branch name from the full reference.
        - The command counts the number of matching heads in the remote repository.
        - If the count is 0, the branch does not exist in the remote repository, and the function returns False.
        - If the count is greater than 0, the branch exists, and the function returns True.
        - For the 'master' branch, the function always returns True without checking.

    """
    if branch != 'master':
        if branch.startswith("origin"):
             cmd = 'git ls-remote --heads {}{} {} | wc -l'.format(repo_path, repo, re.findall("/(.*)",branch )[0])

        else:
            cmd = 'git ls-remote --heads {}{} {} | wc -l'.format(repo_path, repo, branch)

        #print("\tRunning: ", cmd)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
        out = proc.stdout.read()
        if int(out)== 0: # If Feature Branch Exists Does Not exist, revert to master
            return False
        else:
            return True
    else:
        return True

def check_branch_path(branch):
    if branch.startswith("origin"):
        if "release" in branch:
            cmd = ("git checkout {}".format(branch)).split()
        else:
            cmd = ("git checkout origin/{}".format( re.findall("/(.*)",branch)[0])).split()
    else:
        if "release" in branch:
            cmd = ("git checkout origin/{}".format(branch)).split()
        else:
            cmd = ("git checkout origin/{}".format( re.findall("([^/]+$)",branch)[0])).split()

    return cmd

def check_branch_path_merge(branch):
    if branch.startswith("origin"):
        if "release" in branch:
            cmd = ("git merge --no-edit --verbose {}".format(branch)).split()
        else:
            cmd = ("git merge --no-edit --verbose origin/{}".format( re.findall("/(.*)",branch)[0])).split()
    else:
        if "release" in branch:
            cmd = ("git merge --no-edit --verbose origin/{}".format(branch)).split()
        else:
            cmd = ("git merge --no-edit --verbose origin/{}".format( re.findall("([^/]+$)",branch)[0])).split()

    return cmd

def fetch_tests(work_dir, branch, merge_target=None, test_paths=[]):
    """
    Fetches test files from multiple repositories and checks out the specified branch, merging if necessary.

    This function performs the following steps:
    1. Sets up the repository path and source directory.
    2. Checks if HTTPS is restricted and switches to SSH if needed.
    3. Clones the main repository ('casa6') and optionally merges the feature branch into the target branch.
    4. Appends test paths from the 'casa6' repository to the test_paths list.
    5. Clones auxiliary repositories, checks out the specified branch, and appends their test paths to the test_paths list.
    6. If the specified branch does not exist in the remote repository, defaults to 'master' or uses tags from 'build.conf'.

    Args:
        work_dir (str): The working directory where repositories will be cloned.
        branch (str): The branch to be checked out.
        merge_target (str, optional): The target branch to merge the feature branch into. Defaults to None.
        test_paths (list, optional): A list to store the paths to test files. Defaults to an empty list.

    Returns:
        list: A list of paths to the test files from the cloned repositories.

    Example:
        test_paths = fetch_tests('/path/to/workdir', 'feature-branch', 'develop')
        # Output: ['/path/to/workdir/casasources/casa6/casatests/regression/', ...]

    Notes:
        - The function assumes that the remote repositories are hosted on 'https://open-bitbucket.nrao.edu/scm/casa/'.
        - The function relies on `run_shell_command`, `check_branch_path`, `check_branch_path_merge`, `is_in_remote`, and `read_conf` functions.
        - If the merge_target is provided, the function will attempt to merge the feature branch into it.
        - If a branch does not exist in the remote repository, the function defaults to 'master' or tags from 'build.conf' if available.
    """
    if merge_target is not None:
        print("Merge Target Enabled: \n\tTarget Branch: {} \n\tFeature Branch: {}".format(merge_target, branch))

    repo_path = "https://open-bitbucket.nrao.edu/scm/casa/"
    source_dir = work_dir + "/casasources"
    # Test if https is restricted
    p1 = subprocess.Popen(shlex.split("curl -k -X GET https://open-bitbucket.nrao.edu/rest/api/1.0/projects/CASA/repos/casa6"), stdout=subprocess.PIPE)
    p2 = subprocess.Popen(shlex.split('grep "not permitted"'), stdin=p1.stdout)
    p2.communicate()
    if (p2.returncode == 0):
        repo_path = "ssh://git@open-bitbucket.nrao.edu:7999/casa/"
    if not os.path.exists(source_dir):
        os.makedirs(source_dir)

    # All of the repositositories have their tests in different directories
    # so we need a mapping
    def get_repo_test_paths(x):
        return {
            "casa6": ["/casa6/casatests/regression/","/casa6/casatests/stakeholder/","/casa6/casatasks/tests/","/casa6/casatools/tests/"],
            "casampi": ["/casampi/src/casampi/tests"],
            "casaplotms": ["/casaplotms/tests/plotms"],
            "casaviewer": ["/casaviewer/tests/tasks"]
        }[x]

    # Fetch CASA6 Repo First as Base
    repo = "casa6"
    print("Fetching Repository: {}".format(repo))
    cmd = ("git clone " + repo_path + repo).split()
    print("\tRunning: ", " ".join(str(x) for x in cmd))
    run_shell_command(cmd, source_dir)

    if merge_target is not None:

        cmd = check_branch_path(merge_target)
        print("\tRunning: ", " ".join(str(x) for x in cmd))
        run_shell_command(cmd, source_dir + "/" + repo)

        if is_in_remote(branch,repo_path, repo): # Test if the branch is in the remote repository
            print("\tMerging {} into {}".format(branch, merge_target))

            # Locally Checkout Branch
            cmd = check_branch_path(branch)
            out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
            print(out.decode("utf-8"))

            # Locally Checkout Target
            cmd = check_branch_path(merge_target)
            out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
            print(out.decode("utf-8"))

            # Merge Branch into Target
            cmd = check_branch_path_merge(branch)
            print("\tRunning: ", " ".join(str(x) for x in cmd))
            out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
            print(out.decode("utf-8"))

            print("\tRunning: git status")
            out = subprocess.check_output(["git", "status"], cwd=source_dir + "/" + repo)
            print(out.decode("utf-8"))
        else:
            print("\t{} not in Remote Repository {}".format(branch,repo))
    else:
        cmd = check_branch_path(branch)
        if is_in_remote(branch,repo_path, repo):
            print("\tRunning: ", " ".join(str(x) for x in cmd))
        else:
            print("\t{} not in Remote Repository {}. Defaulting to master".format(branch,repo))
        run_shell_command(cmd, source_dir + "/" + repo)

    for x in get_repo_test_paths(repo):
        test_paths.append(source_dir + "/" + x)

    # Clone the auxiliary repositories and checkout branch
    repositories = ["casampi", "casaplotms", "casaviewer"]
    for repo in repositories:
        print("")
        print("Fetching Repository: {}".format(repo))
        cmd = ("git clone " + repo_path + repo).split()
        print("\tRunning: ", " ".join(str(x) for x in cmd))
        run_shell_command(cmd, source_dir)

        if merge_target is not None:

            cmd = check_branch_path(branch)
            print("\tRunning: ", " ".join(str(x) for x in cmd))
            run_shell_command(cmd, source_dir + "/" + repo)

            if is_in_remote(branch,repo_path, repo): # Test if the branch is in the remote repository
                # Locally Checkout Branch
                cmd = check_branch_path(branch)
                out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
                print(out.decode("utf-8"))

                # Locally Checkout Target
                cmd = check_branch_path(merge_target)
                out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
                print(out.decode("utf-8"))
                
                print("\tMerging {} into {}".format(branch, merge_target))
                cmd = check_branch_path_merge(branch)
                print("\tRunning: ", " ".join(str(x) for x in cmd))
                out = subprocess.check_output(cmd, cwd=source_dir + "/" + repo)
                print(out.decode("utf-8"))

                print("\tRunning: git status")
                out = subprocess.check_output(["git", "status"], cwd=source_dir + "/" + repo)
                print(out.decode("utf-8"))
            else:
                print("\t{} not in Remote Repository {}".format(branch,repo))
                if os.path.isfile(source_dir+"/casa6/build.conf"):
                    print("\tCheckout from build.conf")
                    branchtag = "tags/{}".format(read_conf(source_dir+"/casa6/build.conf")[repo])
                    print("\tTag: " + branchtag)
                    cmd = ("git checkout {}".format(branchtag)).split()
                else:
                    print("No casa6/build.conf found. Defaulting to master")
                    cmd = check_branch_path(merge_target)
                print("\tRunning: ", " ".join(str(x) for x in cmd))
                run_shell_command(cmd, source_dir + "/" + repo)

        else:

            # Use Local build.conf to get build tags to git checkout
            if os.path.isfile(source_dir+"/casa6/build.conf"):
                branchtag = "tags/{}".format(read_conf(source_dir+"/casa6/build.conf")[repo])
                print("\tTag: " + branchtag)
                cmd = ("git checkout {}".format(branchtag)).split()
            else:
                # Check If Feature Branch Exists
                if is_in_remote(branch,repo_path, repo):
                    cmd = check_branch_path(branch)
                else:
                    print("\t{} not in Remote Repository {} Defaulting to master.".format(branch,repo))
                    cmd = ("git checkout origin/master").split()

            print("\tRunning: ", " ".join(str(x) for x in cmd))
            run_shell_command(cmd, source_dir + "/" + repo)

        for x in get_repo_test_paths(repo):
            test_paths.append(source_dir + "/" + x)

    return test_paths

def run_cmd(cmd, pytest_args=[]):
    try:
        from casampi.MPIEnvironment import MPIEnvironment
        if MPIEnvironment.is_mpi_enabled:
            pytest.main(cmd)
        else:
            result = subprocess.run([sys.executable,"-m","pytest"] + pytest_args + cmd , env={**os.environ})
    except:
        result = subprocess.run([sys.executable,"-m","pytest"] + pytest_args + cmd, env={**os.environ})

    return result

def setup_and_run(cmd,workdir, workpath, dirname, DRY_RUN, pytest_args ):
    """
    Sets up the environment and runs pytest with specified arguments, generating an XML report.

    This function performs the following steps:
    1. Adds verbosity, report options, and traceback formatting to the pytest command.
    2. Adds a dry run option if DRY_RUN is True.
    3. Creates the necessary directory structure for storing XML reports.
    4. Prepares the pytest command with options for XML output, disabling warnings, and other configurations.
    5. Checks if there are any tests to run in the specified workpath.
    6. Changes to the test directory and executes the pytest command.
    7. Generates pytest.ini and conftest.py files, runs the tests, updates the XML report, and cleans up temporary files.

    Args:
        cmd (list): List of pytest command arguments.
        workdir (str): The base directory where the tests will be run.
        workpath (str): The path where the XML reports and other files will be stored.
        dirname (str): The name of the directory for the current test run.
        DRY_RUN (bool): If True, the pytest command will only collect tests without executing them.
        pytest_args (list): Additional arguments to pass to the pytest command.

    Notes:
        - If no tests are found in the workpath, the function prints a message and exits.
        - The function changes the current working directory to the test directory, runs the tests, and then reverts to the original working directory.
        - Temporary files such as pytest.ini and conftest.py are cleaned up after the test run.
    """
    # https://docs.pytest.org/en/stable/usage.html
    cmd = ["--verbose"] + ["-ra"] + ["--tb=short"] + cmd

    if DRY_RUN:
        cmd = ["--collect-only"] + cmd

    if not os.path.isdir(workpath + '/xml/{}/'.format(dirname)):
        os.makedirs(workpath + '/xml/{}/'.format(dirname))
    xmlfile = workpath + 'xml/{}/nose.xml'.format(dirname)
    cmd = ["--junitxml={}".format(xmlfile)] + ["-s"] + ["--disable-pytest-warnings"] + cmd
    if len(os.listdir(workpath)) < 1: # If only the XML dir was created
        print("No Tests to Run")
        sys.exit()
    else:
        myworkdir = os.getcwd()
        os.chdir(workdir + "{}/".format(dirname))
        print("Test Directory: {}\n".format(os.getcwd()))
        print("Running Command: pytest " + " ".join(str(x) for x in cmd))
        write_pytestini(os.path.join(os.getcwd(),"pytest.ini"),dirname)
        write_conftest(os.path.join(os.getcwd(),"conftest.py"))
        result = run_cmd(cmd, pytest_args)
        update_xml(xmlfile, result, name= os.getcwd().split("/")[-1])
        #os.remove(os.path.join(os.getcwd(),"conftest.py"))
        os.remove(os.path.join(os.getcwd(),"pytest.ini"))
        os.chdir(myworkdir)

########################################################################################################################
##############################################            Run            ###############################################
########################################################################################################################

def run(testnames, branch=None, merge_target=None, DRY_RUN=False, pytest_args=[], test_paths=[]):
    """
    Runs specified tests using pytest, setting up the environment and handling test file management.

    This function performs the following steps:
    1. Checks if pytest is installed; raises an ImportError if not.
    2. Sets up directories for test execution and XML report storage.
    3. Removes duplicates from the test list and prepares test commands.
    4. Fetches and sets up tests from a remote repository if necessary.
    5. Copies test files to the working directory and sets up the environment for each test.
    6. Executes tests using pytest with the specified arguments and options.
    7. Cleans up and resets the working directory after test execution.

    Args:
        testnames (list): List of test names or paths to be executed.
        branch (str, optional): The branch to check out from the remote repository. Defaults to 'master'.
        merge_target (str, optional): The target branch to merge into. Defaults to None.
        DRY_RUN (bool, optional): If True, only collects tests without running them. Defaults to False.
        pytest_args (list, optional): Additional arguments to pass to pytest.
        test_paths (list, optional): List to store paths to test files fetched from the repository.

    Raises:
        ImportError: If pytest is not installed.

    Notes:
        - The function requires the presence of helper functions: clean_working_directory, fetch_tests, gather_all_tests, setup_and_run, getname, gettests, and any necessary imports.
        - It handles both local and remote test paths, setting up directories and copying test files as needed.
        - The function creates necessary directories, sets up pytest configuration, and manages test execution and reporting.
        - The working directory is restored to its original state after execution.
    """

    if not HAVE_PYTEST:
        raise ImportError('No Module Named Pytest. Pytest is Required for runtest.py')

    if HAVE_PYTEST:
        cwd = os.getcwd() + "/"
        workpath = os.getcwd() +"/nosedir/"
        workdir = os.getcwd() +"/nosedir/"

        clean_working_directory(workpath)
        # Copy Tests to Working Directory
        os.makedirs(workdir)

        # Remove Duplicates
        # Since working directory is based on script name
        # We need to remove multiple calls to the same script
        setlist = []
        for duplicate in list(set([ x.split("/")[-1] for x in testnames])):
            inlist = True
            for test in testnames:
                if duplicate in test:
                    if inlist:
                        setlist.append(test)
                        inlist = False
        testnames = setlist
        print("Tests: {}".format(sorted(testnames)))
        gittest = True
        if branch ==None:
            branch = 'master'
        # Only Checkout When Needed
        if any([False if ".py" in x else True for x in testnames ]):
            testpaths = fetch_tests(workdir, branch, merge_target, test_paths)
            os.makedirs(workdir + "tests/")
            for path in testpaths:
                gather_all_tests(path, workdir + "tests/")
            print("Directory Of Tests: ", workdir + "tests/")

        for testname in testnames:
            #print(testname)
            cmd = []

            # Copy Test To nosedir Directory if in cwd
            if testname.startswith("test"):
                test = testname
                # Check if specific tests are requested
                if "[" and "]" in test:
                    testname = getname(test)
                    tests = gettests(test)

                    teststring = ""
                    if len(tests) == 1:
                        teststring = tests[0]
                    elif len(tests) > 1:
                        print(tests)
                        teststring = " or ".join(tests)

                    cmd = ["-k {}".format(teststring)] + cmd
                    test = testname

                # Set up Test Working Directory
                dirname = test if not test.endswith(".py") else test[:-3]
                if not os.path.exists(workdir + "{}/".format(dirname)):
                    print("\nSetting Working Directory: {}".format(workdir + "{}/".format(dirname)))
                    os.makedirs(workdir + "{}/".format(dirname))
                    cmd = [ workdir + "{}/".format(dirname) ] + cmd

                if test.endswith(".py"):
                    try:
                        #print("Copying: {} to {}".format(test, workdir + "{}/".format(dirname)))
                        shutil.copy2(test, workdir + "{}/".format(dirname))
                    except:
                        traceback.print_exc()
                else:
                    try:
                        #print("Copying: {} to {}".format(workdir + "tests/",test), workdir + "{}/".format(dirname))
                        shutil.copy2("{}{}.py".format(workdir + "tests/",test), workdir + "{}/".format(dirname))
                    except:
                        traceback.print_exc()
                setup_and_run(cmd, workdir, workpath, dirname, DRY_RUN , pytest_args)

            ##################################################
            ########## Real Path ##########
            ##################################################
             # Copy Test To nosedir Directory assuming it's in another location
            elif testname.startswith("/"):
                testpath = testname.split("[")[0]
                cmd = []
                dirname = testname.split("/")[-1]
                test = dirname
                if "[" and "]" in test:
                    testname = getname(test)
                    tests = gettests(test)
                    teststring = ""
                    if len(tests) == 1:
                        teststring = tests[0]
                    elif len(tests) > 1:
                        print(tests)
                        teststring = " or ".join(tests)
                    cmd = ["-k {}".format(teststring)] + cmd
                    dirname = testname

                dirname = "{}".format(dirname if not dirname.endswith(".py") else dirname[:-3])

                # Set up Test Working Directory
                if not os.path.exists(workdir + "{}/".format(dirname)):
                    print("\nSetting Working Directory: {}".format(workdir + "{}/".format(dirname)))
                    os.makedirs(workdir + "{}/".format(dirname))
                    cmd = [ workdir + "{}/".format(dirname) ] + cmd
                try:
                    shutil.copy2(testpath, workdir + "{}/".format(dirname))
                except:
                    traceback.print_exc()

                setup_and_run(cmd, workdir, workpath, dirname, DRY_RUN, pytest_args )
        #build_xml(workpath + '/xml/xUnit.xml', workpath + '/xml/')
        os.chdir(cwd)

def run_bamboo_test(r, cmd, timeout, cwd):
    print("Running cmd " + str(cmd) + "in " + cwd)
    if not os.path.exists(cwd):
        os.makedirs(cwd)
    starttime = datetime.datetime.now()
    output = r.runshell(cmd, timeout,cwd)
    endtime = datetime.datetime.now()
    runtime = endtime - starttime

    return output, runtime

########################################################################################################################
#######################################            Run Bamboo Option            ########################################
########################################################################################################################

def run_bamboo(pkg, work_dir, branch = None, test_group = None, test_list= None, test_paths = [], test_config_path=None, ncores=2, verbosity=False, pmode=None, tests_to_ignore=None, merge_target=None):
    """
    Executes a set of tests specified to a package, with configurations and options, similar to using the Bamboo testing framework.

    This function performs the following steps:
    1. Validates input parameters and unpacks the provided package.
    2. Sets up a virtual frame buffer for graphical operations on Linux.
    3. Fetches necessary test paths from a repository if not provided.
    4. Reads and processes a JSON configuration file to identify test scripts and their details.
    5. Filters tests based on provided lists and groups, including handling missing configurations.
    6. Runs the filtered tests, either in parallel or serial mode, depending on the specified parameters.
    7. Generates an XML report of the test results.
    8. Cleans up by stopping the virtual frame buffer.

    Args:
        pkg (str): Path to the package to be tested.
        work_dir (str): Directory where the package will be unpacked and where tests will be executed.
        branch (str, optional): Git branch to fetch test paths from. Defaults to 'master' if not provided.
        test_group (str, optional): Comma-separated string of Jira component groups to filter tests by.
        test_list (str, optional): Comma-separated string of specific test names to include in the test run.
        test_paths (list, optional): List of directories where test scripts are located.
        test_config_path (str, optional): Path to the JSON configuration file that maps test scripts to their details. Defaults to a predefined path if not provided.
        ncores (int, optional): Number of cores to use for parallel tests. Defaults to 2.
        verbosity (bool, optional): If True, prints detailed information about tests and their configurations.
        pmode (str, optional): Mode of execution ('serial', 'parallel', or 'both'). Determines whether tests run in parallel or serially.
        tests_to_ignore (list, optional): List of test names to be ignored during execution.
        merge_target (str, optional): The branch to merge into before fetching tests. Defaults to None.

    Raises:
        Exception: If required parameters `pkg` or `work_dir` are missing or if a specified test is not found.

    Example:
        run_bamboo('path/to/package', 'working/directory', branch='feature-branch', test_group='component1,component2', test_list='test1,test2', verbosity=True, pmode='parallel')

    Notes:
        - Requires helper functions: unpack_pkg, fetch_tests, run_bamboo_test, and necessary imports for JSON handling and file operations.
        - Uses Xvfb for virtual display handling on non-Mac platforms.
        - Generates an XML report of the test results using the Xunit class.
        - The working directory is cleaned and restored after test execution.
    """
    
    if test_list is not None:
        test_list = [x.strip() for x in test_list.split(',')]
    if test_group is not None:
        test_group = [x.strip() for x in test_group.split(',')]

    print ("---------------     run_bamboo     ---------------")
    print ("Test list: {}\nTest group: {}".format(test_list,test_group))

    if pkg is None:
        raise Exception("Missing pkg")
    if work_dir is None:
        raise Exception("Missing work_dir")

    # Unpack the distribution
    exec_path, casatestutils_exec_path = atlassian_helper.unpack_pkg(pkg, work_dir, work_dir + "/pkg", default_timeout)

    print("Executable path: {}\nCasatestutils Path: {}".format(exec_path, casatestutils_exec_path))

    # Start Xvfb on Linux
    xvfb = xvfb_helper.XvfbHelper()
    if sys.platform != "darwin":
        xvfb.start_virtual_frame_buffer()

    if args.branch == None:
        branch = "master"

    print ("run_bamboo fetch_tests branch" + branch)

    # Clone a default set of repositories to if test paths are not provided from command line
    if len(test_paths) == 0 :
        test_paths = fetch_tests(str(work_dir), branch, merge_target, test_paths)

    if test_config_path == None:
       test_config_path = work_dir + "/casasources/casa6/casatestutils/casatestutils/component_to_test_map.json"
    # Read the JSON configuration
    print ("Reading config from: " + test_config_path)
    with open(test_config_path ) as f:
      test_config = json.load(f)

    # Get the actual tests as list
    test_config_elems = test_config['testlist']

    print("Test Paths: ", test_paths)

    # Match the test names provided in the JSON file to the actual test locations.
    tests_to_run = []
    for x in test_config_elems:
        for dir in test_paths:
            for currentpath, folders, files in os.walk(dir):
                for file in files:
                    if file == (x['testScript']+".py"):
                        test_location = os.path.join(currentpath, file)
                        if verbosity:
                            print("Found: " + test_location)
                            print("Script:", x['testScript'], test_location, x['testGroup'])
                        if x['timeout'] < 1:
                            timeout = default_timeout
                        else:
                            timeout = x['timeout']

                        opts = x['testOptions']
                        opts.sort()

                        tests_to_run.append(casa_test(x['testScript'],
                                                      test_location,
                                                      x['testGroup'],
                                                      x['testType'],
                                                      x['maintainer'],
                                                      x['maintainerEmail'],
                                                      tuple(opts),
                                                      x['comment'],
                                                      timeout))

    # Filter tests by test list
    if test_list is not None and len(test_list)>0:
        print ("Test list provided. Filtering tests.")
        tmp_tests_to_run = []
        for test in test_list:
            found = False
            for t1 in tests_to_run:
                if test == t1.name:
                    tmp_tests_to_run.append(t1)
                    found = True
            if not found:
                # If there is a test in the list but no configuration, add it without any options.
                # This can be useful for testing new scripts.
                print ("Test " + test + " configuration not found. Searching for the test...")
                print ("dir: " + dir )
                # Potential test location
                for dir in test_paths:
                    # Search for files (needed for casatools and casatasks)
                    for currentpath, folders, files in os.walk(dir):
                        for file in files:
                            if file == test + ".py":
                                test_location = os.path.join(currentpath, file)
                                print("Found: " + test_location)
                                print("No JSON configuration found. Test will be added to execution list without options.")
                                tmp_tests_to_run.append(casa_test(test, test_location, "","","","","",""))
                                found = True
            if not found:
                raise Exception("Couldn't locate test: " + test)
        tests_to_run = tmp_tests_to_run

    # Filter by Jira components
    if test_group is not None and len(test_group)>0:
        print("Jira component list provided. Filtering tests.")
        tmp_tests_to_run = []
        #print(test_group)
        for jira_component in test_group:
            found = False
            for t1 in tests_to_run:
                if jira_component.lower() in [x.lower() for x in t1.test_group]:
                    tmp_tests_to_run.append(t1)
                    found = True
            # Throw an exception is user provides a component that doesn't exist
            if not found:
                print ("WARNING: No tests found for jira_component " + jira_component + ". Check the contents of " + test_config_path)
        # Remove duplicates
        tests_to_run = set(tmp_tests_to_run)

    print("Subset tests:")
    for t in tests_to_run:
        if verbose:
            print(t.name + " : " +
                  t.path + " : " +
                  t.test_type + " : " +
                  t.maintainer + " : " +
                  t.email + " : " +
                  str(t.options) + " : " +
                  str(t.test_group) + " : " +
                  str(t.timeout) + " : " +
                  t.comment)
        else:
            print(t.name)

    # Run tests
    if tests_to_ignore is not None:
        print("\nTests to Ignore: ",tests_to_ignore )
        indices = []
        for i, t in enumerate(tests_to_run):
            if t.name in tests_to_ignore:
                indices.append(i)
        tests_to_run = [v for i,v in enumerate(tests_to_run) if i not in indices]

    for test in tests_to_run:
        r = ShellRunner()
        xunit = Xunit()
        # Skip MPI on Darwin for now
        if "mpi" in test.options and sys.platform != "darwin" and ( pmode == 'parallel' or pmode == 'both'):
            print("Running test: {} in MPI mode".format(test.name))
            casa_exe = exec_path + "/mpicasa"
            casaopts = "-n " + str(ncores) + " " + exec_path + "/casa" + " --nogui --nologger --log2term --agg " + cachedir + " "
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name
            if pmode == 'both':
                cwd = work_dir + "/" + test.name + '_mpi'
            output, runtime = run_bamboo_test(r, cmd, test.timeout, cwd)
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

            if pmode == 'both':
                if "casampi" in test.name:
                    continue
                print("Running test: {} in Serial mode".format(test.name))
                casaopts = " --nogui --nologger --log2term"
                casa_exe = exec_path + "/casa"
                assert (test != None)
                cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
                cwd = work_dir + "/" + test.name
                output, runtime = run_bamboo_test(r, cmd, test.timeout, cwd)
                xunit.append_result(test.name, str(runtime), len(output), output)
                print("")

        elif "mpi" not in test.options and sys.platform != "darwin" and ( pmode == 'parallel' or pmode == 'both'):
            # Special Case when you have tests in the HPC / parallel list but no mpi test option
            print("Running test: {} in Serial mode".format(test.name))
            casaopts = " --nogui --nologger --log2term"
            casa_exe = exec_path + "/casa"
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name
            output, runtime = run_bamboo_test(r, cmd, test.timeout, cwd)
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

        elif pmode == 'serial':
            if "casampi" in test.name:
                continue
            print("Running test: {} in Serial mode".format(test.name))
            casaopts = " --nogui --nologger --log2term"
            casa_exe = exec_path + "/casa"
            assert (test != None)
            cmd = (casa_exe + " " + casaopts + " -c " + test.path).split()
            cwd = work_dir + "/" + test.name
            output, runtime = run_bamboo_test(r, cmd, test.timeout, cwd)
            xunit.append_result(test.name, str(runtime), len(output), output)
            print("")

    xunit.generateXml("suite")

    # Close Xvfb on exit
    import atexit
    @atexit.register
    def goodbye():
        print("Stopping Xvfb.")
        if sys.platform != "darwin":
            xvfb.signal_stop_virtual_frame_buffer()
        print("Xvfb stopped.")

########################################################################################################################
########################################            Main-Start-Up            ###########################################
########################################################################################################################

if __name__ == "__main__":
    """
    Main entry point for executing tests with configurable options.

    This script is designed to run tests with various configurations, including:
    - Listing available tests and their tags.
    - Executing tests with verbose output.
    - Performing dry runs of tests.
    - Running specific test classes from test scripts.
    - Reading test lists from files.
    - Ignoring specified tests.

    Arguments and Options:
    - `-i` / `--list`: Print the list of tests and tags defined in `component_to_test_map.json`.
    - `-v` / `--verbose`: Enable verbose output for test execution.
    - `-x` / `--dry-run`: Perform a dry run without executing tests.
    - `-s` / `--classes`: Print the classes from a test script.
    - `-f` / `--file`: Run tests listed in an ASCII file, with one test per line.
    - `-e` / `--mapfile`: Specify a component-to-test map file.
    - `-b` / `--branch`: Specify the JIRA branch for test repository checkouts.
    - `--merge_target`: Specify the JIRA branch for test repository merges.
    - `-p` / `--pkg`: Specify the tarball or DMG file.
    - `-w` / `--work_dir`: Specify the working directory.
    - `-n` / `--ncores`: Number of cores for MPI tests (default: 2).
    - `-t` / `--test_paths`: Comma-separated list of paths containing tests.
    - `-l` / `--test_list`: Filter tests by a comma-separated list of test names.
    - `-c` / `--test_config`: Specify the test configuration file.
    - `-j` / `--test_group`: Filter tests by a comma-separated list of components.
    - `-m` / `--pmode`: Set parallelization mode (serial, parallel, both).
    - `--bamboo`: Set the Bamboo flag to True.
    - `-r` / `--cachedir`: Specify the Casa cachedir ( previously --rcdir, which also covered the paths to the startup and config files).
    - `--ignore_list`: Specify a map file or comma-separated list of tests to ignore.

    Execution Flow:
    1. Parse command-line arguments and print them.
    2. Load and process the `ignore_list` if provided.
    3. Handle optional `cachedir ( previously rcdir )` and `test_group` arguments, loading the appropriate component-to-test map if necessary.
    4. List tests and classes, handle dry runs, and process files for test names.
    5. Prepare arguments for pytest and handle unknown arguments, including test case filtering.
    6. Execute tests either through Bamboo integration or locally using the `run()` function.
    7. Print errors and traceback information if exceptions occur.

    The script integrates with Bamboo for CI/CD pipelines when the `--bamboo` flag is set, otherwise it runs tests locally based on the specified configurations.
    """
    
    print("HAVE_PYTEST: {}".format(HAVE_PYTEST))
    print("")

    # List of tests to run
    testnames = []
    test_paths = []

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument("-i", "--list",action='store_true',help='print the list of tests & tags defined in component_to_test_map.json')
    parser.add_argument("-v", "--verbose",action='store_true',help="Verbose Test Execution")
    parser.add_argument("-x", "--dry-run",action='store_true',help="dry run Test Execution")
    parser.add_argument("-s", "--classes",nargs='+',metavar='test',help='print the classes from a test script') # copy of Dry-Run
    parser.add_argument("-f", "--file",nargs='?', type=argparse.FileType('r'),help='run the tests defined in an ASCII file <list>; one test per line')

    # Component Arguments
    parser.add_argument("-e","--mapfile", nargs='?', type=argparse.FileType('r'), help='Component to test map file', required=False)
    parser.add_argument("-b","--branch", help='JIRA Branch for test repository checkouts', required=False)
    parser.add_argument("--merge_target", help='JIRA Branch for test repository merge', required=False)

    # casa-build-utils Arguments
    parser.add_argument('-p','--pkg', help='Tarball or dmg', required=False)
    parser.add_argument('-w','--work_dir', help='Working directory.', required=False)
    parser.add_argument('-n','--ncores', help='Number of cores for MPI tests', default=2)
    parser.add_argument('-t','--test_paths', help='A comma separated list of paths containing tests.', required=False)
    parser.add_argument('-l','--test_list', help='Filter tests by a comma separated list of tests', required=False)
    parser.add_argument('-c','--test_config',  help='Test configuration file', required=False)
    parser.add_argument('-j','--test_group',  help='Filter tests by a comma separated list of components', required=False)
    parser.add_argument('-m','--pmode',  help='Parallelization mode: serial, parallel, both', required=False)
    parser.add_argument('--bamboo', help='Set Bamboo Flag to True',default=False,action='store_true', required=False)
    parser.add_argument('-r','--cachedir',  help='Casa cachedir ( previously --rcdir, which also covered the paths to the startup and config files)', required=False)
    parser.add_argument('--ignore_list',  help='map file of tests to ignore', required=False)

    args, unknownArgs = parser.parse_known_args()

    print(args)
    print("")

    tests_to_ignore = None
    if args.ignore_list is not None:
        if args.ignore_list.endswith(".json"):
            ignore_test_map = json.load(open(args.ignore_list))
            tests_to_ignore = [x["testScript"].strip() for x in ignore_test_map["testlist"]]
        else:
            tests_to_ignore = [x.strip() for x in args.ignore_list.split(",")]

    print("Operating system: " +  platform.system())
    print("")

    cachedir=""
    if args.cachedir is not None:
        cachedir="--cachedir=" + args.cachedir
        print("cachedir: " + cachedir)

    if args.test_group is not None:
        components = args.test_group
        components = [x.strip() for x in components.split(",")]
        if len(components) == 1 and not components[0]:
            print("Component list is empty. Using component 'default'")
            components = ["default"]
        print("Testing Components" + str(components))
        print("")

        if not args.bamboo:
            if args.mapfile is not None:
                component_to_test_map = json.load(args.mapfile)
            else:
                try:
                    import casatestutils as _;
                    with open("{}/{}".format(_.__path__[0], "component_to_test_map.json")) as ctt:
                        component_to_test_map = json.load(ctt)
                except:
                    print("No JSON file to Map")
            no_test_components = []
            for c in components:
                _isComponent = False
                component = c.strip()
                for myDict in component_to_test_map["testlist"]:
                    if component in myDict["testGroup"] or component in myDict["testType"]:
                        _isComponent = True
                        if (myDict["testScript"] not in testnames):
                            if tests_to_ignore is not None:
                                if myDict["testScript"] in tests_to_ignore:
                                    continue
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

    verbose = False
    if args.verbose:
        verbose = True

    if args.list:
        try:
            tmp = {}
            import casatestutils as _;
            with open("{}/{}".format(_.__path__[0], "component_to_test_map.json")) as ctt:
                component_to_test_map = json.load(ctt)
            for myDict in component_to_test_map["testlist"]:
                tmp[myDict["testScript"]] = myDict["testGroup"]
            for key, value in tmp.items():
                print(key,value)
        except:
            list_tests(os.getcwd())
        sys.exit(1)

    ## Dry Run
    DRY_RUN = False
    if args.dry_run or (args.classes is not None):
        DRY_RUN = True

    if args.file is not None:
        for line in args.file:
            try:
                testnames.append(re.sub(r'[\n\r]+', '',line))
            except:
                raise Exception(" The list should contain one test per line.")

    if args.test_paths is not None:
        test_paths = [x.strip() for x in args.test_paths.split(',')]

    temp_storage = []
    pytest_args = []

    for arg in unknownArgs:
        if arg.startswith(("-", "--")):
            pytest_args.append(arg)
        else:
            if '[' in arg:
                tests = [x.strip() for x in arg.split("],")]
                for i in range(len(tests)):
                    test = tests[i]
                    if '[' in test and not test.endswith("]"):
                        tests[i] = tests[i] + "]"
                for i in range(len(tests)):
                    test = tests[i]
                    if test.find(",") < test.find('['):
                        temp_storage = temp_storage + test.split(',',1)
                    else:
                        temp_storage.append(test)
                tests = temp_storage
            else:
                tests = [x.strip() for x in arg.split(",")]
            for test in tests:
                try:
                    testcases = None
                    # Check if testcases are provided
                    if "[" in test:
                        listarray = test.split("[")
                        if not listarray[0].endswith(".py"):
                            testnames.append(test)
                        else:
                            test = listarray[0]
                            testcases = listarray[1]
                            real_path = os.path.realpath(test)
                            if ("test_" not in real_path) or  ("test_" not in real_path) or ( os.path.exists(real_path) ==False):
                                print("{} is not a Test File".format(test))
                                continue
                            else:
                                if testcases is not None: real_path = os.path.realpath(test) + "[" + testcases
                                testnames.append(real_path)

                    # Check if test is real path are provided
                    elif test.endswith(".py"):
                        real_path = os.path.realpath(test)
                        if ("test_" not in real_path) or ( os.path.exists(real_path) ==False):
                            print("{} is not a Test File".format(test))
                            continue
                        else:
                            testnames.append(real_path)

                    # else Assume test exists in bitbucket
                    else:
                        testnames.append(test)
                except:
                    traceback.print_exc()
    print("Arguments Sent Direct To Pytest : ",pytest_args)

    try:
        if args.bamboo:
            from testrunner.shell_runner import ShellRunner
            from testrunner import xvfb_helper
            from testrunner.xunit import Xunit
            from testrunner import atlassian_helper
            if args.pkg:
                print("Package: " + args.pkg)
            print("Test configuration file: " + str(args.test_config))
            print("Number of cores: " + str(args.ncores))
            print("Workdir: " + str(args.work_dir))
            print("branch: " + str(args.branch))
            pmodes = ['serial','parallel','both']
            if args.pmode not in pmodes:
                raise Exception("Invalid pmode: '{}'. Valid modes: '{}'".format(args.pmode ,str(pmodes)))

            run_bamboo(args.pkg, args.work_dir, args.branch, args.test_group, args.test_list, test_paths, args.test_config, args.ncores, args.verbose, args.pmode, tests_to_ignore, args.merge_target)

        else:
            #If no tests are given, no subet tag or --all option
            if args.test_paths is not None:
                tests = []
                test_paths = [x.strip() for x in args.test_paths.split(',')]
                if len(testnames) == 0:
                    for test_path in test_paths:
                        for root, dirs, files in os.walk(test_path):
                            for file in files:
                                if file.endswith(".py") and file.startswith("test_"):
                                     tests.append(os.path.realpath(os.path.join(root, file)))
                else:
                    for test_path in test_paths:
                        for test in testnames:
                            if not test.endswith(".py"):
                                test = test + ".py"
                            for root, dirs, files in os.walk(test_path):
                                for file in files:
                                    if file == test:
                                        tests.append(os.path.realpath(os.path.join(root, file)))
                testnames = tests
            if tests_to_ignore is not None:
                print("\nTests to Ignore: ",tests_to_ignore )
                indices = []
                for i, t in enumerate(testnames):
                     if t.split("/")[-1].replace(".py","") in tests_to_ignore:
                        indices.append(i)
                testnames = [v for i,v in enumerate(testnames) if i not in indices]
            if testnames == [] or len(testnames) == 0:
                print("List of tests is empty")
                parser.print_help(sys.stderr)
                sys.exit(1)
            print("Running {} Test(s)".format(len(testnames)))
            run(testnames, args.branch, args.merge_target, DRY_RUN, pytest_args, test_paths if args.test_paths is not None else [])
    except:
        traceback.print_exc()
