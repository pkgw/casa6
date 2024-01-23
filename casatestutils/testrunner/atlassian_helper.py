import os
import platform

from .shell_runner import ShellRunner

def get_casatestutils_exec_path(pkg_dir):
    # Since runtest is no longer part of casatestutils, this may be removed.
    for currentpath, folders, files in os.walk(pkg_dir):
        for file in files:
            #print(">>>" + os.path.join(currentpath, file))
            if currentpath.endswith('casatestutils') and file == 'runtest.py':
                return(os.path.join(currentpath, file))
    return "/dev/null/"

def unpack_tarball(pkg, outputdir, default_timeout=1800):
    print ("Unpacking tarball: " + pkg + " to " +  outputdir)
    cmd = ("tar -xf " + pkg + " -C " + outputdir).split()
    print(cmd)
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())

    installpath = None

    print("Outputdir contents:" + outputdir)
    for root, dirs, files in os.walk(outputdir):
        for d in dirs:
            print(" " + d)
            if d.startswith("casa-"):
                installpath = d
                print("installpath: " + installpath)
        break

    if installpath is None:
        raise  RuntimeError("Couldn't find a directory that looks like a Casa distribution. Expected directory name to start with 'casa-'")
    return outputdir + "/" + installpath

def unpack_dmg(pkg, work_dir, outputdir, default_timeout=1800):
    mountpoint = work_dir + "/mnt"
    if not os.path.exists(mountpoint):
        os.makedirs(mountpoint)

    print ("Unpacking dmg: " + pkg + " to " +  outputdir)
    cmd = ("hdiutil attach " + pkg + " -mountpoint " + mountpoint).split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    installpath = outputdir + "/CASA.app"
    cmd = ("ditto " + mountpoint + "/CASA.app " + outputdir + "/CASA.app").split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    cmd = ("hdiutil detach " + mountpoint).split()
    r = ShellRunner()
    output = r.runshell(cmd, default_timeout, cwd=os.getcwd())
    return installpath

def unpack_pkg(pkg, work_dir, outputdir, default_timeout=1800):
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if platform.system() == "Linux":
        installpath = unpack_tarball(pkg, outputdir, default_timeout)
        print ("Package root: " + installpath)
        exec_path = installpath + "/bin"
    elif platform.system() == "Darwin":
        installpath = unpack_dmg(pkg,work_dir, outputdir, default_timeout)
        print("Package root: " + installpath)
        exec_path = installpath + "/Contents/MacOS"
    else:
        raise Exception("Unknown operating system")
    if exec_path is None:
        raise Exception ("Couldn't find casa executable path")
    casatestutils_exec_path = get_casatestutils_exec_path(installpath)
    if casatestutils_exec_path == None:
        raise Exception("Couldn't find casatestutils")
    return exec_path, casatestutils_exec_path
