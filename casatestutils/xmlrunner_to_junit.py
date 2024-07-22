import xml.etree.cElementTree as ET
import re
import datetime
import socket
import argparse
import os
def convert(filename, outfile):
    """
    Converts a test result file to an XML format suitable for a test report.

    This function reads a test result file, either creating a default XML file if the
    input file does not exist or converting the contents of an existing test result file
    into a standardized XML format.

    If the input file does not exist:
    - Creates a new XML file with default values indicating that no tests were generated.
    - Sets various attributes such as the hostname, timestamp, and test details.

    If the input file exists:
    - Reads the file contents and extracts test results using regular expressions.
    - Parses and aggregates test results, including the number of tests, failures, errors,
      and skipped tests.
    - Constructs an XML document with the aggregated results and writes it to the specified
      output file.

    Args:
        filename (str): The path to the input test result file. If the file does not exist,
                        a default XML file is generated.
        outfile (str): The path to the output XML file where the converted results will be
                       written.

    Returns:
        None

    Raises:
        Exception: If there is an issue with reading the input file or writing to the output
                   file.

    Example:
        convert('input_results.xml', 'output_results.xml')
    """
    if not os.path.isfile(filename):
        print("File {} does not exist. Generating:...".format(filename))
        fMessage = "{} Not Generated. Check Log".format(filename)
        name = filename.split(".py")[0].split("/")[-1].strip(".xml")
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

        with open(outfile, "wb") as f:
            f.write(b_xml)
        return

    print("Reading: ",filename )
    with open(filename, 'r') as f:
        data = f.read()
    print("Converting... ")
    results = re.findall('(?s)(?<=\<testsuite).*?(?=\<\/testsuite\>)', data)

    start="""<?xml version="1.0" ?>
    <testsuite
    """

    end="""</testsuite>"""
    name = filename.split(".py")[0].split("/")[-1]

    e = datetime.datetime.now()
    timestamp = e.strftime('%Y-%m-%dT%H:%M:%S.%f')
    data = ET.Element('testsuites')
    element1 = ET.SubElement(data, 'testsuite')
    element1.set('name', "'{}'".format(name))

    tests, failures, errors, skipped = 0 , 0 , 0, 0
 
    runtime = 0.0
    for result in results:
        s = "{}{}{}".format(start, result, end)
        root = ET.fromstring(s)
        if "tests" in root.attrib.keys(): tests = tests + int(root.attrib['tests'])
        if "failures" in root.attrib.keys(): failures = failures + int(root.attrib['failures'])
        if "errors" in root.attrib.keys(): errors = errors + int(root.attrib['errors'])
        # if "skipped" in root.attrib.keys(): skipped = skipped + int(root.attrib['skipped'])
        runtime = runtime + float(root.attrib['time'])

        for i in range(len(root)):
            if not root[i].attrib: continue
            #print(root[i].attrib)
            s_elem1 = ET.SubElement(element1, 'testcase')
            s_elem1.set('classname', "{}".format(name))
            s_elem1.set('name', "{}.{}".format(root[i].attrib['classname'],root[i].attrib['name']))
            s_elem1.set('time', root[i].attrib['time'])

            if len(root[i].getchildren()) != 0:
                for child in root[i].getchildren():
                    if child.attrib['type'] == 'skip':
                        ss_elem1 = ET.SubElement(s_elem1, 'skipped')
                        skipped = skipped + 1
                    else:
                        ss_elem1 = ET.SubElement(s_elem1, 'failure')
                    fMessage = child.attrib['message']
                    ss_elem1.set('message', fMessage)
                    ss_elem1.text = fMessage

    element1.set('hostname', socket.gethostname())
    element1.set('timestamp', timestamp)
    element1.set('errors', "0") # for Bamboo, Errors are treated as Failures
    element1.set('failures', str(max(failures, errors))) # for Bamboo, Errors are treated as Failures
    element1.set('skipped', str(skipped))
    element1.set('tests', str(tests))
    element1.set('time', str(round(runtime, 3)))
    b_xml = ET.tostring(data)

    with open(outfile, "wb") as f:
        f.write(b_xml)
    print("File Written To: ",outfile )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", help="XMLRunner module File")
    parser.add_argument("--outfile", help="Converted XMLRunner Output File")
    args = parser.parse_args()

    convert(args.infile, args.outfile)
