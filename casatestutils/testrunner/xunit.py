class Xunit:
    xml_escape_table = {
        "&": "&amp;",
        '"': "&quot;",
        "'": "&apos;",
        ">": "&gt;",
        "<": "&lt;",
    }

    results = []
    fail_total = 0
    def append_result (self,testname, runtime, returncode, testerr):
        result =    {
          "testname": testname,
          "runtime" : runtime,
          "returncode": returncode,
          "testerr": testerr,
        }
        self.results.append(result)

    def xml_escape(self, text):
        return "".join(self.xml_escape_table.get(c,c) for c in text)

    def test_result_to_xml (self,result):
        import signal
        import datetime
        import time
        x = time.strptime(result['runtime'].split(',')[0],'%H:%M:%S.%f')
        runtime = datetime.timedelta(
                hours=x.tm_hour,
                minutes=x.tm_min,
                seconds=x.tm_sec).total_seconds()

        self.fail_total = self.fail_total + len(result['testerr'])
        testxml = '<testcase classname="' + result['testname'].replace(".py", "") + '.{}"'.format(result['testname'].replace(".py", "")) \
              + ' name="'+ result['testname'].replace(".py", "") + '" time="' + str(runtime) + '">'
        if ( result['returncode'] != 0) :
            fMessage = result['testerr']
            try: fMessage = signal.strsignal(abs(result['returncode']))
            except: fMessage = signal.Signals(abs(result['returncode'])).name
            testxml = testxml + '<failure>' + str(fMessage) + '</failure>'
            if self.fail_total == 0:
                self.fail_total=self.fail_total + 1
        testxml = testxml + '</testcase>\n'
        return testxml

    def generateXml(self, testname):
        import datetime, socket
        e = datetime.datetime.now()
        timestamp = e.strftime('%Y-%m-%dT%H:%M:%S.%f')
        xmlResults = list(map(lambda result: self.test_result_to_xml(result), self.results))

        testHeader = '<?xml version="1.0" encoding="UTF-8"?>' + "\n" \
                 + '<testsuites><testsuite name="runtest.bamboo" tests="' \
                 + str(len(self.results)) + '" errors="'+ str(self.fail_total) + '"' \
                 + ' failures="' + str(self.fail_total) + '" skipped="0" timestamp="'+ timestamp \
                 + '" hostname="' + socket.gethostname()+ '">\n'
        print("Results: " + str(self.results))
        testFooter ="\n</testsuite></testsuites>"

        # Write xUnit.xml
        xUnit = open("xUnit-"+testname+".xml", "w+")
        xUnit.write(testHeader + ''.join(xmlResults) + testFooter)
        xUnit.close()
