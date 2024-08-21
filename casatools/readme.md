
## CASAtools
 
If you are looking for the place to check out [CASA](http://casa.nrao.edu/) with [Git](https://en.wikipedia.org/wiki/Git), please refer to the top level directory of [this repository](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse).

CASAtools is a self-contained python module that provides the tools from the [CASA](http://casa.nrao.edu/) project.

## Building casatools

Please refer to the top level [readme.md](../readme.md).

### Testing casatools

A quick test to see if the casatools installation is successful is to try loading the CASAtools module:

```
$ python3
Python 3.8.13 (default, Nov  8 2022, 17:19:05) 
[GCC 8.5.0 20210514 (Red Hat 8.5.0-15)] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from casatools import image
>>> ia = image( )
>>> ia.fromshape("mytest.im",[20,20])
2019-03-25 16:15:08     INFO    ImageFactory::createImage       Created Paged image 'mytest.im' of shape [20, 20] with float valued pixels.
True
>>>
-bash-4.2$
```

## Available Tools

| Tool Name            | Description                                                     |
| -------------------- | --------------------------------------------------------------- |
| agentflagger         | Tool for manual and automated flagging                          |
| atcafiller           | Filler for ATNF/ATCA RPFITS data                                |
| atmosphere           | Atmosphere model                                                |
| calanalysis          | Get and fit data from a calibration table (CASA 3.4 and later). |
| calibrater           | Synthesis calibration (self- and cross-)                        |
| componentlist        | A tool for the manipulation of groups of components             |
| coordsys             | Operations on CoordinateSystems                                 |
| functional           | Functionals handling                                            |
| imagemetadata        | Operations on imagemetadata                                     |
| image                | Operations on images                                            |
| imagepol             | Polarimetric analysis of images                                 |
| imager               | tool for synthesis imaging                                      |
| iterbotsink          | tool for synthesis imaging                                      |
| logsink              | tool for logsink                                                |
| measures             | measures tool                                                   |
| miriadfiller         | Tool for the importmiriad task                                  |
| msmetadata           | Operations to retrieve metadata from a measurment set           |
| ms                   | Operations on measurement sets                                  |
| mstransformer        | Tool to apply spw and frame transformations in MS               |
| quanta               | quanta tool handles units and quantities                        |
| regionmanager        | Create and manipulate regions of interest                       |
| sakura               | New single dish tool interface using sakura                     |
| sdm                  | Manipulate or examine SDM datasets                              |
| sidebandseparator    | Tool for sideband separation                                    |
| simulator            | Tool for simulation                                             |
| singledishms         | New single dish tool interface to process an MS                 |
| spectralline         | spectral line tool                                              |
| synthesisdeconvolver | tool for synthesis imaging                                      |
| synthesisimager      | tool for synthesis imaging                                      |
| synthesisimstore     | tool for synthesis imaging                                      |
| synthesismaskhandler | tool for mask handling in sysnthesis imaging                    |
| synthesisnormalizer  | tool for synthesis imaging                                      |
| synthesisutils       | tool for synthesis imaging                                      |
| table                | Access tables from casapy                                       |
| vlafiller            | tool for VLA filler tasks                                       |
| vpmanager            | Tool for specifying voltage patterns and primary beams          |
| wvr                  | Tool for the analysis of Water Vapour Radiometer data           |

## XML

### Tool Specification

<pre>
    &lt;tool name="<i><font color="red">tool-name</font></i>"&gt;
        <font color="blue">&lt;shortdescription&gt;</font><i>one-line description</i><font color="blue">&lt;/shortdescription&gt;</font>
        <font color="blue">&lt;description&gt;</font><i>paragraph description</i><font color="blue">&lt;/description&gt;</font>

        <font color="#7CFC00">&lt;needs&gt;</font><i><font color="red">tool-name</font></i><font color="#7CFC00">&lt;/needs&gt;</font>

        <font color="blue">&lt;code&gt;</font>
            <font color="blue">&lt;include&gt;</font>path to headerfile<font color="blue">&lt;/include&gt;</font>
            <font color="blue">&lt;private&gt;</font>
                <font color="blue">&lt;include&gt;</font>path to headerfile<font color="blue">&lt;/include&gt;</font>
            <font color="blue">&lt;/private&gt;</font>
        <font color="blue">&lt;/code&gt;</font>

        <font color="#7CFC00">&lt;method</font> name="<i><font color="red">method-name</font></i>"<font color="#7CFC00">&gt;</font>
            <font color="blue">&lt;shortdescription&gt;</font><i>one-line description</i><font color="blue">&lt;/shortdescription&gt;</font>
            <font color="blue">&lt;description&gt;</font><i>paragraph description</i><font color="blue">&lt;/description&gt;</font>
            <font color="blue">&lt;input&gt;</font>
                <font color="#7CFC00">&lt;param</font> name="<font color="red">param-name</font>"<font color="#7CFC00">&gt;</font>
                    <font color="blue">&lt;description&gt;</font><i>short description</i><font color="blue">&lt;/description&gt;</font>
                    <font color="#7CFC00">&lt;type</font> <font color="blue">units="<i><font color="red">unit-name</font></i>" mustexist="true"</font><font color="#7CFC00">&gt;</font><i><font color="red">type-name</font></i><font color="#7CFC00">&lt;/type&gt;</font>
                    <b><font color="blue">&lt;value&gt;value-content&lt;/value&gt;</font></b>
                <font color="#7CFC00">&lt;/param&gt;</font>
            <font color="blue">&lt;/input&gt;</font>
            <font color="blue">&lt;output&gt;</font>
                <font color="#7CFC00">&lt;param</font> name="<font color="red">param-name</font>"<font color="#7CFC00">&gt;</font>
                    <font color="blue">&lt;description&gt;</font><i>short description</i><font color="blue">&lt;/description&gt;</font>
                    <font color="#7CFC00">&lt;type</font> <font color="blue">units="<i><font color="red">unit-name</font></i>" mustexist="true"</font><font color="#7CFC00">&gt;</font><i><font color="red">type-name</font></i><font color="#7CFC00">&lt;/type&gt;</font>
                    <b><font color="blue">&lt;value&gt;value-content&lt;/value&gt;</font></b>
                <font color="#7CFC00">&lt;/param&gt;</font>
            <font color="blue">&lt;/output&gt;</font>

            <font color="blue">&lt;returns</font> type="<font color="red">type-name</font>"<font color="blue">/&gt;</font>
        <font color="#7CFC00">&lt;/method&gt;</font>
    &lt;/tool&gt;
</pre>

* <font color="red">Red</font> text indicate unique names. They may be developer specified (e.g. tool-name), they may be predefined enumerations (e.g. unit-name), or they may be both (e.g. type-name).
* <font color="blue">Blue</font> XML elements or attributes indicates optional items (zero or one).
* <font color="#7CFC00">Green</font> XML elements can be repeated.
* Black XML elements or attributes must be supplied (once)
