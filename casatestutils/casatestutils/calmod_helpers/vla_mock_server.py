#!/usr/bin/env python

# modified version of Josh Marvil's script to work with CASA tests
from flask import Flask,request
import os

APPNAME = 'VLA Calibrator Database'
app = Flask(__name__)

app.url_map.strict_slashes = False
app.config.update(
    APPNAME=APPNAME,
)

@app.route("/", methods = ['GET']) 
def vlacal():

    if 'type' not in request.args:
        return 'type of query required', 400
    if request.args['type'] == 'setjy':
        if (
            (
                ('source' in request.args and request.args['source'] == '3C48')
                or (
                    'position' in request.args 
                    and request.args['position'] == 'J2000 01:37:41.1 33.09.32'
                )
            )
            and request.args['band'] == 'Q' and request.args['obsdate'] == '50000'
        ):
            myfile = os.sep.join(
                [os.path.dirname(os.path.abspath(__file__)), 'query1.json']
            )
            with open(myfile, 'r') as file:
                file_contents = file.read()
            return file_contents, 200
    return request.args, 400


if __name__ == "__main__":
   app.run(debug=True,host='127.0.0.1',port=8080)


# example query
# http://192.168.1.14:5100/?type=setjy&source=3C48&band=C&obsdate=54321
