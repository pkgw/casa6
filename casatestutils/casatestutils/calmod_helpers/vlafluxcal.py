#!/usr/bin/env python

# modified version of Josh Marvil's script to work with CASA tests

from flask import Flask,request
import json
import numpy as np
import os, sys
import sqlite3
from sqlite3 import Error as sqlError
from casatools import quanta, measures, componentlist
from casatasks import casalog

APPNAME = 'VLA Calibrator Database'

# SQLITE = '/Users/jmarvil/vlacal2/mysite/static/vlafluxcal.sqlite'
# SQLITE = '/export/home/eowyn/dmehring/casa_eclipse/projects/casa6-modular-build/vla-flux-cal/static/vlafluxcal.sqlite'
SQLITE = os.sep.join([os.path.dirname(os.path.abspath(__file__)), 'vlafluxcal.sqlite'])

print(SQLITE)

app = Flask(__name__)

app.url_map.strict_slashes = False
app.config.update(
    APPNAME=APPNAME,
    SQLITE=SQLITE,
)




class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)




def create_connection(path):
    connection = None
    try:
        connection = sqlite3.connect(path)
    except sqlError as e:
        raise
    return connection


def execute_read(connection, query):
    cursor = connection.cursor()
    result = None
    try:
        casalog.post(f'Attempt to run {query}', 'INFO')
        cursor.execute(query)
        result = cursor.fetchall()
        return result
    except sqlError as e:
        casalog.post(f'Query failed, Exception was {e}', 'WARN')
        raise


flux_query = '''\
SELECT f2.date_observed, f2.date_added, f2.coefficients \
FROM ( \
    SELECT f.date_observed, f.date_added, f.coefficients, \
           ROW_NUMBER() OVER (PARTITION BY f.date_observed ORDER BY f.date_added DESC) as rn \
    FROM fluxscale f \
    WHERE f.source = '{0}' AND f.date_added <= {1} \
    ) f2 \
WHERE f2.rn = 1;
'''


component_query = '''\
SELECT c2.date_observed, c2.date_added, c2.complist, c2.csys \
FROM ( \
    SELECT c.date_observed, c.date_added, c.complist, c.csys, \
           ROW_NUMBER() OVER (PARTITION BY c.date_observed ORDER BY c.date_added DESC) as rn \
    FROM components c \
    WHERE c.source = '{0}' AND c.band = '{1}' AND c.date_added <= {2} \
    ) c2 \
WHERE c2.rn = 1;
'''


def vlacal_summary( refdate ):
    connection = create_connection(SQLITE)
    connection.row_factory = sqlite3.Row
    fluxscale_summary = execute_read(connection, f"select source,date_observed,date_added from fluxscale where date_added <= {refdate};")
    component_summary = execute_read(connection, f"select source,band,date_observed,date_added from components where date_added <= {refdate};")
    return {'fluxscale':[dict(f) for f in fluxscale_summary], 'components':[dict(c) for c in component_summary]}


def vlacal_fluxscale( source, refdate ):
    connection = create_connection(SQLITE)
    connection.row_factory = sqlite3.Row
    fluxscale_results = execute_read(connection, flux_query.format( source.upper(), refdate ))

    if len(fluxscale_results) == 0:
        return {}
    else:

        dict_items = [dict(f) for f in fluxscale_results]

        dict_keys = [i for i in range(len(dict_items))]

        # have server turn database strings (stringify'd lists) into actual json (like dictionaries)
        for item in dict_items:
            item['coefficients'] = json.loads(item['coefficients'])

        # to give rows a custom name
        return {key:value for (key,value) in zip(dict_keys,dict_items)}


def vlacal_components( source, band, refdate ):
    connection = create_connection(SQLITE)
    connection.row_factory = sqlite3.Row

    try:
        component_results = execute_read(connection, component_query.format( source.upper(), band.upper(), refdate ))
        return_dict = dict(component_results[0])
        for key in ['complist','csys']:
            return_dict[key] = json.loads(return_dict[key])

        return return_dict
    except:
        return {}


def vlacal_calc_flux_from_coeffs( coeff, freq ):
        x = np.log10( freq / 1e9 )
        p1 = np.poly1d( coeff[::-1] )
        return 10**(p1(x))


def vlacal_interpolate_fluxscale( fluxscale_coeffs, freq, target_date, verbose=True ):

    dates = np.array( [ fluxscale_coeffs[key]['date_observed'] for key in fluxscale_coeffs.keys() ], dtype='int')
    coeffs = [ fluxscale_coeffs[key]['coefficients'] for key in fluxscale_coeffs.keys() ]
    to_sort = dates.argsort()

    fluxes = np.array( [ vlacal_calc_flux_from_coeffs( coef, freq ) for coef in coeffs ] )

    if verbose:
        print(f'Available fluxscale data at reference frequency {freq/1e9} GHz:')
        for i,date in enumerate(dates[to_sort]):
            print('   {0}  {1}'.format(date, fluxes[to_sort][i]) )

    new_flux = np.interp( target_date, dates[to_sort], fluxes[to_sort] )

    return new_flux


def vlacal_band_to_freqs( band ):
        
    band_table_txt = '''\
    P, 0.23, 0.47
    L, 1, 2
    S, 2, 4
    C, 4, 8
    X, 8, 12
    U, 12, 18
    K, 18, 26.5
    A, 26.5, 40
    Q, 40, 50
    '''

    band_table = np.genfromtxt( band_table_txt.split('\n'), dtype=None, delimiter=',', encoding='utf-8' )
    f_min, f_max = band_table[band_table['f0']==band][['f1','f2']][0]
    f_min = f_min*1e9
    f_max = f_max*1e9
    f_reference = 0.5*(f_min+f_max)

    return (f_reference, f_min, f_max)



def vlacal_test_direction( position ):
    me = measures()
    default_direction = me.direction()
    try:
        d1 = me.direction( *position.split() )
        assert d1 != {}
        assert d1 != default_direction
        return True
    except:
        return False


def vlacal_position_to_source( position, tolerance ):
    qa = quanta()
    me = measures()
    source_directions = {
        '3C48':  'J2000 01:37:41.299431 33.09.35.13299',
        '3C138': 'J2000 05:21:09.886021 16.38.22.05122',
        '3C147': 'J2000 05:42:36.137916 49.51.7.23356',
        '3C286': 'J2000 13:31:8.287984 30.30.32.95885'
    }
    d1 = me.direction( *position.split() )
    sources = source_directions.keys()
    for source in sources:
        d2 = me.direction( *source_directions[source].split() )
        sep = qa.convert(me.separation(d1,d2),'arcsec')['value']
        if sep < tolerance:
            return source
    return ''



def vlacal_scalefactor( source, f_reference, obsdate, refdate ):

    fluxscale_coeffs = vlacal_fluxscale( source, refdate )

    target_date = obsdate

    new_flux = vlacal_interpolate_fluxscale( fluxscale_coeffs, f_reference, target_date, verbose=False )

    pb17_coeffs = [fluxscale_coeffs[key]['coefficients'] for key in fluxscale_coeffs.keys() if fluxscale_coeffs[key]['date_observed'] == 57754][0]
    pb17_flux = vlacal_calc_flux_from_coeffs( pb17_coeffs, f_reference )

    excess_flux = new_flux - pb17_flux

    return new_flux, excess_flux, pb17_flux, f_reference




def vlacal_cl_for_setjy( source, band, obsdate, refdate ):

    qa = quanta()
    cl = componentlist()

    f_reference, f_min, f_max = vlacal_band_to_freqs( band )

    fluxscale_coeffs = vlacal_fluxscale( source, refdate )

    target_date = obsdate

    new_flux = vlacal_interpolate_fluxscale( fluxscale_coeffs, f_reference, target_date, verbose=True )


    pb17_coeffs = [fluxscale_coeffs[key]['coefficients'] for key in fluxscale_coeffs.keys() if fluxscale_coeffs[key]['date_observed'] == 57754][0]
    pb17_flux = vlacal_calc_flux_from_coeffs( pb17_coeffs, f_reference )

    excess_flux = new_flux - pb17_flux
    pb17_low = vlacal_calc_flux_from_coeffs( pb17_coeffs, f_min )
    pb17_high = vlacal_calc_flux_from_coeffs( pb17_coeffs, f_max )

    flux_low = vlacal_interpolate_fluxscale( fluxscale_coeffs, f_min, target_date)
    flux_high = vlacal_interpolate_fluxscale( fluxscale_coeffs, f_max, target_date)

    excess_low = flux_low - pb17_low
    excess_high = flux_high - pb17_high
    casalog.post(f'excess_high {excess_high}', 'WARN')
    casalog.post(f'excess_low {excess_low}', 'WARN')
    casalog.post(f'excess_high/excess_low {excess_high/excess_low}', 'WARN')
    if excess_high/excess_low > 0 and np.isfinite( np.log(excess_high/excess_low) ):
        core_index = np.log(excess_high/excess_low) / np.log(f_max/f_min)
    else:
        core_index = 0.0


    components = vlacal_components( source, band, refdate )

    direction_rad = components['csys']['direction0']['crval']
    core_direction = 'J2000 {0} {1}'.format( qa.time('{0}rad'.format(direction_rad[0]), prec=12)[0], qa.angle('{0}rad'.format(direction_rad[1]),prec=12 )[0] )

    cl.close(log=False)
    cl.fromrecord(components['complist'])
    cl.addcomponent( flux=excess_flux, fluxunit='Jy', dir=core_direction, shape='point', freq='{0}Hz'.format(f_reference), spectrumtype='spectral index', index=core_index  )
    clrec = cl.torecord()
    cl.done()
    return clrec



@app.route("/", methods = ['GET']) 
def vlacal():

    if 'type' not in request.args:
        return 'type of query required', 400
    else:
        request_type = request.args['type']
        if request_type not in ['summary','fluxscale','components','scalefactor','setjy']:
            return f'request type {request_type} is not supported', 400


    if 'refdate' not in request.args:
        refdate = '99999'
    else:
        refdate = request.args['refdate']


    if request_type == 'summary': 
        return vlacal_summary( refdate )


    if 'source' not in request.args:
        if 'position' not in request.args:
            return 'source or position required', 400  
        else:
            position = request.args['position']
            if not vlacal_test_direction( position ):                
                return f'unable to convert position {position} to a CASA direction measure', 400
            else:
                tolerance = 10 #arcsec
                source = vlacal_position_to_source( position, tolerance )
                if source == '':
                    return f'no flux calibrator found within {tolerance} arcsec of position {position}', 400
    else:
        source = request.args['source'].upper()
        if source not in ['3C48','3C138','3C147','3C286']:
            return 'source name must match 3C48, 3C138, 3C147 or 3C286', 400

    if request_type == 'fluxscale':
        return vlacal_fluxscale( source, refdate ), 200


    if 'band' not in request.args:
        return 'band required', 400 
    else:
        band = request.args['band']
        if band not in list('PLSCXUKAQ'):
            return 'band must be P L S C X U K A or Q', 400


    if request_type == 'components':
        return vlacal_components( source, band, refdate ), 200


    if not 'obsdate' in request.args:
        return 'obsdate required', 400
    else:
        obsdate = request.args['obsdate']


    if not 'reffreq' in request.args:
        f_reference = vlacal_band_to_freqs( band )[0]
    else:
        f_reference = float( request.args['reffreq'] )*1e9


    if request_type == 'scalefactor':
        new_flux, excess_flux, pb17_flux, f_reference = vlacal_scalefactor( source, f_reference, obsdate, refdate )
        return  'Obtained flux of {0:.3f} Jy on mjd {1} for source {2} at reference frequency {3} GHz'.format(new_flux, obsdate, source, f_reference/1e9) + '<br/>' + \
                'This interpolated flux differs from Perley-Butler 2017 by {0:0.2f}%'.format( 100.0*excess_flux/pb17_flux ), 200


    if request_type == 'setjy':
        return json.loads( json.dumps( vlacal_cl_for_setjy( source, band, obsdate, refdate ), cls=NumpyEncoder ) ), 200

    return request.args






if __name__ == "__main__":
   app.run(debug=True,host='127.0.0.1',port=8080)


# example query
# http://192.168.1.14:5100/?type=setjy&source=3C48&band=C&obsdate=54321
