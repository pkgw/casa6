from casatools import componentlist, measures
import json
from urllib import request
from urllib.parse import urlparse


def __is_valid_url_host(url):
    parsed = urlparse(url)
    return bool(parsed.netloc)

def __query(url):
    req = request.Request(url)
    with request.urlopen(req) as response:
        if response.status != 200:
            return None
        components = json.loads(response.read().decode('utf-8'))
    return components


def calmod(
    outfile, source, direction, band, obsdate, refdate, hosts
):
    print(f'outfile is "{outfile}"')
    print(f'source is "{source}"')
    print(f'direction is "{direction}"')
    print(f'band is "{band}"')
    print(f'obsdate is {obsdate}')
    print(f'refdate is {refdate}')
    print(f'hosts is {hosts}')
    if not outfile.strip():
        raise ValueError('outfile must be specified')
    if not (source.strip() or direction.strip()):
        raise ValueError('Exactly one of source or direction must be specified')
    if source and direction:
        raise ValueError('Both source and direction may not be simultaneously specified')
    if source and source.upper() not in ("3C48", "3C286", "3C138", "3C147"):
        raise ValueError(f'Unsupported calibrator {source}')
    if direction:
        dirstr = direction.split(' ')
        if not (len(dirstr) == 3 and measures().direction(dirstr[0], dirstr[1], dirstr[2])):
            raise ValueError(f'Illegal direction specification {direction}')
    src_or_dir = source if source else direction
    if not band:
        raise ValueError('band must be specified')
    if band.upper() not in ["P", "L", "S", "C", "X", "U", "K", "A", "Q"]:
        raise ValueError(f'band {band} not supported')
    if obsdate > 0 and obsdate < 44239:
        raise ValueError('obsdate must be <= 0 or >= 44239')
    if refdate > 0 and refdate < 44239:
        raise ValueError('refdate must be <= 0 or >= 44239')
    if not hosts:
        raise ValueError('hosts must be specified')
    for h in hosts:
        if not __is_valid_url_host(h):
            raise ValueError(f'{h} is not a valid host expressed as a URL')
        url = f'{h}/components/{src_or_dir}/{band}'
        if obsdate > 0:
            url = f'{url}/{obsdate}'
        else:
            print(f'Querying most recent data for {src_or_dir}')
        print (f'Trying {url} ...')
        components = __query(url)
        if components:
            break
    if not components:
        raise RuntimeError('All URLs failed to return a component list')
    cl = componentlist()
    cl.fromrecord(components['complist'])
    cl.rename(outfile)
    print(f'component list {outfile} has been written')
    cl.done()
    
