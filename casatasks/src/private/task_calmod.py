from urllib.parse import urlparse


def __is_valid_url_host(url):
    parsed = urlparse(url)
    return bool(parsed.netloc)


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
        raise RuntimeError('outfile must be specified')
    if not (source.strip() or direction.strip()):
        raise RuntimeError('Exactly one of source or direction must be specified')
    if source and direction:
        raise RuntimeError('Both source and direction may not be simultaneously specified')
    if not band:
        raise RuntimeError('band must be specified')
    if band.upper() not in ["P", "L", "S", "C", "X", "U", "K", "A", "Q"]:
        raise RuntimeError(f'band {band} not supported')
    if obsdate < 44239:
        raise RuntimeError('obsdate must be >= 44239')
    if refdate > 0 and refdate < 44239:
        raise RuntimeError('refdate must be <= 0 or >= 44239')
    if not hosts:
        raise RuntimeError('hosts must be specified')
    for h in hosts:
        if not __is_valid_url_host(h):
            raise RuntimeError(f'{h} is not a valid host expressed as a URL')

    
