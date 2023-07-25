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
