def get_nps_soil_field_name(prefix, startDepth, endDepth):
    '''
    Given a starting and ending depths, generate the nps_int field name
    of a variable (e.g. SM040100 for soil moisture starting 40cm below
    the surface and ending 100cm below the surface)
    '''
    assert startDepth < 1000 and endDepth < 1000
    return "{0}{1:03d}{2:03d}".format(prefix, startDepth, endDepth)

def lis_name_to_nps_prefix(lisname):
    '''
    Given a variable name in LIS, get the corresponding NPS prefix
    '''
    if lisname == 'SoilMoist_tavg':
        return 'SM'
    elif lisname == 'SoilTemp_tavg':
        return 'ST'
    elif lisname == 'SoilWet_inst':
        return 'SW'
    else:
        raise Exception("Do not know NPS prefix for LIS variable {}"
                        .format(lisname))

def get_nps_soil_field_description(prefix, startDepth, endDepth):
    '''
    Given a prefix (of 'sm' or 'st'), a starting depth and an ending
    depth, return a description in the same format the NPS uses,
    e.g. given ('sm', 10, 40), return "Soil Moist 100-200 cm below gr layer"
    '''
    if prefix.lower() == 'sm':
        return "Soil Moist {}-{} cm below gr layer".format(startDepth, endDepth)
    elif prefix.lower() == 'st':
        return "T {}-{} cm below ground layer (Upper)".format(startDepth, endDepth)
    else:
        raise Exception("Prefix must be either 'sm' or 'st'")
