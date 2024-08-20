import os
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplcursors

# beachball function
def beach(strike, dip, rake, linewidth=.5, facecolor='0.75', bgcolor='w', alpha=None,
            edgecolor='k', nofill=False, zorder=0):
    '''
    Return a beach ball as a collection.

    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    from matplotlib import collections

    colors, p = plot_dc(strike,dip,rake)
    col = collections.PatchCollection(p, match_original=False)
    col.set_facecolors([facecolor,bgcolor])
    if alpha is not None:
        col.set_alpha(alpha)

    col.set_edgecolor(edgecolor)
    col.set_linewidth(linewidth)
    col.set_zorder(zorder)

    return col 

#################### 3D Beachball to Cartesian xy ####################
def takeoff_az2xy(takeoff,azimuth,projection='stereographic'):
    '''
    Projects takeoff and azimuths onto focal sphere.
    Supported projections are 'stereographic' and 'lambert'

    Takeoffs:
        >90: downgoing, <90: upgoing
    Azimuths:
        # 0: North, 90: east, etc.
    '''
    takeoff=180-takeoff
    r=np.ones(len(takeoff))
    r[takeoff>90]=-1

    theta=np.deg2rad(takeoff)
    phi=np.deg2rad(90-azimuth)
    xyz=np.empty((3,len(takeoff)),dtype=float)
    xyz[0,:]=r*np.sin(theta)*np.cos(phi)
    xyz[1,:]=r*np.sin(theta)*np.sin(phi)
    xyz[2,:]=r*np.cos(theta)
    if projection=='stereographic':
        xy=xyz[:2,:]/(1+xyz[2,:])
    elif projection=='lambert':
        xy=xyz[:2,:]/np.sqrt(1+xyz[2,:])
    else:
        raise ValueError('Unknown projection: {}'.format(projection))
    return xy.T


#################### STRIKE DIP from NORMAL VECTOR ####################
def strike_dip(n, e, u):
    '''
    Finds strike and dip of plane given normal vector having components n, e,
    and u.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.

    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    r2d = 180 / np.pi
    if u < 0:
        n = -n
        e = -e
        u = -u

    strike = np.arctan2(e, n) * r2d
    strike = strike - 90
    while strike >= 360:
        strike = strike - 360
    while strike < 0:
        strike = strike + 360
    x = np.sqrt(np.power(n, 2) + np.power(e, 2))
    dip = np.arctan2(x, u) * r2d
    return (strike, dip)


# aux_plane function
def aux_plane(s1, d1, r1):
    '''
    Get Strike and dip of second plane.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.

    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    r2d = 180 / np.pi

    z = (s1 + 90) / r2d
    z2 = d1 / r2d
    z3 = r1 / r2d
    # slick vector in plane 1
    sl1 = -np.cos(z3) * np.cos(z) - np.sin(z3) * np.sin(z) * np.cos(z2)
    sl2 = np.cos(z3) * np.sin(z) - np.sin(z3) * np.cos(z) * np.cos(z2)
    sl3 = np.sin(z3) * np.sin(z2)
    (strike, dip) = strike_dip(sl2, sl1, sl3)

    n1 = np.sin(z) * np.sin(z2)  # normal vector to plane 1
    n2 = np.cos(z) * np.sin(z2)
    h1 = -sl2  # strike vector of plane 2
    h2 = sl1

    z = h1 * n1 + h2 * n2
    z = z / np.sqrt(h1 * h1 + h2 * h2)
    # we might get above 1.0 only due to floating point
    # precision. Clip for those cases.
    float64epsilon = 2.2204460492503131e-16
    if 1.0 < abs(z) < 1.0 + 100 * float64epsilon:
        z = np.copysign(1.0, z)
    z = np.arccos(z)
    rake = 0
    if sl3 > 0:
        rake = z * r2d
    if sl3 <= 0:
        rake = -z * r2d
    return (strike, dip, rake)


def pol2cart(th, r):
    '''
    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    x = r * np.cos(th)
    y = r * np.sin(th)
    return (x, y)


def xy2patch(x, y, res, xy):
    '''
    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    # check if one or two resolutions are specified (Circle or Ellipse)
    from matplotlib import path as mplpath
    from matplotlib import patches
    try:
        assert len(res) == 2
    except TypeError:
        res = (res, res)
    # transform into the Path coordinate system
    x = x * res[0] + xy[0]
    y = y * res[1] + xy[1]
    verts = list(zip(x.tolist(), y.tolist()))
    codes = [mplpath.Path.MOVETO]
    codes.extend([mplpath.Path.LINETO] * (len(x) - 2))
    codes.append(mplpath.Path.CLOSEPOLY)
    path = mplpath.Path(verts, codes)
    return patches.PathPatch(path)


def plot_dc(strike,dip,rake, size=100, xy=(0, 0), width=2):

    '''
    Uses one nodal plane of a double couple to draw a beach ball plot.

    :param ax: axis object of a matplotlib figure
    :param np1: :class:`~NodalPlane`

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.

    From ObsPy:
    https://docs.obspy.org/packages/autogen/obspy.imaging.beachball.html
    '''
    # check if one or two widths are specified (Circle or Ellipse)
    try:
        assert len(width) == 2
    except TypeError:
        width = (width, width)
    s_1 = strike
    d_1 = dip
    r_1 = rake
    D2R = np.pi / 180

    # make sure the rake is between -180 and 180
    m = 0
    if r_1 > 180:
        r_1 -= 180
        m = 1
    if r_1 < 0:
        r_1 += 180
        m = 1

    # Get azimuth and dip of second plane
    (s_2, d_2, _r_2) = aux_plane(s_1, d_1, r_1)

    d = size / 2

    # make sure the dip is between 0 and <90 (not exactly 90)
    if d_1 >= 90:
        d_1 = 89.9999
    if d_2 >= 90:
        d_2 = 89.9999

    # arange checked for numerical stability, np.pi is not multiple of 0.1
    phi = np.arange(0, np.pi, .01)
    # 
    l1 = np.sqrt(
        np.power(90 - d_1, 2) / (
            np.power(np.sin(phi), 2) +
            np.power(np.cos(phi), 2) * np.power(90 - d_1, 2) / np.power(90, 2)))
    l2 = np.sqrt(
        np.power(90 - d_2, 2) / (
            np.power(np.sin(phi), 2) + 
            np.power(np.cos(phi), 2) * np.power(90 - d_2, 2) / np.power(90, 2)))

    collect = []
    # plot paths, once for tension areas and once for pressure areas
    for m_ in ((m + 1) % 2, m):
        inc = 1
        (x_1, y_1) = pol2cart(phi + s_1 * D2R, l1)

        if m_ == 1:
            lo = s_1 - 180
            hi = s_2
            if lo > hi:
                inc = -1
            th1 = np.arange(s_1 - 180, s_2, inc)
            (xs_1, ys_1) = pol2cart(th1 * D2R, 90 * np.ones((1, len(th1))))
            (x_2, y_2) = pol2cart(phi + s_2 * D2R, l2)
            th2 = np.arange(s_2 + 180, s_1, -inc)
        else:
            hi = s_1 - 180
            lo = s_2 - 180
            if lo > hi:
                inc = -1
            th1 = np.arange(hi, lo, -inc)
            (xs_1, ys_1) = pol2cart(th1 * D2R, 90 * np.ones((1, len(th1))))
            (x_2, y_2) = pol2cart(phi + s_2 * D2R, l2)
            x_2 = x_2[::-1]
            y_2 = y_2[::-1]
            th2 = np.arange(s_2, s_1, inc)
        (xs_2, ys_2) = pol2cart(th2 * D2R, 90 * np.ones((1, len(th2))))
        x = np.concatenate((x_1, xs_1[0], x_2, xs_2[0]))
        y = np.concatenate((y_1, ys_1[0], y_2, ys_2[0]))

        x = x * d / 90
        y = y * d / 90

        # calculate resolution
        res = [value / float(size) for value in width]

        # construct the patch
        collect.append(xy2patch(y, x, res, xy))
    return ['b', 'w'], collect


#################### PLOT MECHANISM (MAIN FUNCTION) ####################
def plot_mech(
    mech_df,
    pol_df,
    pol_info_df, # this also inclus same information as pol_df + more
    alt_mech_df=None,
    **kwargs):

    '''
    Creates a plot of the mechanism with the P-polarity.
    '''

    # Optional parameters
    acceptable_sdr=kwargs.get('acceptable_sdr', False)
    subplots=kwargs.get('subplots', None) # list of 2 integers
    figsize=kwargs.get('figsize', (6,6))

    # Check if the dataframes are empty, if they are, return None
    if (len(pol_df)==0) | (len(mech_df)==0):
        return
    event_id=str(pol_df.event_id.values[0])
    mech_quality=mech_df['quality'].values[0]
    
    # Convert takeoff and azimuth to x,y coordinates
    xy=takeoff_az2xy(pol_info_df['takeoff'].values,pol_info_df['azimuth'].values)

    # Get the indices of the P-polarity
    up_ind=np.where(pol_info_df['p_polarity']>0)[0]      ## use pol_info_df instead of pol_df if needed
    down_ind=np.where(pol_info_df['p_polarity']<0)[0]
    zero_ind=np.where(pol_info_df['p_polarity']==0)[0]

    # Draw the preferred mechanism
    beach1=beach(mech_df.strike.values[0], mech_df.dip.values[0], mech_df.rake.values[0],
                facecolor='0.75',linewidth=0.5,zorder=0)

    # Add the beachballs to a plot axis
    if subplots:
        fig, axes = plt.subplots(subplots[0],subplots[1],figsize=figsize)
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(1,1,figsize=figsize)
        axes = [axes]

    time = mech_df.time.values[0] if 'time' in mech_df.columns else ''
    ax = axes[0]
    ax.add_collection(beach1)
    ax.set_title(f'{event_id} Quality: {mech_quality}\nTime: {str(time)[:16]}')
    
    # Draw the acceptable mechanisms
    beach_accept=[]
    if acceptable_sdr:      # if True, plot the acceptable mechanisms
        for i in range(len(alt_mech_df)):
            beach_accept.append(
                beach(
                    alt_mech_df.strike.values[i], alt_mech_df.dip.values[i], alt_mech_df.rake.values[i],
                    facecolor='None',bgcolor='None',linewidth=0.25,edgecolor='0.5',zorder=1))
    
    # Add the acceptable mechanisms to the plot axis
    if beach_accept: # not empty
        beach1_2=beach(mech_df.strike.values[0], mech_df.dip.values[0], mech_df.rake.values[0],
                edgecolor='k',linewidth=0.5,facecolor='None',bgcolor='None',zorder=2)
        for beach_ in beach_accept:
            ax.add_collection(copy.deepcopy(beach_))
        ax.add_collection(beach1_2)
        
    # Add the P-polarity (Up/Down) to the plot axis
    sc_up = ax.scatter(xy[up_ind, 0], xy[up_ind, 1], marker='+', s=50, linewidths=.6, c='k', picker=True, zorder=3) if len(up_ind) > 0 else None
    sc_down = ax.scatter(xy[down_ind, 0], xy[down_ind, 1], marker='o', s=50, linewidths=.5, edgecolor='k', facecolor='None', picker=True, zorder=3) if len(down_ind) > 0 else None
    sc_zero = ax.scatter(xy[zero_ind, 0], xy[zero_ind, 1], marker='d', s=50, linewidths=.6, c='k', picker=True, zorder=3) if len(zero_ind) > 0 else None

    # Add Station names to the plot axis
    sta_names = pol_info_df['sta_code'].str.split('.').str[0].values
    for i in range(len(sta_names)):
        ax.text(xy[i, 0], xy[i, 1], sta_names[i],fontsize=6,zorder=4)
    
    # Set the plot axis limits and remove the axis labels
    ax.set_xlim([-1.01,1.01])
    ax.set_ylim([-1.01,1.01])
    ax.set_aspect(1)
    ax.set_xticks([])
    ax.set_yticks([])
    # tighen the layout
    plt.tight_layout()

    return fig, axes, sc_up, sc_down, sc_zero, xy, event_id

def get_nearest_station(xy_clicked, xy_data, pol_info_df):
    '''
    Get the nearest station to the clicked point.
    '''
    distances = np.linalg.norm(xy_clicked - xy_data, axis=1)
    nearest_index = np.argmin(distances)
    # Reset the index of the pol_info_df so that for next events, the index is not messed up
    pol_info_df.reset_index(drop=True, inplace=True)
    return pol_info_df.loc[nearest_index, 'sta_code']

def get_exact_station(xy_clicked, xy_data, pol_info_df, tolarance=1e-3):
    '''
    Get the nearest station to the clicked point.
    '''
    point = np.array(xy_clicked)        # convert click object to numpy array
    index = np.where(np.all(np.isclose(xy_data, point, atol=tolarance), axis=1))[0]     # find the index of the clicked point in xy_data
    if index.size > 0:      # if the index is not empty
        corr_sta_code = pol_info_df['sta_code'].values[index[0]]     # get the station code
        return corr_sta_code
    else:
        return None