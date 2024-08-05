'''
Functions used to compute earthquake focal mechanisms.

Significant portions of the functions in this file are based on the Fortran HASH
code originally written by Jeanne L. Hardebeck & Peter M. Shearer, and all of it
is inspired by their work. Please cite the appropriate references if you use
this code.
'''

# Standard libraries
import os

# External libraries
import numpy as np
import pandas as pd
rng=np.random.default_rng(123) # Used to produce reproducable bootstrapped results.


def lookup_takeoff(table,perturbed_origin_depth_km,sr_dist_km,look_dep,look_del,deptab,delttab,num_velocity_models=1):
	'''
	Given a hypocentral depths and source-receiver distances, queries the takeoff angles
	using the lookup table
	Input:
		table: Lookup table produced by create_takeoff_table(), 3d array.
		perturbed_origin_depth_km: perturbed origin depths in km produced by perturb_eq_locations, 2d array.
		sr_dist_km: source-receiver distances in km produced by perturb_eq_locations, 2d array.
		look_dep: minimum source depth, maximum, and interval for the lookup table, list
		look_del: minimum source-station distance, maximum, and interval for the lookup table, list
		deptab: array of depths (km) produced by create_lookup_table()
		delttab: array of distances (km) produced by create_lookup_table()
		num_velocity_models: the number of velocity models used.
	Output:
		takeoff: the corresponding takeoff angles
	'''
	# Determines the lookup table cell for permiated earthquake depths
	id1=((perturbed_origin_depth_km-look_dep[0])/look_dep[2]).astype(int)
	id2=id1+1

	# If perturbed locations are outside of the modeled range, sets those perturbations to the model edge
	dist_flag=sr_dist_km<delttab[0]
	if np.any(dist_flag):
		print('*WARNING: Perturbed epicenters(s) are closer than lookup table depth range. Setting these perturbed locations to {} km *'.format(delttab[0]))
		sr_dist_km[dist_flag]=delttab[0]
	dist_flag=sr_dist_km>delttab[-1]
	if np.any(dist_flag):
		print('*WARNING: Perturbed epicenters(s) are further than lookup table depth range. Setting these perturbed locations to {} km *'.format(delttab[-1]))
		sr_dist_km[dist_flag]=delttab[-1]

	# Determines the lookup table cell for permiated earthquake epicenters
	ix1=((sr_dist_km-look_del[0])/look_del[2]).astype(int)
	ix2=ix1+1

	# Determines which velocity model to use for each trial
	if num_velocity_models>1:
		randomized_vm_ind=rng.integers(low=0, high=num_velocity_models,size=perturbed_origin_depth_km.shape[-1])
		randomized_vm_ind[0]=0
	else:
		randomized_vm_ind=np.zeros(perturbed_origin_depth_km.shape[-1],dtype=int)

	# Uses the lookup table to determine the travel times / takeoff angles
	xfrac=(sr_dist_km-delttab[ix1])/(delttab[ix2]-delttab[ix1])

	if id1.ndim==1:
		id1=np.repeat(id1[np.newaxis,:],len(sr_dist_km),axis=0)
		id2=np.repeat(id2[np.newaxis,:],len(sr_dist_km),axis=0)

	t1=table[ix1,id1,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)]+\
		xfrac*(table[ix2,id1,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)]-
		table[ix1,id1,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)])

	# ariful islam
	# if id2 > table.shape[1]:
	# 	print(f"Warning: id2 ({id2}) is out of bounds for table with shape {table.shape}. Setting id2 to {table.shape[1] - 1}.")
	# 	id2 = table.shape[1] - 1
	#/ ariful islam

	t2=table[ix1,id2,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)]+\
		xfrac*(table[ix2,id2,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)]-
		table[ix1,id2,np.repeat([randomized_vm_ind],len(sr_dist_km),axis=0)])

	dfrac=(perturbed_origin_depth_km-deptab[id1])/(deptab[id2]-deptab[id1])
	takeoff=t1+dfrac*(t2-t1)

	return takeoff


def perturb_eq_locations(event_pol_df,look_dep,perturb_epicentral_location,nmc=1):
	'''
	Randomly perturbs the hypocentral locations.
	Input:
		event_pol_df: polarity dataframe
		look_dep: minimum source depth, maximum, and interval for the lookup table, list
		perturb_epicentral_location: flag to determine if epicentral locations should be perturbed, boolean
		nmc: number of trials, integer
	Output:
		perturbed_origin_depth_km: Perturbed hypocentral depths
		sr_dist_km: Perturbed source-receiver distances
		sr_azimuth: Perturbed source-receiver azimuths
	'''
	if 'event_id2' in event_pol_df:
		unique_event_id2,unique_ind,unique_inv=np.unique(event_pol_df.event_id2.values,return_index=True,return_inverse=True)
		num_unique_events=len(unique_event_id2)
	else:
		num_unique_events=1
		unique_ind=[0]
		unique_inv=np.zeros(len(event_pol_df),dtype=int)

	rng_z=rng.normal(size=(nmc,num_unique_events))

	rng_z[0]=0
	perturbed_origin_depth_km=event_pol_df.iloc[unique_ind]['origin_depth_km'].values+rng_z*event_pol_df.iloc[unique_ind]['vert_uncert_km'].values
	# if 'event_id2' in event_pol_df:
	perturbed_origin_depth_km=perturbed_origin_depth_km[:,unique_inv]

	if perturb_epicentral_location: # Perturbs horizontal earthquake locations
		perturb_rand_azimuth=rng.integers(low=0, high=360, size=(nmc,num_unique_events))
		rng_h=rng.normal(size=(nmc,num_unique_events))
		rng_h[0]=0

		# Earthquake lon/lat in deg convertered to radians
		eq_lon_r=np.deg2rad(event_pol_df.iloc[unique_ind]['origin_lon'].values)
		eq_lat_r=np.deg2rad(event_pol_df.iloc[unique_ind]['origin_lat'].values)
		perturb_rand_horz_km=event_pol_df.iloc[unique_ind]['horz_uncert_km'].values*rng_h

		earth_radius_km=6371 # approx radius of earth in km
		eq_pert_lat_r=np.arcsin(np.sin(eq_lat_r)*np.cos(perturb_rand_horz_km/earth_radius_km)+
									np.cos(eq_lat_r)*np.sin(perturb_rand_horz_km/earth_radius_km)*np.cos(perturb_rand_azimuth))
		eq_pert_lon_r=eq_lon_r+np.arctan2(np.sin(perturb_rand_azimuth)*np.sin(perturb_rand_horz_km/earth_radius_km)*np.cos(eq_lat_r),
					np.cos(perturb_rand_horz_km/earth_radius_km)-np.sin(eq_lat_r)*np.sin(eq_pert_lat_r))

		aspect=np.cos(np.median(eq_pert_lat_r[:,0]))
		qlon=np.rad2deg(eq_pert_lon_r)
		qlat=np.rad2deg(eq_pert_lat_r)

	else:
		qlon=np.atleast_2d(event_pol_df.iloc[unique_ind]['origin_lon'].values)
		qlat=np.atleast_2d(event_pol_df.iloc[unique_ind]['origin_lat'].values)
		aspect=np.cos(np.deg2rad(np.median(qlat)))

	qlon=qlon[:,unique_inv]
	qlat=qlat[:,unique_inv]

	# Ensures all perturbed depths are within the modeled range
	qdep_v_flag=perturbed_origin_depth_km<look_dep[0]

	if np.any(qdep_v_flag):
		print('*WARNING: Perturbed locations for event_id \'{}\' are shallower than lookup table depth range ({}). Setting these perturbed depths to {} km.'.format(event_pol_df.iloc[unique_ind[0]]['event_id'],look_dep[0],look_dep[0]))
		perturbed_origin_depth_km[qdep_v_flag]=look_dep[0]
	qdep_v_flag=perturbed_origin_depth_km>look_dep[1]
	if np.any(qdep_v_flag):
		print('*WARNING: Perturbed locations for event_id \'{}\' are deeper than lookup table depth range ({}). Setting these perturbed depths to {} km.'.format(event_pol_df.iloc[unique_ind[0]]['event_id'],look_dep[1],look_dep[1]))
		# print('*WARNING: Perturbed event(s) are deeper than lookup table depth range. Setting these perturbed depths to {} km.'.format(look_dep[1]))
		perturbed_origin_depth_km[qdep_v_flag]=look_dep[1]-.01

	# Approximates the source-receiver distances and azimuths
	flon=event_pol_df['station_lon'].values
	flat=event_pol_df['station_lat'].values
	dx=(flon-qlon)*111.2*aspect
	dy=(flat-qlat)*111.2
	sr_dist_km=np.sqrt(dx**2+dy**2)
	sr_azimuth=90.-np.rad2deg(np.arctan2(dy,dx))
	sr_azimuth[sr_azimuth<0]+=360

	return perturbed_origin_depth_km.T,sr_dist_km.T,sr_azimuth.T


def perturb_azimuth_takeoff(event_pol_df,nmc):
	'''
	Perturb predetermined azimuth and takeoff angles using the given uncertainties.
	Input:
		event_pol_df: polarity dataframe
		nmc: number of trials, integer
	Output:
		sr_azimuth: perturbed source-receiver azimuths
		takeoff: perturbed takeoff angles
	'''
	rng_az=rng.normal(size=(len(event_pol_df),nmc))
	rng_az[:,0]=0

	rng_takeoff=rng.normal(size=(len(event_pol_df),nmc))
	rng_takeoff[:,0]=0

	sr_azimuth=event_pol_df['azimuth'].values[:,np.newaxis]+event_pol_df['azimuth_uncertainty'].values[:,np.newaxis]*rng_az
	sr_azimuth[sr_azimuth<0]+=360
	sr_azimuth[sr_azimuth>360]-=360

	takeoff=event_pol_df['takeoff'].values[:,np.newaxis]+event_pol_df['takeoff_uncertainty'].values[:,np.newaxis]*rng_takeoff
	takeoff[takeoff<=0]=1e-6
	takeoff[takeoff>180]=180

	return sr_azimuth,takeoff


def cartesian_transform(the,phi,r):
	'''
	Transforms spherical coordinates to cartesian (theta in degrees)
	'''
	x=r*np.sin(np.deg2rad(the))*np.cos(np.deg2rad(phi))
	y=r*np.sin(np.deg2rad(the))*np.sin(np.deg2rad(phi))
	z=-r*np.cos(np.deg2rad(the))
	return x,y,z


def vector_from_sdr(strike,dip,rake):
	'''
	Gets fault normal vector (fnorm) and slip vector (slip) from [strike,dip,rake].
	Uses (x,y,z) coordinate system with x=north, y=east, z=down
		Reference:  Aki and Richards, p. 115
	Based on code from HASH (Hardebeck & Shearer, 2002).
	'''
	fnorm=np.zeros(3)
	slip=np.zeros(3)

	phi=np.deg2rad(strike)
	delt=np.deg2rad(dip)
	lam=np.deg2rad(rake)

	fnorm[0]=-np.sin(delt)*np.sin(phi)
	fnorm[1]=np.sin(delt)*np.cos(phi)
	fnorm[2]=-np.cos(delt)
	slip[0]=np.cos(lam)*np.cos(phi)+np.cos(delt)*np.sin(lam)*np.sin(phi)
	slip[1]=np.cos(lam)*np.sin(phi)-np.cos(delt)*np.sin(lam)*np.cos(phi)
	slip[2]=-np.sin(lam)*np.sin(delt)

	return fnorm,slip


def sdr_from_vector(faultnorms,slips):
	'''
	Gets [strike,dip,rake] from fault normal vectors (faultnorms) and slip vector (slips).
	faultnorms and slips should each be a 3-by-n array.

	Uses (x,y,z) coordinate system with x=north, y=east, z=down
		Reference:  Aki and Richards, p. 115
	Based on code from HASH (Hardebeck & Shearer, 2002).
	'''
	num_vect=faultnorms.shape[1]
	if faultnorms.shape != slips.shape:
		raise ValueError('***Error: shape of faultnorms and slips must be the same')
	if faultnorms.shape[0]!=3:
		raise ValueError('***Error: faultnorms and slips must be an array of shape 3-by-n')

	phi=np.zeros(num_vect)
	delt=np.zeros(num_vect)
	lam=np.zeros(num_vect)

	undef_flag=(1-np.abs(faultnorms[2,:]))<=(1e-7)
	if np.any(undef_flag):
		undef_ind=np.where(undef_flag)[0]
		print('*sdr_from_vector warning, horz fault, strike undefined')
		phi[undef_ind]=np.arctan2(-slips[0,undef_ind],slips[1,undef_ind])
		clam=np.cos(phi[undef_ind])*slips[0,undef_ind]+np.sin(phi[undef_ind])*slips[1,undef_ind]
		slam=np.sin(phi[undef_ind])*slips[0,undef_ind]-np.cos(phi[undef_ind])*slips[1,undef_ind]
		lam[undef_ind]=np.arctan2(slam,clam)
	if np.any(~undef_flag):
		def_ind=np.where(~undef_flag)[0]

		phi[def_ind]=np.arctan2(-faultnorms[0,def_ind],faultnorms[1,def_ind])
		a=np.sqrt(faultnorms[0,def_ind]*faultnorms[0,def_ind]+faultnorms[1,def_ind]*faultnorms[1,def_ind])
		delt[def_ind]=np.arctan2(a,-faultnorms[2,def_ind])
		clam=np.cos(phi[def_ind])*slips[0,def_ind]+np.sin(phi[def_ind])*slips[1,def_ind]
		slam=-slips[2,def_ind]/np.sin(delt[def_ind])
		lam[def_ind]=np.arctan2(slam,clam)

		tmp_ind=def_ind[np.where(delt[def_ind]>(0.5*np.pi))[0]]
		if len(tmp_ind)>0:
			delt[tmp_ind]=np.pi-delt[tmp_ind]
			phi[tmp_ind]=phi[tmp_ind]+np.pi
			lam[tmp_ind]=-lam[tmp_ind]

	strike=np.rad2deg(phi)
	dip=np.rad2deg(delt)
	rake=np.rad2deg(lam)

	strike[strike<0]+=360
	rake[rake<-180]+=360
	rake[rake>180]-=360

	return strike,dip,rake


def create_takeoff_table(vmodel_depthvp,deptab,delttab,nump,nx0,nd0,takeoff_az_precision):
	'''
	Creates tables of takeoff angles given a 1D velocity model.
	Based on code from HASH (Hardebeck & Shearer, 2002).
	Input:
		vmodel_depthvp: depth (km) and Vp (km/s), 2d array
		deptab: deptab: array of depths (km)
		delttab: array of distances (km)
		nump: number of rays traced
		nx0: maximum source-station distance bins for look-up tables
		nd0: maximum source depth bins for look-up tables
		takeoff_az_precision: number of decimal places to round resulting table
	Output:
		table: produced lookup table of takeoff angles
	'''
	if vmodel_depthvp[:,0][-1]<deptab[-1]:
		vmodel_depthvp=np.vstack((vmodel_depthvp,vmodel_depthvp[-1,:]))
		vmodel_depthvp[-1,0]=deptab[-1]+1
		vmodel_depthvp[-1,1]=vmodel_depthvp[-1,1]+0.001

	table=np.zeros((nx0,len(deptab)))-999

	ndel=len(delttab)
	ndep=len(deptab)
	pmin=0

	z=vmodel_depthvp[:,0]
	alpha=vmodel_depthvp[:,1]

	npts=len(z)
	z=np.hstack((z,z[npts-1]))
	alpha=np.hstack((alpha,alpha[npts-1]))

	for i in range(npts-1,0,-1):
		for idep in range(ndep-1,-1,-1):
			if (z[i-1]<=(deptab[idep]-0.1)) & (z[i]>=(deptab[idep]+.1)):
				z=np.insert(z,i,z[i-1])
				alpha=np.insert(alpha,i,alpha[i-1])
				z[i]=deptab[idep]
				frac=(z[i]-z[i-1])/(z[i+1]-z[i-1])
				alpha[i]=alpha[i-1]+frac*(alpha[i+1]-alpha[i-1])

	slow=1/alpha
	pmax=slow[0]
	pstep=(pmax-pmin)/nump

	# do P-wave ray tracing
	npmax=int((pmax+pstep/2-pmin)/pstep)+1

	depxcor=np.zeros((npmax,nd0))
	depucor=np.zeros((npmax,nd0))
	deptcor=np.zeros((npmax,nd0))

	tmp_ind=np.where(deptab==0)[0]

	tmp_ind=np.where(deptab!=0)[0]
	depxcor[:,tmp_ind]=-999
	deptcor[:,tmp_ind]=-999

	ptab=np.linspace(pmin,pmin+pstep*(npmax-1),num=npmax)

	h_array=z[1:]-z[:-1]
	utop=slow[:-1]
	ubot=slow[1:]

	'''LAYERTRACE equivalent'''
	dx=np.zeros((npmax,len(utop)))
	dt=np.zeros((npmax,len(utop)))
	irtr=np.zeros((npmax,len(utop)),dtype=int)

	qs=np.zeros((npmax,len(utop)));qs[:]=np.nan
	qr=np.zeros((npmax,len(utop)));qr[:]=np.nan
	ytop=utop-ptab[:, np.newaxis]
	ytop_pos_flag=ytop>0
	qs[ytop_pos_flag]=ytop[ytop_pos_flag]*(utop+ptab[:, np.newaxis])[ytop_pos_flag]
	qs[ytop_pos_flag]=np.sqrt(qs[ytop_pos_flag])

	qr=np.arctan2(qs,ptab[:, np.newaxis])

	b=np.ma.divide(-np.log(ubot/utop),h_array)
	b=b.filled(np.nan)

	# integral at upper limit, 1/b factor omitted until end
	etau=qs-qr*ptab[:, np.newaxis]
	ex=qr

	# check lower limit to see if we have turning point
	ybot=ubot-ptab[:, np.newaxis]
	# if turning point, then no contribution from bottom point
	y_subzero_flag=(ybot<=0)
	y_greaterzero_flag=(ybot>0)
	irtr[y_subzero_flag ]=2
	irtr[y_greaterzero_flag]=1

	irtr[~ytop_pos_flag]=0

	dx[y_subzero_flag]=ex[y_subzero_flag]
	dx=dx/b
	dtau=etau/b
	dt[y_subzero_flag]=dtau[y_subzero_flag]+(dx*ptab[:, np.newaxis])[y_subzero_flag] # converts tau to t

	q=np.zeros(ybot.shape);q[:]=np.nan
	q[y_greaterzero_flag]=ybot[y_greaterzero_flag]*((ubot+ptab[:, np.newaxis])[y_greaterzero_flag])
	qs=np.sqrt(q)

	qr=np.arctan2(qs,ptab[:,np.newaxis])
	etau=etau-qs+ptab[:, np.newaxis]*qr
	ex=ex-qr

	exb=ex/b
	dtau=etau/b
	dx[y_greaterzero_flag]=exb[y_greaterzero_flag]
	dt[y_greaterzero_flag]=dtau[y_greaterzero_flag]+(exb*ptab[:, np.newaxis])[y_greaterzero_flag]

	# Ensures values after ray has turned are nan
	x=((irtr==0) | (irtr==2))
	idx=np.arange(npmax), x.argmax(axis=1)
	tmp=(x[idx]==True)

	idx_1=idx[0][tmp],idx[1][tmp]+1
	tmp_1=idx_1[1]<(len(utop)-1)
	idx_1=idx_1[0][tmp_1],idx_1[1][tmp_1]
	for xx in range(len(idx_1[0])):
		row=idx_1[0][xx]
		col=idx_1[1][xx]
		dx[row,col:]=np.nan
		dt[row,col:]=np.nan
	deltab=np.nansum(dx,axis=1)*2
	tttab=np.nansum(dt,axis=1)*2

	idx_2=idx[0][tmp],idx[1][tmp]
	tmp_2=idx_2[1]<(len(utop)-1)
	idx_2=idx_2[0][tmp_2],idx_2[1][tmp_2]
	for xx in range(len(idx_2[0])):
		row=idx_2[0][xx]
		col=idx_2[1][xx]
		dx[row,col:]=np.nan
		dt[row,col:]=np.nan

	depxcor=np.cumsum(dx,axis=1)
	deptcor=np.cumsum(dt,axis=1)
	output_col_ind=np.where(np.isin(z,deptab))[0]-1
	depxcor=depxcor[:,output_col_ind]
	deptcor=deptcor[:,output_col_ind]
	depxcor[:,0]=0
	deptcor[:,0]=0
	depucor[:]=slow[output_col_ind+1]
	depucor[np.isnan(depxcor)]=-999
	depxcor[np.isnan(depxcor)]=-999
	deptcor[np.isnan(deptcor)]=-999

	x=np.diff(depxcor,axis=0)<=0
	idx=x.argmax(axis=0)+1
	tmp=(x[idx,np.arange(nd0)]==False)
	idx[tmp]=npmax-1

	for idep in range(ndep):
		# upgoing rays from source
		xsave_up=depxcor[:(idx[idep]),idep]
		tsave_up=deptcor[:(idx[idep]),idep]
		usave_up=depucor[:(idx[idep]),idep]
		psave_up=-1*ptab[:(idx[idep])]

		# downgoing rays from source
		down_idx=np.where((depxcor[:,idep]!=-999) & (deltab!=-999))[0][::-1]
		xsave_down=(deltab[down_idx]-depxcor[down_idx,idep])
		tsave_down=(tttab[down_idx]-deptcor[down_idx,idep])
		usave_down=depucor[down_idx,idep]
		psave_down=ptab[down_idx]

		# Merges upgoing and downgoing ray arrays
		xsave=np.hstack([xsave_up,xsave_down])
		tsave=np.hstack([tsave_up,tsave_down])
		usave=np.hstack([usave_up,usave_down])
		psave=np.hstack([psave_up,psave_down])

		scr1=np.zeros(ndel)
		for idel in range(1,ndel):
			del_x=delttab[idel]
			ind=np.where( (xsave[:-1]<=del_x) & (xsave[1:]>=del_x) )[0]+1

			frac=(del_x-xsave[ind-1])/(xsave[ind]-xsave[ind-1])
			t1=tsave[ind-1]+frac*(tsave[ind]-tsave[ind-1])

			min_ind=ind[np.argmin(t1)]

			scr1[idel]=psave[min_ind]/usave[min_ind]
		angle=np.rad2deg(np.arcsin(scr1))

		angle_flag=angle>=0
		angle*=-1
		angle[angle_flag]+=180

		table[:,idep]=angle

	if delttab[0]==0:
		table[0,:]=0. # straight up at zero range
	table=np.round(table,takeoff_az_precision)
	return table


def dir_cos_setup(p_dict):
	'''
	Sets up array with direction cosines for all coordinate transformations
	Input:
		p_dict: Parameter values created in SKHASH.py, dictionary
	Output:
		dir_cos_dict: Coordinate transformation variables, dictionary
	'''
	ntab=180
	astep=1/ntab

	dir_cos_dict={
		'ncoor':0,
		'thetable':np.empty(0),
		'phitable':np.empty(0),
		'b1':np.empty(0),
		'b2':np.empty(0),
		'b3':np.empty(0),
		'amptable':np.empty(0)
	}

	num_izeta=int(np.floor(179.9/p_dict['dang']))
	the=np.arange(0,90.001,p_dict['dang']) # angles for grid search

	rthe=np.deg2rad(the)
	costhe=np.cos(rthe)
	sinthe=np.sin(rthe)
	fnumang=360./p_dict['dang']
	dphi=np.round(fnumang*np.sin(rthe))
	dphi[dphi!=0]=360./dphi[dphi!=0]
	dphi[dphi==0]=10000.
	num_iphi=np.floor(359.9/dphi).astype(int)

	dir_cos_dict['ncoor']=np.sum(num_iphi+1)*(num_izeta+1)

	dir_cos_dict['b1']=np.zeros((3,dir_cos_dict['ncoor']))
	dir_cos_dict['b2']=np.zeros((3,dir_cos_dict['ncoor']))
	dir_cos_dict['b3']=np.zeros((3,dir_cos_dict['ncoor']))
	irot=0

	for ithe in range(len(the)):
		bb1=np.zeros(3)
		bb3=np.zeros(3)
		for iphi in range(0,num_iphi[ithe]+1):
			phi=iphi*dphi[ithe]
			rphi=np.deg2rad(phi)
			cosphi=np.cos(rphi)
			sinphi=np.sin(rphi)
			bb3[2]=costhe[ithe]
			bb3[0]=sinthe[ithe]*cosphi
			bb3[1]=sinthe[ithe]*sinphi

			bb1[2]=-sinthe[ithe]
			bb1[0]=costhe[ithe]*cosphi
			bb1[1]=costhe[ithe]*sinphi

			bb2=np.cross(bb1,bb3)*-1

			for izeta in range(0,num_izeta+1):
				if (irot>dir_cos_dict['ncoor']):
					raise ValueError('***FOCAL error: # of rotations too big')
				zeta=izeta*p_dict['dang']
				rzeta=np.deg2rad(zeta)
				coszeta=np.cos(rzeta)
				sinzeta=np.sin(rzeta)

				dir_cos_dict['b3'][2,irot]=bb3[2]
				dir_cos_dict['b3'][0,irot]=bb3[0]
				dir_cos_dict['b3'][1,irot]=bb3[1]
				dir_cos_dict['b1'][0,irot]=bb1[0]*coszeta+bb2[0]*sinzeta
				dir_cos_dict['b1'][1,irot]=bb1[1]*coszeta+bb2[1]*sinzeta
				dir_cos_dict['b1'][2,irot]=bb1[2]*coszeta+bb2[2]*sinzeta
				dir_cos_dict['b2'][0,irot]=bb2[0]*coszeta-bb1[0]*sinzeta
				dir_cos_dict['b2'][1,irot]=bb2[1]*coszeta-bb1[1]*sinzeta
				dir_cos_dict['b2'][2,irot]=bb2[2]*coszeta-bb1[2]*sinzeta
				irot=irot+1

	bbb=-1.+np.arange(0,2*ntab+1)*astep
	_x,_y=np.meshgrid(bbb,bbb)
	dir_cos_dict['thetable']=np.arccos(bbb)
	dir_cos_dict['phitable']=np.arctan2(_x,_y).T
	dir_cos_dict['phitable'][dir_cos_dict['phitable']<0]+=2*np.pi

	# Creates amptable
	if p_dict['ampfile']:
		dir_cos_dict['amptable']=np.zeros((2,ntab,2*ntab))
		phi=np.arange(0,2*ntab)*np.pi*astep
		theta=np.arange(0,ntab)*np.pi*astep

		_x,_y=np.meshgrid(theta,phi)
		dir_cos_dict['amptable'][0,:,:]=np.abs(np.sin(2*_x)*np.cos(_y)).T

		s1=np.cos(2*_x)*np.cos(_y)
		s2=-np.cos(_x)*np.sin(_y)
		dir_cos_dict['amptable'][1,:,:]=np.sqrt(s1*s1+s2*s2).T

		if p_dict['min_amp']>0:
			dir_cos_dict['amptable'][dir_cos_dict['amptable']<p_dict['min_amp']]=0.0

	return dir_cos_dict


def average_mech(norm1in,norm2in):
	'''
	Computes the average mech of the solutions.

	Inputs:
		norm1in: normal to fault plane, array(3,nf)
		norm2in: slip vector, array(3,nf)
	Output:
		norm1_avg: normal to average of plane 1
		norm2_avg: normal to average of plane 2
	'''
	if norm1in.shape != norm2in.shape:
		raise ValueError('***Error in average_mech: shape of norm1in and norm2in must be the same')
	if len(norm1in.shape)!=2:
		raise ValueError('***Error in average_mech: norm1in and norm2in must each be an array of shape n-by-3')
	if norm1in.shape[0]!=3:
		raise ValueError('***Error in average_mech: norm1in and norm2in must each be an array of shape n-by-3')

	norm1=norm1in.copy()
	norm2=norm2in.copy()

	# If there is only one mechanism, return that mechanism
	if norm1.shape[1]==1:
		return norm1[:,0],norm2[:,0]

	norm1_ref=norm1in[:,0].copy()
	norm2_ref=norm2in[:,0].copy()

	rota,temp1,temp2=mech_rotation(norm1_ref,norm1[:,1:],norm2_ref,norm2[:,1:])

	norm1_avg=np.sum(np.hstack((norm1[:,[0]],temp1)),axis=1)
	norm2_avg=np.sum(np.hstack((norm2[:,[0]],temp2)),axis=1)
	ln_norm1=np.sqrt(np.sum(norm1_avg**2))
	ln_norm2=np.sqrt(np.sum(norm2_avg**2))
	norm1_avg=norm1_avg/ln_norm1
	norm2_avg=norm2_avg/ln_norm2

	# Determine the RMS observed angular difference between the average
	# Normal vectors and the normal vectors of each mechanism
	rota,temp1,temp2=mech_rotation(norm1_avg,norm1,norm2_avg,norm2)
	d11=temp1[0,:]*norm1_avg[0]+temp1[1,:]*norm1_avg[1]+temp1[2,:]*norm1_avg[2]
	d22=temp2[0,:]*norm2_avg[0]+temp2[1,:]*norm2_avg[1]+temp2[2,:]*norm2_avg[2]

	d11[d11>1]=1
	d11[d11<-1]=-1
	d22[d22>1]=1
	d22[d22<-1]=-1

	a11=np.arccos(d11)
	a22=np.arccos(d22)

	avang1=np.sqrt(np.sum(a11**2)/len(a11))
	avang2=np.sqrt(np.sum(a22**2)/len(a22))

	# the average normal vectors may not be exactly orthogonal (although
	# usually they are very close) - find the misfit from orthogonal and
	# adjust the vectors to make them orthogonal - adjust the more poorly
	# constrained plane more
	if (avang1+avang2)>=0.0001:
		maxmisf=0.01
		fract1=avang1/(avang1+avang2)
		for icount in range(100):
			dot1=norm1_avg[0]*norm2_avg[0]+norm1_avg[1]*norm2_avg[1]+norm1_avg[2]*norm2_avg[2]
			misf=90-np.rad2deg(np.arccos(dot1))
			if abs(misf)<=maxmisf:
				break
			else:
				theta1=np.deg2rad(misf*fract1)
				theta2=np.deg2rad(misf*(1-fract1))
				temp=norm1_avg
				norm1_avg=norm1_avg-norm2_avg*np.sin(theta1)
				norm2_avg=norm2_avg-temp*np.sin(theta2)
				ln_norm1=np.sqrt(np.sum(norm1_avg*norm1_avg))
				ln_norm2=np.sqrt(np.sum(norm2_avg*norm2_avg))
				norm1_avg=norm1_avg/ln_norm1
				norm2_avg=norm2_avg/ln_norm2
	return norm1_avg,norm2_avg


def mech_rotation(norm1_in,norm2_in,slip1_in,slip2_in):
	'''
	Finds the minimum rotation angle between two mechanisms.
	Does not assume that the normal and slip vectors are matched.
	Input:
		norm1_in: normal to fault plane 1
		norm2_in: normal to fault plane 2
		slip1_in: slip vector 1
		slip2_in: slip vector 2
	Output:
		rota: rotation angle
		norm2: normal to fault plane, best combination
		slip2: slip vector, best combination
	'''
	if norm1_in.shape != slip1_in.shape:
		raise ValueError('***Error in mech_rotation: shape of norm1_in and slip1_in must be the same')
	if norm2_in.shape != slip2_in.shape:
		raise ValueError('***Error in mech_rotation: shape of norm2_in and slip2_in must be the same')
	if len(norm1_in.shape)!=1:
		raise ValueError('***Error in mech_rotation: norm1_in and slip1_in must each be an array of length 3')
	if norm1_in.shape[0]!=3:
		raise ValueError('***Error in mech_rotation: norm1_in and slip1_in must each be an array of length 3')
	if len(norm2_in.shape)!=2:
		raise ValueError('***Error in mech_rotation: norm2_in and slip2_in must each be an array of shape 3-by-n')
	if norm2_in.shape[0]!=3:
		raise ValueError('***Error in mech_rotation: norm2_in and slip2_in must each be an array of shape 3-by-n')

	norm1=norm1_in.copy()
	norm2=norm2_in.copy().T
	slip1=slip1_in.copy()
	slip2=slip2_in.copy().T

	num_vect=norm2.shape[0]
	rotemp=np.zeros((num_vect,4))
	for iter_x in range(0,4): # Iteration over the 4 possibilities
		if iter_x<2:
			norm2_temp=norm2.copy()
			slip2_temp=slip2.copy()
		else:
			norm2_temp=slip2.copy()
			slip2_temp=norm2.copy()
		if (iter_x==1) | (iter_x==3):
			norm2_temp=-norm2_temp
			slip2_temp=-slip2_temp

		B1=np.cross(slip1,norm1)*-1
		B2=np.cross(slip2_temp,norm2_temp)*-1

		phi=np.zeros((num_vect,3))
		phi[:,0]=norm1[0]*norm2_temp[:,0]+norm1[1]*norm2_temp[:,1]+norm1[2]*norm2_temp[:,2]
		phi[:,1]=slip1[0]*slip2_temp[:,0]+slip1[1]*slip2_temp[:,1]+slip1[2]*slip2_temp[:,2]
		phi[:,2]=B1[0]*B2[:,0]+B1[1]*B2[:,1]+B1[2]*B2[:,2]
		phi[phi>1]=1
		phi[phi<-1]=-1
		phi=np.arccos(phi)

		phi_flag=(phi<(1e-3))
		# if the mechanisms are very close, rotation = 0. Otherwise, calculate the rotation
		rot_ind=np.where(np.any(~phi_flag,axis=1))[0]

		# if one vector is the same, it is the rotation axis
		tmp_ind=rot_ind[np.where(phi_flag[rot_ind,2])[0]]
		rotemp[tmp_ind,iter_x]=(phi[tmp_ind,0])
		tmp_ind=rot_ind[np.where(phi_flag[rot_ind,0])[0]]
		rotemp[tmp_ind,iter_x]=(phi[tmp_ind,1]) 
		tmp_ind=rot_ind[np.where(phi_flag[rot_ind,1])[0]]
		rotemp[tmp_ind,iter_x]=(phi[tmp_ind,2])

		# find difference vectors - the rotation axis must be orthogonal to all three vectors
		rot_ind=np.where(np.all(~phi_flag,axis=1))[0]

		if len(rot_ind)==0:
			continue

		n=np.zeros((len(rot_ind),3,3))
		n[:,:,0]=norm1-norm2_temp[rot_ind,:]
		n[:,:,1]=slip1-slip2_temp[rot_ind,:]
		n[:,:,2]=B1-B2[rot_ind,:]
		scale=np.sqrt(n[:,0,:]**2+n[:,1,:]**2+n[:,2,:]**2)
		n=n/scale[:,np.newaxis,:]

		qdot=np.zeros((len(rot_ind),3))
		qdot[:,2]=n[:,0,0]*n[:,0,1]+n[:,1,0]*n[:,1,1]+n[:,2,0]*n[:,2,1]
		qdot[:,1]=n[:,0,0]*n[:,0,2]+n[:,1,0]*n[:,1,2]+n[:,2,0]*n[:,2,2]
		qdot[:,0]=n[:,0,1]*n[:,0,2]+n[:,1,1]*n[:,1,2]+n[:,2,1]*n[:,2,2]

		# use the two largest difference vectors, as long as they aren't orthogonal
		iout=np.zeros(len(rot_ind),dtype=int)-1
		qdot_flag=np.any(qdot>0.9999,axis=1)
		tmp_row=np.where(qdot_flag)[0]
		if len(tmp_row)>0:
			iout[tmp_row]=np.argmax(qdot[tmp_row,:],axis=1)
		tmp_row=np.where(~qdot_flag)[0]
		if len(tmp_row)>0:
			iout[tmp_row]=np.argmin(scale[tmp_row,:],axis=1)

		n1=np.zeros((len(rot_ind),3))
		n2=np.zeros((len(rot_ind),3))
		k=np.ones(len(rot_ind),dtype=bool)
		for j in range(3):
			tmp_ind=np.where(j!=iout)[0]
			tmp_ind_1=tmp_ind[k[tmp_ind]==True]
			tmp_ind_2=tmp_ind[k[tmp_ind]==False]

			if len(tmp_ind_1)>0:
				n1[tmp_ind_1,:]=n[tmp_ind_1,:,j]
				k[tmp_ind_1]=False
			if len(tmp_ind_2)>0:
				n2[tmp_ind_2,:]=n[tmp_ind_2,:,j]

		#  find rotation axis by taking cross product
		R=np.cross(n2,n1)*-1
		scaleR=np.sqrt(np.sum(R**2,axis=1))

		if np.any(scaleR==0):
			tmp_ind=np.where(scaleR==0)[0]
			rotemp[rot_ind[tmp_ind],iter_x]=9999
			rot_ind=np.delete(rot_ind,tmp_ind)
			scaleR=np.delete(scaleR,tmp_ind)
			R=np.delete(R,tmp_ind,axis=0)

		R=R/scaleR[:,np.newaxis]
		theta=np.zeros((len(rot_ind),3))
		theta[:,0]=norm1[0]*R[:,0]+norm1[1]*R[:,1]+norm1[2]*R[:,2]
		theta[:,1]=slip1[0]*R[:,0]+slip1[1]*R[:,1]+slip1[2]*R[:,2]
		theta[:,2]=B1[0]*R[:,0]+B1[1]*R[:,1]+B1[2]*R[:,2]
		theta[theta>1]=1
		theta[theta<-1]=-1
		theta=np.arccos(theta)

		iuse=np.argmin(np.abs(theta-(np.pi/2)),axis=1)
		tmp_ind=np.arange(len(iuse))

		tmp_rotemp=(np.cos(phi[rot_ind,iuse])-np.cos(theta[tmp_ind,iuse])**2)/(np.sin(theta[tmp_ind,iuse])**2)
		tmp_rotemp[tmp_rotemp>1]=1
		tmp_rotemp[tmp_rotemp<-1]=-1
		tmp_rotemp=np.arccos(tmp_rotemp)
		rotemp[rot_ind,iter_x]=tmp_rotemp

	rotemp=np.rad2deg(rotemp)
	rotemp=np.abs(rotemp)
	irot=np.argmin(rotemp,axis=1)

	tmp_ind=np.arange(len(irot))
	rota=rotemp[tmp_ind,irot]

	tmp_ind=np.where(irot>=2)[0]
	qtemp=slip2[tmp_ind,:]
	slip2[tmp_ind,:]=norm2[tmp_ind,:]
	norm2[tmp_ind,:]=qtemp

	tmp_ind=np.where( (irot==1) | (irot==3) )[0]
	norm2[tmp_ind,:]*=-1
	slip2[tmp_ind,:]*=-1

	return rota,norm2.T,slip2.T


def mech_probability(norm1in,norm2in,cangle,prob_max,iterative_avg=False):
	'''
	Determines the probability of mechanism solutions.

	Inputs:
		norm1in: normal to fault plane
		norm2in: slip vector
		cangle: cutoff angle
		prob_max: cutoff percent for mechanism multiples
		iterative_avg:
			if True: Compute average by removing one mechanism far from the average at a time following HASH
			if False: Compute average by considering the solutions within $cangle of the average.
	Outputs:
		mech_df: DataFrame of mechanism probabilities
	'''
	if norm1in.shape != norm2in.shape:
		raise ValueError('***Error in mech_probability: shape of norm1in and norm2in must be the same')
	if len(norm1in.shape)!=2:
		raise ValueError('***Error in mech_probability: norm1in and norm2in must each be an array of shape 3-by-n')
	if norm1in.shape[0]!=3:
		raise ValueError('***Error in mech_probability: norm1in and norm2in must each be an array of shape 3-by-n')

	nf=norm1in.shape[1]
	if nf==0:
		print('*mech_probability: norm1 is empty')
		return -1

	str_avg=np.zeros(5)-999
	dip_avg=np.zeros(5)-999
	rak_avg=np.zeros(5)-999
	prob=np.zeros(5)-999
	rms_diff=np.zeros((2,5))

	if nf==1: # If there's only one mech, return that mech
		str_avg[0],dip_avg[0],rak_avg[0]=sdr_from_vector(norm1in[:,[0]],norm2in[:,[0]])
		prob[0]=1
		mech_df=pd.DataFrame(columns=['str_avg','dip_avg','rak_avg','prob','rms_diff','rms_diff_aux'],index=[0])
		mech_df['str_avg']=str_avg[0]
		mech_df['dip_avg']=dip_avg[0]
		mech_df['rak_avg']=rak_avg[0]
		mech_df['prob']=prob[0]
		mech_df['rms_diff']=rms_diff[0][0]
		mech_df['rms_diff_aux']=rms_diff[1][0]
		return mech_df

	norm1=norm1in.copy()
	norm2=norm2in.copy()
	norm_ind=np.arange(norm1.shape[1])
	rota=np.zeros(nf)
	nsltn=0

	for imult in range(5):
		unused_norm_ind=[]
		if iterative_avg:
			for icount in range(len(norm_ind)):
				norm1_avg,norm2_avg=average_mech(norm1[:,norm_ind],norm2[:,norm_ind])
				temp_rota,temp1,temp2=mech_rotation(norm1_avg,norm1[:,norm_ind],norm2_avg,norm2[:,norm_ind])

				temp_rota=np.abs(temp_rota)

				imax=np.argmax(temp_rota)
				imax_ind=norm_ind[imax]
				maxrot=temp_rota[imax]

				if maxrot<=cangle:
					break
				else:
					unused_norm_ind.append(imax_ind)
					norm_ind=np.delete(norm_ind,imax)
			prob[imult]=len(norm_ind)/nf
		else:
			# Compute the average mech of the solutions and then determine the rotation angle for each solution
			norm1_avg,norm2_avg=average_mech(norm1[:,norm_ind],norm2[:,norm_ind])
			temp_rota,temp1,temp2=mech_rotation(norm1_avg,norm1[:,norm_ind],norm2_avg,norm2[:,norm_ind])
			tmp_ind=norm_ind[np.abs(temp_rota)<=cangle]

			if len(tmp_ind)==0:
				break

			# Considering the solutions within $cangle of the average mech, recompute the average mech and find similar solutions
			norm1_avg,norm2_avg=average_mech(norm1[:,tmp_ind],norm2[:,tmp_ind])
			temp_rota,temp1,temp2=mech_rotation(norm1_avg,norm1[:,norm_ind],norm2_avg,norm2[:,norm_ind])
			tmp2_flag=(np.abs(temp_rota)<=cangle)
			norm1_avg,norm2_avg=average_mech(norm1[:,norm_ind[tmp2_flag]],norm2[:,norm_ind[tmp2_flag]])

			unused_norm_ind=norm_ind[~tmp2_flag]
			norm_ind=norm_ind[tmp2_flag]

			prob[imult]=np.sum(tmp2_flag)/nf

		if (imult>0) & (prob[imult]<prob_max):
			break

		# determine the RMS observed angular difference between the average
		# normal vectors and the normal vectors of each mechanism
		rota,temp1,temp2=mech_rotation(norm1_avg,norm1in,norm2_avg,norm2in)
		d11=temp1[0]*norm1_avg[0]+temp1[1]*norm1_avg[1]+temp1[2]*norm1_avg[2]
		d22=temp2[0]*norm2_avg[0]+temp2[1]*norm2_avg[1]+temp2[2]*norm2_avg[2]

		d11[d11>1]=1
		d11[d11<-1]=-1
		d22[d22>1]=1
		d22[d22<-1]=-1

		a11=np.arccos(d11)
		a22=np.arccos(d22)
		rms_diff[0,imult]=np.rad2deg(np.sqrt(np.sum(a11**2)/nf))
		rms_diff[1,imult]=np.rad2deg(np.sqrt(np.sum(a22**2)/nf))

		str_avg[imult],dip_avg[imult],rak_avg[imult]=sdr_from_vector(norm1_avg[np.newaxis].T,norm2_avg[np.newaxis].T)
		nsltn+=1

		if len(unused_norm_ind)>0:
			norm_ind=np.asarray(unused_norm_ind)
		else:
			break
	str_avg=str_avg[:nsltn]
	dip_avg=dip_avg[:nsltn]
	rak_avg=rak_avg[:nsltn]
	prob=prob[:nsltn]
	rms_diff=rms_diff[:,:nsltn]

	sort_ind=np.argsort(prob)[::-1]

	# Returns result as a DataFrame
	mech_df=pd.DataFrame(columns=['str_avg','dip_avg','rak_avg','prob','rms_diff','rms_diff_aux'],index=np.arange(len(sort_ind)))
	mech_df['str_avg']=str_avg[sort_ind]
	mech_df['dip_avg']=dip_avg[sort_ind]
	mech_df['rak_avg']=rak_avg[sort_ind]
	mech_df['prob']=prob[sort_ind]
	mech_df['rms_diff']=rms_diff[0][sort_ind]
	mech_df['rms_diff_aux']=rms_diff[1][sort_ind]

	return mech_df


def focal_gridsearch(sr_azimuth,takeoff,p_pol,sp_amp,dir_cos_dict,nextra,ntotal,qextra,qtotal,maxout,ncoor,min_ratio_trial_solutions=0.5,min_num_sp_solutions=10):
	'''
	Performs a grid search to find focal mechanisms using P-polarity and S/P ratio information using the python routine.
	Input:
		sr_azimuth: source-receiver azimuths, 2d array
		takeoff: takeoff angles, 2d array
		p_pol: polarity weights, 1d array
		sp_amp: polarity weights, 1d array
		dir_cos_dict: coordinate transformation dictionary, produced by dir_cos_setup()
		nextra: number of polarity additional misfits allowed above minimum
		ntotal: total number of allowed polarity misfits
		qextra: additional amplitude misfit allowed above minimum
		qtotal: total allowed amplitude misfit
		maxout: maximum number of fault planes to return
		ncoor: number of test mechanisms
		min_ratio_trial_solutions: minimum ratio of trial solutions from polarities before criteria loosened
		min_num_sp_solutions: minimum ratio of trial solutions from S/P ratios before criteria loosened
	Output:
		faultnorms_all: fault normal vectors
		faultslips_all: fault slip vectors
	'''
	ntab=180
	astep=1/ntab

	'''
	Transforms the spherical coordinates to cartesian
	'''
	pol_ind=np.where(p_pol!=0)[0]

	takeoff_r=np.deg2rad(takeoff[pol_ind])
	sr_azimuth_r=np.deg2rad(sr_azimuth[pol_ind])

	xyz=np.zeros((3,takeoff_r.shape[0],takeoff_r.shape[1]))
	xyz[0,:]=np.sin(takeoff_r)*np.cos(sr_azimuth_r)
	xyz[1,:]=np.sin(takeoff_r)*np.sin(sr_azimuth_r)
	xyz[2,:]=-np.cos(takeoff_r)

	p_b1=np.tensordot(xyz,dir_cos_dict['b1'],axes=[[0],[0]])
	p_b3=np.tensordot(xyz,dir_cos_dict['b3'],axes=[[0],[0]])

	# # Slower (but probably more comprehensible) implementation of the above
	# p_b1=b1[0,:]*xyz[0,:,:,np.newaxis]+\
	#         b1[1,:]*xyz[1,:,:,np.newaxis]+\
	#         b1[2,:]*xyz[2,:,:,np.newaxis]
	# p_b3=b3[0,:]*xyz[0,:,:,np.newaxis]+\
	#         b3[1,:]*xyz[1,:,:,np.newaxis]+\
	#         b3[2,:]*xyz[2,:,:,np.newaxis]

	prod=((p_b1<0)!=(p_b3<0)) # If True, predicted sign is negative. If False, predicted sign is positive

	qmiss=(prod != (p_pol[pol_ind]<0)[:,np.newaxis,np.newaxis])*np.abs(p_pol[pol_ind])[:,np.newaxis,np.newaxis]

	fit=np.sum(qmiss,axis=0)

	# Calculates max misfit for each trial
	qmissmax=fit.min(axis=1)+nextra
	qmissmax[qmissmax<ntotal]=ntotal

	sp_finite_ind=np.where(np.isfinite(sp_amp))[0]
	qacount=len(sp_finite_ind)
	if qacount>0:
		takeoff_r=np.deg2rad(takeoff[sp_finite_ind])
		sr_azimuth_r=np.deg2rad(sr_azimuth[sp_finite_ind])

		xyz=np.zeros((3,takeoff_r.shape[0],takeoff_r.shape[1]))
		xyz[0,:]=np.sin(takeoff_r)*np.cos(sr_azimuth_r)
		xyz[1,:]=np.sin(takeoff_r)*np.sin(sr_azimuth_r)
		xyz[2,:]=-np.cos(takeoff_r)

		p_b3=np.tensordot(xyz,dir_cos_dict['b3'],axes=[[0],[0]])

		p_proj1=xyz[0,:,:,np.newaxis]-(p_b3*dir_cos_dict['b3'][0,:])
		p_proj2=xyz[1,:,:,np.newaxis]-(p_b3*dir_cos_dict['b3'][1,:])
		p_proj3=xyz[2,:,:,np.newaxis]-(p_b3*dir_cos_dict['b3'][2,:])

		plen=np.sqrt(p_proj1**2+p_proj2**2+p_proj3**2)
		p_proj1=p_proj1/plen
		p_proj2=p_proj2/plen
		p_proj3=p_proj3/plen

		pp_b1=dir_cos_dict['b1'][0,:]*p_proj1+dir_cos_dict['b1'][1,:]*p_proj2+dir_cos_dict['b1'][2,:]*p_proj3
		pp_b2=dir_cos_dict['b2'][0,:]*p_proj1+dir_cos_dict['b2'][1,:]*p_proj2+dir_cos_dict['b2'][2,:]*p_proj3
		i=np.round((p_b3+1.)/astep).astype(int)

		theta=dir_cos_dict['thetable'][i]
		i=np.round((pp_b2+1.)/astep).astype(int)
		j=np.round((pp_b1+1.)/astep).astype(int)
		phi=dir_cos_dict['phitable'][i,j]

		i=np.round(phi/(np.pi*astep)).astype(int)
		i[i>(2*ntab-1)]=0
		j=np.round(theta/(np.pi*astep)).astype(int)
		j[j>(ntab-1)]=0

		p_amp=dir_cos_dict['amptable'][0,j,i]
		s_amp=dir_cos_dict['amptable'][1,j,i]

		sp_rat=np.zeros(p_amp.shape)
		sp_rat[p_amp==0]=4.0
		sp_rat[s_amp==0]=-2.0
		nonzero_flag=((p_amp!=0) & (s_amp!=0))
		sp_rat[nonzero_flag]=np.log10(4.9*s_amp[nonzero_flag]/p_amp[nonzero_flag])

		qamiss=np.abs(sp_amp[sp_finite_ind][:,np.newaxis,np.newaxis]-sp_rat)
		afit=np.sum(qamiss,axis=0)

		# Calculates max misfit for each trial
		qamissmax=afit.min(axis=1)+qextra
		qamissmax[qamissmax<qtotal]=qtotal

		goodmech_flag=( (fit<=qmissmax[:,np.newaxis]) & (afit<=qamissmax[:,np.newaxis]) )

		# Ratio of the number of trials with a solution to trials with no solution
		ratio_trials_with_solutions=np.sum(np.any(goodmech_flag,axis=1))/goodmech_flag.shape[0]
		# Number of unique mech solutions across the trials
		good_fp_ind=np.where(np.any(goodmech_flag,axis=0))[0]

		# If there are no solutions that meet the criteria, loosen the amplitude criteria
		if (len(good_fp_ind)<min_num_sp_solutions) | (ratio_trials_with_solutions<min_ratio_trial_solutions):
			masked_afit=np.ma.masked_array(afit,mask=(fit>qmissmax[:,np.newaxis]))
			qmis0min=masked_afit.min(axis=1).data
			qamissmax=qmis0min+nextra
			qamissmax[qamissmax<qtotal]=qtotal
			goodmech_flag=( (fit<=qmissmax[:,np.newaxis]) & (afit<=qamissmax[:,np.newaxis]) )
			good_fp_ind=np.where(np.any(goodmech_flag,axis=0))[0]

	else:
		goodmech_flag=( (fit<=qmissmax[:,np.newaxis]) )
		good_fp_ind=np.where(np.any(goodmech_flag,axis=0))[0]

	if len(good_fp_ind)>maxout: # If more than maxout solutions meet criteria, randomly select maxout solutions
		good_fp_ind=rng.choice(good_fp_ind,maxout,replace=False)

	faultnorms_all=np.vstack((dir_cos_dict['b3'][0,good_fp_ind],dir_cos_dict['b3'][1,good_fp_ind],dir_cos_dict['b3'][2,good_fp_ind]))
	faultslips_all=np.vstack((dir_cos_dict['b1'][0,good_fp_ind],dir_cos_dict['b1'][1,good_fp_ind],dir_cos_dict['b1'][2,good_fp_ind]))

	return faultnorms_all,faultslips_all


def determine_max_gap(in_azimuth_deg,in_takeoff_deg):
	'''
	Given an array of azimuths and takeoff angles, calculates the maximum
	azimuthal and takeoff gaps.
	Input:
		in_azimuth_deg: source-receiver azimuths
		in_takeoff_deg: takeoff angles
	Output:
		max_azimuthal_gap: maximum azimuthal gap between measurements
		max_takeoff_gap: maximum takeoff angle gap between measurements
	'''
	azimuth_deg=in_azimuth_deg.copy()
	takeoff_deg=in_takeoff_deg.copy()
	azimuth_deg[takeoff_deg>90]-=180
	takeoff_deg[takeoff_deg>90]=180-takeoff_deg[takeoff_deg>90]
	azimuth_deg[azimuth_deg<0]+=360

	azimuth_deg=np.sort(azimuth_deg)
	takeoff_deg=np.sort(takeoff_deg)

	max_azimuthal_gap=np.max([np.max(np.diff(azimuth_deg)),azimuth_deg[0]+360-azimuth_deg[-1]])
	max_takeoff_gap=np.max([np.max(np.diff(takeoff_deg)),takeoff_deg[0],90-takeoff_deg[-1]])

	return max_azimuthal_gap,max_takeoff_gap


def mech_misfit(mech_df,p_azi_mc,p_the_mc,p_pol,sp_amp):
	'''
	Loops over the potential mech solution dataframe and calcultes the polarity
	and S/P misfits. The polarity and S/P agreements for the first mech solutions
	are returned.
	Input:
		mech_df: mechanism solution dataframe, produced by mech_probability()
		p_azi_mc: source-receiver azimuth for the measurement, 1d array
		p_the_mc: takeoff angles for the measurement, 1d array
		p_pol: polarity weights, 1d array
		sp_amp: polarity weights, 1d array
	Output:
		mech_df: mechanism solution dataframe, with additional columns
		pol_agreement_out: Flag for if the polarity measurement (dis)agrees with the solution
		sp_diff_out: Flag for if the S/P measurement (dis)agrees with the solution
	'''
	mech_df['mfrac']=np.nan
	mech_df['mavg']=np.nan
	mech_df['stdr']=np.nan

	for imult in range(len(mech_df)):
		mech_df.loc[imult,'mfrac'],mech_df.loc[imult,'mavg'],mech_df.loc[imult,'stdr'],tmp_pol_agreement_out,tmp_sp_diff_out=calculate_misfit(p_azi_mc,p_the_mc,p_pol,sp_amp,mech_df.loc[imult,'str_avg'],mech_df.loc[imult,'dip_avg'],mech_df.loc[imult,'rak_avg'])
		if imult==0:
			pol_agreement_out=tmp_pol_agreement_out
			sp_diff_out=tmp_sp_diff_out
	return mech_df,pol_agreement_out,sp_diff_out


def mech_quality(mech_df,qual_criteria_dict):
	'''
	Given the mech solutions and the quality criteria, adds the mech qualities
	to the mech solution dataframe.
	Input:
		mech_df: mechanism solution dataframe, produced by mech_probability()
		qual_criteria_dict: dictionary of quality criteria, created in SKHASH.py
	Output:
		mech_df: mechanism solution dataframe, with additional quality column
	'''
	qual_ind=np.zeros((len(mech_df),4),dtype=int)-999
	qual_ind[:,0]=np.argmax(mech_df.prob.values[:,np.newaxis]>=qual_criteria_dict['probs'],axis=1)
	qual_ind[:,1]=np.argmax(mech_df[['rms_diff','rms_diff_aux']].mean(axis=1).values[:,np.newaxis]<=qual_criteria_dict['var_avg'],axis=1)
	qual_ind[:,2]=np.argmax(mech_df.mfrac.values[:,np.newaxis]<=qual_criteria_dict['mfrac'],axis=1)
	qual_ind[:,3]=np.argmax(mech_df.stdr.values[:,np.newaxis]>=qual_criteria_dict['stdr'],axis=1)
	mech_df['qual']=qual_criteria_dict['qual_letter'][np.max(qual_ind,axis=1)]

	return mech_df


def calculate_misfit(p_azi_mc,p_the_mc,p_pol,sp_amp,str_avg,dip_avg,rak_avg):
	'''
	Calculates the polarity misfit percent and S/P misfit for a given mechanism solution.

	Inputs:
		p_azi_mc: azimuths
		p_the_mc: takeoff angles
		p_pol: polarity measurements
		sp_amp: S/P ratios
		str_avg,dip_avg,rak_avg: the given mechanism solution

	Outputs:
		mfrac: weighted fraction misfit polarities
		mavg: average S/P misfit (log10)
		stdr: station distribution ratio
		pol_agreement_out: boolean of whether the polarity agrees with the solution
		sp_diff_out: difference between the measured and expected S/P ratio given the solution
	'''

	strike=np.deg2rad(str_avg)
	dip=np.deg2rad(dip_avg)
	rake=np.deg2rad(rak_avg)

	M=np.zeros((3,3))
	M[0,0]=-np.sin(dip)*np.cos(rake)*np.sin(2*strike)-np.sin(2*dip)*np.sin(rake)*np.sin(strike)*np.sin(strike)
	M[1,1]=np.sin(dip)*np.cos(rake)*np.sin(2*strike)-np.sin(2*dip)*np.sin(rake)*np.cos(strike)*np.cos(strike)
	M[2,2]=np.sin(2*dip)*np.sin(rake)
	M[0,1]=np.sin(dip)*np.cos(rake)*np.cos(2*strike)+0.5*np.sin(2*dip)*np.sin(rake)*np.sin(2*strike)
	M[1,0]=M[0,1]
	M[0,2]=-np.cos(dip)*np.cos(rake)*np.cos(strike)-np.cos(2*dip)*np.sin(rake)*np.sin(strike)
	M[2,0]=M[0,2]
	M[1,2]=-np.cos(dip)*np.cos(rake)*np.sin(strike)+np.cos(2*dip)*np.sin(rake)*np.cos(strike)
	M[2,1]=M[1,2]

	bb3,bb1=vector_from_sdr(strike,dip,rake)
	bb2=np.cross(bb1,bb3)*-1

	mfrac=0.
	qcount=0.
	stdr=0.
	scount=0.
	mavg=0.
	acount=0.

	p_a1,p_a2,p_a3=cartesian_transform(p_the_mc,p_azi_mc,1)
	p_b1=bb1[0]*p_a1+bb1[1]*p_a2+bb1[2]*p_a3
	p_b3=bb3[0]*p_a1+bb3[1]*p_a2+bb3[2]*p_a3
	p_proj1=p_a1-p_b3*bb3[0]
	p_proj2=p_a2-p_b3*bb3[1]
	p_proj3=p_a3-p_b3*bb3[2]
	plen=np.sqrt(p_proj1*p_proj1+p_proj2*p_proj2+p_proj3*p_proj3)
	p_proj1=p_proj1/plen
	p_proj2=p_proj2/plen
	p_proj3=p_proj3/plen
	pp_b1=bb1[0]*p_proj1+bb1[1]*p_proj2+bb1[2]*p_proj3
	pp_b2=bb2[0]*p_proj1+bb2[1]*p_proj2+bb2[2]*p_proj3
	phi=np.arctan2(pp_b2,pp_b1)
	theta=np.arccos(p_b3)
	p_amp=np.abs(np.sin(2*theta)*np.cos(phi))

	pol_ind=np.where(p_pol!=0)[0]
	wt=np.sqrt(p_amp)

	scount=np.sum(np.abs(p_pol[pol_ind]))
	a=np.zeros((3,len(pol_ind)))
	b=np.zeros((3,len(pol_ind)))

	azi=np.deg2rad(p_azi_mc[pol_ind])
	toff=np.deg2rad(p_the_mc[pol_ind])
	a[0,:]=np.sin(toff)*np.cos(azi)
	a[1,:]=np.sin(toff)*np.sin(azi)
	a[2,:]=-np.cos(toff)
	b=np.sum(M[:,:,np.newaxis]*a,axis=1)

	neg_pol_v=(a[0,:]*b[0,:]+a[1,:]*b[1,:]+a[2,:]*b[2,:])*p_pol[pol_ind]
	mfrac=np.sum(wt[pol_ind][neg_pol_v<0])

	pol_agreement_out=np.empty(len(p_pol))
	pol_agreement_out[:]=np.nan
	pol_agreement_out[pol_ind[neg_pol_v>=0]]=True
	pol_agreement_out[pol_ind[neg_pol_v<0]]=False

	qcount=np.sum(wt[pol_ind])
	stdr=np.sum(wt[pol_ind])

	sp_ind=np.where(~np.isnan(sp_amp))[0]
	acount=len(sp_ind)
	sp_diff_out=np.empty(len(sp_amp),dtype=float)
	sp_diff_out[:]=np.nan
	if acount>0:
		s1=np.cos(2*theta[sp_ind])*np.cos(phi[sp_ind])
		s2=-np.cos(theta[sp_ind])*np.sin(phi[sp_ind])
		s_amp=np.sqrt(s1*s1+s2*s2)
		sp_rat=np.log10(4.9*s_amp/p_amp[sp_ind])

		#Difference between S/P ratio and expected given the solution
		sp_diff=sp_amp[sp_ind]-sp_rat
		sp_diff_out[sp_ind]=sp_diff

		mavg=np.sum(sp_diff)
		stdr+=np.sum(wt[sp_ind])
		scount+=acount

	if qcount==0:
		mfrac=0
	else:
		mfrac=mfrac/qcount

	if acount==0:
		mavg=0
	else:
		mavg=mavg/acount

	if scount==0:
		stdr=0
	else:
		stdr=stdr/scount

	return mfrac,mavg,stdr,pol_agreement_out,sp_diff_out


def create_lookup_table(p_dict):
	'''
	Reads in the velocity model inputs, does minor QCing, and creates lookup table (if necessary).
	If the lookup table already has been computed, it can read it from the disk.

	Inputs:
		p_dict['vmodel_paths']: a list of filepaths to the velocity models. Each file should be a whitespace delimited file following the format:
			depth(km) vp(km/s)
		p_dict['look_dep']: list of: minimum source depth, maximum, and interval for the lookup table
		p_dict['look_del']: list of: minimum source-station distance, maximum, and interval for the lookup table
		p_dict['nump']: number of rays traced
		p_dict['nx0']: maximum source-station distance (km) bins for look-up tables
		p_dict['nd0']: maximum source depth (km) bins for look-up tables
		p_dict['output_angle_precision']: Number of decimal places to output for take off angles
		p_dict['recompute_lookup_table']: a boolean. If True, the lookup tables will be recomputed, even if a matching lookup table exists.
		p_dict['write_lookup_table']: a boolean. If True, the lookup table will be written to the disk. A '.npz' file suffix will be added to the path.
	Output:
		deptab: array of depths (km)
		delttab: array of distances (km)
		table: Lookup table (array)
	'''

	deptab=np.arange(p_dict['look_dep'][0],p_dict['look_dep'][1]+p_dict['look_dep'][2],p_dict['look_dep'][2]) # array of depths (km)
	delttab=np.arange(p_dict['look_del'][0],p_dict['look_del'][1]+p_dict['look_del'][2],p_dict['look_del'][2]) # array of distances (km)

	if (p_dict['recompute_lookup_table']==False) or (p_dict['write_lookup_table']):
		lookup_table_params=np.asarray([p_dict['nump'],p_dict['nx0'],p_dict['nd0'],p_dict['output_angle_precision']])

	table_list=[]
	for vmodel_ind,vmodel_path in enumerate(p_dict['vmodel_paths']):
		lookup_vmodel_path=vmodel_path+'.lookup.npz'
		vmodel_depthvp=pd.read_csv(vmodel_path,names=['depth','vp_km_s'],sep=',',comment='#').values
		if len(vmodel_depthvp)==1:
			raise ValueError('***Velocity model ({}) has only a single velocity. There must be at least two points.***'.format(vmodel_depthvp))
		else:
			if not(all(np.diff(vmodel_depthvp[:,0])>0)):
				raise ValueError('***Velocity model ({}) is expected to be ordered in terms of increasing depth.'.format(vmodel_depthvp))
		if vmodel_depthvp[:,0].max()>6400:
			print('*WARNING: The velocity model ({}) has depths as large as {} km. Is this intentional?'.format(vmodel_path,vmodel_depthvp[:,0].max()))

		# If there are constant velocity layers, merge those rows
		drop_constant_vel_ind=np.where(np.diff(vmodel_depthvp[:,1])<=0)[0]+1
		if len(drop_constant_vel_ind)>0:
			vmodel_depthvp=np.delete(vmodel_depthvp,drop_constant_vel_ind,axis=0)

		compute_lookup_table=True
		if os.path.exists(lookup_vmodel_path) & (p_dict['recompute_lookup_table']==False):
			print('Loading precomputed lookup table ({}/{}): {}'.format(vmodel_ind,len(p_dict['vmodel_paths'])-1,lookup_vmodel_path))
			try:
				tmp_lookup=np.load(lookup_vmodel_path)
			except:
				print('\tError loading lookup table. Recomputing.')

			# Ensures the lookup table .npz file contains the expected variable names
			if tmp_lookup.files==['table', 'lookup_table_params', deptab, delttab]:
				# Ensures that the parameters used to create the lookup table are identical to the current parameters.
				if np.array_equal(tmp_lookup['lookup_table_params'],lookup_table_params) &\
					np.array_equal(tmp_lookup[deptab],deptab) &\
					np.array_equal(tmp_lookup[delttab],delttab):
					table=tmp_lookup['table']
					compute_lookup_table=False
					print('\tLoad successful.')
				else:
					print('\tThe parameters used to create the lookup table differ from your current set parameters. Recomputing lookup table.')
			else:
				print('\tIssue with interpreting previously saved lookup table. Recomputing lookup table.')

		if compute_lookup_table:
			print('Creating lookup table ({}/{}): {}'.format(vmodel_ind,len(p_dict['vmodel_paths'])-1,vmodel_path))
			table=create_takeoff_table(vmodel_depthvp,deptab,delttab,p_dict['nump'],p_dict['nx0'],p_dict['nd0'],p_dict['output_angle_precision'])
			print('\tCreated table.')

			if p_dict['write_lookup_table']:
				np.savez(lookup_vmodel_path,table=table,lookup_table_params=lookup_table_params,deptab=deptab,delttab=delttab)
				print('\tSaved lookup table: {}'.format(lookup_vmodel_path))
		table_list.append(table)
	table=np.dstack(table_list)

	return {'deptab':deptab,'delttab':delttab,'table':table}
