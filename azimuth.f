	call azimth(0,slat,slon,rlat,rlon,delta,azim,bazim)
c
c
	subroutine azimth(ellips,slat,slon,rlat,rlon,delta,azim,bazim)

c   This routine uses Euler angles to find the geocentric distance,
c   azimuth, and back azimuth for a source-reciever pair.
c
c   Input
c
c     slat  - source geographic latitude in decimal degrees
c     slon  - source longitude in decimal degrees
c     rlat  - reciever geographic latitude in decimal degrees
c     rlon  - reciever longitude in decimal degrees
c
c   Output
c
c     delta - geocentric source-reciever distance in decimal degrees of arc
c     azim  - geocentric azimuth from the source to the reciever
c     bazim - geocentric back azimuth from the reciever to the source
c
c   The distance calculated here delta is always between 0 and 180 degrees. 
c   Accordingly, the azimuth and back azimuth are defined for the minor 
c   arc between (slat,slon) and (rlat,rlon).

        implicit real*8 (a-h,o-z)
	integer ellips

        data flt/298.25d0/

	pi=4.d0*datan(1.d0)
	dtor=pi/180.d0
	
	if(ellips.ne.0) then
          e=1.d0/flt
	else
	  e=0.d0
	endif
	
c   Convert to geocentric coordinates and from latitude to colatitude.

        slatra=dtor*slat
        w=dsin(slatra)
        s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(slatra)
        scolat=pi/2.d0-slatra+s
        rlatra=dtor*rlat
        w=dsin(rlatra)
        s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(rlatra)
        rcolat=pi/2.d0-rlatra+s

        slonra=slon*dtor
        rlonra=rlon*dtor
        c2=dcos(scolat)
        s2=dsin(scolat)
        c1=dcos(slonra)
        s1=dsin(slonra)
        slatrc=dsin(rcolat)

c   Find the azimuth and distance by rotating the source to the north pole.

        x0=slatrc*dcos(rlonra)
        y0=slatrc*dsin(rlonra)
        z0=dcos(rcolat)
        x1=c1*x0+s1*y0

        z0=dcos(rcolat)
        x1=c1*x0+s1*y0
        y1=-s1*x0+c1*y0
        z1=z0
        x2=c2*x1-s2*z1
        y2=y1
        z2=c2*z1+s2*x1
        call angles(x2,y2,z2,delta,azim)
        azim=180.d0-azim

c   Find the back azimuth by rotating the receiver to the north pole.

        c2=dcos(rcolat)
        s2=dsin(rcolat)
        c1=dcos(rlonra)
        s1=dsin(rlonra)
        slatrc=dsin(scolat)
        x0=slatrc*dcos(slonra)
        y0=slatrc*dsin(slonra)
        z0=dcos(scolat)
        x1=c1*x0+s1*y0
        y1=-s1*x0+c1*y0
        z1=z0
        x2=c2*x1-s2*z1
        y2=y1
        z2=c2*z1+s2*x1
        call angles(x2,y2,z2,delta,bazim)
        bazim=180.d0-bazim

        return
        end
