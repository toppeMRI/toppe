function [gxr,gyr,t,npts] = makevdspiral2(FOV,npix,shots,dsamp,gmax,smax)

%{
FOV:cm
npix: num of pix
shots: num of shots for spiral
gmax: max grad (G/cm)
smax: max slew (mT/m/ms) 150 should be good enough
%}

[gvd,kvd,t,slvd,npts] = genspivd2(FOV, npix, shots, gmax, smax, dsamp);
gx=real(gvd);
gy=imag(gvd);

rotangle = 360/shots;%*180/length(zlist);
rotlist = 0:rotangle:(shots-1)*rotangle;
for n =1:shots
    ga=deg2rad(rotlist(n));
    gxr(:,n) = cos(ga) * gx - sin(ga) * gy;
    gyr(:,n) = sin(ga) * gx + cos(ga) * gy;
end

