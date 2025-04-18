clc;
clear;
close all;

N = 128;
nModes = 200;
L = 9.0/100 *2.0*pi;

km0 = 2.0*pi/L;
  
lx = L; ly = L; lz = L;
nx = N; ny = N; nz = N;
dx = lx/nx; dy = ly/ny; dz = lz/nz;
hdx = dx/2.0; hdy = dy/2.0; hdz = dz/2.0;

phi = zeros(nModes+1, 1);
nu = zeros(nModes+1, 1);
theta = zeros(nModes+1, 1);
psi = zeros(nModes+1, 1);

for ii = 1:nModes + 1
    phi(ii) = 2.0*pi*rand();
    nu(ii) = rand();
    theta(ii) = acos(2*nu(ii)-1.0);
    psi(ii) = pi*rand() - pi/2.0;
end

kmmax = pi / dx;
dk = (kmmax - km0) / nModes;
km = linspace(km0, kmmax, nModes+1);

kx = zeros(nModes+1, 1);
ky = zeros(nModes+1, 1);
kz = zeros(nModes+1, 1);

for ii = 1 : nModes+1
    kx(ii) = sin(theta(ii))*cos(phi(ii))*km(ii);
    ky(ii) = sin(theta(ii))*sin(phi(ii))*km(ii);
    kz(ii) = cos(theta(ii))*km(ii);
end

ktx = zeros(nModes+1, 1);
kty = zeros(nModes+1, 1);
ktz = zeros(nModes+1, 1);

for ii = 1 : nModes+1
    ktx(ii) = sin(kx(ii)*hdx)/dx;
    kty(ii) = sin(ky(ii)*hdy)/dy;
    ktz(ii) = sin(kz(ii)*hdz)/dz;
end

sxm = zeros(nModes+1, 1);
sym = zeros(nModes+1, 1);
szm = zeros(nModes+1, 1);

for ii = 1:nModes + 1
    phi(ii) = 2.0*pi*rand();
    nu(ii) = rand();
    theta(ii) = acos(2*nu(ii)-1.0);
end

for ii = 1:nModes + 1
    zetax = sin(theta(ii)) * cos(phi(ii));
    zetay = sin(theta(ii)) * sin(phi(ii));
    zetaz = cos(theta(ii));

    sxm(ii) = (zetay * ktz(ii) - zetaz * kty(ii));
    sym(ii) = (zetax * ktz(ii) - zetaz * ktx(ii));
    szm(ii) = (zetax * kty(ii) - zetay * ktx(ii));

    smag = 1.0/sqrt(sxm(ii)*sxm(ii) + sym(ii)*sym(ii) + szm(ii)*szm(ii));
    sxm(ii) = sxm(ii) * smag;
    sym(ii) = sym(ii) * smag;
    szm(ii) = szm(ii) * smag;
end

for ii =1:nModes+1
    espec = 2.0 * sqrt(karman_spec(km(ii))*dk);
    sxm(ii) = sxm(ii) * espec;
    sym(ii) = sym(ii) * espec;
    szm(ii) = szm(ii) * espec;
end

%% visulization
x_1d = linspace(0, lx, nx);
y_1d = linspace(0, ly, ny);
z_1d = linspace(0, lz, nz);

uu = zeros(nx, ny, nz);
vv = zeros(nx, ny, nz);
ww = zeros(nx, ny, nz);

for ii = 1 : nx
for jj = 1 : ny
for kk = 1 : nz
    for mm = 1 : nModes+1
        arg = kx(mm) * x_1d(ii) + ky(mm) * y_1d(jj) + kz(mm) * z_1d(kk) - psi(mm);
        uu(ii, jj, kk) = uu(ii, jj, kk) + cos(arg)*sxm(mm);
        vv(ii, jj, kk) = vv(ii, jj, kk) + cos(arg)*sym(mm);
        ww(ii, jj, kk) = ww(ii, jj, kk) + cos(arg)*szm(mm);
    end
end
end
end

[knyquist, wave_numbers, tke_spectrum] = cal_tke_spectrum(uu,vv,ww,lx,ly,lz);
loglog(wave_numbers, tke_spectrum);
grid on;

karman_spectrum = zeros(length(wave_numbers), 1);
for ii = 1:length(wave_numbers)
    karman_spectrum(ii) = karman_spec(wave_numbers(ii));
end

hold on;
loglog(wave_numbers, karman_spectrum, 'ro');

%% function
function res = pow(aa, bb)
    res = aa ^ bb;
end

function tke = karman_spec( kk )
    ke_ = 40.0;
    nu = 1.0e-5;
    urms = 0.25;
    
    ke = sqrt(5.0/12)*ke_;
    L = 0.746834/ke;
    alpha = 1.452762113;
    epsilone = urms*urms*urms/L;
    keta = pow(epsilone, 0.25)*pow(nu,-3.0/4.0);

    r1 = kk/ke;
    r2 = kk/keta;
    tke = alpha*(urms*urms/ke)*(r1*r1*r1*r1 / pow(1.0 + r1*r1, 17.0/6.0))*exp(-2.0*r2*r2);

end

function [knyquist, wave_numbers, tke_spectrum] = cal_tke_spectrum(uu, vv, ww, lx, ly, lz)
    dimension = size(uu);

    nx = dimension(1);
    ny = dimension(2);
    nz = dimension(3);

    nt = nx * ny * nz;
    n = nx;

    uh = fftn(uu)/nt;
    vh = fftn(vv)/nt;
    wh = fftn(ww)/nt;

    tkeh = real(0.5*(uh.*conj(uh) + vh.*conj(vh) + wh.*conj(wh)));

    k0x = 2.0*pi/lx;
    k0y = 2.0*pi/ly;
    k0z = 2.0*pi/lz;

    knorm = (k0x + k0y + k0z)/3.0;

    kxmax = nx/2.0;
    kymax = ny/2.0;
    kzmax = nz/2.0;

    wave_numbers = knorm*linspace(0,n-1,n);

    tke_spectrum = zeros(n, 1);

    for kx = 0:nx-1
        rkx = kx;
        if(kx > kxmax)
            rkx = rkx - nx;
        end
        for ky = 0:ny-1
            rky = ky;
            if(ky > kymax)
                rky = rky - ny;
            end
            for kz = 0:nz-1
                rkz = kz;
                if(kz > kzmax)
                    rkz = rkz - nz;
                end
                rk = sqrt(rkx*rkx+rky*rky+rkz*rkz);
                k = round(rk);
                tke_spectrum(k+1) = tke_spectrum(k+1) + tkeh(kx+1,ky+1,kz+1);
            end
        end
    end

    tke_spectrum = tke_spectrum/knorm;

    knyquist = knorm * min([nx,ny,nz])/2;
end
