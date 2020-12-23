% Purpose: Convert TRDI Sentinel V velocities from beam to Earth coordinates
% and apply some QC steps.
clear all; close all; clc

addpath('/home/andre/code/ADCPtools/', '/home/andre/code/ADCPtools/qc', '/home/andre/code/ADCPtools/utils');

tl = '2017-01-24 00:00:00';
tr = '2017-12-24 06:00:00';
theta = 25;
threshsd = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select mooring and deployment number %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mooring = 'OC40S-A';
deployment = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deployment = num2str(deployment);
disp(['Mooring: ' mooring ', deployment #', deployment]);
disp('===============================');
disp(' ');

fin = ['/home/andre/phd/data/innershelfdri_moorings/' mooring '/deployment' deployment '/' mooring(1:end-2) '-Ad' deployment 'beam.mat'];

load(fin);
tl = datenum(tl);
tr = datenum(tr);
ft = A.t >= tl & A.t < tr;
hdng = A.hdng(ft);
ptch = A.ptch(ft);
roll = A.roll(ft);
b1 = A.b1(:,ft);
b2 = A.b2(:,ft);
b3 = A.b3(:,ft);
b4 = A.b4(:,ft);
b5 = A.b5(:,ft);
r = A.r;
z = A.z;
t = A.t(ft);
p = A.p;
clear A

Gimbaled = true;
uvwBeam5 = true;
Use3beamSol = true;

[unbmp vnbmp wnbmp w5nbmp] = janus5beam2earth(hdng, ptch, roll, theta, b1, b2, b3, b4, b5, r, Gimbaled, 'none', uvwBeam5, Use3beamSol);

[u v w w5] = janus5beam2earth(hdng, ptch, roll, theta, b1, b2, b3, b4, b5, r, Gimbaled, 'linear', uvwBeam5, Use3beamSol);

clear b1 b2 b3 b4 b5

warning off
disp('Despiking u.');
u = despike_filt(u, threshsd);
u = despike_filt(u', threshsd)';
disp('Despiking v.');
v = despike_filt(v, threshsd);
v = despike_filt(v', threshsd)';
disp('Despiking w.');
w = despike_filt(w, threshsd);
w = despike_filt(w', threshsd)';
disp('Despiking w5.');
w5 = despike_filt(w5, threshsd);
w5 = despike_filt(w5', threshsd)';

A = struct();
A.t = t;
A.r = r;
A.z = z;
A.hdng = hdng;
A.ptch = ptch;
A.roll = roll;
A.p = p;
A.u = u;
A.v = v;
A.w = w;
A.w5 = w5;

fout = ['/home/andre/phd/data/innershelfdri_moorings/' mooring '/deployment' deployment '/' mooring(1:end-2) '-Ad' deployment 'enu'];

save(fout, 'A', '-v7.3');

% Time-averaged ENU velocities.
Ntavg = 30;

favg = @(A) nanmean(A);
A.u = blkproc(A.u, [1 Ntavg], favg);
A.v = blkproc(A.v, [1 Ntavg], favg);
A.w = blkproc(A.w, [1 Ntavg], favg);
A.w5 = blkproc(A.w5, [1 Ntavg], favg);
A.p = blkproc(A.p, [1 Ntavg], favg);
A.hdng = blkproc(A.hdng, [1 Ntavg], favg);
A.ptch = blkproc(A.ptch, [1 Ntavg], favg);
A.roll = blkproc(A.roll, [1 Ntavg], favg);
A.t = blkproc(A.t, [1 Ntavg], favg);

fout = ['/home/andre/phd/data/innershelfdri_moorings/' mooring '/deployment' deployment '/' mooring(1:end-2) '-Ad' deployment 'enu30s'];

save(fout, 'A', '-v7.3');

% % u
% figure;
% p1 = subplot(311);
% pcolor(t, r, unbmp); shading flat; caxis([-0.5 0.5]); colorbar; datetick('x');
% p2 = subplot(312);
% pcolor(t, r, u); shading flat; caxis([-0.5 0.5]); colorbar; datetick('x');
% p3 = subplot(313);
% pcolor(t, r, u-unbmp); shading flat; caxis([-0.05 0.05]); colorbar; datetick('x');
%
% linkaxes([p1 p2 p3], 'x'); linkaxes([p1 p2 p3], 'y');
% axis tight;
% colormap(bluewhitered);
%
% % v
% figure;
% colormap(bluewhitered);
% p1 = subplot(311);
% pcolor(t, r, vnbmp); shading flat; caxis([-0.5 0.5]); colorbar; datetick('x')
% p2 = subplot(312);
% pcolor(t, r, v); shading flat; caxis([-0.5 0.5]); colorbar; datetick('x')
% p3 = subplot(313);
% pcolor(t, r, v-vnbmp); shading flat; caxis([-0.05 0.05]); colorbar; datetick('x')
%
% linkaxes([p1 p2 p3], 'x'); linkaxes([p1 p2 p3], 'y');
% axis tight;
% colormap(bluewhitered);

% u
figure(13);
plot(nanmean(unbmp, 2),z, 'bo'); hold on;
plot(nanmean(u, 2),z, 'ro');

% v
figure(14);
plot(nanmean(vnbmp, 2),z, 'bo'); hold on;
plot(nanmean(v, 2),z, 'ro');

% w
figure(15);
plot(nanmean(wnbmp, 2),z, 'bo'); hold on;
plot(nanmean(w, 2),z, 'ro');

% w5
figure(16);
plot(nanmean(w5nbmp, 2),z, 'bo'); hold on;
plot(nanmean(w5, 2),z, 'ro');
