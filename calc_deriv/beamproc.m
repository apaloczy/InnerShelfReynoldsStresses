% Purpose: Apply some QC steps on raw along-beam TRDI Sentinel V velocities.
clear all; close all; clc

addpath('/home/andre/code/ADCPtools/', '/home/andre/code/ADCPtools/qc', '/home/andre/code/ADCPtools/utils');

PRES_THRESH = 10;
AMP_THRESH = 30;
COR_THRESH = 80;
SPD_THRESH = 1.0;
ECHODIFF_FISH = 30; % 40
threshsd = 2;

theta = 25;
instrument_height = 0.7; % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select mooring and deployment number %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mooring = 'OC40S-A';
deployment = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deployment = num2str(deployment);
disp(['Mooring: ' mooring ', deployment #', deployment]);
disp('===============================');
disp(' ');
fdir = ['/home/andre/phd/data/innershelfdri_moorings/' mooring '/deployment' deployment '_raw/'];

mm2m = 1e-3;
fnames = dir([fdir,'*_000*']);
% fnames = fnames(18);
nfiles = length(fnames);

i = 1;
A = struct();
A.t = []; A.p = []; A.hdng = []; A.ptch = []; A.roll = [];
A.b1 = []; A.b2 = []; A.b3 = []; A.b4 = []; A.b5 = [];
A.amp1 = []; A.amp2 = []; A.amp3 = []; A.amp4 = []; A.amp5 = [];
% A.cor1 = []; A.cor2 = []; A.cor3 = []; A.cor4 = [];

while i<=nfiles
	fname=fnames(i).name;
	disp(['File number ' int2str(i) ' of ' int2str(nfiles)])
	[adcp] = rdradcpJmkFast_5beam([fdir fname]);

	adcp.vel = adcp.vel*mm2m;
	adcp.vel5 = adcp.vel5*mm2m;

	if i==1
		r = squeeze(adcp.depths);
		A.r = r;
		A.z = r*cosd(theta) + instrument_height;
	end

  if strcmp(mooring, 'OC25SB-A')
    adcp.time = adcp.time + datenum(2000,1,1,0,0,0) - 1;
  end

  if strcmp(mooring, 'OC25SB-A') && strcmp(deployment, '2')
    if adcp.time(end)<datenum('2017-10-20 08:00:00')
      disp('Skipping file - part of the record with excessive tilt (OC25SB-A, deployment 2).')
      disp(' ')

      i = i + 1;
      continue
    end
  end

	% Threshold despiking.
  if ~strcmp(mooring, 'OC40N-A') % Pressure sensor in OC40N-A seems to have gone bad after few days.
  	ib=find(adcp.pr<PRES_THRESH); adcp.vel(:,:,ib)=NaN; adcp.vel5(:,ib)=NaN;
  end
	ib=find(adcp.cor<COR_THRESH); adcp.vel(ib)=NaN;
	ib=find(adcp.int<AMP_THRESH); adcp.vel(ib)=NaN;
	ib=find(adcp.int5<AMP_THRESH); adcp.vel5(ib)=NaN;

	ib=find(abs(adcp.vel)>SPD_THRESH); adcp.vel(ib)=NaN;
	ib=find(abs(adcp.vel5)>SPD_THRESH); adcp.vel5(ib)=NaN;

	A.b1 = [A.b1, squeeze(adcp.vel(1,:,:))];
	A.b2 = [A.b2, squeeze(adcp.vel(2,:,:))];
	A.b3 = [A.b3, squeeze(adcp.vel(3,:,:))];
	A.b4 = [A.b4, squeeze(adcp.vel(4,:,:))];
	A.b5 = [A.b5, adcp.vel5];

	A.amp1 = [A.amp1, squeeze(adcp.int(1,:,:))];
	A.amp2 = [A.amp2, squeeze(adcp.int(2,:,:))];
	A.amp3 = [A.amp3, squeeze(adcp.int(3,:,:))];
	A.amp4 = [A.amp4, squeeze(adcp.int(4,:,:))];
	A.amp5 = [A.amp5, adcp.int5];

  % A.cor1 = [A.cor1, squeeze(adcp.cor(1,:,:))];
	% A.cor2 = [A.cor2, squeeze(adcp.cor(2,:,:))];
	% A.cor3 = [A.cor3, squeeze(adcp.cor(3,:,:))];
	% A.cor4 = [A.cor4, squeeze(adcp.cor(4,:,:))];

	offset_pitchroll = 655.36;
	ib=find(adcp.pitch<-100);
	adcp.pitch(ib) = adcp.pitch(ib) + offset_pitchroll;
	ib=find(adcp.roll<-100);
	adcp.roll(ib) = adcp.roll(ib) + offset_pitchroll;

  A.t = [A.t adcp.time];
	A.hdng = [A.hdng adcp.heading];
	A.ptch = [A.ptch adcp.pitch];
	A.roll = [A.roll adcp.roll];
	A.p = [A.p adcp.pr];

	i = i + 1;
end % while

% Rescale beam amplitudes for fish detection algorithm to work right on all beams.
% amp1 = A.amp1 - nanmedian(A.amp1(:));
% amp2 = A.amp2 - nanmedian(A.amp2(:));
% amp3 = A.amp3 - nanmedian(A.amp3(:));
% amp4 = A.amp4 - nanmedian(A.amp4(:));
% amp5 = A.amp5 - nanmedian(A.amp5(:));
% amp1 = A.amp1;
% amp2 = A.amp2;
% amp3 = A.amp3;
% amp4 = A.amp4;
% amp5 = A.amp5;

% Cap off near-surface sidelobe contamination (plus an offset) following pressure.
if ~strcmp(mooring, 'OC40N-A') % Pressure sensor in OC40N-A seems to have gone bad after few days.
  Roffset = 2.0; % 0; % Extra cap.
  costh = cosd(theta);
  for k=1:numel(A.t)
  	H = A.p(k);
  	h = H*costh - Roffset;
  	ib = find(A.r>=h);
  	A.b1(ib,k) = NaN;
  	A.b2(ib,k) = NaN;
  	A.b3(ib,k) = NaN;
  	A.b4(ib,k) = NaN;
  	ib = find(A.r>=(H - Roffset));
  	A.b5(ib,k) = NaN;
  end
else
  H = 39;
  Roffset = 2.0; % 0;
  costh = cosd(theta);
  h = H*costh - Roffset;
  ib = find(A.r>=h);
  A.b1(ib,:) = NaN;
  A.b2(ib,:) = NaN;
  A.b3(ib,:) = NaN;
  A.b4(ib,:) = NaN;
  ib = find(A.r>=(H - Roffset));
  A.b5(ib,:) = NaN;
end

% [A.b1, A.b2, A.b3, A.b4, A.b5] = mskfish5(A.b1, A.b2, A.b3, A.b4, A.b5, amp1, amp2, amp3, amp4, amp5, 'Threshold', ECHODIFF_FISH);

if ~strcmp(mooring, 'OC25SB-A')
  A.t = A.t + datenum(2000,1,1,0,0,0) - 1; % Fix time vector.
end

warning off
disp('Despiking beam 1.');
A.b1 = despike_filt(A.b1, threshsd);
A.b1 = despike_filt(A.b1', threshsd)';
disp('Despiking beam 2.');
A.b2 = despike_filt(A.b2, threshsd);
A.b2 = despike_filt(A.b2', threshsd)';
disp('Despiking beam 3.');
A.b3 = despike_filt(A.b3, threshsd);
A.b3 = despike_filt(A.b3', threshsd)';
disp('Despiking beam 4.');
A.b4 = despike_filt(A.b4, threshsd);
A.b4 = despike_filt(A.b4', threshsd)';
disp('Despiking beam 5.');
A.b5 = despike_filt(A.b5, threshsd);
A.b5 = despike_filt(A.b5', threshsd)';

% Correct heading for magnetic declination.
% https://www.ngdc.noaa.gov/geomag-web/#declination
% Queried magnetic declination at 35.0429,-120.6680 on 2017-sep-28. Model: WMM.
% Result: 12.70 +- 0.33 degrees E.
mag_decl = 12.70;
A.hdng = A.hdng + mag_decl;

% save file in beam coordinates.
fout = ['/home/andre/phd/data/innershelfdri_moorings/' mooring '/deployment' deployment '/' mooring(1:end-2) '-Ad' deployment 'beam'];

save(fout, 'A', '-v7.3');

% r = A.r;
% t = A.t;
% b1 = A.b1;
% b2 = A.b2;
% b3 = A.b3;
% b4 = A.b4;
% b5 = A.b5;
% amp1 = A.amp1;
% amp2 = A.amp2;
% amp3 = A.amp3;
% amp4 = A.amp4;
% amp5 = A.amp5;
% cor1 = A.cor1;
% cor2 = A.cor2;
% cor3 = A.cor3;
% cor4 = A.cor4;
% clear A;
%
% figure;
% pltbeamsvel(t, r, b1, b2, b3, b4, b5);
% figure;
% pltbeamsamp(t, r, amp1, amp2, amp3, amp4, amp5);
% figure;
% pltbeamscor(t, r, cor1, cor2, cor3, cor4, cor4*NaN);
