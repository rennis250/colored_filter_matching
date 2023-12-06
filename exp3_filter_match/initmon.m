function initmon(INITDATA)
%INITMON  Initializes DKL<-->RGB conversion matrices
% OBJECTIVE
%   INITMON(INITDATA) initializes DKL<-->RGB conversion matrices for the
%   CIE xyY coordinates of monitor phosphors given in INITDATA. 
% PARTICULARITY
%   Function INITMON must be called once prior to subsequent calls of
%   conversion routines DKL2RGB or RGB2DKL.
%   See also DKL2RGB, RGB2DKL.
% INPUT
%   INITDATA inputs the 3x3 matrix of xyY values for each phosphor. 
%   format:
%       rx ry rY
%       gx gy gY
%       bx by bY
%   May be specified in 2 ways:
%   a.) As the 3x3 matrix.
%   b.) As a filename for a file with the 3x3 matrix. 
%   DEFAULT (nargin==0): INITMON without arguments uses default xyY values
%   to initialize the conversion matrices.
% EXAMPLE
%   Run without input.
% GENEALOGY 
%   2004 Thorsten Hansen.
%   2009mar17 enable INPUT of xyYmon as 3x3 matrix; tidied up; removed
%   colormap for LUTs (no use?) [cw]. 
%   2009jun23 added default for "isempty(INITDATA)". [cw]

% define xyY monitor coordinates
if nargin == 0 | isempty(INITDATA) % use default calibration files from a sample monitor
  moncie = samplemonitor;
  disp('Initialize conversion matrices from default values.');
else
  if ischar(INITDATA)
      moncie = textread(INITDATA); % Original (th) method.      
      disp(['Initialize conversion matrices from ''' INITDATA '''.']);
  elseif isnumeric(INITDATA) % New (cw) Method.
      moncie = INITDATA;
      disp(['Initialize conversion matrices from INPUT MATRIX.']);
  end
  if sum(size(moncie) ~= [3 3])
    error('INPUT ERROR: INITDATA must be 3x3.')
  end
end


% initialize conversion matrices M_dkl2rgb and M_rgb2dkl 
% from monitor coordinates moncie
global M_dkl2rgb M_rgb2dkl


M_dkl2rgb = getdkl(moncie);
M_rgb2dkl = inv(M_dkl2rgb);

% initialize conversion matrices M_rgb2lms and M_lms2rgb
%
global M_rgb2lms M_lms2rgb

% TEST: new vectorized implementation
monxyY = moncie; 
x = monxyY(:,1); y = monxyY(:,2); Y = monxyY(:,3);

% xyY -> xyz conversion
if prod(y) == 0, error('y column contains zero value.'), end
z = 1-x-y;
monxyz = [x y z];

white = Y/2;

% xyY -> XYZ conversion
X = x./y.*Y;
Z = z./y.*Y;
monXYZ = [X Y Z]; % this should be monCIE
% end TEST

monCIE = zeros(3,3);
monCIE(:,2) = moncie(:,3);

for i=1:3
  moncie(i,3) = 1.0 - moncie(i,1) - moncie(i,2);
  monCIE(i,1) = (moncie(i,1)/moncie(i,2))*monCIE(i,2);
  monCIE(i,3) = (moncie(i,3)/moncie(i,2))*monCIE(i,2);
  monRGB(i,1) = 0.15514 * monCIE(i,1) + ...
      0.54313 * monCIE(i,2) - 0.03386 * monCIE(i,3);
  monRGB(i,2) = -0.15514 * monCIE(i,1) + ...
      0.45684 * monCIE(i,2) + 0.03386 * monCIE(i,3);
  monRGB(i,3) = 0.01608 * monCIE(i,3);
  tsum = monRGB(i,1) + monRGB(i,2);
  monrgb(i,1) = monRGB(i,1) / tsum;
  monrgb(i,2) = monRGB(i,2) / tsum;
  monrgb(i,3) = monRGB(i,3) / tsum;
end

Xmon = monCIE; % who needs Xmon?

Xmat= inv(Xmon); % who needs Xmat?

M_rgb2lms = monRGB; % M_rgb2lms used in mon2cones
                    % why not directly compute on M_rgb2lms ?

M_lms2rgb = inv(M_rgb2lms);

% white point
%w = mon2cones(0.5, 0.5, 0.5);
w = [0.5 0.5 0.5] * M_rgb2lms; % who needs this?

%% *************************** SUBFUNCTIONS *******************************

%% M_dkl2rgb (determine global variable
function M_dkl2rgb = getdkl(monxyY) 
% compute dkl2rgb conversion matrix from moncie coordinates
% (compare function "getdkl" in color.c)

x = monxyY(:,1); y = monxyY(:,2); Y = monxyY(:,3);
if prod(y) == 0, error('y column contains zero value.'), end
xyz = [x y 1-x-y];
white = Y/2;

% Smith & Pokorny cone fundamentals 
% V. C. Smith & J. Pokorny (1975), Vision Res. 15, 161-172.
%        X          Y       Z [cw] 
M = [ 0.15514  0.54312  -0.03286    % L alias R 
     -0.15514  0.45684   0.03286    % M alias G
      0.0      0.0       0.01608];  % S alias B

RGB = xyz*M'; % R, G  and B cones (i.e, long, middle and short wavelength)
% More precisely: RGB = excitation of the L, M and S-cones by the monitor
% primaries [cw]
RG_sum = RGB(:,1) + RGB(:,2); % R G sum (L+M)
R = RGB(:,1)./RG_sum;
B = RGB(:,3)./RG_sum; 
G = 1 - R;

% constant blue axis
a = white(1)*B(1);
b = white(1)*(R(1)+G(1));
c = B(2);
d = B(3);
e = R(2)+G(2);
f = R(3)+G(3);
dGcb = (a*f/d - b)/(c*f/d - e); % solve x
dBcb = (a*e/c - b)/(d*e/c - f); % solve y

% tritanopic confusion axis
a = white(3)*R(3); % L-Cone excitation of the luminance component of the blue primary
b = white(3)*G(3); % M-Cone excitation ~
c = R(1);
d = R(2);
e = G(1);
f = G(2);
dRtc = (a*f/d - b)/(c*f/d - e); % solve x
dGtc = (a*e/c - b)/(d*e/c - f); % solve y

IMAX = 1;
M_dkl2rgb = IMAX * [1        1         dRtc/white(1)
                    1  -dGcb/white(2)  dGtc/white(2)
                    1  -dBcb/white(3)     -1]; 

%% samplemonitor
function xyYmon = samplemonitor
% INITDATA 'sony-th.xyY' from:
% xyY = textread('/home/hansen/bib/data/lut/sony-th.xyY'); 
xyYmon = [...
    0.6130 0.3489 20.2888;
    0.2829 0.6054 64.0547;
    0.1565 0.0709 8.6309;];