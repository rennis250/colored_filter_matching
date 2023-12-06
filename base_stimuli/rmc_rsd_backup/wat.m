clear all; close all; clc;

fnsx = dir('*.xml');
fnse = dir('*.exr');

for fnc = 1:length(fnse)
 fne = fnse(fnc).name;

 for fxc = 1:length(fnsx)
  fnx = fnsx(fxc).name;

  fne_oxt = fne(1:end-4);
  fnx_oxt = fnx(1:end-4);
  if strcmp(fne_oxt, fnx_oxt)
       fne_oxt
       fnx_oxt
   movefile(fne, ['./finished/' fne])
   movefile(fnx, ['./finished/' fnx])
  end
 end
end
