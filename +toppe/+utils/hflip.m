% newarr = hflip(arr [,cp]);
%
% Flips the array about the horizontal axis.
%
% If cp is given, the array is flipped about the point 
% cp, where 2*cp must be an integer.
%
% $Id: hflip.m,v 1.1 2018/10/26 01:38:45 jfnielse Exp $


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: hflip.m,v $
%	Revision 1.1  2018/10/26 01:38:45  jfnielse
%
%	: Committing in .
%	:
%	: Modified Files:
%	: 	+utils/makebalanced.m +utils/trapwave.m
%	: Added Files:
%	: 	+utils/hflip.m
%
%	Revision 1.1  2016/03/22 18:47:26  jfnielse
%
%	: Committing in .
%	:
%	: Modified Files:
%	: 	bipolar3.m
%	: Added Files:
%	: 	bipolarpar.m hflip.m
%	: ----------------------------------------------------------------------
%
%	Revision 1.3  2013/09/04 16:31:24  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	hflip.m makecrusher.m makeevenlength.m mat2wav.m mybridged.m
%	: 	myslrrf.m plotrf.m plotwav.m readwav.m trapwave.m
%
%	Revision 1.2  2013/09/04 16:22:04  jfnielse
%
%	: Committing in .
%	:
%	: Modified Files:
%	: 	hflip.m
%
%	Revision 1.1  2013/09/04 16:12:03  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	hflip.m makecrusher.m makeevenlength.m mybridged.m trapwave.m
%
%	Revision 1.1  2013/08/23 17:46:03  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	hflip.m makedess.m makeevenlength.m myrfstat.m readwav.m
%	: 	trapwave.m writewav.m
%
%	Revision 1.1  2010/03/15 02:18:35  jfnielse
%	
%	: Committing in .
%	:
%	: Added Files:
%	: 	README hflip.m make2dft.m makeEPI.m makecrusher.m makegz.m
%	: 	readjfnwav.m trapwave.m writejfnwav.m
%	
%	Revision 1.1  2008/09/21 14:37:51  jfnielse
%	
%	: Modified Files:
%	: 	makeEPI.m
%	: Added Files:
%	: 	hflip.m trapwave.m
%	
%	Revision 1.1  2004/04/05 16:20:35  knayak
%	initial Matlab Directory
%	
%	Revision 1.1  2002/03/28 01:35:30  bah
%	Added to CVS
%	
%
% ===========================================================


function newarr = hflip(arr, cp);

import toppe.utils.*
import toppe.utils.spiral.*

s = size(arr);
s2 = s(2);

if (nargin < 2)
	cp = (s2+1)/2;
end;

cp = cp+s2;
temp = [arr, arr, arr];

newarr = 0*arr;

for k = 1:s2
  newarr(:,k) = temp(:,2*cp-k-s2);
end;






