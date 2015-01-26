%  This function was adapted from the following software package:
%
% -----rbMIT Software Copyright MIT 2006-07----
%  DBP  Huynh,NC Nguyen, AT Patera, G Rozza
% Henceforth, Software shall refer to the rbMIT_System software package: 
% the contents of the rbMIT_Library and rbMIT_Aux folders.
% 
% This License governs use of all accompanying Software, and your download and/or use of the Software constitutes acceptance 
% of this License.
% 
% You may use this Software for any non-commercial purpose, subject to the restrictions in this license. 
% Some purposes which can be non-commercial are teaching and academic research.
% 
% You may not use or distribute this Software or any derivative works in any form for commercial purposes. 
% Examples of commercial purposes would be licensing, leasing, or selling the Software, or distributing the 
% Software for use with commercial products.
% 
% You may modify this Software and distribute the modified Software for non-commercial purposes,
%  however, you may not grant rights to the Software or derivative works that are broader than those 
% provided by this License. For example, you may not distribute modifications of the Software under terms 
% that would permit commercial use.
% 
% In return, we require that you agree:
% 
% not to remove any copyright or other notices from the Software;
% 
% that if you distribute the Software in source or object form, you will include a verbatim copy of this license;
% 
% that if you distribute derivative works of the Software in source code form you do so only under a license
% that includes all of the provisions of this License, and if you distribute derivative works of the Software solely 
% in object form you do so only under a license that complies with this License;
% 
% that if you have modified the Software or created derivative works, and distribute such modifications or derivative
% works, you will cause the modified files to carry prominent notices so that recipients know that they are not receiving the original Software;
% 
% THAT THE SOFTWARE COMES "AS IS", WITH NO WARRANTIES: THIS MEANS NO EXPRESS, IMPLIED OR STATUTORY WARRANTY, 
% INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY WARRANTY 
% OF TITLE OR NON-INFRINGEMENT; YOU MUST PASS THIS DISCLAIMER ON WHENEVER YOU DISTRIBUTE THE SOFTWARE OR DERIVATIVE WORKS;
% 
% THAT NEITHER MIT NOR THE AUTHORS WILL BE LIABLE FOR ANY DAMAGES RELATED TO THE SOFTWARE OR THIS LICENSE, 
% INCLUDING DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL OR INCIDENTAL DAMAGES; YOU MUST PASS THIS LIMITATION OF LIABILITY 
% ON WHENEVER YOU DISTRIBUTE THE SOFTWARE OR DERIVATIVE WORKS;
% 
% that your rights under the License end automatically if you breach the License in any way;
% 
% that MIT reserves all rights not expressly granted to you in this License.
% 
% (Note: The wording of this license agreement is derived from a Microsoft Shared Source license agreement.)
% 
% License at http://augustine.mit.edu/rbMIT_SystemLicense.htm.
%
% rbMIT Software Part_II_and_IV
%
% ---------- 
% -----rbMIT_System Software Copyright MIT 2006-07-----
% ----- Geometry file conv_diff_pdedom.m -----

function [x,y] = conv_diff_pdedom(bs,s)

nbs = 17;

if nargin==0
    x=nbs;
    return
end

mu1 = 0.100000;

psym = [-1, 0; 0, -1; 0, -8; -8, -8; -8, 8; 0, 8; 0, 1; -1825444424271651/2251799813685248*2^(1/2), -1825444424271651/2251799813685248*2^(1/2); -1825444424271651/2251799813685248*2^(1/2), 1825444424271651/2251799813685248*2^(1/2)];

dl(:,1) = [0.000000; 1.000000; 0; 7];
dl(:,2) = [0.000000; 1.000000; 3; 0];
dl(:,3) = [0.000000; 1.000000; 7; 2];
dl(:,4) = [0.000000; 1.000000; 2; 3];
dl(:,5) = [0.000000; 1.000000; 0; 6];
dl(:,6) = [0.000000; 1.000000; 6; 7];
dl(:,7) = [0.000000; 1.000000; 0; 1];
dl(:,8) = [0.000000; 1.000000; 1; 6];
dl(:,9) = [0.000000; 1.000000; 0; 9];
dl(:,10) = [0.000000; 1.000000; 8; 1];
dl(:,11) = [0.000000; 1.000000; 9; 8];
dl(:,12) = [0.000000; 1.000000; 0; 4];
dl(:,13) = [0.000000; 1.000000; 4; 9];
dl(:,14) = [0.000000; 1.000000; 0; 5];
dl(:,15) = [0.000000; 1.000000; 5; 4];
dl(:,16) = [0.000000; 1.000000; 3; 5];
dl(:,17) = [0.000000; 1.000000; 8; 2];

bs1=bs(:)';

if nargin==1
    x=dl(:,bs1);
    return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1,
    bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2),
    error('PDE:conv_diff_pdedom:SizeBs', 'bs must be scalar or of same size as s.');
end

if ~isempty(s)
   % boundary segment 1
   ii=find(bs==1);
   if length(ii)
       for iii = 1:length(ii)
           t = s(ii(iii))*(1/2*pi)+(pi);
           x(ii(iii)) = cos(t);
           y(ii(iii)) = sin(t);
            if (s(ii(iii)) == 0)
                x(ii(iii)) = psym(1,1);
                y(ii(iii)) = psym(1,2);
            end
            if (s(ii(iii)) == 1)
                x(ii(iii)) = psym(2,1);
                y(ii(iii)) = psym(2,2);
            end
       end
   end
   % boundary segment 2
   ii=find(bs==2);
   if length(ii)
       for iii = 1:length(ii)
           t = s(ii(iii))*(-1/2*pi)+(pi);
           x(ii(iii)) = cos(t);
           y(ii(iii)) = sin(t);
            if (s(ii(iii)) == 0)
                x(ii(iii)) = psym(1,1);
                y(ii(iii)) = psym(1,2);
            end
            if (s(ii(iii)) == 1)
                x(ii(iii)) = psym(7,1);
                y(ii(iii)) = psym(7,2);
            end
       end
   end
   % boundary segment 3
   ii=find(bs==3);
   if length(ii)
       t = s(ii);
       x(ii) = psym(1,1) + t*(psym(8,1)-psym(1,1));
       y(ii) = psym(1,2) + t*(psym(8,2)-psym(1,2));
   end
   % boundary segment 4
   ii=find(bs==4);
   if length(ii)
       t = s(ii);
       x(ii) = psym(1,1) + t*(psym(9,1)-psym(1,1));
       y(ii) = psym(1,2) + t*(psym(9,2)-psym(1,2));
   end
   % boundary segment 5
   ii=find(bs==5);
   if length(ii)
       t = s(ii);
       x(ii) = psym(2,1) + t*(psym(3,1)-psym(2,1));
       y(ii) = psym(2,2) + t*(psym(3,2)-psym(2,2));
   end
   % boundary segment 6
   ii=find(bs==6);
   if length(ii)
       t = s(ii);
       x(ii) = psym(2,1) + t*(psym(8,1)-psym(2,1));
       y(ii) = psym(2,2) + t*(psym(8,2)-psym(2,2));
   end
   % boundary segment 7
   ii=find(bs==7);
   if length(ii)
       t = s(ii);
       x(ii) = psym(3,1) + t*(psym(4,1)-psym(3,1));
       y(ii) = psym(3,2) + t*(psym(4,2)-psym(3,2));
   end
   % boundary segment 8
   ii=find(bs==8);
   if length(ii)
       t = s(ii);
       x(ii) = psym(3,1) + t*(psym(8,1)-psym(3,1));
       y(ii) = psym(3,2) + t*(psym(8,2)-psym(3,2));
   end
   % boundary segment 9
   ii=find(bs==9);
   if length(ii)
       t = s(ii);
       x(ii) = psym(4,1) + t*(psym(5,1)-psym(4,1));
       y(ii) = psym(4,2) + t*(psym(5,2)-psym(4,2));
   end
   % boundary segment 10
   ii=find(bs==10);
   if length(ii)
       t = s(ii);
       x(ii) = psym(4,1) + t*(psym(8,1)-psym(4,1));
       y(ii) = psym(4,2) + t*(psym(8,2)-psym(4,2));
   end
   % boundary segment 11
   ii=find(bs==11);
   if length(ii)
       t = s(ii);
       x(ii) = psym(4,1) + t*(psym(9,1)-psym(4,1));
       y(ii) = psym(4,2) + t*(psym(9,2)-psym(4,2));
   end
   % boundary segment 12
   ii=find(bs==12);
   if length(ii)
       t = s(ii);
       x(ii) = psym(5,1) + t*(psym(6,1)-psym(5,1));
       y(ii) = psym(5,2) + t*(psym(6,2)-psym(5,2));
   end
   % boundary segment 13
   ii=find(bs==13);
   if length(ii)
       t = s(ii);
       x(ii) = psym(5,1) + t*(psym(9,1)-psym(5,1));
       y(ii) = psym(5,2) + t*(psym(9,2)-psym(5,2));
   end
   % boundary segment 14
   ii=find(bs==14);
   if length(ii)
       t = s(ii);
       x(ii) = psym(6,1) + t*(psym(7,1)-psym(6,1));
       y(ii) = psym(6,2) + t*(psym(7,2)-psym(6,2));
   end
   % boundary segment 15
   ii=find(bs==15);
   if length(ii)
       t = s(ii);
       x(ii) = psym(6,1) + t*(psym(9,1)-psym(6,1));
       y(ii) = psym(6,2) + t*(psym(9,2)-psym(6,2));
   end
   % boundary segment 16
   ii=find(bs==16);
   if length(ii)
       t = s(ii);
       x(ii) = psym(7,1) + t*(psym(9,1)-psym(7,1));
       y(ii) = psym(7,2) + t*(psym(9,2)-psym(7,2));
   end
   % boundary segment 17
   ii=find(bs==17);
   if length(ii)
       t = s(ii);
       x(ii) = psym(8,1) + t*(psym(9,1)-psym(8,1));
       y(ii) = psym(8,2) + t*(psym(9,2)-psym(8,2));
   end
end

